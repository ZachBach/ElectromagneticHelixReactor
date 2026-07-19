/* Phase 3C — sustained-discharge simulator.

   The decay-mode studies (run-ambipolar.js) answered "does the pitch effect
   survive ambipolar fields?" — mostly no, for bulk confinement. This runner
   answers the question those results point to:

     If pitch doesn't strongly change confinement, does it change the
     EQUILIBRIUM density and potential structure?

   A constant volumetric source injects electron-ion PAIRS (charge-neutral by
   construction) at --rate super-pairs/us with Maxwellian velocities (Te for
   electrons, Tgas for ions). Both species evolve as in run-ambipolar.js
   (Boris + MCC + quasi-static Poisson feedback, mobile ions, reduced ion
   mass). The plasma fills until losses balance the source, then statistics
   are taken over the second half of the run:

     - steady inventory <Ne>, <Ni>  →  tau_eff = <N>/S (no censoring —
       confinement expressed as equilibrium inventory)
     - time-averaged phi(r,z), ne(r,z), ni(r,z) → the equilibrium structure
       (written as a fields CSV; analyze with analyze-fields.js)
     - wall-flux balance Gamma_i/Gamma_e, which steady state must drive → 1;
       reported as the equilibrium check, with per-channel splits.

   deposit()/interpE() are duplicated from run-ambipolar.js deliberately (the
   two runners stay independently runnable; consolidate into a shared module
   if a third PIC-style runner appears).

   Usage: node sim/run-sustained.js [--rate=200] [--ttotal=6e-5] [--n0seed=800]
     [--mion=1] [--ratios=0,1,3] [--model=screw|beltrami|powerlaw] [--nexp=1]
     [--bwall=0.01] [--pressure=10] [--te=3] [--r=0.1] [--l=0.4] [--seed=1]
     [--nr=48] [--nz=96] [--fieldevery=<steps>] [--cap=24000]               */
'use strict';

const fs = require('fs');
const path = require('path');
const C = require('./constants');
const { makeField } = require('./field');
const { borisRotate } = require('./boris');
const coll = require('./collisions');
const { KIND } = require('./diagnostics');
const { PoissonRZ } = require('./poisson');

const args = {};
for (const a of process.argv.slice(2)) {
  const m = /^--([^=]+)(?:=(.*))?$/.exec(a);
  if (m) args[m[1]] = m[2] === undefined ? true : m[2];
}
const cfg = {
  rate: parseFloat(args.rate || '200'),        // super-pairs per microsecond
  Ttotal: parseFloat(args.ttotal || '6e-5'),
  n0seed: parseInt(args.n0seed || '800', 10),  // pre-seeded pairs (shortens transient)
  mIonAmu: parseFloat(args.mion || '1'),
  ratios: (args.ratios || '0,1,3').split(',').map(Number),
  model: args.model || 'screw',
  nExp: parseFloat(args.nexp || '1'),
  Bwall: parseFloat(args.bwall || '0.01'),
  pmTorr: parseFloat(args.pressure || '10'),
  TeV: parseFloat(args.te || '3'),
  Tgas: 300,
  R: parseFloat(args.r || '0.1'),
  L: parseFloat(args.l || '0.4'),
  seed: parseInt(args.seed || '1', 10),
  Nr: parseInt(args.nr || '48', 10),
  Nz: parseInt(args.nz || '96', 10),
  cap: parseInt(args.cap || '24000', 10),
};

const mI = cfg.mIonAmu * C.AMU;
const vth = Math.sqrt((cfg.TeV * C.EV) / C.ME);
const vthI = Math.sqrt((C.KB * cfg.Tgas) / mI);
const nn = C.neutralDensity(cfg.pmTorr, cfg.Tgas);
const qOverM = -C.QE / C.ME;
const qOverMI = C.QE / mI;
const cs = Math.sqrt((cfg.TeV * C.EV) / mI);
const wcMax = Math.abs(qOverM) * cfg.Bwall * (cfg.model === 'beltrami' ? 1.5 : 1);
const dt = Math.min(0.2 / wcMax, 0.05 / coll.electronMaxRate(nn), 0.05 / coll.ionMaxRate(nn, 5 * cs));
const h = (qOverM * dt) / 2, hI = (qOverMI * dt) / 2;
const HiStar = (qOverMI * cfg.Bwall) / (nn * (C.SIGMA_CX + C.SIGMA_I_EL) * cs);
cfg.fieldEvery = args.fieldevery ? parseInt(args.fieldevery, 10) : Math.max(1, Math.round(5.7e-10 / dt));

// super-particle weight fixed at the Phase 3 space-charge scale so the
// quasi-static validity constraints carry over unchanged
const Vseed = Math.PI * (0.9 * cfg.R) ** 2 * (0.8 * cfg.L);
const wSuper = (3e11 * Vseed) / 2000;
const Ssuper = cfg.rate * 1e6; // super-pairs per second

const Nr = cfg.Nr, Nz = cfg.Nz, dr = cfg.R / Nr, dz = cfg.L / Nz;

console.log(`EHR Phase 3C sustained discharge — ${cfg.model}` +
  `${cfg.model === 'powerlaw' ? `(n=${cfg.nExp})` : ''} field, ` +
  `m_ion = ${cfg.mIonAmu} amu, H_i* = ${HiStar.toFixed(2)} at c_s`);
console.log(`  |B(R)| = ${(cfg.Bwall * 1e4).toFixed(0)} G, p = ${cfg.pmTorr} mTorr, ` +
  `Te = ${cfg.TeV} eV, source = ${cfg.rate}/us into the seed region, grid ${Nr}x${Nz}`);
console.log(`  Ttotal = ${cfg.Ttotal} s (measure over 2nd half), dt = ${dt.toExponential(2)} s, ` +
  `field update every ${cfg.fieldEvery} steps, cap = ${cfg.cap}`);
console.log('');
console.log('  ratio | <Ne>  <Ni>  | tau_eff us | phiMax V | Gi/Ge | e-loss end%/rad% |    s');
console.log('  ' + '-'.repeat(84));

function depositInto(grid, pos, alive, nAlive) {
  grid.fill(0);
  for (let k = 0; k < nAlive; k++) {
    const i3 = 3 * alive[k];
    const x = pos[i3], y = pos[i3 + 1], z = pos[i3 + 2];
    const r = Math.sqrt(x * x + y * y);
    const rf = r / dr - 0.5, zf = z / dz - 0.5;
    const i0 = Math.floor(rf), j0 = Math.floor(zf);
    const fr = rf - i0, fz = zf - j0;
    const iA = i0 < 0 ? 0 : i0;
    const iB = i0 + 1 >= Nr ? Nr - 1 : i0 + 1;
    const jA = j0 < 0 ? 0 : j0 >= Nz ? Nz - 1 : j0;
    const jB = j0 + 1 < 0 ? 0 : j0 + 1 >= Nz ? Nz - 1 : j0 + 1;
    grid[jA * Nr + iA] += (1 - fr) * (1 - fz);
    grid[jB * Nr + iA] += (1 - fr) * fz;
    grid[jA * Nr + iB] += fr * (1 - fz);
    grid[jB * Nr + iB] += fr * fz;
  }
  for (let j = 0; j < Nz; j++) {
    for (let i = 0; i < Nr; i++) {
      const vol = 2 * Math.PI * (i + 0.5) * dr * dr * dz;
      grid[j * Nr + i] *= wSuper / vol;
    }
  }
}

function interpE(Er, Ez, x, y, z, out) {
  const r = Math.sqrt(x * x + y * y);
  const rf = r / dr - 0.5, zf = z / dz - 0.5;
  const i0 = Math.floor(rf), j0 = Math.floor(zf);
  const fr = rf - i0, fz = zf - j0;
  let iA = i0, sA = 1;
  if (iA < 0) { iA = 0; sA = -1; } else if (iA >= Nr) iA = Nr - 1;
  let iB = i0 + 1, sB = 1;
  if (iB < 0) { iB = 0; sB = -1; } else if (iB >= Nr) iB = Nr - 1;
  const jA = j0 < 0 ? 0 : j0 >= Nz ? Nz - 1 : j0;
  const jB = j0 + 1 < 0 ? 0 : j0 + 1 >= Nz ? Nz - 1 : j0 + 1;
  const wAA = (1 - fr) * (1 - fz), wAB = (1 - fr) * fz;
  const wBA = fr * (1 - fz), wBB = fr * fz;
  const cAA = jA * Nr + iA, cAB = jB * Nr + iA, cBA = jA * Nr + iB, cBB = jB * Nr + iB;
  const er = wAA * sA * Er[cAA] + wAB * sA * Er[cAB] + wBA * sB * Er[cBA] + wBB * sB * Er[cBB];
  const ez = wAA * Ez[cAA] + wAB * Ez[cAB] + wBA * Ez[cBA] + wBB * Ez[cBB];
  if (r > 1e-12) { out[0] = (er * x) / r; out[1] = (er * y) / r; } else { out[0] = 0; out[1] = 0; }
  out[2] = ez;
}

const outDir = path.join(__dirname, 'output');
fs.mkdirSync(outDir, { recursive: true });
const stamp = new Date().toISOString().replace(/[:T]/g, '-').slice(0, 19);
const runTag = `${cfg.model}${cfg.model === 'powerlaw' ? '-n' + cfg.nExp : ''}-${(cfg.Bwall * 1e4).toFixed(0)}G-${cfg.pmTorr}mT`;
const summary = [
  `# EHR Phase 3C sustained discharge  ${new Date().toISOString()}`,
  `# model=${cfg.model} nexp=${cfg.nExp} mion_amu=${cfg.mIonAmu} HiStar=${HiStar} Bwall_T=${cfg.Bwall} p_mTorr=${cfg.pmTorr}`,
  `# rate_per_us=${cfg.rate} Ttotal_s=${cfg.Ttotal} Tseed_eV=${cfg.TeV} R_m=${cfg.R} L_m=${cfg.L} dt_s=${dt} grid=${Nr}x${Nz} seed=${cfg.seed}`,
  'ratio,NeSS,NiSS,tauEff_s,phiMax_V,GiOverGe,fracEndE,fracRadialE,capHits',
];

const R2 = cfg.R * cfg.R;
const B = new Float64Array(3), Evec = new Float64Array(3);

for (let ri = 0; ri < cfg.ratios.length; ri++) {
  const ratio = cfg.ratios[ri];
  const t0Wall = Date.now();
  const rng = coll.makeRng(cfg.seed * 7919 + ri * 104729 + 1);
  const gauss = coll.makeGauss(rng);
  const field = makeField({ model: cfg.model, Bwall: cfg.Bwall, ratio, R: cfg.R, nExp: cfg.nExp });
  const eColl = coll.makeElectronCollider(nn, dt, rng);
  const iColl = coll.makeIonCollider(nn, cfg.Tgas, dt, rng, gauss, mI);
  const sor = new PoissonRZ(Nr, Nz, cfg.R, cfg.L);
  const ne = new Float64Array(Nr * Nz), ni = new Float64Array(Nr * Nz), S = new Float64Array(Nr * Nz);
  const phiAvg = new Float64Array(Nr * Nz), neAvg = new Float64Array(Nr * Nz), niAvg = new Float64Array(Nr * Nz);

  const cap = cfg.cap;
  const pos = new Float64Array(3 * cap), vel = new Float64Array(3 * cap);
  const posI = new Float64Array(3 * cap), velI = new Float64Array(3 * cap);
  const alive = new Int32Array(cap), aliveI = new Int32Array(cap);
  for (let i = 0; i < cap; i++) { alive[i] = i; aliveI[i] = i; }
  let nAlive = 0, nAliveI = 0;

  function samplePos(out3, o) {
    const r = 0.9 * cfg.R * Math.sqrt(rng()), th = rng() * 2 * Math.PI;
    out3[o] = r * Math.cos(th);
    out3[o + 1] = r * Math.sin(th);
    out3[o + 2] = cfg.L * (0.1 + 0.8 * rng());
  }
  function injectPair() {
    if (nAlive >= cap || nAliveI >= cap) return false;
    const se = alive[nAlive], i3 = 3 * se;
    samplePos(pos, i3);
    vel[i3] = gauss() * vth; vel[i3 + 1] = gauss() * vth; vel[i3 + 2] = gauss() * vth;
    // ion born at the same spot — the pair is charge-neutral at injection
    const si = aliveI[nAliveI], j3 = 3 * si;
    posI[j3] = pos[i3]; posI[j3 + 1] = pos[i3 + 1]; posI[j3 + 2] = pos[i3 + 2];
    velI[j3] = gauss() * vthI; velI[j3 + 1] = gauss() * vthI; velI[j3 + 2] = gauss() * vthI;
    nAlive++; nAliveI++;
    return true;
  }
  for (let i = 0; i < cfg.n0seed; i++) injectPair();

  const nSteps = Math.ceil(cfg.Ttotal / dt);
  const measStart = nSteps >> 1;
  let injAcc = 0, capHits = 0;
  let neSum = 0, niSum = 0, nSamp = 0, fieldSamp = 0;
  let lostEEnd = 0, lostERad = 0, lostIEnd = 0, lostIRad = 0; // measurement window only

  for (let step = 0; step < nSteps; step++) {
    const measuring = step >= measStart;
    injAcc += Ssuper * dt;
    while (injAcc >= 1) { if (!injectPair()) capHits++; injAcc -= 1; }

    if (step % cfg.fieldEvery === 0) {
      depositInto(ne, pos, alive, nAlive);
      depositInto(ni, posI, aliveI, nAliveI);
      for (let c = 0; c < S.length; c++) S[c] = (C.QE * (ni[c] - ne[c])) / C.EPS0;
      sor.solve(S, { maxIter: step === 0 ? 3000 : 60, tol: 1e-4 });
      sor.gradients();
      if (measuring) {
        for (let c = 0; c < S.length; c++) {
          phiAvg[c] += sor.phi[c]; neAvg[c] += ne[c]; niAvg[c] += ni[c];
        }
        fieldSamp++;
      }
    }

    for (let k = 0; k < nAlive; k++) {
      const idx = alive[k], i3 = 3 * idx;
      let x = pos[i3], y = pos[i3 + 1], z = pos[i3 + 2];
      interpE(sor.Er, sor.Ez, x, y, z, Evec);
      vel[i3] += h * Evec[0]; vel[i3 + 1] += h * Evec[1]; vel[i3 + 2] += h * Evec[2];
      field.at(x, y, z, B);
      borisRotate(vel, i3, B[0], B[1], B[2], h);
      vel[i3] += h * Evec[0]; vel[i3 + 1] += h * Evec[1]; vel[i3 + 2] += h * Evec[2];
      x += vel[i3] * dt; y += vel[i3 + 1] * dt; z += vel[i3 + 2] * dt;
      pos[i3] = x; pos[i3 + 1] = y; pos[i3 + 2] = z;
      let lost = 0;
      if (x * x + y * y >= R2) lost = KIND.RADIAL;
      else if (z <= 0 || z >= cfg.L) lost = KIND.END_LOW;
      if (lost) {
        if (measuring) { if (lost === KIND.RADIAL) lostERad++; else lostEEnd++; }
        alive[k] = alive[--nAlive]; alive[nAlive] = idx;
        k--;
        continue;
      }
      eColl.collide(vel, i3);
    }
    for (let k = 0; k < nAliveI; k++) {
      const idx = aliveI[k], i3 = 3 * idx;
      let x = posI[i3], y = posI[i3 + 1], z = posI[i3 + 2];
      interpE(sor.Er, sor.Ez, x, y, z, Evec);
      velI[i3] += hI * Evec[0]; velI[i3 + 1] += hI * Evec[1]; velI[i3 + 2] += hI * Evec[2];
      field.at(x, y, z, B);
      borisRotate(velI, i3, B[0], B[1], B[2], hI);
      velI[i3] += hI * Evec[0]; velI[i3 + 1] += hI * Evec[1]; velI[i3 + 2] += hI * Evec[2];
      x += velI[i3] * dt; y += velI[i3 + 1] * dt; z += velI[i3 + 2] * dt;
      posI[i3] = x; posI[i3 + 1] = y; posI[i3 + 2] = z;
      let lost = 0;
      if (x * x + y * y >= R2) lost = KIND.RADIAL;
      else if (z <= 0 || z >= cfg.L) lost = KIND.END_LOW;
      if (lost) {
        if (measuring) { if (lost === KIND.RADIAL) lostIRad++; else lostIEnd++; }
        aliveI[k] = aliveI[--nAliveI]; aliveI[nAliveI] = idx;
        k--;
        continue;
      }
      iColl.collide(velI, i3);
    }
    if (measuring) { neSum += nAlive; niSum += nAliveI; nSamp++; }
  }

  const NeSS = neSum / nSamp, NiSS = niSum / nSamp;
  const tauEff = NeSS / Ssuper;
  for (let c = 0; c < phiAvg.length; c++) {
    phiAvg[c] /= fieldSamp; neAvg[c] /= fieldSamp; niAvg[c] /= fieldSamp;
  }
  let phiMax = 0;
  for (let c = 0; c < phiAvg.length; c++) if (phiAvg[c] > phiMax) phiMax = phiAvg[c];
  const lostE = lostEEnd + lostERad, lostI = lostIEnd + lostIRad;
  const gRatio = lostE > 0 ? lostI / lostE : Infinity;
  const wallSec = ((Date.now() - t0Wall) / 1000).toFixed(1);
  console.log(
    `  ${ratio.toFixed(2).padStart(5)} | ${NeSS.toFixed(0).padStart(5)} ${NiSS.toFixed(0).padStart(5)} | ` +
    `${(tauEff * 1e6).toFixed(2).padStart(10)} | ${phiMax.toFixed(2).padStart(8)} | ` +
    `${gRatio.toFixed(2).padStart(5)} | ` +
    `${lostE ? ((100 * lostEEnd) / lostE).toFixed(1).padStart(7) : '    n/a'}/${lostE ? ((100 * lostERad) / lostE).toFixed(1) : 'n/a'} | ${wallSec.padStart(4)}`);
  summary.push([ratio, NeSS, NiSS, tauEff, phiMax, gRatio,
    lostE ? lostEEnd / lostE : 0, lostE ? lostERad / lostE : 0, capHits].join(','));
  if (capHits > 0) console.log(`        WARNING: source capped ${capHits} times — raise --cap or lower --rate`);

  const rows = ['r_m,z_m,phi_V,ne_m3,ni_m3'];
  for (let j = 0; j < Nz; j++) {
    for (let i = 0; i < Nr; i++) {
      const c = j * Nr + i;
      rows.push([((i + 0.5) * dr).toFixed(5), ((j + 0.5) * dz).toFixed(5),
        phiAvg[c].toExponential(4), neAvg[c].toExponential(4), niAvg[c].toExponential(4)].join(','));
    }
  }
  fs.writeFileSync(path.join(outDir, `sustained-fields-${runTag}-ratio${ratio}-${stamp}.csv`), rows.join('\n'));
}

const outFile = path.join(outDir, `sustained-${runTag}-${stamp}.csv`);
fs.writeFileSync(outFile, summary.join('\n'));
console.log(`\n  tau_eff = <Ne>/S (uncensored, from steady inventory); fields CSVs are time-averaged`);
console.log(`  CSV written: ${path.relative(process.cwd(), outFile)}`);
