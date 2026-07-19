/* Phase 3 — ambipolar feedback study (Poisson Milestones 2 + 3):

     Does the pitch-dependent confinement survive once the plasma is allowed
     to generate its own electric field?

   Model: kinetic electrons (Boris + MCC as in Phase 1) over a FROZEN ion
   background equal to the initial electron deposition, so rho = 0 at t = 0
   and space charge appears only as electrons leave or redistribute. This is
   the standard immobile-ion approximation on electron timescales: electron
   escape charges the plasma positive, the resulting potential well pulls
   electrons back — ambipolar confinement without ion dynamics.

   Loop (every --fieldevery steps): CIC-deposit n_e → rho = e(n_i - n_e) →
   SOR Poisson (warm start) → E = -grad phi → Boris with half-E kicks.

   Numerical domain of validity (printed at startup): omega_pe·dt ≲ 0.2 and
   lambda_D ≳ grid spacing set the usable density; the default n0 = 1e14 m^-3
   satisfies both on the default 48×96 grid.

   Mobile ions (--ions=mobile): kinetic ions replace the frozen background —
   both species are pushed, both deposit, both are lost to walls, and the
   system relaxes toward true ambipolar equilibrium (the wall-flux ratio
   Gamma_i/Gamma_e → 1 is reported as the equilibrium check). `--mion` sets
   the ion AND neutral mass in amu: the default 1 (hydrogen-like model gas,
   electron cross sections kept from argon) pulls the H_i = 1 ion
   magnetization boundary down into an affordable field range so a B-scan can
   cross it. H_i is printed per run. Regime conclusions depend on the Hall
   parameter, not species identity.

   Usage: node sim/run-ambipolar.js [--feedback=on|off] [--ions=frozen|mobile]
     [--mion=1] [--n0=3e11] [--n=2000] [--ratios=0,1,3]
     [--model=screw|beltrami|powerlaw] [--nexp=1] [--bwall=0.01]
     [--pressure=10] [--te=3] [--r=0.1] [--l=0.4] [--tmax=2e-5] [--seed=1]
     [--nr=48] [--nz=96] [--fieldevery=<steps>]                             */
'use strict';

const fs = require('fs');
const path = require('path');
const C = require('./constants');
const { makeField } = require('./field');
const { borisRotate } = require('./boris');
const coll = require('./collisions');
const { KIND, LossStats } = require('./diagnostics');
const { PoissonRZ } = require('./poisson');

const args = {};
for (const a of process.argv.slice(2)) {
  const m = /^--([^=]+)(?:=(.*))?$/.exec(a);
  if (m) args[m[1]] = m[2] === undefined ? true : m[2];
}
const cfg = {
  feedback: (args.feedback || 'on') !== 'off',
  /* n0 sets the space-charge reservoir. The full-depletion potential
     e·n0·R²/(4·eps0) must stay at a few Te for a well-behaved quasi-static
     demonstration — at helicon densities (1e17+) this model would need an
     implicit or quasi-neutral formulation instead. 3e11 → ~13 V ≈ 4.5·Te. */
  n0: parseFloat(args.n0 || '3e11'),
  n: parseInt(args.n || '2000', 10),
  ratios: (args.ratios || '0,1,3').split(',').map(Number),
  model: args.model || 'screw',
  nExp: parseFloat(args.nexp || '1'),
  Bwall: parseFloat(args.bwall || '0.01'),
  pmTorr: parseFloat(args.pressure || '10'),
  TeV: parseFloat(args.te || '3'),
  Tgas: 300,
  R: parseFloat(args.r || '0.1'),
  L: parseFloat(args.l || '0.4'),
  Tmax: parseFloat(args.tmax || '2e-5'),
  seed: parseInt(args.seed || '1', 10),
  Nr: parseInt(args.nr || '48', 10),
  Nz: parseInt(args.nz || '96', 10),
  ions: args.ions || 'frozen',
  mIonAmu: parseFloat(args.mion || '1'),
};

const mob = cfg.ions === 'mobile';
const mI = cfg.mIonAmu * C.AMU;
const vth = Math.sqrt((cfg.TeV * C.EV) / C.ME);
const nn = C.neutralDensity(cfg.pmTorr, cfg.Tgas);
const qOverM = -C.QE / C.ME;
const qOverMI = C.QE / mI;
const wpe = Math.sqrt((cfg.n0 * C.QE * C.QE) / (C.EPS0 * C.ME));
const lamD = Math.sqrt((C.EPS0 * cfg.TeV * C.EV) / (cfg.n0 * C.QE * C.QE));
const wcMax = Math.abs(qOverM) * cfg.Bwall * (cfg.model === 'beltrami' ? 1.5 : 1);
const nuMax = coll.electronMaxRate(nn);
const cs = Math.sqrt((cfg.TeV * C.EV) / mI);
const nuMaxI = coll.ionMaxRate(nn, 5 * cs);
const dt = Math.min(0.2 / wcMax, 0.05 / nuMax, 0.2 / wpe, mob ? 0.05 / nuMaxI : Infinity);
const h = (qOverM * dt) / 2;
const eKick = (qOverM * dt) / 2; // half-step E acceleration factor
const hI = (qOverMI * dt) / 2;
// ion Hall parameter, two conventions: thermal (the paper's), and evaluated
// at the ambipolar exit speed c_s — escaping ions are accelerated to ~c_s by
// the well, where CX collisionality is higher, so Hi_star governs whether
// the LOSS channel is magnetized. First B-scan (10 mTorr) showed the gate
// sits at Hi_star, not thermal H_i.
const vBarI = Math.sqrt((8 * C.KB * cfg.Tgas) / (Math.PI * mI));
const Hi = (qOverMI * cfg.Bwall) / (nn * (C.SIGMA_CX + C.SIGMA_I_EL) * vBarI);
const HiStar = (qOverMI * cfg.Bwall) / (nn * (C.SIGMA_CX + C.SIGMA_I_EL) * cs);
// field update every ~0.57 ns of physical time regardless of dt
cfg.fieldEvery = args.fieldevery
  ? parseInt(args.fieldevery, 10)
  : Math.max(1, Math.round(5.7e-10 / dt));

const Nr = cfg.Nr, Nz = cfg.Nz, dr = cfg.R / Nr, dz = cfg.L / Nz;
const Vseed = Math.PI * (0.9 * cfg.R) ** 2 * (0.8 * cfg.L);
const wSuper = (cfg.n0 * Vseed) / cfg.n; // real electrons per super-particle

console.log(`EHR Phase 3 ambipolar study — feedback ${cfg.feedback ? 'ON' : 'OFF'}, ` +
  `ions ${mob ? `MOBILE (m=${cfg.mIonAmu} amu, H_i=${Hi.toFixed(2)}, H_i*=${HiStar.toFixed(2)} at c_s)` : 'frozen'}, ` +
  `${cfg.model}${cfg.model === 'powerlaw' ? `(n=${cfg.nExp})` : ''} field`);
console.log(`  n0 = ${cfg.n0.toExponential(1)} m^-3, |B(R)| = ${(cfg.Bwall * 1e4).toFixed(0)} G, ` +
  `p = ${cfg.pmTorr} mTorr, Te = ${cfg.TeV} eV, grid ${Nr}x${Nz}`);
const phiCap = (C.QE * cfg.n0 * cfg.R * cfg.R) / (4 * C.EPS0);
console.log(`  N = ${cfg.n}, Tmax = ${cfg.Tmax} s, dt = ${dt.toExponential(2)} s | ` +
  `wpe*dt = ${(wpe * dt).toFixed(4)}, lambda_D/dr = ${(lamD / dr).toFixed(2)}, ` +
  `phi_cap = ${phiCap.toFixed(1)} V (${(phiCap / cfg.TeV).toFixed(1)} Te), ` +
  `field update every ${cfg.fieldEvery} steps`);
console.log('');
console.log(mob
  ? '  ratio | survE% survI% | tau_mean us | tau_leak us | phiMax V | Gi/Ge |    s'
  : '  ratio | surv%  | tau_mean us | tau_leak us | phiMax V | radial%  end% |    s');
console.log('  ' + '-'.repeat(82));

// CIC deposit of super-particles into a number-density grid (per m^3)
function deposit(grid, pos, alive, nAlive) {
  grid.fill(0);
  for (let k = 0; k < nAlive; k++) {
    const i3 = 3 * alive[k];
    const x = pos[i3], y = pos[i3 + 1], z = pos[i3 + 2];
    const r = Math.sqrt(x * x + y * y);
    const rf = r / dr - 0.5, zf = z / dz - 0.5;
    const i0 = Math.floor(rf), j0 = Math.floor(zf);
    const fr = rf - i0, fz = zf - j0;
    const iA = i0 < 0 ? 0 : i0;                       // axis fold
    const iB = i0 + 1 >= Nr ? Nr - 1 : i0 + 1;
    const jA = j0 < 0 ? 0 : j0 >= Nz ? Nz - 1 : j0;
    const jB = j0 + 1 < 0 ? 0 : j0 + 1 >= Nz ? Nz - 1 : j0 + 1;
    grid[jA * Nr + iA] += (1 - fr) * (1 - fz);
    grid[jB * Nr + iA] += (1 - fr) * fz;
    grid[jA * Nr + iB] += fr * (1 - fz);
    grid[jB * Nr + iB] += fr * fz;
  }
  // counts → density: divide by cell volume 2·pi·rc·dr·dz, times wSuper
  for (let j = 0; j < Nz; j++) {
    for (let i = 0; i < Nr; i++) {
      const vol = 2 * Math.PI * (i + 0.5) * dr * dr * dz;
      grid[j * Nr + i] *= wSuper / vol;
    }
  }
}

// bilinear E interpolation with axis parity (Er odd, Ez even) and clamping;
// allocation-free — runs once per particle per step
function interpE(Er, Ez, x, y, z, out) {
  const r = Math.sqrt(x * x + y * y);
  const rf = r / dr - 0.5, zf = z / dz - 0.5;
  const i0 = Math.floor(rf), j0 = Math.floor(zf);
  const fr = rf - i0, fz = zf - j0;
  let iA = i0, sA = 1;
  if (iA < 0) { iA = 0; sA = -1; }                    // mirror: Er(-r) = -Er(r)
  else if (iA >= Nr) iA = Nr - 1;
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
const runTag = `${cfg.model}${cfg.model === 'powerlaw' ? '-n' + cfg.nExp : ''}-fb${cfg.feedback ? 'on' : 'off'}` +
  (mob ? `-mob${(cfg.Bwall * 1e4).toFixed(0)}G` : '');
const summary = [
  `# EHR Phase 3 ambipolar study  ${new Date().toISOString()}`,
  `# feedback=${cfg.feedback} ions=${cfg.ions} mion_amu=${cfg.mIonAmu} Hi=${Hi} HiStar=${HiStar} n0_m3=${cfg.n0} model=${cfg.model} nexp=${cfg.nExp} n=${cfg.n} Bwall_T=${cfg.Bwall} p_mTorr=${cfg.pmTorr}`,
  `# Tseed_eV=${cfg.TeV} R_m=${cfg.R} L_m=${cfg.L} Tmax_s=${cfg.Tmax} dt_s=${dt} grid=${Nr}x${Nz} fieldevery=${cfg.fieldEvery} seed=${cfg.seed}`,
  'ratio,surviveFrac,tauMean_s,tauLeak_s,phiMax_V,fracRadialLoss,fracEndLoss' +
    (mob ? ',surviveFracIon,GiOverGe,Hi' : ''),
];

const R2 = cfg.R * cfg.R;
const B = new Float64Array(3), Evec = new Float64Array(3);

for (let ri = 0; ri < cfg.ratios.length; ri++) {
  const ratio = cfg.ratios[ri];
  const t0Wall = Date.now();
  const rng = coll.makeRng(cfg.seed * 7919 + ri * 104729 + 1);
  const gauss = coll.makeGauss(rng);
  const field = makeField({ model: cfg.model, Bwall: cfg.Bwall, ratio, R: cfg.R, nExp: cfg.nExp });
  const collider = coll.makeElectronCollider(nn, dt, rng);
  const sor = new PoissonRZ(Nr, Nz, cfg.R, cfg.L);
  const ne = new Float64Array(Nr * Nz), ni = new Float64Array(Nr * Nz);
  const S = new Float64Array(Nr * Nz);

  const N = cfg.n;
  const pos = new Float64Array(3 * N), vel = new Float64Array(3 * N);
  const alive = new Int32Array(N);
  let nAlive = N;
  for (let i = 0; i < N; i++) {
    alive[i] = i;
    const r = 0.9 * cfg.R * Math.sqrt(rng()), th = rng() * 2 * Math.PI;
    pos[3 * i] = r * Math.cos(th);
    pos[3 * i + 1] = r * Math.sin(th);
    pos[3 * i + 2] = cfg.L * (0.1 + 0.8 * rng());
    vel[3 * i] = gauss() * vth;
    vel[3 * i + 1] = gauss() * vth;
    vel[3 * i + 2] = gauss() * vth;
  }
  deposit(ni, pos, alive, nAlive); // frozen background, or t=0 state if mobile

  // mobile ions: same initial positions as the electrons (exact initial
  // neutrality — space charge appears only dynamically), thermal velocities
  const vthI = Math.sqrt((C.KB * cfg.Tgas) / mI);
  let posI = null, velI = null, aliveI = null, nAliveI = 0, statsI = null, ionColl = null;
  if (mob) {
    posI = pos.slice();
    velI = new Float64Array(3 * N);
    aliveI = new Int32Array(N);
    nAliveI = N;
    for (let i = 0; i < N; i++) {
      aliveI[i] = i;
      velI[3 * i] = gauss() * vthI;
      velI[3 * i + 1] = gauss() * vthI;
      velI[3 * i + 2] = gauss() * vthI;
    }
    statsI = new LossStats(N);
    ionColl = coll.makeIonCollider(nn, cfg.Tgas, dt, rng, gauss, mI);
  }

  const stats = new LossStats(N);
  const nSteps = Math.ceil(cfg.Tmax / dt);
  // alive-count trace for the quasi-steady leak rate
  const traceEvery = Math.max(1, (nSteps / 200) | 0);
  const traceT = [], traceA = [];
  let phiMax = 0;

  for (let step = 0; step < nSteps && nAlive > 0; step++) {
    if (cfg.feedback && step % cfg.fieldEvery === 0) {
      deposit(ne, pos, alive, nAlive);
      if (mob) deposit(ni, posI, aliveI, nAliveI);
      for (let c = 0; c < S.length; c++) S[c] = (C.QE * (ni[c] - ne[c])) / C.EPS0;
      sor.solve(S, { maxIter: step === 0 ? 3000 : 60, tol: 1e-4 });
      sor.gradients();
    }
    const t = step * dt;
    for (let k = 0; k < nAlive; k++) {
      const idx = alive[k], i3 = 3 * idx;
      let x = pos[i3], y = pos[i3 + 1], z = pos[i3 + 2];
      if (cfg.feedback) {
        interpE(sor.Er, sor.Ez, x, y, z, Evec);
        vel[i3] += eKick * Evec[0]; vel[i3 + 1] += eKick * Evec[1]; vel[i3 + 2] += eKick * Evec[2];
      }
      field.at(x, y, z, B);
      borisRotate(vel, i3, B[0], B[1], B[2], h);
      if (cfg.feedback) {
        vel[i3] += eKick * Evec[0]; vel[i3 + 1] += eKick * Evec[1]; vel[i3 + 2] += eKick * Evec[2];
      }
      x += vel[i3] * dt; y += vel[i3 + 1] * dt; z += vel[i3 + 2] * dt;
      pos[i3] = x; pos[i3 + 1] = y; pos[i3 + 2] = z;
      let lost = 0;
      if (x * x + y * y >= R2) lost = KIND.RADIAL;
      else if (z <= 0) lost = KIND.END_LOW;
      else if (z >= cfg.L) lost = KIND.END_HIGH;
      if (lost) {
        stats.record(t + dt, lost);
        alive[k] = alive[--nAlive];
        k--;
        continue;
      }
      collider.collide(vel, i3);
    }
    // ---- ion push (mobile mode): same dt, opposite charge ----
    if (mob) {
      for (let k = 0; k < nAliveI; k++) {
        const idx = aliveI[k], i3 = 3 * idx;
        let x = posI[i3], y = posI[i3 + 1], z = posI[i3 + 2];
        if (cfg.feedback) {
          interpE(sor.Er, sor.Ez, x, y, z, Evec);
          velI[i3] += hI * Evec[0]; velI[i3 + 1] += hI * Evec[1]; velI[i3 + 2] += hI * Evec[2];
        }
        field.at(x, y, z, B);
        borisRotate(velI, i3, B[0], B[1], B[2], hI);
        if (cfg.feedback) {
          velI[i3] += hI * Evec[0]; velI[i3 + 1] += hI * Evec[1]; velI[i3 + 2] += hI * Evec[2];
        }
        x += velI[i3] * dt; y += velI[i3 + 1] * dt; z += velI[i3 + 2] * dt;
        posI[i3] = x; posI[i3 + 1] = y; posI[i3 + 2] = z;
        let lost = 0;
        if (x * x + y * y >= R2) lost = KIND.RADIAL;
        else if (z <= 0) lost = KIND.END_LOW;
        else if (z >= cfg.L) lost = KIND.END_HIGH;
        if (lost) {
          statsI.record(t + dt, lost);
          aliveI[k] = aliveI[--nAliveI];
          k--;
          continue;
        }
        ionColl.collide(velI, i3);
      }
    }
    if (step % traceEvery === 0) { traceT.push(t); traceA.push(nAlive); }
    if (cfg.feedback && step % (cfg.fieldEvery * 10) === 0) {
      for (let c = 0; c < sor.phi.length; c++) if (sor.phi[c] > phiMax) phiMax = sor.phi[c];
    }
  }

  // quasi-steady leak rate: fit ln(alive) over the last half of the trace
  let tauLeak = Infinity;
  const half = traceT.length >> 1;
  if (traceA[traceT.length - 1] > 20 && traceA[half] > traceA[traceT.length - 1]) {
    tauLeak = (traceT[traceT.length - 1] - traceT[half]) /
      Math.log(traceA[half] / traceA[traceT.length - 1]);
  }
  const s = stats.summary(N, cfg.Tmax);
  const wallSec = ((Date.now() - t0Wall) / 1000).toFixed(1);
  if (mob) {
    const sI = statsI.summary(N, cfg.Tmax);
    // ambipolar equilibrium check: wall-flux ratio over the second half
    let lateE = 0, lateI = 0;
    for (let q = 0; q < stats.count; q++) if (stats.times[q] > cfg.Tmax / 2) lateE++;
    for (let q = 0; q < statsI.count; q++) if (statsI.times[q] > cfg.Tmax / 2) lateI++;
    const gRatio = lateE > 0 ? lateI / lateE : Infinity;
    console.log(
      `  ${ratio.toFixed(2).padStart(5)} | ${(100 * s.surviveFrac).toFixed(1).padStart(6)} ` +
      `${(100 * sI.surviveFrac).toFixed(1).padStart(6)} | ` +
      `${s.censored ? '>' : ' '}${(s.tauMean * 1e6).toFixed(2).padStart(9)} | ` +
      `${(tauLeak === Infinity ? Infinity : tauLeak * 1e6).toFixed(1).padStart(11)} | ` +
      `${phiMax.toFixed(2).padStart(8)} | ${gRatio.toFixed(2).padStart(5)} | ${wallSec.padStart(4)}`);
    summary.push([ratio, s.surviveFrac, s.tauMean, tauLeak, phiMax, s.fracRadial, s.fracEnd,
      sI.surviveFrac, gRatio, Hi].join(','));
  } else {
    console.log(
      `  ${ratio.toFixed(2).padStart(5)} | ${(100 * s.surviveFrac).toFixed(1).padStart(5)} | ` +
      `${s.censored ? '>' : ' '}${(s.tauMean * 1e6).toFixed(2).padStart(9)} | ` +
      `${(tauLeak === Infinity ? Infinity : tauLeak * 1e6).toFixed(1).padStart(11)} | ` +
      `${phiMax.toFixed(2).padStart(8)} | ${(100 * s.fracRadial).toFixed(1).padStart(6)} ` +
      `${(100 * s.fracEnd).toFixed(1).padStart(5)} | ${wallSec.padStart(4)}`);
    summary.push([ratio, s.surviveFrac, s.tauMean, tauLeak, phiMax, s.fracRadial, s.fracEnd].join(','));
  }

  // fields snapshot (phi, ne, ni) for this ratio
  if (cfg.feedback) {
    deposit(ne, pos, alive, nAlive);
    const rows = ['r_m,z_m,phi_V,ne_m3,ni_m3'];
    for (let j = 0; j < Nz; j++) {
      for (let i = 0; i < Nr; i++) {
        const c = j * Nr + i;
        rows.push([((i + 0.5) * dr).toFixed(5), ((j + 0.5) * dz).toFixed(5),
          sor.phi[c].toExponential(4), ne[c].toExponential(4), ni[c].toExponential(4)].join(','));
      }
    }
    fs.writeFileSync(path.join(outDir, `ambipolar-fields-${runTag}-ratio${ratio}-${stamp}.csv`), rows.join('\n'));
  }
}

const outFile = path.join(outDir, `ambipolar-${runTag}-${stamp}.csv`);
fs.writeFileSync(outFile, summary.join('\n'));
console.log(`\n  ('>' = censored at Tmax; tau_leak from ln(alive) slope over the last half)`);
console.log(`  CSV written: ${path.relative(process.cwd(), outFile)}`);
