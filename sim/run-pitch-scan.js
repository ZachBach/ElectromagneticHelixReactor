/* Phase 1 headline study — the simulator side of Prediction 4:

     Does scanning the field-line pitch Bθ/Bz, at fixed |B(R)|, pressure, and
     injection energy, change the confinement time and loss channels of a
     test-particle ensemble?

   For each pitch ratio, an ensemble is seeded in the chamber with a Maxwellian
   velocity distribution and integrated (Boris + MCC vs neutral argon, E = 0,
   no heating) until it hits a wall or Tmax. Outputs a console table and a CSV
   in sim/output/.

   Usage:
     node sim/run-pitch-scan.js [--n=4000] [--species=electron|ion]
       [--model=screw|beltrami] [--ratios=0,0.25,0.5,1,1.5,2,3]
       [--bwall=0.01] [--pressure=10] [--te=3] [--ti=0.026] [--tgas=300]
       [--r=0.1] [--l=0.4] [--tmax=<s>] [--seed=1] [--density]                 */
'use strict';

const fs = require('fs');
const path = require('path');
const C = require('./constants');
const { makeField } = require('./field');
const { borisRotate } = require('./boris');
const coll = require('./collisions');
const { KIND, LossStats, DensityRZ } = require('./diagnostics');

// ---- configuration ----
const args = {};
for (const a of process.argv.slice(2)) {
  const m = /^--([^=]+)(?:=(.*))?$/.exec(a);
  if (m) args[m[1]] = m[2] === undefined ? true : m[2];
}
const cfg = {
  species: args.species || 'electron',
  model: args.model || 'screw',
  n: parseInt(args.n || '4000', 10),
  ratios: (args.ratios || '0,0.25,0.5,0.75,1,1.5,2,3').split(',').map(Number),
  Bwall: parseFloat(args.bwall || '0.01'),   // T (0.01 T = 100 G, paper baseline)
  pmTorr: parseFloat(args.pressure || '10'), // paper baseline
  TeV: parseFloat(args.te || '3'),           // electron seed temperature, eV
  TiV: parseFloat(args.ti || '0.026'),       // ion seed temperature, eV (~300 K)
  Tgas: parseFloat(args.tgas || '300'),      // neutral gas temperature, K
  R: parseFloat(args.r || '0.1'),            // chamber radius, m
  L: parseFloat(args.l || '0.4'),            // chamber length, m
  seed: parseInt(args.seed || '1', 10),
  density: !!args.density,
  nocoll: !!args.nocoll,                     // disable MCC: geometry + drifts only
  dtscale: parseFloat(args.dtscale || '1'),  // timestep multiplier for convergence checks
};
cfg.Tmax = parseFloat(args.tmax || (cfg.species === 'ion' ? '0.05' : '5e-5'));

const isIon = cfg.species === 'ion';
const mass = isIon ? C.MAR : C.ME;
const qOverM = (isIon ? C.QE : -C.QE) / mass;
const TseedV = isIon ? cfg.TiV : cfg.TeV;
const vth = Math.sqrt((TseedV * C.EV) / mass); // per-component thermal speed
const nn = C.neutralDensity(cfg.pmTorr, cfg.Tgas);

// ---- timestep: resolve gyration and keep ν·dt small ----
const wcMax = Math.abs(qOverM) * cfg.Bwall * (cfg.model === 'beltrami' ? 1.5 : 1);
const nuMax = isIon ? coll.ionMaxRate(nn, 5 * Math.sqrt(3) * vth) : coll.electronMaxRate(nn);
const dt = Math.min(0.2 / wcMax, 0.05 / nuMax) * cfg.dtscale;
const h = (qOverM * dt) / 2;

console.log(`EHR Phase 1 pitch scan — ${cfg.species}s, ${cfg.model} field` +
  `${cfg.nocoll ? ' [collisions OFF]' : ''}${cfg.dtscale !== 1 ? ` [dt x${cfg.dtscale}]` : ''}`);
console.log(
  `  |B(R)| = ${(cfg.Bwall * 1e4).toFixed(0)} G, p = ${cfg.pmTorr} mTorr ` +
  `(n_n = ${nn.toExponential(2)} m^-3), T_seed = ${TseedV} eV, ` +
  `R = ${cfg.R} m, L = ${cfg.L} m`);
console.log(
  `  N = ${cfg.n}, Tmax = ${cfg.Tmax} s, dt = ${dt.toExponential(2)} s ` +
  `(wc*dt = ${(wcMax * dt).toFixed(3)}, nu_max*dt = ${(nuMax * dt).toExponential(1)}), seed = ${cfg.seed}`);
console.log('');

const tUnit = isIon ? 1e3 : 1e6, tLabel = isIon ? 'ms' : 'us';
const header =
  `  ratio |  Bz(G) BthR(G) | Lp(R/2) cm | tau_mean ${tLabel}   tau_med | ` +
  `radial%  end%  surv% | ion/pt |    s`;
console.log(header);
console.log('  ' + '-'.repeat(header.length - 2));

const csv = [
  `# EHR Phase 1 pitch scan  ${new Date().toISOString()}`,
  `# species=${cfg.species} model=${cfg.model} n=${cfg.n} Bwall_T=${cfg.Bwall} p_mTorr=${cfg.pmTorr}`,
  `# Tseed_eV=${TseedV} Tgas_K=${cfg.Tgas} R_m=${cfg.R} L_m=${cfg.L} Tmax_s=${cfg.Tmax} dt_s=${dt} seed=${cfg.seed} nocoll=${cfg.nocoll} dtscale=${cfg.dtscale}`,
  'ratio,Bz_T,BthWall_T,pitch_halfR_m,nTotal,nLost,surviveFrac,tauMean_s,tauSE_s,tauMedian_s,fracRadialLoss,fracEndLoss,fracEndLow,fracEndHigh,ionEventsPerParticle',
];

const outDir = path.join(__dirname, 'output');
fs.mkdirSync(outDir, { recursive: true });
const stamp = new Date().toISOString().replace(/[:T]/g, '-').slice(0, 19);

const R2 = cfg.R * cfg.R;
const B = new Float64Array(3);

for (let ri = 0; ri < cfg.ratios.length; ri++) {
  const ratio = cfg.ratios[ri];
  const t0 = Date.now();
  const rng = coll.makeRng(cfg.seed * 7919 + ri * 104729 + 1);
  const gauss = coll.makeGauss(rng);
  const field = makeField({ model: cfg.model, Bwall: cfg.Bwall, ratio, R: cfg.R });
  const collider = cfg.nocoll ? null : (isIon
    ? coll.makeIonCollider(nn, cfg.Tgas, dt, rng, gauss)
    : coll.makeElectronCollider(nn, dt, rng));

  const N = cfg.n;
  const pos = new Float64Array(3 * N), vel = new Float64Array(3 * N);
  for (let i = 0; i < N; i++) {
    const r = 0.9 * cfg.R * Math.sqrt(rng()), th = rng() * 2 * Math.PI;
    pos[3 * i] = r * Math.cos(th);
    pos[3 * i + 1] = r * Math.sin(th);
    pos[3 * i + 2] = cfg.L * (0.1 + 0.8 * rng());
    vel[3 * i] = gauss() * vth;
    vel[3 * i + 1] = gauss() * vth;
    vel[3 * i + 2] = gauss() * vth;
  }

  const stats = new LossStats(N);
  const dens = cfg.density ? new DensityRZ(40, 80, cfg.R, cfg.L) : null;
  let ionEvents = 0;

  for (let i = 0; i < N; i++) {
    const i3 = 3 * i;
    let t = 0, step = 0;
    while (t < cfg.Tmax) {
      field.at(pos[i3], pos[i3 + 1], pos[i3 + 2], B);
      borisRotate(vel, i3, B[0], B[1], B[2], h);
      const x = pos[i3] + vel[i3] * dt;
      const y = pos[i3 + 1] + vel[i3 + 1] * dt;
      const z = pos[i3 + 2] + vel[i3 + 2] * dt;
      pos[i3] = x; pos[i3 + 1] = y; pos[i3 + 2] = z;
      t += dt;
      if (x * x + y * y >= R2) { stats.record(t, KIND.RADIAL); break; }
      if (z <= 0) { stats.record(t, KIND.END_LOW); break; }
      if (z >= cfg.L) { stats.record(t, KIND.END_HIGH); break; }
      if (collider && collider.collide(vel, i3) === 2) ionEvents++;
      if (dens && (step & 7) === 0) dens.deposit(x, y, z, dt * 8);
      step++;
    }
  }

  const s = stats.summary(N, cfg.Tmax);
  const pitchHalf = field.pitchAt(cfg.R / 2);
  const meta = field.meta;
  const BzG = (meta.Bz !== undefined ? meta.Bz : meta.B0) * 1e4;
  const BthG = (meta.BthWall !== undefined
    ? meta.BthWall
    : Math.sqrt(cfg.Bwall * cfg.Bwall - (cfg.Bwall / Math.sqrt(1 + ratio * ratio)) ** 2)) * 1e4;
  const wallSec = ((Date.now() - t0) / 1000).toFixed(1);
  const cens = s.censored ? '>' : ' ';
  console.log(
    `  ${ratio.toFixed(2).padStart(5)} | ${BzG.toFixed(1).padStart(6)} ${BthG.toFixed(1).padStart(7)} | ` +
    `${(pitchHalf === Infinity ? Infinity : pitchHalf * 100).toFixed(1).padStart(10)} | ` +
    `${cens}${(s.tauMean * tUnit).toFixed(2).padStart(8)} ±${(s.tauSE * tUnit).toFixed(2)}` +
    ` ${(s.tauMedian * tUnit).toFixed(2).padStart(8)} | ` +
    `${(100 * s.fracRadial).toFixed(1).padStart(6)} ${(100 * s.fracEnd).toFixed(1).padStart(5)} ` +
    `${(100 * s.surviveFrac).toFixed(1).padStart(5)} | ${(ionEvents / N).toFixed(2).padStart(6)} | ${wallSec.padStart(4)}`);

  csv.push([
    ratio, meta.Bz !== undefined ? meta.Bz : meta.B0, BthG / 1e4, pitchHalf, N, s.nLost,
    s.surviveFrac, s.tauMean, s.tauSE, s.tauMedian, s.fracRadial, s.fracEnd,
    s.fracEndLow, s.fracEndHigh, ionEvents / N,
  ].join(','));

  if (dens) {
    const dcsv = ['r_m,z_m,density_arb'];
    for (const row of dens.rows()) dcsv.push(row.join(','));
    fs.writeFileSync(path.join(outDir, `density-${stamp}-ratio${ratio}.csv`), dcsv.join('\n'));
  }
}

const outFile = path.join(outDir, `pitch-scan-${stamp}.csv`);
fs.writeFileSync(outFile, csv.join('\n'));
console.log(`\n  ('>' = ensemble censored at Tmax; tau_mean is then a lower bound)`);
console.log(`  CSV written: ${path.relative(process.cwd(), outFile)}`);
