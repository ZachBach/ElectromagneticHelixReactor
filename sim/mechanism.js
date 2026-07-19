/* Mechanism decomposition, layer 0 — the free-streaming geometric model.

   Hypothesis: the pitch effect on confinement is geometric. Field lines lie on
   cylinders (Br = 0), so a collisionless particle keeps its parallel speed
   (|B| constant along each line → no mirror force) and streams axially at
     v_z = v_par · Bz(r)/|B(r)|,
   i.e. tilting the field slows end loss by the local factor |B|/Bz. This model
   integrates nothing: each particle's loss time is distance/|v_z|, censored at
   Tmax exactly like the Boris runs. No gyration, no drifts, no collisions.

   Compare tau(ratio) from this model against `run-pitch-scan.js --nocoll`
   (adds gyration + grad-B/curvature drifts) and the full run (adds MCC):
   agreement pins the mechanism; discrepancies attribute the remainder.

   Usage: node sim/mechanism.js [--model=screw|beltrami] [--n=20000]
     [--ratios=...] [--bwall=0.01] [--te=3] [--r=0.1] [--l=0.4]
     [--tmax=5e-5] [--seed=1]                                                */
'use strict';

const C = require('./constants');
const { makeField } = require('./field');
const coll = require('./collisions');
const { KIND, LossStats } = require('./diagnostics');

const args = {};
for (const a of process.argv.slice(2)) {
  const m = /^--([^=]+)(?:=(.*))?$/.exec(a);
  if (m) args[m[1]] = m[2] === undefined ? true : m[2];
}
const cfg = {
  model: args.model || 'screw',
  n: parseInt(args.n || '20000', 10),
  ratios: (args.ratios || '0,0.25,0.5,0.75,1,1.5,2,3').split(',').map(Number),
  Bwall: parseFloat(args.bwall || '0.01'),
  TeV: parseFloat(args.te || '3'),
  R: parseFloat(args.r || '0.1'),
  L: parseFloat(args.l || '0.4'),
  Tmax: parseFloat(args.tmax || '5e-5'),
  seed: parseInt(args.seed || '1', 10),
};
const vth = Math.sqrt((cfg.TeV * C.EV) / C.ME);

console.log(`Free-streaming geometric model — ${cfg.model} field, electrons at ${cfg.TeV} eV`);
console.log(`  |B(R)| = ${(cfg.Bwall * 1e4).toFixed(0)} G, R = ${cfg.R} m, L = ${cfg.L} m, ` +
  `N = ${cfg.n}, Tmax = ${cfg.Tmax} s, seed = ${cfg.seed}\n`);
console.log('  ratio | tau_mean us   tau_med | surv%');
console.log('  ------------------------------------');

const B = new Float64Array(3);
for (let ri = 0; ri < cfg.ratios.length; ri++) {
  const ratio = cfg.ratios[ri];
  const rng = coll.makeRng(cfg.seed * 7919 + ri * 104729 + 1);
  const gauss = coll.makeGauss(rng);
  const field = makeField({ model: cfg.model, Bwall: cfg.Bwall, ratio, R: cfg.R });
  const stats = new LossStats(cfg.n);
  for (let i = 0; i < cfg.n; i++) {
    // same seeding distribution as the Boris runner
    const r = 0.9 * cfg.R * Math.sqrt(rng()), th = rng() * 2 * Math.PI;
    const x = r * Math.cos(th), y = r * Math.sin(th);
    const z = cfg.L * (0.1 + 0.8 * rng());
    const vx = gauss() * vth, vy = gauss() * vth, vz = gauss() * vth;
    field.at(x, y, z, B);
    const bm = Math.hypot(B[0], B[1], B[2]);
    const vpar = (vx * B[0] + vy * B[1] + vz * B[2]) / bm;
    const w = vpar * (B[2] / bm); // axial streaming speed
    const t = w > 0 ? (cfg.L - z) / w : w < 0 ? z / -w : Infinity;
    if (t < cfg.Tmax) stats.record(t, w > 0 ? KIND.END_HIGH : KIND.END_LOW);
  }
  const s = stats.summary(cfg.n, cfg.Tmax);
  console.log(
    `  ${ratio.toFixed(2).padStart(5)} | ${s.censored ? '>' : ' '}${(s.tauMean * 1e6).toFixed(2).padStart(8)} ` +
    `±${(s.tauSE * 1e6).toFixed(2)} ${(s.tauMedian * 1e6).toFixed(2).padStart(8)} | ` +
    `${(100 * s.surviveFrac).toFixed(1).padStart(5)}`);
}
console.log(`\n  ('>' = censored at Tmax, as in the Boris runs)`);
