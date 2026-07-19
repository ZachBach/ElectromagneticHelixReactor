/* Quantitative analysis of ambipolar field snapshots (Phase 3 follow-up).

   Reads ambipolar-fields-*.csv snapshots (r,z grids of phi, ne, ni) and
   prints, per file:
     - plasma potential: phi_max, its location, the on-axis plateau mean,
       radial and axial half-max widths (well *shape*, not just depth)
     - electron density profile at the axial midsection: peak radius,
       core/edge ratio n(0)/n(0.7R), and uniformity sigma/mu over r < 0.7R
       (the chuck footprint used by the web app)
   Convention: electrons are trapped by a POSITIVE plasma potential, so the
   "well depth" for electrons is phi_max relative to grounded walls.

   Usage: node sim/analyze-fields.js sim/output/ambipolar-fields-*.csv      */
'use strict';

const fs = require('fs');

function stats(file) {
  const lines = fs.readFileSync(file, 'utf8').trim().split('\n').slice(1);
  const rs = new Set(), zs = new Set();
  const rows = lines.map((l) => l.split(',').map(Number));
  for (const [r, z] of rows) { rs.add(r); zs.add(z); }
  const rArr = [...rs].sort((a, b) => a - b), zArr = [...zs].sort((a, b) => a - b);
  const Nr = rArr.length, Nz = zArr.length;
  const idx = new Map(rArr.map((v, i) => [v, i]));
  const jdx = new Map(zArr.map((v, j) => [v, j]));
  const phi = new Float64Array(Nr * Nz), ne = new Float64Array(Nr * Nz);
  for (const [r, z, p, n] of rows) { const c = jdx.get(z) * Nr + idx.get(r); phi[c] = p; ne[c] = n; }

  // --- potential structure ---
  let phiMax = -Infinity, rAt = 0, zAt = 0;
  for (let j = 0; j < Nz; j++) for (let i = 0; i < Nr; i++) {
    const v = phi[j * Nr + i];
    if (v > phiMax) { phiMax = v; rAt = rArr[i]; zAt = zArr[j]; }
  }
  // on-axis plateau: mean of phi(r0, z) over the central 50% of z
  let plat = 0, nPlat = 0;
  for (let j = (Nz * 0.25) | 0; j < (Nz * 0.75) | 0; j++) { plat += phi[j * Nr]; nPlat++; }
  plat /= nPlat;
  // radial half-max width at the axial midplane
  const jm = (Nz / 2) | 0;
  let rHalf = rArr[Nr - 1];
  for (let i = 0; i < Nr; i++) if (phi[jm * Nr + i] < phiMax / 2) { rHalf = rArr[i]; break; }
  // axial half-max width on axis
  let zLo = zArr[0], zHi = zArr[Nz - 1];
  for (let j = 0; j < Nz; j++) if (phi[j * Nr] >= phiMax / 2) { zLo = zArr[j]; break; }
  for (let j = Nz - 1; j >= 0; j--) if (phi[j * Nr] >= phiMax / 2) { zHi = zArr[j]; break; }

  // --- electron density profile, averaged over the central 50% of z ---
  const prof = new Float64Array(Nr);
  for (let i = 0; i < Nr; i++) {
    let s = 0, c = 0;
    for (let j = (Nz * 0.25) | 0; j < (Nz * 0.75) | 0; j++) { s += ne[j * Nr + i]; c++; }
    prof[i] = s / c;
  }
  let peakN = 0, rPeak = 0;
  for (let i = 0; i < Nr; i++) if (prof[i] > peakN) { peakN = prof[i]; rPeak = rArr[i]; }
  const R = rArr[Nr - 1] + (rArr[1] - rArr[0]) / 2;
  const i07 = Math.min(Nr - 1, Math.round((0.7 * R) / (rArr[1] - rArr[0]) - 0.5));
  const coreEdge = prof[i07] > 0 ? prof[0] / prof[i07] : Infinity;
  // uniformity over the chuck footprint r < 0.7R, area-weighted
  let mu = 0, w = 0;
  for (let i = 0; i <= i07; i++) { mu += prof[i] * rArr[i]; w += rArr[i]; }
  mu /= w;
  let va = 0;
  for (let i = 0; i <= i07; i++) { const d = prof[i] - mu; va += d * d * rArr[i]; }
  const unif = Math.sqrt(va / w) / mu;

  return { phiMax, rAt, zAt, plat, rHalf, zLo, zHi, rPeak, coreEdge, unif };
}

console.log('file | phiMax V (plateau) | well r_half m, z FWHM m | ne: r_peak m, core/edge, sigma/mu(r<0.7R)');
for (const file of process.argv.slice(2)) {
  const s = stats(file);
  const short = file.replace(/^.*ambipolar-fields-/, '').replace(/-2026.*$/, '');
  console.log(
    `${short.padEnd(22)} | ${s.phiMax.toFixed(2)} (${s.plat.toFixed(2)}) | ` +
    `${s.rHalf.toFixed(3)}, ${(s.zHi - s.zLo).toFixed(3)} | ` +
    `${s.rPeak.toFixed(3)}, ${s.coreEdge.toFixed(2)}, ${(100 * s.unif).toFixed(1)}%`);
}
