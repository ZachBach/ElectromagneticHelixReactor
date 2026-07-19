/* Study A — fit the confinement exponent alpha in

     tau_c ∝ <|B|/Bz>^alpha,

   where <|B|/Bz> is the field-line path-length factor averaged over the
   seeding distribution (computed numerically from the field model, so any
   model/profile works). Free streaming predicts alpha ≈ 1 (residence time ∝
   path length); fully diffusive parallel transport predicts alpha → 2
   (residence ∝ path²). The measured alpha locates each run on that continuum
   and is the core of the proposed transport law.

   Usage: node sim/fit-alpha.js sim/output/scan1.csv [scan2.csv ...]
   Fits both tau_mean and tau_median (medians are less sensitive to the
   censored slow-v_z tail).                                                  */
'use strict';

const fs = require('fs');
const { makeField } = require('./field');

// seeding-averaged path factor: <|B|/Bz> over r ∈ [0, 0.9R], weight 2r dr
function pathFactor(field, R) {
  let s = 0, w = 0;
  const M = 400, B = new Float64Array(3);
  for (let i = 0; i < M; i++) {
    const r = (0.9 * R * (i + 0.5)) / M;
    field.at(r, 0, 0, B);
    const mag = Math.hypot(B[0], B[1], B[2]);
    s += (mag / B[2]) * r;
    w += r;
  }
  return s / w;
}

// least squares y = a + alpha·x with slope standard error
function fitSlope(xs, ys) {
  const n = xs.length;
  let sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; }
  const den = n * sxx - sx * sx;
  const slope = (n * sxy - sx * sy) / den;
  const a = (sy - slope * sx) / n;
  let sse = 0;
  for (let i = 0; i < n; i++) { const e = ys[i] - a - slope * xs[i]; sse += e * e; }
  const se = Math.sqrt((sse / Math.max(n - 2, 1)) * (n / den));
  return { slope, se };
}

for (const file of process.argv.slice(2)) {
  const lines = fs.readFileSync(file, 'utf8').trim().split('\n');
  const meta = {};
  for (const l of lines.filter((l) => l.startsWith('#'))) {
    for (const m of l.matchAll(/([A-Za-z_]+)=([^\s]+)/g)) meta[m[1]] = m[2];
  }
  const header = lines.find((l) => l.startsWith('ratio,')).split(',');
  const iRatio = header.indexOf('ratio');
  const iMean = header.indexOf('tauMean_s');
  const iMed = header.indexOf('tauMedian_s');
  const model = meta.model, R = parseFloat(meta.R_m), Bwall = parseFloat(meta.Bwall_T);
  const nExp = meta.nexp !== undefined ? parseFloat(meta.nexp) : 1;

  const xs = [], yMean = [], yMed = [];
  for (const l of lines) {
    if (l.startsWith('#') || l.startsWith('ratio,')) continue;
    const c = l.split(',');
    const ratio = parseFloat(c[iRatio]);
    const field = makeField({ model, Bwall, ratio, R, nExp });
    xs.push(Math.log(pathFactor(field, R)));
    yMean.push(Math.log(parseFloat(c[iMean])));
    yMed.push(Math.log(parseFloat(c[iMed])));
  }
  const fMean = fitSlope(xs, yMean), fMed = fitSlope(xs, yMed);
  const tag = [
    model + (model === 'powerlaw' ? `(n=${nExp})` : ''),
    meta.p_mTorr ? `${meta.p_mTorr} mTorr` : 'collisionless',
    meta.nocoll === 'true' ? 'nocoll' : '',
    meta.Tseed_eV ? `${meta.Tseed_eV} eV` : '',
  ].filter(Boolean).join(', ');
  console.log(`${file}`);
  console.log(`  [${tag}]  alpha_mean = ${fMean.slope.toFixed(2)} ± ${fMean.se.toFixed(2)}` +
    `   alpha_median = ${fMed.slope.toFixed(2)} ± ${fMed.se.toFixed(2)}   (${xs.length} pts)`);
}
