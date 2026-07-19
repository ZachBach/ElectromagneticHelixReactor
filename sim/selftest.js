/* Physics self-test for the Phase 1 simulator: run `node sim/selftest.js`.
   Checks the Boris pusher against exact gyromotion, both field models against
   ∇·B = 0, the Beltrami model against ∇×B = αB, and the MCC sampling against
   the analytic collision rate. Exits nonzero on failure. */
'use strict';

const C = require('./constants');
const { makeField } = require('./field');
const { borisRotate } = require('./boris');
const coll = require('./collisions');

let failures = 0;
function check(name, ok, detail) {
  console.log(`  ${ok ? 'PASS' : 'FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
  if (!ok) failures++;
}

// ---- 1. Boris pusher: energy conservation and gyroradius in uniform Bz ----
{
  const Bz = 0.01; // T
  const E0 = 3; // eV, all perpendicular
  const v0 = Math.sqrt((2 * E0 * C.EV) / C.ME);
  const wc = (C.QE / C.ME) * Bz;
  const dt = 0.2 / wc;
  const h = (-C.QE / C.ME) * dt / 2;
  const vel = new Float64Array([v0, 0, 0]);
  let x = 0, y = 0;
  let xMin = 0, xMax = 0, yMin = 0, yMax = 0;
  const steps = 200000;
  for (let s = 0; s < steps; s++) {
    borisRotate(vel, 0, 0, 0, Bz, h);
    x += vel[0] * dt; y += vel[1] * dt;
    if (x < xMin) xMin = x; if (x > xMax) xMax = x;
    if (y < yMin) yMin = y; if (y > yMax) yMax = y;
  }
  const vf = Math.sqrt(vel[0] ** 2 + vel[1] ** 2 + vel[2] ** 2);
  const rL = (C.ME * v0) / (C.QE * Bz);
  const dia = ((xMax - xMin) + (yMax - yMin)) / 2; // should be 2 rL
  check('Boris energy conservation', Math.abs(vf / v0 - 1) < 1e-10,
    `|v| drift ${(vf / v0 - 1).toExponential(1)} over ${steps} steps`);
  check('Boris gyroradius', Math.abs(dia / (2 * rL) - 1) < 0.02,
    `orbit diameter/2rL = ${(dia / (2 * rL)).toFixed(4)}, rL = ${(rL * 1e3).toFixed(2)} mm`);
}

// ---- 2. divergence and curl of the field models ----
{
  const R = 0.1, eps = 1e-6;
  const Bp = new Float64Array(3), Bm = new Float64Array(3);
  function divB(f, x, y, z) {
    let d = 0;
    for (const [ax, i] of [[0, 0], [1, 1], [2, 2]]) {
      const p = [x, y, z]; p[ax] += eps;
      const m = [x, y, z]; m[ax] -= eps;
      f.at(p[0], p[1], p[2], Bp); f.at(m[0], m[1], m[2], Bm);
      d += (Bp[i] - Bm[i]) / (2 * eps);
    }
    return d;
  }
  function curlB(f, x, y, z) {
    const J = [];
    for (const ax of [0, 1, 2]) {
      const p = [x, y, z]; p[ax] += eps;
      const m = [x, y, z]; m[ax] -= eps;
      f.at(p[0], p[1], p[2], Bp); f.at(m[0], m[1], m[2], Bm);
      J.push([(Bp[0] - Bm[0]) / (2 * eps), (Bp[1] - Bm[1]) / (2 * eps), (Bp[2] - Bm[2]) / (2 * eps)]);
    }
    return [J[1][2] - J[2][1], J[2][0] - J[0][2], J[0][1] - J[1][0]];
  }
  const rng = coll.makeRng(42);
  for (const model of ['screw', 'beltrami']) {
    const f = makeField({ model, Bwall: 0.01, ratio: 1.2, R });
    let maxDiv = 0, maxCurlErr = 0;
    for (let k = 0; k < 200; k++) {
      const r = 0.95 * R * Math.sqrt(rng()), th = rng() * 2 * Math.PI;
      const x = r * Math.cos(th), y = r * Math.sin(th), z = rng() * 0.4;
      maxDiv = Math.max(maxDiv, Math.abs(divB(f, x, y, z)));
      if (model === 'beltrami') {
        const cb = curlB(f, x, y, z);
        f.at(x, y, z, Bp);
        const a = f.meta.alpha;
        const err = Math.hypot(cb[0] - a * Bp[0], cb[1] - a * Bp[1], cb[2] - a * Bp[2]) /
          (a * Math.hypot(Bp[0], Bp[1], Bp[2]));
        maxCurlErr = Math.max(maxCurlErr, err);
      }
    }
    check(`div B = 0 (${model})`, maxDiv < 1e-6 * 0.01 / R,
      `max |divB| = ${maxDiv.toExponential(1)} T/m`);
    if (model === 'beltrami') {
      check('Beltrami force-free (curl B = aB)', maxCurlErr < 1e-3,
        `max rel err ${maxCurlErr.toExponential(1)}`);
    }
  }
}

// ---- 3. |B(R)| held fixed across the pitch scan ----
{
  const R = 0.1, Bwall = 0.01, out = new Float64Array(3);
  let worst = 0;
  for (const model of ['screw', 'beltrami']) {
    for (const ratio of [0.01, 0.5, 1, 2, 3.5]) {
      const f = makeField({ model, Bwall, ratio, R });
      f.at(R, 0, 0.1, out);
      const mag = Math.hypot(out[0], out[1], out[2]);
      worst = Math.max(worst, Math.abs(mag / Bwall - 1));
      const rTest = ratio > 0.01 ? out[1] / out[2] : 0; // Bθ/Bz at (R,0)
      worst = Math.max(worst, ratio > 0.01 ? Math.abs(rTest / ratio - 1) : 0);
    }
  }
  check('|B(R)| and Bth/Bz(R) match requested values', worst < 1e-6,
    `worst rel err ${worst.toExponential(1)}`);
}

// ---- 4. MCC electron collision rate vs analytic nu(E) ----
{
  const nn = C.neutralDensity(10, 300);
  const E = 3; // eV — below ionization threshold, elastic only
  const v = Math.sqrt((2 * E * C.EV) / C.ME);
  const nuTheory = nn * C.sigmaElectronElastic(E) * v;
  const dt = 0.02 / nuTheory;
  const rng = coll.makeRng(7);
  const collider = coll.makeElectronCollider(nn, dt, rng);
  const vel = new Float64Array(3);
  const M = 400000;
  let hits = 0;
  for (let s = 0; s < M; s++) {
    // reset speed each step so E stays fixed despite scattering losses
    vel[0] = v; vel[1] = 0; vel[2] = 0;
    if (collider.collide(vel, 0) !== 0) hits++;
  }
  const nuMeas = hits / (M * dt);
  check('MCC elastic rate matches nu(E)', Math.abs(nuMeas / nuTheory - 1) < 0.05,
    `measured/theory = ${(nuMeas / nuTheory).toFixed(3)}, nu = ${nuTheory.toExponential(2)} /s`);
}

console.log(failures ? `\n${failures} check(s) FAILED` : '\nAll checks passed.');
process.exit(failures ? 1 : 0);
