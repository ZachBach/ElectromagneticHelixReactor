/* Module 2 — Boris velocity rotation for a static magnetic field (E = 0).
   Exactly energy-conserving; phase error O((ω_c·dt)²), so the runner keeps
   ω_c·dt ≤ 0.2. Operates in place on a flat [vx,vy,vz,...] typed array. */
'use strict';

// h = (q/m)·dt/2, SI. B components in tesla.
function borisRotate(vel, i3, Bx, By, Bz, h) {
  const tx = h * Bx, ty = h * By, tz = h * Bz;
  const t2 = tx * tx + ty * ty + tz * tz;
  const sf = 2 / (1 + t2);
  const sx = tx * sf, sy = ty * sf, sz = tz * sf;
  const vx = vel[i3], vy = vel[i3 + 1], vz = vel[i3 + 2];
  const ux = vx + (vy * tz - vz * ty);
  const uy = vy + (vz * tx - vx * tz);
  const uz = vz + (vx * ty - vy * tx);
  vel[i3] = vx + (uy * sz - uz * sy);
  vel[i3 + 1] = vy + (uz * sx - ux * sz);
  vel[i3 + 2] = vz + (ux * sy - uy * sx);
}

module.exports = { borisRotate };
