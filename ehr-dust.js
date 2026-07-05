/* EHR dusty plasma module — charged dust in the EHR plasma.
   Forces: gravity, wall/end sheath confinement (∝ plasma density), ion drag from −∇n_i
   (sampled from the live ion density grid) + axial flow toward the chuck, acoustic
   radiation force from the neutral modulation, Epstein neutral drag, and pairwise
   Yukawa (screened Coulomb) repulsion. Diagnostics: kinetic temperature T_d, mean
   nearest-neighbor spacing a, κ=a/λ_D, coupling parameter Γ, state (gas/liquid/crystal),
   and pair-correlation g(r). */
(function () {
  'use strict';
  const TWO_PI = Math.PI * 2;
  let _g = null;
  function gauss() {
    if (_g !== null) { const v = _g; _g = null; return v; }
    let u = 0, v = 0, s = 0;
    do { u = Math.random() * 2 - 1; v = Math.random() * 2 - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
    const m = Math.sqrt(-2 * Math.log(s) / s);
    _g = v * m; return u * m;
  }

  class EHRDust {
    constructor(engine) {
      this.e = engine;
      this.cap = 1200;
      this.pos = new Float32Array(this.cap * 3);
      this.vel = new Float32Array(this.cap * 3);
      this.n = 0;
      this.gofr = new Float32Array(48);
      this._gAcc = new Float32Array(48);
      this._nn = new Float32Array(this.cap);
      this._F = new Float32Array(this.cap * 3);
      this.gamma = 0; this.kappa = 0; this.Td = 0; this.spacing = 0;
      this.state = 'OFF';
    }

    setN(n) {
      n = Math.min(this.cap, Math.max(0, n | 0));
      for (let i = this.n; i < n; i++) this.spawn(i);
      this.n = n;
      if (!n) this.state = 'OFF';
    }

    reseed() {
      for (let i = 0; i < this.n; i++) this.spawn(i);
      this.gofr.fill(0);
      this.gamma = 0; this.kappa = 0; this.Td = 0; this.spacing = 0;
      this.state = this.n ? 'GAS' : 'OFF';
    }
    spawn(i) {
      this.pos[3 * i] = (Math.random() * 2 - 1) * 0.55;
      this.pos[3 * i + 1] = -0.05 - Math.random() * 0.6;
      this.pos[3 * i + 2] = (Math.random() * 2 - 1) * 1.5;
      this.vel[3 * i] = gauss() * 0.02;
      this.vel[3 * i + 1] = gauss() * 0.02;
      this.vel[3 * i + 2] = gauss() * 0.02;
    }

    update(dtIn) {
      const e = this.e, p = e.params;
      if ((p.dustN | 0) !== this.n) this.setN(p.dustN);
      const n = this.n;
      if (!n) { this.state = 'OFF'; this.gamma = 0; return; }
      const dt = Math.min(Math.max(dtIn, 0.001), 0.04);
      const pos = this.pos, vel = this.vel, F = this._F, nnD = this._nn;
      const q = p.dustQ / 4;
      const lam = Math.max(p.dustLam, 0.02);
      const cut = lam * 6, cut2 = cut * cut;
      const nu = 0.9 * (p.pressure / 3);
      const nRel = e._nCal ? Math.min(e.diag.nMean / Math.max(e._nCal, 1e-4), 2) : 0.5;
      const S0 = 2.4 * q * (0.25 + nRel);
      const grav = 0.4;
      const kac = TWO_PI / p.acLambda;
      const ph = TWO_PI * p.acFreq * e.tReal;
      const acF = 0.7 * p.acA;
      const standing = p.acMode === 'standing';
      const cosPh = Math.cos(ph);
      const g = e.grid;
      const gs = e.gridMax > 0.002 ? 0.9 / e.gridMax : 0;
      const KY = 0.02 * q * q;

      // ---- pair pass: Yukawa forces, nearest-neighbor, g(r) histogram ----
      nnD.fill(9, 0, n);
      const gAcc = this._gAcc; gAcc.fill(0);
      F.fill(0, 0, n * 3);
      for (let i = 0; i < n; i++) {
        const i3 = 3 * i;
        const xi = pos[i3], yi = pos[i3 + 1], zi = pos[i3 + 2];
        for (let j = i + 1; j < n; j++) {
          const j3 = 3 * j;
          const dx = xi - pos[j3], dy = yi - pos[j3 + 1], dz = zi - pos[j3 + 2];
          const d2 = dx * dx + dy * dy + dz * dz;
          if (d2 < 0.25) {
            const b = (Math.sqrt(d2) * 96) | 0;
            if (b < 48) gAcc[b]++;
          }
          if (d2 > cut2 || d2 < 1e-8) continue;
          const d = Math.sqrt(d2);
          if (d < nnD[i]) nnD[i] = d;
          if (d < nnD[j]) nnD[j] = d;
          const f = KY * (1 / d2) * (1 + d / lam) * Math.exp(-d / lam) / d;
          const fx = f * dx, fy = f * dy, fz = f * dz;
          F[i3] += fx; F[i3 + 1] += fy; F[i3 + 2] += fz;
          F[j3] -= fx; F[j3 + 1] -= fy; F[j3 + 2] -= fz;
        }
      }

      // ---- single pass: body forces + integrate ----
      let vxm = 0, vym = 0, vzm = 0;
      for (let i = 0; i < n; i++) {
        const i3 = 3 * i;
        let x = pos[i3], y = pos[i3 + 1], z = pos[i3 + 2];
        let fx = F[i3], fy = F[i3 + 1] - grav, fz = F[i3 + 2];
        const r = Math.sqrt(x * x + y * y) + 1e-9;
        // sheath confinement (radial wall + ends), strength tied to plasma density
        const sw = S0 * Math.exp(-(1 - Math.min(r, 1)) / 0.09);
        fx -= sw * x / r; fy -= sw * y / r;
        fz -= S0 * Math.exp(-(2 - Math.min(z, 2)) / 0.09);
        fz += S0 * Math.exp(-(Math.max(z, -2) + 2) / 0.09);
        // ion drag: −∇n from live density grid + net flow toward chuck
        if (gs) {
          const ix = Math.min(62, Math.max(1, ((x + 1) * 32) | 0));
          const iy = Math.min(62, Math.max(1, ((y + 1) * 32) | 0));
          const iz = Math.min(126, Math.max(1, ((z + 2) * 32) | 0));
          const idx = (iz * 64 + iy) * 64 + ix;
          const gd = 0.55 * q * gs * 16;
          fx -= gd * (g[idx + 1] - g[idx - 1]);
          fy -= gd * (g[idx + 64] - g[idx - 64]);
          fz -= gd * (g[idx + 4096] - g[idx - 4096]);
          fz += 0.12 * q * nRel;
        }
        // acoustic radiation force (toward pressure nodes)
        const dnn = standing ? kac * Math.cos(kac * z) * cosPh : kac * Math.cos(kac * z - ph);
        fz -= acF * dnn;
        // semi-implicit integrate with Epstein drag
        let vx = vel[i3], vy = vel[i3 + 1], vz = vel[i3 + 2];
        const dd = 1 + nu * dt;
        vx = (vx + fx * dt) / dd;
        vy = (vy + fy * dt) / dd;
        vz = (vz + fz * dt) / dd;
        x += vx * dt; y += vy * dt; z += vz * dt;
        const r2b = x * x + y * y;
        if (r2b > 0.9216) { const rb = Math.sqrt(r2b); x *= 0.955 / rb; y *= 0.955 / rb; vx *= -0.3; vy *= -0.3; }
        if (z > 1.94) { z = 1.94; vz *= -0.3; }
        if (z < -1.94) { z = -1.94; vz *= -0.3; }
        pos[i3] = x; pos[i3 + 1] = y; pos[i3 + 2] = z;
        vel[i3] = vx; vel[i3 + 1] = vy; vel[i3 + 2] = vz;
        vxm += vx; vym += vy; vzm += vz;
      }
      vxm /= n; vym /= n; vzm /= n;

      // ---- diagnostics ----
      let ke = 0;
      for (let i = 0; i < n; i++) {
        const i3 = 3 * i;
        const ax = vel[i3] - vxm, ay = vel[i3 + 1] - vym, az = vel[i3 + 2] - vzm;
        ke += ax * ax + ay * ay + az * az;
      }
      this.Td = 0.5 * ke / n + 1e-6;
      let sa = 0, ca = 0;
      for (let i = 0; i < n; i++) if (nnD[i] < 9) { sa += nnD[i]; ca++; }
      this.spacing = ca ? sa / ca : 0.1;
      this.kappa = this.spacing / lam;
      this.gamma = Math.min(999, 0.012 * q * q * Math.exp(-this.kappa) / (Math.max(this.spacing, 1e-3) * this.Td));
      this.state = this.gamma > 170 ? 'CRYSTAL' : this.gamma > 30 ? 'LIQUID' : 'GAS';
      for (let b = 0; b < 48; b++) {
        const rr = (b + 0.5) / 96;
        this.gofr[b] = gAcc[b] / (n * rr);
      }
    }
  }
  window.EHRDust = EHRDust;
})();
