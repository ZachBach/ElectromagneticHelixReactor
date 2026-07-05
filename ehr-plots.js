/* EHR Simulator — tiny canvas plot lib */
(function () {
  'use strict';

  const FONT = '10px "IBM Plex Mono", monospace';
  const COL = {
    grid: '#101c28', axis: '#22384c', label: '#5b7186',
    accent: '#3fd9ff', amber: '#ffb454', green: '#5cf2a6', red: '#ff5d6c',
  };

  function fit(canvas) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = canvas.clientWidth | 0, h = canvas.clientHeight | 0;
    if (w < 4 || h < 4) return null;
    if (canvas.width !== w * dpr || canvas.height !== h * dpr) {
      canvas.width = w * dpr; canvas.height = h * dpr;
    }
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, w, h };
  }

  function niceTicks(min, max, n) {
    const span = max - min || 1;
    const step0 = span / Math.max(1, n);
    const mag = Math.pow(10, Math.floor(Math.log10(step0)));
    let step = mag;
    for (const m of [1, 2, 2.5, 5, 10]) { if (mag * m >= step0) { step = mag * m; break; } }
    const ticks = [];
    for (let v = Math.ceil(min / step) * step; v <= max + step * 1e-6; v += step) ticks.push(v);
    return ticks;
  }

  function fmt(v) {
    const a = Math.abs(v);
    if (a >= 1000) return v.toExponential(1);
    if (a >= 100) return v.toFixed(0);
    if (a >= 1) return +v.toFixed(2) + '';
    if (a >= 0.01) return +v.toFixed(3) + '';
    if (a === 0) return '0';
    return v.toExponential(1);
  }

  // opts: {series:[{x,y,color,width,dots,dotR,alpha,dash}], xLog, xLabel, yLabel,
  //        xMin,xMax,yMin,yMax, vlines:[{x,color,label,dash}], curX, yPad}
  function chart(canvas, opts) {
    const f = fit(canvas); if (!f) return;
    const { ctx, w, h } = f;
    ctx.clearRect(0, 0, w, h);
    const padL = 44, padR = 10, padT = 8, padB = opts.xLabel ? 28 : 16;
    const pw = w - padL - padR, ph = h - padT - padB;
    if (pw < 10 || ph < 10) return;

    const series = (opts.series || []).filter(s => s && s.x && s.x.length > 0);
    let xMin = opts.xMin, xMax = opts.xMax, yMin = opts.yMin, yMax = opts.yMax;
    if (xMin == null || xMax == null || yMin == null || yMax == null) {
      let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
      for (const s of series) for (let i = 0; i < s.x.length; i++) {
        const xv = s.x[i], yv = s.y[i];
        if (!isFinite(xv) || !isFinite(yv)) continue;
        if (xv < x0) x0 = xv; if (xv > x1) x1 = xv;
        if (yv < y0) y0 = yv; if (yv > y1) y1 = yv;
      }
      if (!isFinite(x0)) { x0 = 0; x1 = 1; y0 = 0; y1 = 1; }
      if (xMin == null) xMin = x0; if (xMax == null) xMax = x1;
      if (yMin == null) yMin = y0; if (yMax == null) yMax = y1;
    }
    if (xMax - xMin < 1e-12) { xMax = xMin + 1; }
    if (yMax - yMin < 1e-12) { yMax = yMin + (Math.abs(yMin) || 1) * 0.5; yMin -= (Math.abs(yMin) || 1) * 0.1; }
    const padY = (opts.yPad != null ? opts.yPad : 0.12) * (yMax - yMin);
    yMax += padY; if (opts.yMin == null) yMin -= padY * 0.4;

    const lx = opts.xLog;
    const X = v => padL + (lx
      ? (Math.log(Math.max(v, 1e-9)) - Math.log(xMin)) / (Math.log(xMax) - Math.log(xMin))
      : (v - xMin) / (xMax - xMin)) * pw;
    const Y = v => padT + (1 - (v - yMin) / (yMax - yMin)) * ph;

    // grid + ticks
    ctx.font = FONT; ctx.textBaseline = 'middle';
    ctx.strokeStyle = COL.grid; ctx.fillStyle = COL.label; ctx.lineWidth = 1;
    const yt = niceTicks(yMin, yMax, 4);
    for (const v of yt) {
      const y = Y(v);
      if (y < padT - 1 || y > padT + ph + 1) continue;
      ctx.beginPath(); ctx.moveTo(padL, y); ctx.lineTo(padL + pw, y); ctx.stroke();
      ctx.textAlign = 'right'; ctx.fillText(fmt(v), padL - 5, y);
    }
    const xt = lx ? [0.25, 0.5, 1, 2, 4].filter(v => v >= xMin && v <= xMax) : niceTicks(xMin, xMax, 5);
    ctx.textAlign = 'center';
    for (const v of xt) {
      const x = X(v);
      ctx.strokeStyle = COL.grid;
      ctx.beginPath(); ctx.moveTo(x, padT); ctx.lineTo(x, padT + ph); ctx.stroke();
      ctx.fillText(fmt(v), x, padT + ph + 9);
    }
    ctx.strokeStyle = COL.axis;
    ctx.strokeRect(padL + 0.5, padT + 0.5, pw - 1, ph - 1);

    // vlines
    for (const vl of (opts.vlines || [])) {
      if (vl.x < xMin || vl.x > xMax) continue;
      const x = X(vl.x);
      ctx.strokeStyle = vl.color || COL.amber; ctx.lineWidth = 1;
      ctx.setLineDash(vl.dash || [4, 4]);
      ctx.beginPath(); ctx.moveTo(x, padT); ctx.lineTo(x, padT + ph); ctx.stroke();
      ctx.setLineDash([]);
      if (vl.label) { ctx.fillStyle = vl.color || COL.amber; ctx.textAlign = 'left'; ctx.fillText(vl.label, x + 4, padT + 8); }
    }

    // series
    ctx.save();
    ctx.beginPath(); ctx.rect(padL, padT, pw, ph); ctx.clip();
    for (const s of series) {
      ctx.globalAlpha = s.alpha != null ? s.alpha : 1;
      ctx.strokeStyle = s.color || COL.accent;
      ctx.fillStyle = s.color || COL.accent;
      ctx.lineWidth = s.width || 1.5;
      if (s.dash) ctx.setLineDash(s.dash); else ctx.setLineDash([]);
      if (!s.dotsOnly) {
        ctx.beginPath();
        let started = false;
        for (let i = 0; i < s.x.length; i++) {
          if (!isFinite(s.y[i])) { started = false; continue; }
          const x = X(s.x[i]), y = Y(s.y[i]);
          if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
        }
        ctx.stroke();
      }
      if (s.dots || s.dotsOnly) {
        const r = s.dotR || 2;
        for (let i = 0; i < s.x.length; i++) {
          if (!isFinite(s.y[i])) continue;
          ctx.beginPath(); ctx.arc(X(s.x[i]), Y(s.y[i]), r, 0, 6.2832); ctx.fill();
        }
      }
      ctx.setLineDash([]);
    }
    ctx.globalAlpha = 1;
    // current-x marker
    if (opts.curX != null && isFinite(opts.curX) && opts.curX >= xMin && opts.curX <= xMax) {
      const x = X(opts.curX);
      ctx.strokeStyle = 'rgba(92,242,166,0.7)'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x, padT); ctx.lineTo(x, padT + ph); ctx.stroke();
    }
    ctx.restore();

    ctx.fillStyle = COL.label; ctx.font = FONT;
    if (opts.xLabel) { ctx.textAlign = 'center'; ctx.fillText(opts.xLabel, padL + pw / 2, h - 7); }
    if (opts.yLabel) {
      ctx.save(); ctx.translate(10, padT + ph / 2); ctx.rotate(-Math.PI / 2);
      ctx.textAlign = 'center'; ctx.fillText(opts.yLabel, 0, 0); ctx.restore();
    }
  }

  // Palettes: t in [0,1] -> [r,g,b]
  function makeRamp(stops) {
    return function (t) {
      t = Math.max(0, Math.min(1, t));
      for (let i = 1; i < stops.length; i++) {
        if (t <= stops[i][0]) {
          const [t0, c0] = stops[i - 1], [t1, c1] = stops[i];
          const u = (t - t0) / (t1 - t0 || 1);
          return [c0[0] + (c1[0] - c0[0]) * u, c0[1] + (c1[1] - c0[1]) * u, c0[2] + (c1[2] - c0[2]) * u];
        }
      }
      return stops[stops.length - 1][1];
    };
  }
  const inferno = makeRamp([
    [0.0, [4, 2, 18]], [0.18, [44, 16, 92]], [0.38, [113, 31, 129]],
    [0.58, [182, 54, 121]], [0.75, [238, 96, 94]], [0.9, [253, 174, 50]], [1.0, [252, 255, 164]],
  ]);
  const plasma5 = makeRamp([
    [0.0, [6, 16, 60]], [0.28, [30, 109, 255]], [0.5, [25, 213, 255]],
    [0.7, [61, 255, 136]], [0.86, [255, 233, 74]], [1.0, [255, 255, 255]],
  ]);
  const ice = makeRamp([
    [0.0, [5, 12, 40]], [0.35, [24, 74, 190]], [0.65, [66, 165, 245]],
    [0.85, [160, 226, 255]], [1.0, [255, 255, 255]],
  ]);
  const ember = makeRamp([
    [0.0, [22, 6, 6]], [0.3, [140, 32, 18]], [0.55, [235, 100, 24]],
    [0.8, [255, 190, 60]], [1.0, [255, 250, 220]],
  ]);
  const PALETTES = { inferno, plasma: plasma5, ice, ember };

  // Heatmap: data Float32Array (w*h), row-major y*w+x
  function heat(canvas, data, dw, dh, opts) {
    const f = fit(canvas); if (!f) return;
    const { ctx, w, h } = f;
    opts = opts || {};
    const pal = PALETTES[opts.palette || 'inferno'] || inferno;
    let max = opts.max;
    if (!max) { max = 1e-9; for (let i = 0; i < data.length; i++) if (data[i] > max) max = data[i]; }
    const off = heat._off || (heat._off = document.createElement('canvas'));
    off.width = dw; off.height = dh;
    const octx = off.getContext('2d');
    const img = octx.createImageData(dw, dh);
    const px = img.data;
    for (let i = 0; i < dw * dh; i++) {
      const t = Math.pow(Math.min(1, data[i] / max), opts.gamma || 0.6);
      const c = pal(t);
      px[i * 4] = c[0]; px[i * 4 + 1] = c[1]; px[i * 4 + 2] = c[2]; px[i * 4 + 3] = 255;
    }
    octx.putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = 'high';
    ctx.clearRect(0, 0, w, h);
    // keep square aspect, centered
    const s = Math.min(w / dw, h / dh);
    const rw = dw * s, rh = dh * s;
    ctx.drawImage(off, (w - rw) / 2, (h - rh) / 2, rw, rh);
    if (opts.circleMask) {
      ctx.save();
      ctx.globalCompositeOperation = 'destination-in';
      ctx.beginPath(); ctx.arc(w / 2, h / 2, Math.min(rw, rh) / 2, 0, 6.2832); ctx.fill();
      ctx.restore();
      ctx.strokeStyle = '#22384c'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(w / 2, h / 2, Math.min(rw, rh) / 2 - 0.5, 0, 6.2832); ctx.stroke();
    }
  }

  window.EHRPlots = { chart, heat, palettes: PALETTES, COL };
})();
