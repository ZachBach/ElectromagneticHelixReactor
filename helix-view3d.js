/* EHR Simulator — Three.js WebGL view: chamber, antenna, field lines, particles,
   ray-marched density volume + nested isosurfaces. THREE is injected (dynamic import). */
(function () {
  'use strict';

  const VERT = `
    out vec3 vPos;
    void main() {
      vPos = position;
      gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    }`;

  const FRAG = `
    precision highp float;
    precision highp sampler3D;
    in vec3 vPos;
    out vec4 fragColor;
    uniform sampler3D uTex;
    uniform vec3 uCamObj;
    uniform float uMode;      // 0 volume, 1 iso, 2 both
    uniform float uAlpha;
    uniform float uGain;
    uniform int uCmap;

    vec3 tf(float t) {
      vec3 c;
      if (uCmap == 1) { // ice
        c = mix(vec3(0.02,0.05,0.16), vec3(0.09,0.29,0.75), smoothstep(0.0,0.35,t));
        c = mix(c, vec3(0.26,0.65,0.96), smoothstep(0.35,0.65,t));
        c = mix(c, vec3(0.63,0.89,1.0), smoothstep(0.65,0.85,t));
        c = mix(c, vec3(1.0), smoothstep(0.85,1.0,t));
      } else if (uCmap == 2) { // ember
        c = mix(vec3(0.09,0.02,0.02), vec3(0.55,0.13,0.07), smoothstep(0.0,0.3,t));
        c = mix(c, vec3(0.92,0.39,0.09), smoothstep(0.3,0.55,t));
        c = mix(c, vec3(1.0,0.75,0.24), smoothstep(0.55,0.8,t));
        c = mix(c, vec3(1.0,0.98,0.86), smoothstep(0.8,1.0,t));
      } else { // plasma: blue→cyan→green→yellow→white
        c = mix(vec3(0.02,0.06,0.24), vec3(0.12,0.43,1.0), smoothstep(0.0,0.28,t));
        c = mix(c, vec3(0.10,0.84,1.0), smoothstep(0.28,0.5,t));
        c = mix(c, vec3(0.24,1.0,0.53), smoothstep(0.5,0.7,t));
        c = mix(c, vec3(1.0,0.91,0.29), smoothstep(0.7,0.86,t));
        c = mix(c, vec3(1.0), smoothstep(0.86,1.0,t));
      }
      return c;
    }

    float dens(vec3 p) {
      if (p.x*p.x + p.y*p.y > 1.0) return 0.0;
      vec3 tc = vec3((p.x+1.0)*0.5, (p.y+1.0)*0.5, (p.z+2.0)*0.25);
      return texture(uTex, tc).r * uGain;
    }

    vec3 grad(vec3 p) {
      const float e = 0.03;
      return vec3(
        dens(p+vec3(e,0,0)) - dens(p-vec3(e,0,0)),
        dens(p+vec3(0,e,0)) - dens(p-vec3(0,e,0)),
        dens(p+vec3(0,0,e)) - dens(p-vec3(0,0,e)));
    }

    void main() {
      vec3 dir = normalize(vPos - uCamObj);
      // slab intersect box [-1,1]x[-1,1]x[-2,2]
      vec3 bmin = vec3(-1.0,-1.0,-2.0), bmax = vec3(1.0,1.0,2.0);
      vec3 inv = 1.0 / dir;
      vec3 t0s = (bmin - uCamObj) * inv, t1s = (bmax - uCamObj) * inv;
      vec3 tsm = min(t0s, t1s), tbg = max(t0s, t1s);
      float t0 = max(max(tsm.x, tsm.y), tsm.z);
      float t1 = min(min(tbg.x, tbg.y), tbg.z);
      t0 = max(t0, 0.0);
      if (t1 <= t0) discard;

      const int STEPS = 88;
      float dt = (t1 - t0) / float(STEPS);
      vec4 acc = vec4(0.0);
      float prev = -1.0;
      float th[4]; th[0]=0.25; th[1]=0.5; th[2]=0.75; th[3]=0.9;
      vec3 shellCol[4];
      shellCol[0]=vec3(0.15,0.4,1.0); shellCol[1]=vec3(0.2,1.0,0.5);
      shellCol[2]=vec3(1.0,0.9,0.25); shellCol[3]=vec3(1.0);
      float shellA[4]; shellA[0]=0.22; shellA[1]=0.32; shellA[2]=0.45; shellA[3]=0.85;

      for (int i = 0; i < STEPS; i++) {
        float t = t0 + (float(i) + 0.5) * dt;
        vec3 p = uCamObj + dir * t;
        float d = dens(p);
        if (uMode < 0.5 || uMode > 1.5) { // volume
          float a = d * d * uAlpha * dt * 30.0;
          a = clamp(a, 0.0, 0.9);
          vec3 col = tf(d);
          acc.rgb += (1.0 - acc.a) * col * a;
          acc.a   += (1.0 - acc.a) * a;
        }
        if (uMode > 0.5 && prev >= 0.0) { // iso shells
          for (int s = 0; s < 4; s++) {
            if ((prev - th[s]) * (d - th[s]) < 0.0) {
              vec3 n = normalize(grad(p) + vec3(1e-5));
              if (dot(n, dir) > 0.0) n = -n;
              vec3 l = -dir;
              float dif = max(dot(n, l), 0.0);
              float spec = pow(max(dot(reflect(-l, n), -dir), 0.0), 28.0);
              vec3 col = shellCol[s] * (0.22 + 0.78 * dif) + vec3(spec * 0.5);
              float a = shellA[s];
              acc.rgb += (1.0 - acc.a) * col * a;
              acc.a   += (1.0 - acc.a) * a;
            }
          }
        }
        prev = d;
        if (acc.a > 0.96) break;
      }
      if (acc.a < 0.004) discard;
      fragColor = acc;
    }`;

  class EHRView3D {
    constructor(THREE, canvas, engine) {
      this.T = THREE; this.canvas = canvas; this.engine = engine;
      this.orbit = { theta: 1.05, phi: 1.12, dist: 5.4 };
      this.layers = { lines: true, particles: true, volume: true, iso: false, antenna: true, dust: true, chuck: true };
      this.frame = 0;

      const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: false, powerPreference: 'high-performance' });
      renderer.setClearColor(0x04070b, 1);
      renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 1.6));
      this.renderer = renderer;

      this.scene = new THREE.Scene();
      this.camera = new THREE.PerspectiveCamera(42, 1, 0.05, 80);
      this.group = new THREE.Group();
      this.group.rotation.y = Math.PI / 2; // chamber axis along screen-X
      this.scene.add(this.group);

      const NG = window.EHREngine.GRID;
      this.texData = new Uint8Array(NG.NX * NG.NY * NG.NZ);
      const tex = new THREE.Data3DTexture(this.texData, NG.NX, NG.NY, NG.NZ);
      tex.format = THREE.RedFormat; tex.type = THREE.UnsignedByteType;
      tex.minFilter = THREE.LinearFilter; tex.magFilter = THREE.LinearFilter;
      tex.unpackAlignment = 1; tex.needsUpdate = true;
      this.tex3d = tex;

      this.volMat = new THREE.ShaderMaterial({
        glslVersion: THREE.GLSL3,
        vertexShader: VERT, fragmentShader: FRAG,
        uniforms: {
          uTex: { value: tex }, uCamObj: { value: new THREE.Vector3() },
          uMode: { value: 0 }, uAlpha: { value: 1.0 }, uGain: { value: 1.05 }, uCmap: { value: 0 },
        },
        transparent: true, depthWrite: false, side: THREE.BackSide,
      });
      this.volMesh = new THREE.Mesh(new THREE.BoxGeometry(2, 2, 4), this.volMat);
      this.volMesh.renderOrder = 10;
      this.group.add(this.volMesh);

      this._buildChamber();
      this._buildParticles();
      this._buildDust();
      this.rebuildGeometry();
      this._bindOrbit();
      this.setColormap(0);
      this.resize();
    }

    _buildChamber() {
      const T = this.T;
      const pts = [];
      const SEG = 72;
      for (const z of [-2, 2]) {
        for (let i = 0; i < SEG; i++) {
          const a0 = i / SEG * Math.PI * 2, a1 = (i + 1) / SEG * Math.PI * 2;
          pts.push(1.02 * Math.cos(a0), 1.02 * Math.sin(a0), z, 1.02 * Math.cos(a1), 1.02 * Math.sin(a1), z);
        }
      }
      for (let i = 0; i < 12; i++) {
        const a = i / 12 * Math.PI * 2;
        pts.push(1.02 * Math.cos(a), 1.02 * Math.sin(a), -2, 1.02 * Math.cos(a), 1.02 * Math.sin(a), 2);
      }
      const g = new T.BufferGeometry();
      g.setAttribute('position', new T.Float32BufferAttribute(pts, 3));
      const mat = new T.LineBasicMaterial({ color: 0x1c2e40, transparent: true, opacity: 0.85 });
      this.chamber = new T.LineSegments(g, mat);
      this.group.add(this.chamber);
      // wafer chuck at z=+2: disk + concentric rings
      const chuck = new T.Group();
      const disk = new T.Mesh(
        new T.CircleGeometry(0.78, 48),
        new T.MeshBasicMaterial({ color: 0x101a24, transparent: true, opacity: 0.92, side: T.DoubleSide })
      );
      disk.position.z = 1.985;
      chuck.add(disk);
      for (const rr of [0.2, 0.42, 0.64, 0.78]) {
        const rp = [];
        for (let i = 0; i <= 48; i++) { const a = i / 48 * Math.PI * 2; rp.push(new T.Vector3(rr * Math.cos(a), rr * Math.sin(a), 1.98)); }
        chuck.add(new T.Line(new T.BufferGeometry().setFromPoints(rp), new T.LineBasicMaterial({ color: 0x2c4258, transparent: true, opacity: 0.9 })));
      }
      this.chuckObj = chuck;
      this.group.add(chuck);
    }

    _buildDust() {
      const T = this.T;
      this.dustCap = 1200;
      this.dPos = new Float32Array(this.dustCap * 3);
      const g = new T.BufferGeometry();
      this.dPosAttr = new T.BufferAttribute(this.dPos, 3); this.dPosAttr.setUsage(T.DynamicDrawUsage);
      g.setAttribute('position', this.dPosAttr);
      g.setDrawRange(0, 0);
      const mat = new T.PointsMaterial({ color: 0xffd9a0, size: 0.032, sizeAttenuation: true, transparent: true, opacity: 0.95, depthWrite: false });
      this.dustObj = new T.Points(g, mat);
      this.dustObj.renderOrder = 6;
      this.dustObj.frustumCulled = false;
      this.group.add(this.dustObj);
    }

    _buildParticles() {
      const T = this.T, N = this.engine.params.N;
      this.pCount = Math.min(N, 22000);
      this.pPos = new Float32Array(this.pCount * 3);
      this.pCol = new Float32Array(this.pCount * 3);
      const g = new T.BufferGeometry();
      this.pPosAttr = new T.BufferAttribute(this.pPos, 3); this.pPosAttr.setUsage(T.DynamicDrawUsage);
      this.pColAttr = new T.BufferAttribute(this.pCol, 3); this.pColAttr.setUsage(T.DynamicDrawUsage);
      g.setAttribute('position', this.pPosAttr);
      g.setAttribute('color', this.pColAttr);
      const mat = new T.PointsMaterial({
        size: 0.014, sizeAttenuation: true, vertexColors: true,
        transparent: true, opacity: 0.8, blending: T.AdditiveBlending, depthWrite: false,
      });
      this.points = new T.Points(g, mat);
      this.points.renderOrder = 5;
      this.points.frustumCulled = false;
      this.group.add(this.points);
    }

    rebuildParticleBuffers() {
      this.group.remove(this.points);
      this.points.geometry.dispose(); this.points.material.dispose();
      this._buildParticles();
    }

    rebuildGeometry() {
      const T = this.T, e = this.engine;
      // field lines
      if (this.linesObj) { this.group.remove(this.linesObj); this.linesObj.geometry.dispose(); this.linesObj.material.dispose(); }
      const lines = e.fieldLines(72);
      let segCount = 0;
      for (const l of lines) segCount += (l.pts.length / 3 - 1);
      const posArr = new Float32Array(segCount * 6), colArr = new Float32Array(segCount * 6);
      let minB = Infinity, maxB = -Infinity;
      for (const l of lines) for (let i = 0; i < l.mag.length; i++) { const m = l.mag[i]; if (m < minB) minB = m; if (m > maxB) maxB = m; }
      const spanB = Math.max(maxB - minB, 1e-6);
      let o = 0;
      const c0 = [0.08, 0.22, 0.36], c1 = [0.31, 0.88, 1.0];
      for (const l of lines) {
        const n = l.pts.length / 3;
        for (let i = 0; i < n - 1; i++) {
          for (let k = 0; k < 2; k++) {
            const idx = (i + k) * 3;
            posArr[o] = l.pts[idx]; posArr[o + 1] = l.pts[idx + 1]; posArr[o + 2] = l.pts[idx + 2];
            const t = (l.mag[i + k] - minB) / spanB;
            colArr[o] = c0[0] + (c1[0] - c0[0]) * t;
            colArr[o + 1] = c0[1] + (c1[1] - c0[1]) * t;
            colArr[o + 2] = c0[2] + (c1[2] - c0[2]) * t;
            o += 3;
          }
        }
      }
      const g = new T.BufferGeometry();
      g.setAttribute('position', new T.BufferAttribute(posArr, 3));
      g.setAttribute('color', new T.BufferAttribute(colArr, 3));
      const mat = new T.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.5 });
      this.linesObj = new T.LineSegments(g, mat);
      this.linesObj.visible = this.layers.lines;
      this.group.add(this.linesObj);

      // antenna helix straps: θ(z) = (kz/m)·z + offset, at r=1.09
      if (this.antObj) {
        this.group.remove(this.antObj);
        this.antObj.traverse(ch => { if (ch.geometry) ch.geometry.dispose(); if (ch.material) ch.material.dispose(); });
      }
      this.antObj = new T.Group();
      const m = Math.max(1, e.params.mMode), kz = e.params.kzA;
      for (let j = 0; j < m; j++) {
        const pts = [];
        for (let i = 0; i <= 60; i++) {
          const z = -1.25 + 2.5 * i / 60;
          const th = (kz / m) * z + j / m * Math.PI * 2;
          pts.push(new T.Vector3(1.09 * Math.cos(th), 1.09 * Math.sin(th), z));
        }
        const curve = new T.CatmullRomCurve3(pts);
        const tube = new T.Mesh(
          new T.TubeGeometry(curve, 60, 0.028, 6, false),
          new T.MeshBasicMaterial({ color: 0xffb454, transparent: true, opacity: 0.9 })
        );
        const glow = new T.Mesh(
          new T.TubeGeometry(curve, 60, 0.06, 6, false),
          new T.MeshBasicMaterial({ color: 0xff9d40, transparent: true, opacity: 0.16, depthWrite: false, blending: T.AdditiveBlending })
        );
        this.antObj.add(tube); this.antObj.add(glow);
      }
      // feed rings
      for (const z of [-1.25, 1.25]) {
        const ringPts = [];
        for (let i = 0; i <= 48; i++) { const a = i / 48 * Math.PI * 2; ringPts.push(new T.Vector3(1.09 * Math.cos(a), 1.09 * Math.sin(a), z)); }
        const rg = new T.BufferGeometry().setFromPoints(ringPts);
        this.antObj.add(new T.Line(rg, new T.LineBasicMaterial({ color: 0x8a5a2a, transparent: true, opacity: 0.7 })));
      }
      this.antObj.visible = this.layers.antenna;
      this.group.add(this.antObj);
    }

    setLayers(l) {
      Object.assign(this.layers, l);
      if (this.linesObj) this.linesObj.visible = this.layers.lines;
      if (this.points) this.points.visible = this.layers.particles;
      if (this.antObj) this.antObj.visible = this.layers.antenna;
      if (this.dustObj) this.dustObj.visible = this.layers.dust;
      if (this.chuckObj) this.chuckObj.visible = this.layers.chuck;
      const vol = this.layers.volume, iso = this.layers.iso;
      this.volMesh.visible = vol || iso;
      this.volMat.uniforms.uMode.value = vol && iso ? 2 : (iso ? 1 : 0);
    }

    setColormap(i) { this.volMat.uniforms.uCmap.value = i; this.cmapIdx = i; }

    _bindOrbit() {
      const c = this.canvas, o = this.orbit;
      let drag = false, lx = 0, ly = 0, pinchDist = 0;
      const pts = new Map(); // pointerId -> {x,y}, tracks active touches for pinch-zoom
      const pinchLen = () => { const [a, b] = [...pts.values()]; return Math.hypot(a.x - b.x, a.y - b.y); };
      this._onPointerDown = e => {
        pts.set(e.pointerId, { x: e.clientX, y: e.clientY });
        c.setPointerCapture(e.pointerId);
        if (pts.size === 1) { drag = true; lx = e.clientX; ly = e.clientY; }
        else if (pts.size === 2) { drag = false; pinchDist = pinchLen(); }
      };
      this._onPointerMove = e => {
        if (!pts.has(e.pointerId)) return;
        pts.set(e.pointerId, { x: e.clientX, y: e.clientY });
        if (pts.size >= 2) {
          const d = pinchLen();
          if (pinchDist > 0) o.dist = Math.min(14, Math.max(2.4, o.dist * (pinchDist / d)));
          pinchDist = d;
          return;
        }
        if (!drag) return;
        o.theta -= (e.clientX - lx) * 0.0055;
        o.phi = Math.min(Math.PI - 0.12, Math.max(0.12, o.phi - (e.clientY - ly) * 0.0055));
        lx = e.clientX; ly = e.clientY;
      };
      this._onPointerUp = e => {
        pts.delete(e.pointerId);
        if (pts.size === 1) { const [p] = [...pts.values()]; drag = true; lx = p.x; ly = p.y; pinchDist = 0; }
        else { drag = false; pinchDist = 0; }
      };
      this._onWheel = e => {
        e.preventDefault();
        o.dist = Math.min(14, Math.max(2.4, o.dist * Math.exp(e.deltaY * 0.0011)));
      };
      c.addEventListener('pointerdown', this._onPointerDown);
      c.addEventListener('pointermove', this._onPointerMove);
      c.addEventListener('pointerup', this._onPointerUp);
      c.addEventListener('pointercancel', this._onPointerUp);
      c.addEventListener('wheel', this._onWheel, { passive: false });
    }

    _disposeObj(obj) {
      if (!obj) return;
      obj.traverse(ch => {
        if (ch.geometry) ch.geometry.dispose();
        if (ch.material) (Array.isArray(ch.material) ? ch.material : [ch.material]).forEach(m => m.dispose());
      });
    }

    resize() {
      const w = this.canvas.clientWidth | 0, h = this.canvas.clientHeight | 0;
      if (w < 4 || h < 4) return;
      if (this._w !== w || this._h !== h) {
        this._w = w; this._h = h;
        this.renderer.setSize(w, h, false);
        this.camera.aspect = w / h;
        this.camera.updateProjectionMatrix();
      }
    }

    _energyColor(E, out, o) {
      // matches volume transfer roughly; E in eV
      const t = Math.min(1, E / 34);
      let r, g, b;
      if (t < 0.3) { const u = t / 0.3; r = 0.06 + 0.06 * u; g = 0.16 + 0.5 * u; b = 0.7 + 0.3 * u; }
      else if (t < 0.6) { const u = (t - 0.3) / 0.3; r = 0.12 + 0.2 * u; g = 0.66 + 0.34 * u; b = 1.0 - 0.5 * u; }
      else if (t < 0.85) { const u = (t - 0.6) / 0.25; r = 0.32 + 0.68 * u; g = 1.0 - 0.09 * u; b = 0.5 - 0.25 * u; }
      else { const u = (t - 0.85) / 0.15; r = 1.0; g = 0.91 + 0.09 * u; b = 0.25 + 0.75 * u; }
      out[o] = r; out[o + 1] = g; out[o + 2] = b;
    }

    update() {
      const e = this.engine, T = this.T;
      this.resize();
      this.frame++;

      if (this.pCount !== Math.min(e.params.N, 22000)) this.rebuildParticleBuffers();

      // particle positions every frame
      if (this.layers.particles) {
        this.pPos.set(e.pos.subarray(0, this.pCount * 3));
        this.pPosAttr.needsUpdate = true;
        if (this.frame % 3 === 0) {
          const v = e.vel;
          for (let i = 0; i < this.pCount; i++) {
            const i3 = 3 * i;
            const E = 5 * (v[i3] * v[i3] + v[i3 + 1] * v[i3 + 1] + v[i3 + 2] * v[i3 + 2]);
            this._energyColor(E, this.pCol, i3);
          }
          this.pColAttr.needsUpdate = true;
        }
      }

      // dust positions
      const dust = this.engine.dust;
      if (dust && this.layers.dust) {
        const dn = Math.min(dust.n, this.dustCap);
        this.dPos.set(dust.pos.subarray(0, dn * 3));
        this.dPosAttr.needsUpdate = true;
        this.dustObj.geometry.setDrawRange(0, dn);
      }

      // volume texture upload every 6 frames
      if ((this.layers.volume || this.layers.iso) && this.frame % 6 === 0) {
        const g = e.grid, td = this.texData;
        const scale = 255 / Math.max(e.gridMax, 1e-3);
        for (let i = 0; i < g.length; i++) {
          const v = g[i] * scale;
          td[i] = v > 255 ? 255 : v;
        }
        this.tex3d.needsUpdate = true;
      }

      // camera
      const o = this.orbit;
      const sp = Math.sin(o.phi), cp = Math.cos(o.phi);
      this.camera.position.set(
        o.dist * sp * Math.cos(o.theta),
        o.dist * cp,
        o.dist * sp * Math.sin(o.theta));
      this.camera.lookAt(0, 0, 0);
      this.camera.updateMatrixWorld();

      // camera in volume object space
      this.volMesh.updateMatrixWorld();
      const inv = new T.Matrix4().copy(this.volMesh.matrixWorld).invert();
      this.volMat.uniforms.uCamObj.value.copy(this.camera.position).applyMatrix4(inv);

      this.renderer.render(this.scene, this.camera);
    }

    dispose() {
      const c = this.canvas;
      if (this._onPointerDown) {
        c.removeEventListener('pointerdown', this._onPointerDown);
        c.removeEventListener('pointermove', this._onPointerMove);
        c.removeEventListener('pointerup', this._onPointerUp);
        c.removeEventListener('pointercancel', this._onPointerUp);
        c.removeEventListener('wheel', this._onWheel);
      }
      this._disposeObj(this.chamber);
      this._disposeObj(this.chuckObj);
      this._disposeObj(this.dustObj);
      this._disposeObj(this.points);
      this._disposeObj(this.linesObj);
      this._disposeObj(this.antObj);
      this._disposeObj(this.volMesh);
      if (this.tex3d) this.tex3d.dispose();
      this.renderer.dispose();
    }
  }

  window.EHRView3D = EHRView3D;
})();
