// Canvas scatter-plot viewer for lattice frames: a rotatable 3D view and
// 2D slices through one lattice layer. No dependencies.
//
// Data (all in ångström): positions is a Float32Array of frames x particles x 3.

const MIN_ZOOM = 0.2;
const MAX_ZOOM = 40;

export function niceTicks(lo, hi, target = 5) {
  const span = hi - lo;
  if (!(span > 0) || !Number.isFinite(span)) return { ticks: [lo], step: 1 };
  const raw = span / Math.max(1, target);
  const mag = 10 ** Math.floor(Math.log10(raw));
  const norm = raw / mag;
  const step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
  const ticks = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-6; v += step) {
    ticks.push(Math.abs(v) < step * 1e-6 ? 0 : v);
  }
  return { ticks, step };
}

function formatTick(value, step) {
  const decimals = Math.max(0, Math.min(6, -Math.floor(Math.log10(step) + 1e-9)));
  return value.toFixed(decimals);
}

export class ScatterViewer {
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.data = null;
    this.frame = 0;
    this.mode = '3d';
    this.layer = 0;
    this.autoFit = false;
    this.pending = false;
    this.drag = null;
    this.resetView(false);

    new ResizeObserver(() => this.requestDraw()).observe(canvas);
    matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => this.requestDraw());

    canvas.addEventListener('pointerdown', (e) => {
      if (e.button !== 0) return;
      canvas.setPointerCapture(e.pointerId);
      this.drag = { x: e.clientX, y: e.clientY };
      canvas.classList.add('dragging');
    });
    canvas.addEventListener('pointermove', (e) => {
      if (!this.drag) return;
      const dx = e.clientX - this.drag.x;
      const dy = e.clientY - this.drag.y;
      this.drag = { x: e.clientX, y: e.clientY };
      if (this.mode === '3d') {
        this.yaw += dx * 0.01;
        this.pitch = Math.max(-1.55, Math.min(1.55, this.pitch + dy * 0.01));
      } else {
        this.pan[0] += dx;
        this.pan[1] += dy;
      }
      this.requestDraw();
    });
    const endDrag = () => {
      this.drag = null;
      canvas.classList.remove('dragging');
    };
    canvas.addEventListener('pointerup', endDrag);
    canvas.addEventListener('pointercancel', endDrag);
    canvas.addEventListener('wheel', (e) => {
      if (!this.data) return;
      e.preventDefault();
      const factor = Math.exp(-e.deltaY * (e.deltaMode === 1 ? 0.05 : 0.0015));
      this.zoom = Math.max(MIN_ZOOM, Math.min(MAX_ZOOM, this.zoom * factor));
      this.requestDraw();
    }, { passive: false });
    canvas.addEventListener('dblclick', () => this.resetView());
  }

  resetView(redraw = true) {
    this.yaw = -0.62;
    this.pitch = 0.42;
    this.zoom = 1;
    this.pan = [0, 0];
    if (redraw) this.requestDraw();
  }

  setData(data) {
    this.data = data;
    this.frame = 0;
    this.layer = 0;
    this.resetView();
  }

  clear() {
    this.data = null;
    this.requestDraw();
  }

  setFrame(index) {
    this.frame = index;
    this.requestDraw();
  }

  setMode(mode) {
    this.mode = mode;
    this.zoom = 1;
    this.pan = [0, 0];
    this.requestDraw();
  }

  setLayer(layer) {
    this.layer = layer;
    this.requestDraw();
  }

  setAutoFit(on) {
    this.autoFit = on;
    this.requestDraw();
  }

  // Bounds of the data (all frames), or of the current frame when auto-fitting.
  bounds(onlyLayer = false) {
    const { positions, particles, layer, lo, hi } = this.data;
    if (!this.autoFit) return { lo, hi };
    const base = this.frame * particles * 3;
    const flo = [Infinity, Infinity, Infinity];
    const fhi = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < particles; i++) {
      if (onlyLayer && layer[i] !== this.layer) continue;
      for (let k = 0; k < 3; k++) {
        const v = positions[base + i * 3 + k];
        if (v < flo[k]) flo[k] = v;
        if (v > fhi[k]) fhi[k] = v;
      }
    }
    if (!Number.isFinite(flo[0])) return { lo, hi };
    const span = Math.max(fhi[0] - flo[0], fhi[1] - flo[1], fhi[2] - flo[2]);
    const pad = 0.04 * span + 1;
    return { lo: flo.map((v) => v - pad), hi: fhi.map((v) => v + pad) };
  }

  requestDraw() {
    if (this.pending) return;
    this.pending = true;
    requestAnimationFrame(() => {
      this.pending = false;
      this.draw();
    });
  }

  colors() {
    const style = getComputedStyle(this.canvas);
    const get = (name) => style.getPropertyValue(name).trim();
    return {
      bg: get('--plot-bg'), axis: get('--axis'), grid: get('--grid'), text: get('--text'),
      muted: get('--muted'), deut: get('--deut'), li: get('--li'),
      font: `12px ${style.fontFamily}`, titleFont: `600 14px ${style.fontFamily}`,
    };
  }

  draw() {
    const { canvas, ctx } = this;
    const dpr = window.devicePixelRatio || 1;
    const w = canvas.clientWidth;
    const h = canvas.clientHeight;
    if (!w || !h) return;
    const pw = Math.round(w * dpr);
    const ph = Math.round(h * dpr);
    if (canvas.width !== pw || canvas.height !== ph) {
      canvas.width = pw;
      canvas.height = ph;
    }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const c = this.colors();
    ctx.fillStyle = c.bg;
    ctx.fillRect(0, 0, w, h);
    if (!this.data) return;
    if (this.mode === '3d') this.draw3d(w, h, c);
    else this.draw2d(w, h, c);
  }

  legend(w, c, counts, marker) {
    const { ctx } = this;
    ctx.font = c.font;
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    const items = [['Deuterium', c.deut, counts[0]], ['Lithium', c.li, counts[1]]];
    const width = Math.max(...items.map(([label, , n]) => ctx.measureText(`${label} (${n})`).width)) + 30;
    items.forEach(([label, color, n], i) => {
      const y = 18 + i * 18;
      const x = w - width - 8;
      ctx.fillStyle = color;
      marker(x + 6, y, 4.5);
      ctx.fillStyle = c.text;
      ctx.fillText(`${label} (${n})`, x + 18, y);
    });
  }

  title(c, text) {
    const { ctx } = this;
    ctx.font = c.titleFont;
    ctx.fillStyle = c.text;
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.fillText(text, 12, 18);
  }

  draw3d(w, h, c) {
    const { ctx, data } = this;
    const { positions, particles, isLi } = data;
    const { lo, hi } = this.bounds();
    const centre = [0, 1, 2].map((k) => (lo[k] + hi[k]) / 2);
    const radius = 0.5 * Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1;
    const scale = (Math.min(w, h) * 0.4 / radius) * this.zoom;
    const ox = w / 2 + this.pan[0];
    const oy = h / 2 + this.pan[1] + 8;
    const cy = Math.cos(this.yaw), sy = Math.sin(this.yaw);
    const cp = Math.cos(this.pitch), sp = Math.sin(this.pitch);
    // Orthographic camera: rotate about z (yaw), then tilt (pitch). Right-handed, z up.
    const project = (x, y, z) => {
      const dx = x - centre[0], dy = y - centre[1], dz = z - centre[2];
      const rx = dx * cy - dy * sy;
      const ry = dx * sy + dy * cy;
      return [ox + rx * scale, oy - (ry * sp + dz * cp) * scale, ry * cp - dz * sp];
    };

    // Bounding box.
    const corner = (i) => [i & 1 ? hi[0] : lo[0], i & 2 ? hi[1] : lo[1], i & 4 ? hi[2] : lo[2]];
    const pc = Array.from({ length: 8 }, (_, i) => project(...corner(i)));
    const edges = [];
    for (let i = 0; i < 8; i++) for (const bit of [1, 2, 4]) if (!(i & bit)) edges.push([i, i | bit, bit]);
    ctx.lineWidth = 1;
    ctx.strokeStyle = c.grid;
    ctx.beginPath();
    for (const [a, b] of edges) {
      ctx.moveTo(pc[a][0], pc[a][1]);
      ctx.lineTo(pc[b][0], pc[b][1]);
    }
    ctx.stroke();

    // Axis ticks on the front-bottom edges (x, y) and the leftmost vertical edge (z).
    const centreScreen = project(...centre);
    const axisEdge = (bit) => {
      const candidates = edges.filter(([, , b]) => b === bit);
      const score = ([a, b]) => bit === 4
        ? Math.min(pc[a][0], pc[b][0]) * 1000 + (pc[a][2] + pc[b][2])
        : -(pc[a][1] + pc[b][1]) * 1000 + (pc[a][2] + pc[b][2]);
      return candidates.reduce((best, e) => (score(e) < score(best) ? e : best));
    };
    ctx.font = c.font;
    ctx.fillStyle = c.muted;
    ctx.strokeStyle = c.axis;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ['x', 'y', 'z'].forEach((name, k) => {
      const bit = 1 << k;
      const [a] = axisEdge(bit);
      const start = corner(a);
      const mid = [...start];
      mid[k] = centre[k];
      const pm = project(...mid);
      let ux = pm[0] - centreScreen[0], uy = pm[1] - centreScreen[1];
      const len = Math.hypot(ux, uy) || 1;
      ux /= len; uy /= len;
      const { ticks, step } = niceTicks(lo[k], hi[k], 4);
      ctx.beginPath();
      for (const t of ticks) {
        const p = [...start];
        p[k] = t;
        const [px, py] = project(...p);
        ctx.moveTo(px, py);
        ctx.lineTo(px + ux * 5, py + uy * 5);
        ctx.fillText(formatTick(t, step), px + ux * 16, py + uy * 16);
      }
      ctx.stroke();
      ctx.fillStyle = c.text;
      ctx.fillText(`${name} (Å)`, pm[0] + ux * 38, pm[1] + uy * 38);
      ctx.fillStyle = c.muted;
    });

    // Particles, far to near, fading with depth.
    const base = this.frame * particles * 3;
    const sx = new Float32Array(particles);
    const syArr = new Float32Array(particles);
    const depth = new Float32Array(particles);
    let dmin = Infinity, dmax = -Infinity;
    for (let i = 0; i < particles; i++) {
      const j = base + i * 3;
      const p = project(positions[j], positions[j + 1], positions[j + 2]);
      sx[i] = p[0]; syArr[i] = p[1]; depth[i] = p[2];
      if (p[2] < dmin) dmin = p[2];
      if (p[2] > dmax) dmax = p[2];
    }
    const order = Array.from({ length: particles }, (_, i) => i).sort((a, b) => depth[b] - depth[a]);
    const spacing = 2.04;
    const r = Math.max(particles > 5000 ? 1.2 : 1.8, Math.min(9, spacing * scale * 0.16));
    const range = dmax - dmin || 1;
    let nLi = 0;
    for (const i of order) {
      if (isLi[i]) nLi++;
      ctx.globalAlpha = 1 - 0.5 * ((depth[i] - dmin) / range);
      ctx.fillStyle = isLi[i] ? c.li : c.deut;
      ctx.beginPath();
      ctx.arc(sx[i], syArr[i], r, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.globalAlpha = 1;
    this.title(c, data.name);
    this.legend(w, c, [particles - nLi, nLi], (x, y, rr) => {
      ctx.beginPath();
      ctx.arc(x, y, rr, 0, Math.PI * 2);
      ctx.fill();
    });
  }

  draw2d(w, h, c) {
    const { ctx, data } = this;
    const { positions, particles, isLi, layer } = data;
    const { lo, hi } = this.bounds(true);
    const m = { left: 64, right: 18, top: 40, bottom: 48 };
    const pw = Math.max(10, w - m.left - m.right);
    const ph = Math.max(10, h - m.top - m.bottom);
    const xr = hi[0] - lo[0] || 1;
    const yr = hi[1] - lo[1] || 1;
    const scale = Math.min(pw / xr, ph / yr) * this.zoom;
    const cx = m.left + pw / 2 + this.pan[0];
    const cy = m.top + ph / 2 + this.pan[1];
    const mx = (lo[0] + hi[0]) / 2;
    const my = (lo[1] + hi[1]) / 2;
    const X = (x) => cx + (x - mx) * scale;
    const Y = (y) => cy - (y - my) * scale;
    const xmin = mx + (m.left - cx) / scale, xmax = mx + (m.left + pw - cx) / scale;
    const ymax = my + (cy - m.top) / scale, ymin = my + (cy - m.top - ph) / scale;

    // Grid, ticks and frame.
    const xt = niceTicks(xmin, xmax, Math.max(3, pw / 90));
    const yt = niceTicks(ymin, ymax, Math.max(3, ph / 70));
    ctx.lineWidth = 1;
    ctx.strokeStyle = c.grid;
    ctx.beginPath();
    for (const t of xt.ticks) { ctx.moveTo(X(t), m.top); ctx.lineTo(X(t), m.top + ph); }
    for (const t of yt.ticks) { ctx.moveTo(m.left, Y(t)); ctx.lineTo(m.left + pw, Y(t)); }
    ctx.stroke();
    ctx.strokeStyle = c.axis;
    ctx.strokeRect(m.left, m.top, pw, ph);
    ctx.font = c.font;
    ctx.fillStyle = c.muted;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (const t of xt.ticks) ctx.fillText(formatTick(t, xt.step), X(t), m.top + ph + 6);
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const t of yt.ticks) ctx.fillText(formatTick(t, yt.step), m.left - 6, Y(t));
    ctx.fillStyle = c.text;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'bottom';
    ctx.fillText('x (Å)', m.left + pw / 2, h - 6);
    ctx.save();
    ctx.translate(16, m.top + ph / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textBaseline = 'middle';
    ctx.fillText('y (Å)', 0, 0);
    ctx.restore();

    // Particles in the chosen layer, as diamonds.
    const diamond = (x, y, s) => {
      ctx.beginPath();
      ctx.moveTo(x, y - s);
      ctx.lineTo(x + s, y);
      ctx.lineTo(x, y + s);
      ctx.lineTo(x - s, y);
      ctx.closePath();
      ctx.fill();
    };
    const s = Math.max(2.5, Math.min(10, 2.04 * scale * 0.14));
    ctx.save();
    ctx.beginPath();
    ctx.rect(m.left, m.top, pw, ph);
    ctx.clip();
    const base = this.frame * particles * 3;
    let nD = 0, nLi = 0;
    for (let i = 0; i < particles; i++) {
      if (layer[i] !== this.layer) continue;
      const j = base + i * 3;
      if (isLi[i]) nLi++; else nD++;
      ctx.fillStyle = isLi[i] ? c.li : c.deut;
      diamond(X(positions[j]), Y(positions[j + 1]), s);
    }
    ctx.restore();
    this.title(c, `${data.name}: layer z = ${this.layer}`);
    this.legend(w, c, [nD, nLi], (x, y, rr) => diamond(x, y, rr + 0.5));
  }
}
