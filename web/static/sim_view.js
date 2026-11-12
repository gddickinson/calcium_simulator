// Drawing for the simulator page: the field image with the desktop app's overlays (ER green, mitochondria blue,
// plasma membrane yellow, IP3R clusters red when closed / green when open), the probe marker, the colour legend and
// the dF/F0 + [Ca2+] trace at the probe point.
import { lut, fluorescence } from './ca_readout.js';

// field key -> label, unit, colour map, default display max (the desktop GUI levels for ER and mitochondria)
export const FIELDS = {
  dff: { label: 'ΔF/F0 (indicator)', unit: '', map: 'fire', max: 1 },
  calcium: { label: 'Cytosolic [Ca2+]', unit: 'uM', map: 'fire', max: 0.5 },
  er: { label: 'ER [Ca2+]', unit: 'uM', map: 'green', max: 1000 },
  mito: { label: 'Mitochondrial [Ca2+]', unit: 'uM', map: 'violet', max: 10 },
  ip3: { label: '[IP3]', unit: 'uM', map: 'ice', max: 2 },
};
const LUTS = Object.fromEntries(['fire', 'ice', 'green', 'violet'].map(n => [n, lut(n)]));

// Values of the chosen field for every grid point
function values(m, key, ind) {
  const N = m.calcium.length, out = new Float64Array(N);
  if (key === 'calcium') out.set(m.calcium);
  else if (key === 'ip3') out.set(m.ip3_conc);
  else if (key === 'er') for (let k = 0; k < N; k++) out[k] = m.er_calcium[k] * m.er[k];
  else if (key === 'mito') for (let k = 0; k < N; k++) out[k] = m.mito_calcium[k] * m.mitochondria[k];
  else {
    const f0 = fluorescence(m.eq_calcium, ind.kd, ind.rf);
    for (let k = 0; k < N; k++) out[k] = fluorescence(m.calcium[k], ind.kd, ind.rf) / f0 - 1;
  }
  return out;
}

const off = document.createElement('canvas');

export function drawField(canvas, m, st) {
  const n = m.grid_size, f = FIELDS[st.field], L = LUTS[f.map], v = values(m, st.field, st.indicator);
  off.width = off.height = n;
  const octx = off.getContext('2d'), img = octx.createImageData(n, n), d = img.data, max = st.max[st.field];
  const blend = (k, r, g, b) => {   // alpha 100/255, as the desktop overlay
    const a = 100 / 255;
    d[k * 4] = d[k * 4] * (1 - a) + r * a; d[k * 4 + 1] = d[k * 4 + 1] * (1 - a) + g * a; d[k * 4 + 2] = d[k * 4 + 2] * (1 - a) + b * a;
  };
  for (let k = 0; k < n * n; k++) {
    const i = Math.max(0, Math.min(255, Math.round(v[k] / max * 255)));
    d[k * 4] = L[i * 3]; d[k * 4 + 1] = L[i * 3 + 1]; d[k * 4 + 2] = L[i * 3 + 2]; d[k * 4 + 3] = 255;
    if (st.overlay.er && m.er[k]) blend(k, 0, 255, 0);
    if (st.overlay.mito && m.mitochondria[k]) blend(k, 0, 0, 255);
    if (st.overlay.pm && m.pm[k]) blend(k, 255, 255, 0);
  }
  octx.putImageData(img, 0, 0);
  const ctx = canvas.getContext('2d'), s = canvas.width / n;
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(off, 0, 0, canvas.width, canvas.height);
  if (st.overlay.ip3r) {
    for (let k = 0; k < n * n; k++) {
      if (!m.ip3r_clusters[k]) continue;
      const open = m.ip3r_open[k] > 0;
      ctx.strokeStyle = open ? '#3dff6a' : 'rgba(255,60,60,.85)'; ctx.lineWidth = open ? 2 : 1;
      ctx.beginPath(); ctx.arc(((k % n) + 0.5) * s, (Math.floor(k / n) + 0.5) * s, open ? 3 * s : 1.6 * s, 0, 7); ctx.stroke();
    }
  }
  const px = (st.probe.x + 0.5) * s, py = (st.probe.y + 0.5) * s, r = 1.5 * s;
  ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5;
  ctx.strokeRect(px - r, py - r, 2 * r, 2 * r);
  ctx.beginPath(); ctx.moveTo(px - 3 * r, py); ctx.lineTo(px - r, py); ctx.moveTo(px + r, py); ctx.lineTo(px + 3 * r, py);
  ctx.moveTo(px, py - 3 * r); ctx.lineTo(px, py - r); ctx.moveTo(px, py + r); ctx.lineTo(px, py + 3 * r); ctx.stroke();
  const bar = 5 / (n * m.dx) * canvas.width;   // 5 um scale bar
  ctx.fillStyle = 'rgba(0,0,0,.55)'; ctx.fillRect(6, canvas.height - 26, bar + 56, 20);
  ctx.fillStyle = '#fff'; ctx.font = '12px ui-monospace,Menlo,monospace';
  ctx.fillRect(12, canvas.height - 14, bar, 3); ctx.fillText('5 µm', 18 + bar, canvas.height - 10);
}

export function drawLegend(host, st) {
  const f = FIELDS[st.field], L = LUTS[f.map];
  host.innerHTML = `<span>0</span><canvas width="256" height="1" aria-hidden="true"></canvas><span>${st.max[st.field]} ${f.unit.replace('uM', 'µM')}</span>`;
  const c = host.querySelector('canvas').getContext('2d'), img = c.createImageData(256, 1);
  for (let i = 0; i < 256; i++) { img.data.set([L[i * 3], L[i * 3 + 1], L[i * 3 + 2], 255], i * 4); }
  c.putImageData(img, 0, 0);
}

// Trace: top half dF/F0, bottom half [Ca2+] (uM), against time (s)
export function drawTrace(canvas, trace) {
  const dpr = window.devicePixelRatio || 1, W = canvas.clientWidth, H = canvas.clientHeight;
  if (canvas.width !== Math.round(W * dpr)) { canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr); }
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
  ctx.font = '11px ui-monospace,Menlo,monospace';
  const left = 52, right = 10, half = (H - 26) / 2;
  const t = trace.t, t0 = t.length ? t[0] : 0, t1 = t.length ? Math.max(t.at(-1), t0 + 1e-3) : 1;
  const X = s => left + (s - t0) / (t1 - t0) * (W - left - right);
  const panel = (vals, top, color, name, floor) => {
    let lo = floor, hi = floor + 1e-3;
    for (const v of vals) { if (v > hi) hi = v; if (v < lo) lo = v; }
    hi += (hi - lo) * 0.08;
    const Y = v => top + half - 6 - (v - lo) / (hi - lo) * (half - 12);
    ctx.strokeStyle = '#212c4c'; ctx.strokeRect(left, top, W - left - right, half);
    ctx.fillStyle = '#94a3c6';
    const fmt = v => (Math.abs(v) < 1e-3 * (hi - lo) ? '0' : v.toPrecision(3));
    ctx.fillText(fmt(hi), 4, top + 10); ctx.fillText(fmt(lo), 4, top + half - 2);
    ctx.fillStyle = color; ctx.fillText(name, left + 6, top + 13);
    ctx.strokeStyle = color; ctx.lineWidth = 1.4; ctx.beginPath();
    vals.forEach((v, i) => { const x = X(t[i]), y = Y(v); if (i) ctx.lineTo(x, y); else ctx.moveTo(x, y); });
    ctx.stroke();
  };
  panel(trace.dff, 4, '#ffb020', 'ΔF/F0 at probe', 0);
  panel(trace.ca, 8 + half, '#2fd3c0', '[Ca2+] at probe (µM)', trace.ca.length ? Math.min(...trace.ca) : 0);
  ctx.fillStyle = '#94a3c6';
  ctx.fillText(`${t0.toFixed(3)} s`, left, H - 4);
  const lab = `${t1.toFixed(3)} s`; ctx.fillText(lab, W - right - ctx.measureText(lab).width, H - 4);
}
