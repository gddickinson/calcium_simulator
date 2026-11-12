// Control panel of the simulator page. Model parameters come from /api/specs (param_specs.py, the desktop GUI's
// ranges and units) and the presets from /api/states (configs/*.json, the Cell States menu). "Apply" rebuilds the
// cell like the desktop app's Apply Settings; the display controls only change the picture.
import { esc } from './core.js';
import { FIELDS } from './sim_view.js';

const unitHtml = u => esc(u).replace(/uM/g, 'µM').replace(/um\^2/g, 'µm²').replace(/\^-1/g, '⁻¹');

function numRow(key, s, value) {
  return `<div class="num"><label for="p-${key}">${esc(s.label)} <span class="mute">(${unitHtml(s.unit)})</span></label>
    <input type="number" id="p-${key}" data-key="${key}" min="${s.min}" max="${s.max}" step="any" value="${value}"></div>`;
}

export function renderPanel(host, data, st, h) {
  const { specs, states } = data, groups = {};
  for (const [key, s] of Object.entries(specs.params)) (groups[s.group] ||= []).push([key, s]);
  const u = specs.uncaging;
  host.innerHTML = `
    <div class="btns">
      <button id="run">${st.running ? 'Pause' : 'Run'}</button><button class="alt" id="reset">Reset</button>
      <button class="alt" id="png">PNG</button><button class="alt" id="csv">Trace CSV</button>
    </div>
    <label for="preset">Cell state (configs/)</label>
    <select id="preset">${states.map(s => `<option value="${esc(s.slug)}" ${s.slug === st.preset ? 'selected' : ''}>${esc(s.name)} (ER ${s.params.er_calcium_init} µM)</option>`).join('')}</select>
    <details open><summary>View</summary>
      <label for="field-sel">Field</label>
      <select id="field-sel">${Object.entries(FIELDS).map(([k, f]) => `<option value="${k}" ${k === st.field ? 'selected' : ''}>${esc(f.label)}</option>`).join('')}</select>
      <div class="num"><label for="vmax">Display max <span class="mute" id="vmax-unit"></span></label><input type="number" id="vmax" min="0" step="any" value="${st.max[st.field]}"></div>
      <div class="chk" role="group" aria-label="Overlays">
        ${[['er', 'ER'], ['mito', 'Mitochondria'], ['pm', 'Plasma membrane'], ['ip3r', 'IP3R clusters']].map(([k, l]) =>
          `<label><input type="checkbox" data-ov="${k}" ${st.overlay[k] ? 'checked' : ''}>${l}</label>`).join('')}
      </div>
      <label for="click-sel">Clicking the cell</label>
      <select id="click-sel"><option value="probe">moves the dF/F0 probe</option><option value="uncage">uncages IP3 there (local)</option></select>
      <div class="num"><label for="budget">Compute per frame (ms)</label><input type="number" id="budget" min="1" max="200" value="${st.budget}"></div>
    </details>
    <details open><summary>IP3 uncaging</summary>
      ${['amount', 'duration', 'radius'].map(k => numRow('u-' + k, u[k], st.uncage[k])).join('')}
      <button id="uncage-global">Global</button><button class="alt" id="uncage-probe">At probe</button>
    </details>
    <details><summary>Indicator (dF/F0 readout)</summary>
      <p class="small mute">Equilibrium single-wavelength dye: F ∝ 1 + (Rf − 1)·c/(c + Kd). A display assumption of this page, not part of the Python model; set it to your dye.</p>
      <div class="num"><label for="kd">Indicator Kd (µM)</label><input type="number" id="kd" min="0.001" step="any" value="${st.indicator.kd}"></div>
      <div class="num"><label for="rf">Rf = Fmax/Fmin</label><input type="number" id="rf" min="1" step="any" value="${st.indicator.rf}"></div>
    </details>
    <details><summary>Grid, time step and seed</summary>
      <label for="grid-sel">Grid</label>
      <select id="grid-sel">
        <option value="200" ${st.grid === 200 ? 'selected' : ''}>200 × 200, dx 0.1 µm (as the desktop app)</option>
        <option value="100" ${st.grid === 100 ? 'selected' : ''}>100 × 100, dx 0.2 µm (faster; structures are drawn in grid points, so twice as large)</option>
      </select>
      <div class="num"><label for="p-dt">Time step dt <span class="mute">(s)</span></label><input type="number" id="p-dt" data-key="dt" min="0.00001" max="0.01" step="any" value="${st.params.dt}"></div>
      <div class="num"><label for="seed">Random seed (np.random.seed)</label><input type="number" id="seed" min="0" step="1" value="${st.seed}"></div>
    </details>
    ${Object.entries(groups).map(([g, rows]) => `<details ${g === 'IP3R' ? 'open' : ''}><summary>${esc(g)}</summary>
      ${rows.map(([key, s]) => numRow(key, s, st.params[key])).join('')}</details>`).join('')}
    <button id="apply">Apply (rebuild cell)</button>
    <p class="small mute">Space: run/pause · R: reset (canvas focused).</p>`;

  const q = sel => host.querySelector(sel);
  const setUnit = () => { q('#vmax-unit').textContent = FIELDS[st.field].unit ? `(${FIELDS[st.field].unit.replace('uM', 'µM')})` : ''; };
  setUnit();
  q('#run').onclick = h.toggle; q('#reset').onclick = h.reset; q('#png').onclick = h.png; q('#csv').onclick = h.csv;
  q('#preset').onchange = e => h.preset(e.target.value);
  q('#field-sel').onchange = e => { st.field = e.target.value; q('#vmax').value = st.max[st.field]; setUnit(); h.redraw(); };
  q('#vmax').onchange = e => { const v = +e.target.value; if (v > 0) { st.max[st.field] = v; h.redraw(); } };
  host.querySelectorAll('[data-ov]').forEach(c => { c.onchange = () => { st.overlay[c.dataset.ov] = c.checked; h.redraw(); }; });
  q('#click-sel').onchange = e => { st.click = e.target.value; };
  q('#budget').onchange = e => { st.budget = Math.min(200, Math.max(1, +e.target.value || 25)); };
  for (const k of ['amount', 'duration', 'radius']) q('#p-u-' + k).onchange = e => { st.uncage[k] = +e.target.value; };
  q('#uncage-global').onclick = () => h.uncage(null);
  q('#uncage-probe').onclick = () => h.uncage(st.probe);
  q('#kd').onchange = e => { st.indicator.kd = Math.max(0.001, +e.target.value); h.indicator(); };
  q('#rf').onchange = e => { st.indicator.rf = Math.max(1, +e.target.value); h.indicator(); };
  q('#apply').onclick = () => {
    const params = { ...st.params };
    host.querySelectorAll('input[data-key]').forEach(i => { if (i.value !== '' && Number.isFinite(+i.value)) params[i.dataset.key] = +i.value; });
    params.ip3r_per_cluster = Math.round(params.ip3r_per_cluster);
    h.apply(params, +q('#grid-sel').value, Math.max(0, Math.round(+q('#seed').value || 0)));
  };
  return { setRunning: r => { q('#run').textContent = r ? 'Pause' : 'Run'; } };
}
