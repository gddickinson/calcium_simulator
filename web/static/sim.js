// Simulator page (#/sim, #/sim?state=<slug>): the full model running in the browser. CalciumModel (ca_engine.js, a
// bit-exact port of calcium_model.py) advances as many 1 ms steps per animation frame as the compute budget allows;
// the field view mirrors the desktop app's tabs and overlays, plus a dF/F0 trace at a probe point.
import { route, app, api, esc, animate, download } from './core.js';
import { CalciumModel } from './ca_engine.js';
import { fluorescence } from './ca_readout.js';
import { drawField, drawLegend, drawTrace, FIELDS } from './sim_view.js';
import { renderPanel } from './sim_panel.js';

const TRACE_MAX = 4000;

route('sim', async (_args, query) => {
  const [specs, states] = await Promise.all([api('/specs'), api('/states')]);
  const first = states.find(s => s.slug === query.state) || states[0];
  const reduced = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  const U = specs.uncaging;
  const st = {
    preset: first.slug, params: { ...first.params }, seed: 0, grid: first.params.grid_size === 100 ? 100 : 200,
    field: 'dff', max: Object.fromEntries(Object.entries(FIELDS).map(([k, f]) => [k, f.max])),
    overlay: { er: false, mito: false, pm: false, ip3r: true }, click: 'probe', budget: 25, running: !reduced,
    uncage: { amount: U.amount.default, duration: U.duration.default, radius: U.radius.default },
    indicator: { kd: 0.35, rf: 100 }, probe: { x: 0, y: 0 },
  };
  let m = null, panel = null, dirty = true, trace = null, f0 = 1, speed = 0;

  app().innerHTML = `
    <h1 class="glow">Simulator</h1>
    <p class="lead" id="note"></p>
    <section class="sim">
      <div class="controls" id="panel" aria-label="Simulation controls"></div>
      <div class="view">
        <div class="fields">
          <canvas id="field" width="600" height="600" tabindex="0" role="img"
            aria-label="Simulated cell: the chosen field as an image, IP3R clusters as circles (red closed, green open), the probe as a white square. Space runs or pauses, R resets."></canvas>
          <div class="legend" id="legend"></div>
          <canvas class="trace" id="trace" role="img" aria-label="dF/F0 and [Ca2+] at the probe point against time"></canvas>
        </div>
        <div class="readout" id="readout" aria-live="off"></div>
      </div>
    </section>`;
  const canvas = document.getElementById('field'), traceCanvas = document.getElementById('trace');

  function probeSample() {
    const n = m.grid_size, { x, y } = st.probe, { kd, rf } = st.indicator;
    let f = 0, c = 0, k = 0;
    for (let i = Math.max(0, y - 1); i <= Math.min(n - 1, y + 1); i++)
      for (let j = Math.max(0, x - 1); j <= Math.min(n - 1, x + 1); j++) { const v = m.calcium[i * n + j]; c += v; f += fluorescence(v, kd, rf); k++; }
    return { dff: f / k / f0 - 1, ca: c / k };
  }
  const resetTrace = () => { trace = { t: [], dff: [], ca: [] }; record(); };
  function record() {
    const s = probeSample();
    trace.t.push(m.time); trace.dff.push(s.dff); trace.ca.push(s.ca);
    if (trace.t.length > TRACE_MAX) for (const a of Object.values(trace)) a.shift();
  }
  const indicatorChanged = () => {
    f0 = fluorescence(m.eq_calcium, st.indicator.kd, st.indicator.rf);
    trace.dff = trace.ca.map(c => fluorescence(c, st.indicator.kd, st.indicator.rf) / f0 - 1);
    dirty = true;
  };

  function build(startup) {
    const params = { ...st.params, grid_size: st.grid, dx: 20 / st.grid };   // a 20 um square cell either way
    m = new CalciumModel(params, st.seed);
    let best = -1;   // probe on the largest IP3R cluster
    m.ip3r_clusters.forEach((v, k) => { if (v > best) { best = v; st.probe = { x: k % st.grid, y: Math.floor(k / st.grid) }; } });
    f0 = fluorescence(m.eq_calcium, st.indicator.kd, st.indicator.rf);
    if (startup) m.addIp3Global(st.uncage.amount, st.uncage.duration);
    resetTrace();
    const nCl = m.ip3r_clusters.reduce((a, v) => a + (v > 0), 0), nCh = m.ip3r_clusters.reduce((a, v) => a + v, 0);
    document.getElementById('note').innerHTML = `${esc(states.find(s => s.slug === st.preset)?.name || 'Custom')} cell, seed ${st.seed}:
      ${nCl} IP3R clusters (${nCh} channels) on the ER, resting [Ca2+] ${(m.eq_calcium * 1000).toFixed(0)} nM.
      ${startup ? `Started with a global IP3 uncaging (${st.uncage.amount} µM over ${st.uncage.duration} s, the desktop defaults).` : 'No IP3 yet: uncage some to open channels.'}`;
    dirty = true;
  }

  const handlers = {
    toggle: () => { st.running = !st.running; panel.setRunning(st.running); },
    reset: () => { build(false); },
    redraw: () => { dirty = true; },
    indicator: indicatorChanged,
    preset: slug => {
      const s = states.find(x => x.slug === slug); if (!s) return;
      st.preset = slug; st.params = { ...s.params };
      panel = renderPanel(document.getElementById('panel'), { specs, states }, st, handlers);
      build(true);
    },
    apply: (params, grid, seed) => { st.params = params; st.grid = grid; st.seed = seed; st.preset = 'custom'; build(false); },
    uncage: at => {
      const { amount, duration, radius } = st.uncage;
      if (!(duration > 0)) return;
      if (at) m.addIp3Local(at.x, at.y, radius, amount, duration); else m.addIp3Global(amount, duration);
      dirty = true;
    },
    png: () => canvas.toBlob(b => download(`calcium_${st.field}_t${m.time.toFixed(3)}s.png`, b)),
    csv: () => download(`calcium_trace_x${st.probe.x}_y${st.probe.y}.csv`,
      'time_s,dF_over_F0,calcium_uM\n' + trace.t.map((t, i) => `${t.toFixed(6)},${trace.dff[i]},${trace.ca[i]}`).join('\n'), 'text/csv'),
  };

  canvas.addEventListener('click', e => {
    const r = canvas.getBoundingClientRect(), n = m.grid_size;
    const x = Math.min(n - 1, Math.max(0, Math.floor((e.clientX - r.left) / r.width * n)));
    const y = Math.min(n - 1, Math.max(0, Math.floor((e.clientY - r.top) / r.height * n)));
    if (st.click === 'uncage') handlers.uncage({ x, y });
    else { st.probe = { x, y }; resetTrace(); }
    dirty = true;
  });
  canvas.addEventListener('keydown', e => {
    if (e.key === ' ') { e.preventDefault(); handlers.toggle(); }
    else if (e.key === 'r' || e.key === 'R') handlers.reset();
  });

  const readout = document.getElementById('readout');
  function show() {
    let open = 0, ch = 0, mean = 0, ip3 = 0;
    for (let k = 0; k < m.calcium.length; k++) { open += m.ip3r_open[k]; ch += m.ip3r_clusters[k]; mean += m.calcium[k]; ip3 += m.ip3_conc[k]; }
    const s = probeSample(), N = m.calcium.length;
    const rows = [['time', `${m.time.toFixed(3)} s`], ['dt · Ca / IP3 sub-steps', `${m.dt} s · ${m.lastSubsteps ?? m.calciumSubsteps()} / ${m.ip3Substeps()}`],
      ['open IP3Rs', `${open} of ${ch}`], ['mean [Ca2+]', `${(mean / N * 1000).toFixed(1)} nM`],
      ['probe (x, y)', `${(st.probe.x * m.dx).toFixed(1)}, ${(st.probe.y * m.dx).toFixed(1)} µm`],
      ['probe [Ca2+] · ΔF/F0', `${(s.ca * 1000).toFixed(1)} nM · ${s.dff.toFixed(3)}`],
      ['mean [IP3]', `${(ip3 / N).toFixed(3)} µM`], ['speed', `${speed.toFixed(0)} ms simulated / s`]];
    readout.innerHTML = rows.map(([k, v]) => `<div><span>${esc(k)}</span><b>${esc(v)}</b></div>`).join('');
  }

  panel = renderPanel(document.getElementById('panel'), { specs, states }, st, handlers);
  build(true);
  if (reduced) { for (let i = 0; i < 30; i++) { m.step(); record(); } }
  let lastShow = 0, simAcc = 0, wallAcc = 0;
  animate(t => {
    if (st.running) {
      const t0 = performance.now(), tm = m.time;
      do { m.step(); record(); } while (performance.now() - t0 < st.budget);
      simAcc += (m.time - tm) * 1000; wallAcc += performance.now() - t0;
      dirty = true;
    }
    if (dirty) { drawField(canvas, m, st); drawLegend(document.getElementById('legend'), st); drawTrace(traceCanvas, trace); dirty = false; }
    if (t - lastShow > 250) {
      if (wallAcc) speed = simAcc / ((t - lastShow) / 1000);
      simAcc = 0; wallAcc = 0; lastShow = t; show();
    }
  });
});
