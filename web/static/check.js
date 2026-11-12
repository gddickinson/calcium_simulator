// Cross-check page (#/check): repeats the seeded Python run served by /api/check (web/reference.py, produced by
// calcium_model.py) with the JavaScript engine in this browser and compares the results.
import { route, app, api, esc } from './core.js';
import { CalciumModel } from './ca_engine.js';

route('check', async () => {
  const R = await api('/check');
  const c = R.case;
  app().innerHTML = `<h1 class="glow">JS = Python</h1>
    <p class="lead">The Python model was run once when this site was built: <code>np.random.seed(${c.seed})</code>, a
      ${c.params.grid_size} × ${c.params.grid_size} cell, local and global IP3 uncaging, ${c.steps} steps of 1 ms. Your browser now
      runs the same case with <code>ca_engine.js</code>.</p><div id="out"><p class="mute">Running...</p></div>`;
  const t0 = performance.now();
  const m = new CalciumModel(c.params, c.seed);
  for (const call of c.uncage) { if (call[0] === 'local') m.addIp3Local(...call.slice(1)); else m.addIp3Global(...call.slice(1)); }
  const clusters = [...m.ip3r_clusters].flatMap((v, k) => (v ? [[k, v]] : []));
  const opens = [];
  for (let s = 0; s < c.steps; s++) { m.step(); opens.push(m.ip3r_open.reduce((a, v) => a + v, 0)); }
  const ms = performance.now() - t0;
  const rel = (js, py) => { const sc = Math.max(...py.map(Math.abs)) || 1; let w = 0; py.forEach((v, k) => { w = Math.max(w, Math.abs(js[k] - v) / sc); }); return w; };
  const rows = [
    ['IP3R clusters (position, channels)', JSON.stringify(clusters) === JSON.stringify(R.clusters), `${clusters.length} clusters`],
    ['open channels after every step', JSON.stringify(opens) === JSON.stringify(R.open_trace), `${opens.reduce((a, v) => a + v, 0)} channel-steps open`],
    ['resting [Ca²⁺]', m.eq_calcium === R.eq_calcium, `${m.eq_calcium} µM`],
    ['cytosolic [Ca²⁺] field', rel(m.calcium, R.calcium) < 1e-12, `max relative difference ${rel(m.calcium, R.calcium).toExponential(1)}`],
    ['IP3 field', rel(m.ip3_conc, R.ip3) < 1e-12, `max relative difference ${rel(m.ip3_conc, R.ip3).toExponential(1)}`],
  ];
  const all = rows.every(r => r[1]);
  document.getElementById('out').innerHTML = `<p class="${all ? 'ok-big' : 'bad-big'}">${all ? 'All match' : 'Mismatch'}</p>
    <div class="scroll"><table><tr><th>Quantity</th><th>Match</th><th>Detail</th></tr>
    ${rows.map(([q, okv, d]) => `<tr><td>${esc(q)}</td><td class="${okv ? 'ok' : 'err'}">${okv ? 'yes' : 'no'}</td><td class="mono">${esc(d)}</td></tr>`).join('')}
    </table></div><p class="small mute">JS run took ${ms.toFixed(0)} ms. The repository's node test (web/static/ca_engine.test.mjs) checks two such
    cases and the default 200 × 200 cell on every build.</p>`;
});
