// Home page (#/): hero, counts taken from the Python model (/api/info, /api/states) and links to the pages.
import { route, app, api, esc } from './core.js';

route('home', async () => {
  const [info, states] = await Promise.all([api('/info'), api('/states')]);
  const tiles = [[`${info.grid} × ${info.grid}`, 'grid points'], [`${info.cell_um} µm`, `square cell, dx ${info.dx_um} µm`],
    [`${info.dt_s * 1000} ms`, 'time step dt'], [`${(info.resting_calcium_uM * 1000).toFixed(0)} nM`, 'resting [Ca²⁺] (SERCA = leak)'],
    [info.parameters, 'model parameters'], [info.cell_states, 'cell states']];
  app().innerHTML = `
    <section class="hero">
      <p class="kicker">Intracellular Ca²⁺ signalling</p>
      <h1 class="glow">Calcium Simulator</h1>
      <p class="lead">${esc(info.description)} The simulator on this site is the repository's Python model
        (<code>calcium_model.py</code>) ported to JavaScript, reproducing a seeded Python run bit for bit.</p>
      <a class="btn" href="#/sim">Open the simulator</a> <a class="btn alt" href="#/model">How the model works</a>
    </section>
    <div class="tiles">${tiles.map(([b, s]) => `<div class="tile"><b>${esc(b)}</b><span>${esc(s)}</span></div>`).join('')}</div>
    <h2>Explore</h2>
    <div class="grid">
      <a class="card" href="#/sim"><b>Simulator</b><small>runs in your browser</small><span>Cytosolic, ER and mitochondrial
        [Ca²⁺] and IP3 as live images with the desktop app's overlays, IP3 uncaging (global or at a point) and a ΔF/F0 trace at a
        probe you place.</span></a>
      <a class="card" href="#/model"><b>Model</b><small>equations, parameters, units</small><span>The reaction-diffusion
        equations, every parameter with its unit and default, and how the explicit scheme is kept stable.</span></a>
      <a class="card" href="#/check"><b>JS = Python</b><small>live cross-check</small><span>Repeats a seeded Python run in
        your browser and compares every IP3R gating decision and concentration.</span></a>
      <a class="card" href="#/about"><b>About</b><small>run it yourself</small><span>Desktop app (PyQt5), local web app and
        credits.</span></a>
    </div>
    <h2>Cell states</h2>
    <div class="grid">${states.map(s => `<a class="card" href="#/sim?state=${esc(s.slug)}"><b>${esc(s.name)}</b>
      <small>${esc(s.file)}</small><span>Initial ER [Ca²⁺] ${s.params.er_calcium_init} µM; all other parameters at their
      defaults.</span></a>`).join('')}</div>`;
});
