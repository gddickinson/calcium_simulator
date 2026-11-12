// Model page (#/model): the equations of calcium_model.py, every parameter with unit and default (/api/specs,
// param_specs.py), the numerical scheme and the corrections made in this version.
import { route, app, api, esc } from './core.js';

const unit = u => esc(u).replace(/uM/g, 'µM').replace(/um\^2/g, 'µm²').replace(/\^-1/g, '⁻¹');

route('model', async () => {
  const [S, states] = await Promise.all([api('/specs'), api('/states')]), D = S.defaults;
  app().innerHTML = `
    <h1 class="glow">The model</h1>
    <p class="lead">A 2D square of cytosol (${D.grid_size} × ${D.grid_size} grid points, dx = ${D.dx} µm) containing a random
      reticulated ER, mitochondria and a plasma-membrane edge. Concentrations are in µM, time in s, distance in µm.
      Every grid point has its own cytosolic, ER and mitochondrial [Ca²⁺], bound buffer and [IP3].</p>
    <h2>Equations</h2>
    <div class="eqs"><pre>d[Ca]cyt/dt = D_Ca ∇²[Ca]cyt + J_IP3R + J_leak − J_SERCA − J_PMCA − J_MCU − J_buffer
d[Ca]ER/dt  = −(J_IP3R + J_leak − J_SERCA)            (ER pixels)
d[Ca]mito/dt = J_MCU                                  (mitochondrial pixels)
d[CaB]/dt   = J_buffer                                (immobile buffer)
d[IP3]/dt   = D_IP3 ∇²[IP3] − k_deg [IP3] + uncaging

J_IP3R   = 5 s⁻¹ × n_open × ([Ca]ER − [Ca]cyt)        on ER pixels
J_leak   = leak_rate × ([Ca]ER − [Ca]cyt)             on ER pixels
J_SERCA  = serca_rate × c² / (c² + K²)                on ER pixels
J_PMCA   = pmca_rate × c                              on the plasma-membrane edge
J_MCU    = mcu_rate × (c − [Ca]mito)                  on mitochondrial pixels
J_buffer = k_on × (c × (B_total − [CaB]) − K_d × [CaB])</pre></div>
    <h3>Stochastic IP3 receptors</h3>
    <p>Each IP3R cluster has a Poisson number of channels (mean ${D.ip3r_per_cluster}) and sits on a random
      ${D.ip3r_cluster_density * 100} % of ER pixels. Every time step each cluster opens one more channel with probability
      k_open · dt · n_closed and closes one with probability k_close · dt · n_open, with</p>
    <pre>k_open  = ip3r_open_rate  × c² [IP3]² / ((c + 0.3)³ ([IP3] + 0.2)²)
k_close = ip3r_close_rate × c / (c + 0.3)            (c, [IP3] in µM)</pre>
    <h3>Resting state</h3>
    <p>The cell starts at the [Ca²⁺] where SERCA uptake balances the ER leak at the initial ER [Ca²⁺]
      (serca_rate · c²/(c² + K²) = leak_rate · ([Ca]ER − c)), with the buffer in equilibrium: about 115 nM for the default
      state, 200 nM for the high-calcium state (ER 1000 µM) and 76 nM for the low-calcium state (ER 250 µM).</p>
    <h2>Parameters</h2>
    <div class="scroll"><table><tr><th>Parameter</th><th>Default</th><th>Unit</th><th>Range in the apps</th></tr>
      ${Object.entries(S.params).map(([k, s]) => `<tr><td>${esc(s.label)} <span class="mute mono">${esc(k)}</span></td>
        <td class="mono">${esc(D[k])}</td><td>${unit(s.unit)}</td><td class="mono">${s.min} – ${s.max}</td></tr>`).join('')}
      <tr><td>Time step <span class="mute mono">dt</span></td><td class="mono">${D.dt}</td><td>s</td><td></td></tr>
      <tr><td>Grid spacing <span class="mute mono">dx</span></td><td class="mono">${D.dx}</td><td>µm</td><td></td></tr>
      <tr><td>Initial mitochondrial [Ca²⁺] <span class="mute mono">mito_calcium_init</span></td><td class="mono">${D.mito_calcium_init}</td><td>µM</td><td></td></tr>
    </table></div>
    <h2>Cell states and sources</h2>
    ${states.map(st => `<h3>${esc(st.name)}</h3><p>${esc(st.description)}</p>${st.references.length ? `<ul>${st.references.map(r =>
      `<li>${esc(r.cite)}, <a href="https://doi.org/${esc(r.doi)}" target="_blank" rel="noopener">doi:${esc(r.doi)}</a></li>`).join('')}</ul>` : ''}`).join('')}
    <p>Tested event classes (seed 0, 0.5 µM IP3 over 50 ms, 400 ms): Blips never have more than one channel of a cluster
      open and raise [Ca²⁺] above twice rest over &lt; 0.5 µm²; Puffs open ≥ 3 channels together (7 in the test) over 1–20 µm²
      (about 3 µm²).</p>
    <p><b>No Waves state.</b> The opening rate c²/(c + 0.3)³ peaks at c = 0.6 µM only ~2.6-fold above its resting value,
      so Ca²⁺ from one cluster barely raises the opening rate of the next: with literature cluster sizes, open times,
      diffusion coefficients, ER load and IP3 levels no propagating wave appears. Waves would need a different IP3R gating
      model, which would change the model rather than its parameters.</p>
    <h2>Numerics</h2>
    <p>Diffusion uses the isotropic 9-point Laplacian (1/6dx²)·[[1,4,1],[4,−20,4],[1,4,1]] with zero-flux (mirror)
      boundaries. Each dt is split into n forward-Euler sub-steps with n = ⌈dt × (sum of first-order loss rates)⌉, where
      the diffusion term contributes 20D/(6dx²) and buffer binding k_on(B_free + c + K_d); this keeps the update stable and
      non-negative. At the defaults that is ${'~'}15 sub-steps for Ca²⁺ and 67 for IP3 (D_IP3 = ${D.D_ip3} µm²/s) per 1 ms step;
      IP3 is not sub-stepped while it is uniform. Channel gating happens once per dt.</p>
    <h2>Corrected in this version</h2>
    <ul>
      <li>The old diffusion kernel was 0.3 × the Laplacian (diffusion ran at 30 % of the stated D), and dt = 1 ms was far
        beyond the explicit stability limit: the cell turned into a 0 / 1000 µM checkerboard within 50 ms.</li>
      <li>The resting [Ca²⁺] formula dropped the ER [Ca²⁺] and gave 0.5 nM; it now solves the SERCA = leak balance.</li>
      <li>IP3R rates were used as per-step probabilities; they are now per second (defaults restated as 10 s⁻¹ and
        10000 s⁻¹, which gives exactly the old behaviour at dt = 1 ms).</li>
      <li>IP3 uncaging delivered amount × dt / duration once; it now delivers the amount over the duration.</li>
      <li>Desktop GUI: Apply Settings ignored the initial conditions and set the ER leak to 0 (two-decimal spin box).</li>
    </ul>
    <p class="small mute">Not in the model: indicator dye buffering, optical blur, mobile buffers, Ca²⁺ entry across the
      plasma membrane. The ΔF/F0 trace in the simulator assumes an equilibrium indicator and is computed by the page.</p>`;
});
