// About page (#/about): what the project is, how to run the desktop and web apps, credits (README, CITATION.cff).
import { route, app } from './core.js';

route('about', () => {
  app().innerHTML = `
    <h1 class="glow">About</h1>
    <p class="lead">Calcium Simulator is a tool for modelling and visualising intracellular Ca²⁺ signals: an interactive
      2D simulator of cytosolic, ER and mitochondrial Ca²⁺ with stochastic IP3 receptor clusters, written in Python with a
      PyQt5 / pyqtgraph desktop interface.</p>
    <h2>Run the desktop app</h2>
    <pre>git clone https://github.com/gddickinson/calcium_simulator
cd calcium_simulator
pip install -r requirements.txt
python main.py</pre>
    <p>Four tabs (cytoplasm, ER, mitochondria, IP3), parameter docks, IP3 uncaging and the Cell States menu.</p>
    <h2>Run the web app</h2>
    <pre>pip install -r requirements-web.txt
python run_web.py            # http://127.0.0.1:8205</pre>
    <p>The local web app serves this same front end; the parameter specs, cell states and the reference run come from the
      Python code through its API. <code>python scripts/build_static.py</code> exports the static site published on GitHub
      Pages.</p>
    <h2>How the web simulator relates to the Python</h2>
    <p><code>ca_engine.js</code> is a line-by-line port of <code>CalciumModel</code>, including numpy's random number
      generator, so a run started with the same seed gives the same cell, the same channel openings and the same
      concentrations as Python (checked by a node test in the repository and live on the <a href="#/check">JS = Python</a>
      page). The ΔF/F0 readout and the colour maps are additions of the web page.</p>
    <h2>Credits</h2>
    <p>George D. Dickinson, University of California, Irvine (<a href="https://orcid.org/0000-0001-7711-0388">ORCID</a>).
      MIT licence. Built with AI assistance from Claude (Anthropic).</p>`;
});
