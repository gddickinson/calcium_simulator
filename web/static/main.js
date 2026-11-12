// Entry point: registers the lazily loaded pages, sets the footer mode text, starts the router and the background.
import { start, route, lazy, app, isStatic, $ } from './core.js';
import { startBackground } from './bg.js';

const PAGES = { 'home.js': ['home'], 'sim.js': ['sim'], 'model.js': ['model'], 'check.js': ['check'], 'about.js': ['about'] };
lazy(Object.fromEntries(Object.entries(PAGES).flatMap(([file, names]) => names.map(n => [n, () => import('./' + file)]))));

route('notfound', () => { app().innerHTML = '<h1>Not found</h1><p><a href="#/">Back to the start</a></p>'; });

$('#mode').textContent = isStatic()
  ? 'static site: the simulator runs in your browser (JS port of calcium_model.py)'
  : 'local app: Python server on this machine';
startBackground($('#bg'));
start();
