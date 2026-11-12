// node ca_engine.test.mjs: the JavaScript CalciumModel against seeded runs of the Python one
// (tests/fixtures/ca_reference.json, written by scripts/make_ca_fixture.py from web/reference.py), plus the
// physics checks of tests/test_physics.py repeated on the JS engine.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { CalciumModel, DEFAULTS, equilibriumCalcium, laplacian9 } from './ca_engine.js';
import { dfOverF } from './ca_readout.js';

const F = JSON.parse(readFileSync(new URL('../../tests/fixtures/ca_reference.json', import.meta.url)));
let n = 0;
const ok = (name, fn) => { fn(); n++; console.log('ok  ', name); };
const clusters = m => [...m.ip3r_clusters].flatMap((v, k) => (v ? [[k, v]] : []));

ok('default 200 x 200 cell for seed 0 equals Python (ER, mitochondria, IP3R clusters)', () => {
  const G = F.geometry, m = new CalciumModel({}, G.seed);
  let er = 0, erSum = 0, mito = 0;
  m.er.forEach((v, k) => { if (v) { er++; erSum += k; } });
  m.mitochondria.forEach(v => { mito += v; });
  assert.equal(er, G.er_pixels); assert.equal(erSum, G.er_index_sum); assert.equal(mito, G.mito_pixels);
  assert.deepEqual(clusters(m), G.clusters);
  assert.equal(m.eq_calcium, G.eq_calcium);
  assert.equal(m.calciumSubsteps(), G.calcium_substeps); assert.equal(m.ip3Substeps(), G.ip3_substeps);
});

for (const R of F.cases) {
  const c = R.case;
  ok(`case "${c.slug}": ${c.steps} seeded steps match Python (gating decisions exact, fields to 1e-12)`, () => {
    const m = new CalciumModel(c.params, c.seed);
    for (const call of c.uncage) {
      if (call[0] === 'local') m.addIp3Local(...call.slice(1)); else m.addIp3Global(...call.slice(1));
    }
    assert.deepEqual(clusters(m), R.clusters);
    const opens = [], subs = [];
    for (let s = 0; s < c.steps; s++) {
      subs.push(m.calciumSubsteps());
      m.step();
      opens.push(m.ip3r_open.reduce((a, v) => a + v, 0));
    }
    assert.deepEqual(subs, R.substep_trace);
    assert.deepEqual(opens, R.open_trace);
    assert.deepEqual([...m.ip3r_open], R.open);
    assert.equal(m.ip3Substeps(), R.ip3_substeps);
    let worst = 0;
    for (const [js, py] of [[m.calcium, R.calcium], [m.ip3_conc, R.ip3]]) {
      const scale = Math.max(...py.map(Math.abs));
      py.forEach((v, k) => { worst = Math.max(worst, Math.abs(js[k] - v) / scale); });
    }
    assert.ok(worst < 1e-12, `max relative difference ${worst}`);
    const sum = a => a.reduce((s, v) => s + v, 0);
    for (const [js, key] of [[m.er_calcium, 'er_calcium_sum'], [m.mito_calcium, 'mito_calcium_sum'], [m.buffer_bound, 'buffer_bound_sum']])
      assert.ok(Math.abs(sum(js) - R[key]) <= 1e-12 * Math.abs(R[key]), `${key} ${sum(js)} vs ${R[key]}`);
    console.log(`      max relative difference ${worst.toExponential(2)}, open channels at end ${opens.at(-1)}`);
  });
}

ok('Laplacian of c = x^2 is 2 in the interior (diffusion coefficient means what it says)', () => {
  const N = 21, dx = 0.1, u = new Float64Array(N * N), out = new Float64Array(N * N);
  for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) u[i * N + j] = (j * dx) ** 2;
  laplacian9(u, N, dx, out);
  for (let i = 2; i < N - 2; i++) for (let j = 2; j < N - 2; j++) assert.ok(Math.abs(out[i * N + j] - 2) < 1e-9);
});

ok('resting Ca2+ balances SERCA and leak (0.115 uM for the defaults)', () => {
  const D = DEFAULTS, c = equilibriumCalcium(D.serca_rate, D.serca_k, D.leak_rate, D.er_calcium_init);
  const serca = D.serca_rate * c * c / (c * c + D.serca_k ** 2), leak = D.leak_rate * (D.er_calcium_init - c);
  assert.ok(Math.abs(serca - leak) < 1e-12 * leak); assert.ok(c > 0.11 && c < 0.12);
});

ok('default cell stays bounded and at rest without IP3 (no checkerboard blow-up)', () => {
  const m = new CalciumModel({ grid_size: 60 }, 0);
  for (let s = 0; s < 100; s++) m.step();
  const mx = Math.max(...m.calcium), mn = Math.min(...m.calcium);
  assert.ok(mx < 5 && mn > 0, `${mn}..${mx}`);
});

ok('uncaging delivers amount over duration; overlapping events finish', () => {
  const m = new CalciumModel({ grid_size: 30, D_ip3: 0, ip3_degradation_rate: 0 }, 0);
  m.addIp3Local(10, 10, 2, 0.4, 0.002); m.addIp3Local(20, 20, 2, 0.6, 0.005); m.addIp3Global(0.1, 0.003);
  for (let s = 0; s < 8; s++) m.step();
  assert.equal(m.uncaging.length, 0);
  assert.ok(Math.abs(m.ip3_conc[10 * 30 + 10] - 0.5) < 1e-9 && Math.abs(m.ip3_conc[20 * 30 + 20] - 0.7) < 1e-9);
});

ok('dF/F0 of an equilibrium indicator: 0 at rest, (Rf - 1) c/(c + Kd) scaling, saturates below Rf - 1', () => {
  assert.equal(dfOverF(0.1, 0.1, 0.35, 100), 0);
  const big = dfOverF(1000, 0.1, 0.35, 100);
  const f0 = 1 + 99 * 0.1 / 0.45;
  assert.ok(Math.abs(big - (100 - f0) / f0) < 2e-3 && big < (100 - f0) / f0);
  assert.ok(dfOverF(0.5, 0.1, 0.35, 100) > dfOverF(0.3, 0.1, 0.35, 100));
});

// Event classes of the cell states (tests/test_presets.py): same seeded run as web/reference.event_stats
for (const [name, ev] of Object.entries(F.events)) {
  ok(`"${name}" state: same event statistics as Python, and the expected event class`, () => {
    const m = new CalciumModel(ev.params, 0);
    m.addIp3Global(0.5, 0.05);
    let maxOpen = 0, maxArea = 0;
    for (let s = 0; s < 400; s++) {
      m.step();
      let o = 0, a = 0;
      for (let k = 0; k < m.calcium.length; k++) { if (m.ip3r_open[k] > o) o = m.ip3r_open[k]; if (m.calcium[k] > 2 * m.eq_calcium) a++; }
      maxOpen = Math.max(maxOpen, o); maxArea = Math.max(maxArea, a * m.dx * m.dx);
    }
    assert.equal(maxOpen, ev.max_open); assert.ok(Math.abs(maxArea - ev.max_area_um2) < 1e-9);
    if (name.startsWith('Blips')) assert.ok(maxOpen === 1 && maxArea < 0.5);
    else assert.ok(maxOpen >= 3 && maxArea >= 1 && maxArea < 20);
    console.log(`      max channels open in a cluster ${maxOpen}, max area above 2x rest ${maxArea.toFixed(2)} um^2`);
  });
}

console.log(`${n} tests passed`);
