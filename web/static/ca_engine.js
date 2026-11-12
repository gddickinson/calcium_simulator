// CalciumModel (calcium_model.py + model_numerics.py) in JavaScript: 2-D reaction-diffusion of cytosolic Ca2+ with
// stochastic IP3R clusters on the ER, SERCA uptake, ER leak, PMCA on the plasma membrane, MCU into mitochondria, an
// immobile buffer and diffusing, degrading IP3 with uncaging. Units: uM, s, um. Every operation is written in the
// same order as the Python so that a seeded run matches it to the last bit (ca_engine.test.mjs, against
// tests/fixtures/ca_reference.json written by scripts/make_ca_fixture.py).
import { NumpyRandom } from './np_random.js';
import { createCellStructure } from './ca_geometry.js';

// CalciumModel.__init__ defaults (units in PARAM_UNITS)
export const DEFAULTS = {
  grid_size: 200, dx: 0.1, dt: 0.001, ip3r_cluster_density: 0.01, ip3r_per_cluster: 10,
  ip3r_open_rate: 10.0, ip3r_close_rate: 10000.0, D_ca: 20, D_ip3: 200, leak_rate: 0.0002,
  serca_rate: 0.4, serca_k: 0.2, ip3_degradation_rate: 0.1, pmca_rate: 0.1, mcu_rate: 0.05,
  buffer_total: 100, buffer_kd: 0.5, buffer_kon: 100, er_calcium_init: 500, mito_calcium_init: 0.1,
};

// model_numerics.diffusion_rate / substeps
export const diffusionRate = (D, dx) => (20.0 * D) / (6.0 * dx * dx);
export const substeps = (dt, rate) => Math.max(1, Math.ceil(dt * rate));

// model_numerics.equilibrium_calcium: SERCA uptake = ER leak, by bisection
export function equilibriumCalcium(serca, k, leak, er) {
  if (er <= 0 || leak <= 0) return 0.0;
  let lo = 0.0, hi = er;
  const k2 = k * k;
  for (let i = 0; i < 200; i++) {
    const mid = 0.5 * (lo + hi);
    if (mid === lo || mid === hi) break;
    if (serca * mid * mid / (mid * mid + k2) - leak * (er - mid) > 0) hi = mid; else lo = mid;
  }
  return 0.5 * (lo + hi);
}

// model_numerics.laplacian9 into `out`: 9-point Laplacian, mirror (zero-flux) boundaries
export function laplacian9(u, n, dx, out) {
  const den = 6.0 * dx * dx;
  for (let i = 0; i < n; i++) {
    const r = i * n, rn = (i > 0 ? i - 1 : 0) * n, rs = (i < n - 1 ? i + 1 : n - 1) * n;
    for (let j = 0; j < n; j++) {
      const jw = j > 0 ? j - 1 : 0, je = j < n - 1 ? j + 1 : n - 1;
      out[r + j] = (4.0 * (((u[rn + j] + u[rs + j]) + u[r + jw]) + u[r + je])
        + (((u[rn + jw] + u[rn + je]) + u[rs + jw]) + u[rs + je]) - 20.0 * u[r + j]) / den;
    }
  }
  return out;
}

const maxOf = a => { let m = -Infinity; for (let i = 0; i < a.length; i++) if (a[i] > m) m = a[i]; return m; };
const minOf = a => { let m = Infinity; for (let i = 0; i < a.length; i++) if (a[i] < m) m = a[i]; return m; };
const clipAll = (a, lo, hi) => { for (let i = 0; i < a.length; i++) a[i] = Math.min(Math.max(a[i], lo), hi); };

export class CalciumModel {
  // params: any subset of DEFAULTS; seed: np.random.seed value used before building the cell (Python:
  // np.random.seed(seed); CalciumModel(**params)), or an existing NumpyRandom to share its stream.
  constructor(params = {}, seed = 0) {
    Object.assign(this, DEFAULTS, params);
    this.rng = seed instanceof NumpyRandom ? seed : new NumpyRandom(seed);
    this.ip3r_flux_rate = 5.0;
    this.reset();
    this.createCellStructure();
  }

  get n() { return this.grid_size; }

  calculateEquilibriumCalcium() {
    return equilibriumCalcium(this.serca_rate, this.serca_k, this.leak_rate, this.er_calcium_init);
  }

  reset() {
    const N = this.grid_size * this.grid_size;
    this.eq_calcium = this.calculateEquilibriumCalcium();
    this.calcium = new Float64Array(N).fill(this.eq_calcium);
    this.er_calcium = new Float64Array(N).fill(this.er_calcium_init);
    this.mito_calcium = new Float64Array(N).fill(this.mito_calcium_init);
    this.ip3r_open = new Int32Array(N);
    if (this.ip3r_clusters) this.ip3r_clusters.fill(0);   // as Python: reset() clears the clusters
    this.ip3_conc = new Float64Array(N);
    this.uncaging = [];
    this.buffer_bound = new Float64Array(N).fill(this.buffer_total * this.eq_calcium / (this.eq_calcium + this.buffer_kd));
    this.time = 0;
    this._lap = new Float64Array(N);
    this._u1 = new Float64Array(N);
    this._u2 = new Float64Array(N);
  }

  createCellStructure() {
    const s = createCellStructure(this.rng, this.grid_size, this.ip3r_cluster_density, this.ip3r_per_cluster);
    this.er = s.er; this.mitochondria = s.mitochondria; this.pm = s.pm; this.ip3r_clusters = s.ip3rClusters;
  }

  step() {
    const N = this.calcium.length, dt = this.dt;
    const cl = this.ip3r_clusters, op = this.ip3r_open;
    // IP3R gating: rate (s^-1) x dt per channel; at most one opening and one closing per cluster per step
    const uOpen = this.rng.fill(this._u1), uClose = this.rng.fill(this._u2);
    for (let k = 0; k < N; k++) {
      const c = Math.min(this.calcium[k], 1000), ip3 = this.ip3_conc[k];
      const a = c + 0.3, b = ip3 + 0.2;
      const openRate = this.ip3r_open_rate * (c * c) * (ip3 * ip3) / ((a * a * a) * (b * b));
      const closeRate = this.ip3r_close_rate * c / (c + 0.3);
      const opening = uOpen[k] < openRate * dt * (cl[k] - op[k]) ? 1 : 0;
      const closing = uClose[k] < closeRate * dt * op[k] ? 1 : 0;
      op[k] = Math.min(Math.max(op[k] + (opening - closing), 0), cl[k]);
    }
    const n = this.calciumSubsteps();
    this.lastSubsteps = n;
    for (let s = 0; s < n; s++) this._calciumSubstep(dt / n);

    for (const ev of this.uncaging) this._applyUncaging(ev);
    this.uncaging = this.uncaging.filter(ev => ev.remaining > 0);
    this._ip3Update();

    clipAll(this.calcium, 0, 1000);
    clipAll(this.er_calcium, 0, 10000);
    clipAll(this.mito_calcium, 0, 1000);
    clipAll(this.ip3_conc, 0, 10);
    clipAll(this.buffer_bound, 0, this.buffer_total);
    this.time += dt;
  }

  calciumSubsteps() {
    let maxFree = -Infinity;
    for (let k = 0; k < this.buffer_bound.length; k++) {
      const f = this.buffer_total - this.buffer_bound[k];
      if (f > maxFree) maxFree = f;
    }
    const rate = diffusionRate(this.D_ca, this.dx)
      + this.buffer_kon * (maxFree + maxOf(this.calcium) + this.buffer_kd)
      + this.ip3r_flux_rate * maxOf(this.ip3r_open)
      + this.leak_rate + this.pmca_rate + this.mcu_rate;
    return substeps(this.dt, rate);
  }

  ip3Substeps() { return substeps(this.dt, diffusionRate(this.D_ip3, this.dx) + this.ip3_degradation_rate); }

  _calciumSubstep(h) {
    const n = this.grid_size, N = n * n;
    const c = this.calcium, erc = this.er_calcium, er = this.er, mito = this.mito_calcium, mm = this.mitochondria;
    const pm = this.pm, op = this.ip3r_open, bb = this.buffer_bound;
    const lap = laplacian9(c, n, this.dx, this._lap);
    const k2 = this.serca_k * this.serca_k;
    for (let k = 0; k < N; k++) {
      const ck = c[k];
      const jIp3r = this.ip3r_flux_rate * op[k] * (erc[k] - ck) * er[k];
      const jLeak = this.leak_rate * (erc[k] - ck) * er[k];
      const jSerca = this.serca_rate * (ck * ck / (ck * ck + k2)) * er[k];
      const jPmca = this.pmca_rate * ck * pm[k];
      const jMcu = this.mcu_rate * (ck - mito[k]) * mm[k];
      const bFree = this.buffer_total - bb[k];
      const jBuffer = this.buffer_kon * (ck * bFree - this.buffer_kd * bb[k]);
      const diff = this.D_ca * lap[k];
      c[k] = ck + (jIp3r + jLeak - jSerca + diff - jPmca - jMcu - jBuffer) * h;
      erc[k] = erc[k] - (jIp3r + jLeak - jSerca) * h * er[k];
      mito[k] = mito[k] + jMcu * h;
      bb[k] = bb[k] + jBuffer * h;
    }
  }

  _ip3Update() {
    const u = this.ip3_conc, k = this.ip3_degradation_rate, N = u.length;
    if (maxOf(u) === minOf(u)) {   // uniform field: no diffusion, degradation only
      for (let i = 0; i < N; i++) u[i] = u[i] - k * u[i] * this.dt;
      return;
    }
    const n = this.ip3Substeps(), h = this.dt / n, lap = this._lap;
    for (let s = 0; s < n; s++) {
      laplacian9(u, this.grid_size, this.dx, lap);
      for (let i = 0; i < N; i++) u[i] = u[i] + (this.D_ip3 * lap[i] - k * u[i]) * h;
    }
  }

  // IP3 uncaging: `amount` uM at a constant rate over `duration` s, everywhere or in a disc (grid points)
  addIp3Global(amount, duration) { this._startUncaging(null, amount, duration); }

  addIp3Local(x, y, radius, amount, duration) {
    const n = this.grid_size, mask = new Uint8Array(n * n);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) mask[i * n + j] = (j - x) ** 2 + (i - y) ** 2 <= radius ** 2 ? 1 : 0;
    this._startUncaging(mask, amount, duration);
  }

  _startUncaging(mask, amount, duration) {
    const ev = { mask, rate: amount / duration, remaining: duration };
    this._applyUncaging(ev);
    if (ev.remaining > 0) this.uncaging.push(ev);
  }

  _applyUncaging(ev) {
    const amt = ev.rate * Math.min(this.dt, ev.remaining), u = this.ip3_conc;
    for (let i = 0; i < u.length; i++) if (!ev.mask || ev.mask[i]) u[i] += amt;
    clipAll(u, 0, 10);
    ev.remaining -= this.dt;
  }

  // Mean of a field over a (2r+1)^2 square around column x, row y (clipped to the grid)
  probe(field, x, y, r = 1) {
    const n = this.grid_size;
    let s = 0, m = 0;
    for (let i = Math.max(0, y - r); i <= Math.min(n - 1, y + r); i++)
      for (let j = Math.max(0, x - r); j <= Math.min(n - 1, x + r); j++) { s += field[i * n + j]; m++; }
    return s / m;
  }
}
