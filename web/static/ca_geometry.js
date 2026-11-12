// Cell structure of CalciumModel.create_cell_structure() (calcium_model.py) in JavaScript: 50 random ER tubules
// (skimage.draw.line, clipped, thickened by two scipy binary_dilation passes with the 4-neighbour cross),
// 20 mitochondria (10 x 4 px blocks), the plasma membrane on the grid edge, and IP3R clusters on a random
// fraction of ER pixels with Poisson-distributed channel counts. Random numbers come from NumpyRandom in the same
// order as numpy's global RNG, so a seeded run gives the same cell as Python.

// skimage.draw.line (Bresenham): rows and columns of the pixels from (r0, c0) to (r1, c1)
export function line(r0, c0, r1, c1) {
  let r = r0, c = c0;
  let dr = Math.abs(r1 - r0), dc = Math.abs(c1 - c0);
  let sc = (c1 - c) > 0 ? 1 : -1, sr = (r1 - r) > 0 ? 1 : -1;
  let steep = false;
  if (dr > dc) { steep = true; [c, r] = [r, c]; [dc, dr] = [dr, dc]; [sc, sr] = [sr, sc]; }
  let d = 2 * dr - dc;
  const rr = new Array(dc + 1), cc = new Array(dc + 1);
  for (let i = 0; i < dc; i++) {
    if (steep) { rr[i] = c; cc[i] = r; } else { rr[i] = r; cc[i] = c; }
    while (d >= 0) { r += sr; d -= 2 * dc; }
    c += sc;
    d += 2 * dr;
  }
  rr[dc] = r1; cc[dc] = c1;
  return [rr, cc];
}

// scipy.ndimage.binary_dilation with the 4-neighbour cross, border_value 0
export function dilateCross(mask, n) {
  const out = new Uint8Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const k = i * n + j;
      out[k] = mask[k] || (i > 0 && mask[k - n]) || (i < n - 1 && mask[k + n])
        || (j > 0 && mask[k - 1]) || (j < n - 1 && mask[k + 1]) ? 1 : 0;
    }
  }
  return out;
}

const clip = (v, lo, hi) => Math.min(Math.max(v, lo), hi);

// Returns { er, mitochondria, pm } as Float64Array 0/1 masks and ip3rClusters as Int32Array channel counts.
export function createCellStructure(rng, n, clusterDensity, perCluster) {
  let er = new Uint8Array(n * n);
  for (let t = 0; t < 50; t++) {
    const x = rng.randint(0, n), y = rng.randint(0, n);
    const length = rng.randint(20, 50);
    const angle = rng.random() * 2 * Math.PI;
    const dx = Math.trunc(length * Math.cos(angle)), dy = Math.trunc(length * Math.sin(angle));
    const [rr, cc] = line(x, y, x + dx, y + dy);
    for (let i = 0; i < rr.length; i++) er[clip(rr[i], 0, n - 1) * n + clip(cc[i], 0, n - 1)] = 1;
  }
  er = dilateCross(dilateCross(er, n), n);

  const mitochondria = new Float64Array(n * n);
  for (let t = 0; t < 20; t++) {
    const x = rng.randint(0, n), y = rng.randint(0, n);
    for (let i = Math.max(0, x - 5); i < Math.min(n, x + 5); i++)
      for (let j = Math.max(0, y - 2); j < Math.min(n, y + 2); j++) mitochondria[i * n + j] = 1;
  }

  const pm = new Float64Array(n * n);
  for (let k = 0; k < n; k++) { pm[k] = pm[(n - 1) * n + k] = pm[k * n] = pm[k * n + n - 1] = 1; }

  const sites = [];
  for (let k = 0; k < n * n; k++) if (er[k]) sites.push(k);
  const nClusters = Math.trunc(sites.length * clusterDensity);
  const ip3rClusters = new Int32Array(n * n);
  for (const idx of rng.choice(sites.length, nClusters)) ip3rClusters[sites[idx]] = rng.poisson(perCluster);

  return { er: Float64Array.from(er), mitochondria, pm, ip3rClusters };
}
