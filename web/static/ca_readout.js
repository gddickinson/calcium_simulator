// Fluorescence readout for the web simulator (not part of the Python model): an ideal single-wavelength Ca2+
// indicator in equilibrium with the local free [Ca2+], F = Fmin (1 + (Rf - 1) c / (c + Kd)), Rf = Fmax / Fmin.
// dF/F0 = (F - F0) / F0 with F0 at the resting [Ca2+]. No optical blur, no indicator buffering of the model's Ca2+.

export function fluorescence(c, kd, rf) { return 1 + (rf - 1) * c / (c + kd); }

export function dfOverF(c, c0, kd, rf) {
  const f0 = fluorescence(c0, kd, rf);
  return (fluorescence(c, kd, rf) - f0) / f0;
}

// Colour maps (0..1 -> [r, g, b] 0..255): 'fire' for Ca2+ (black-red-yellow-white), 'ice' for IP3, 'ocean' for ER
const STOPS = {
  fire: [[0, 0, 0, 0], [0.35, 150, 20, 10], [0.65, 255, 120, 0], [0.9, 255, 230, 80], [1, 255, 255, 255]],
  ice: [[0, 0, 0, 0], [0.4, 20, 40, 160], [0.75, 40, 180, 255], [1, 220, 250, 255]],
  green: [[0, 0, 0, 0], [0.5, 0, 140, 60], [1, 160, 255, 140]],
  violet: [[0, 0, 0, 0], [0.5, 110, 40, 170], [1, 240, 180, 255]],
};

export function lut(name) {
  const s = STOPS[name] || STOPS.fire, out = new Uint8ClampedArray(256 * 3);
  for (let i = 0; i < 256; i++) {
    const t = i / 255;
    let k = 0;
    while (k < s.length - 2 && t > s[k + 1][0]) k++;
    const [t0, ...a] = s[k], [t1, ...b] = s[k + 1], f = (t - t0) / (t1 - t0);
    for (let ch = 0; ch < 3; ch++) out[i * 3 + ch] = a[ch] + (b[ch] - a[ch]) * f;
  }
  return out;
}
