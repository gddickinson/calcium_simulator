# Session log

## 2026-11-11 -- Web interface, static site and physics fixes (branch `web-interface`)

Physics / numerics bugs found and fixed (each its own commit with a regression test in tests/test_physics.py
or tests/test_gui.py):
- Diffusion kernel was 0.3 x the 9-point Laplacian (D acted at 30 %), and dt = 1 ms with dx = 0.1 um was beyond
  the explicit limit (D dt/dx^2 = 2 for Ca2+, 20 for IP3; buffer k_on B ~ 1e4 s^-1 stiff): default cell became a
  0/1000 uM checkerboard in 50 steps. Now normalised stencil + stable sub-stepping (model_numerics.py).
- Resting Ca2+ formula dropped ER [Ca2+]: 0.5 nM -> now SERCA = leak balance (115 nM default); recomputed in reset().
- IP3R open/close rates used as per-step probabilities: now rate x dt; defaults restated as 10 / 10000 s^-1
  (identical behaviour at dt = 1 ms).
- IP3 uncaging delivered amount*dt/duration once: now delivers amount over duration; overlapping events crashed.
- GUI: Apply Settings ignored initial conditions (ER box capped at 100 uM), loaded states did not update widgets,
  and the 2-decimal spin box set the ER leak to 0. gui.py split (gui_docks.py) to stay < 500 lines.
- README units (PMCA/MCU are s^-1), kernel and layout corrected.

Web: FastAPI app (web/app.py, port 8205), static GitHub Pages site; JS engine ca_engine.js + np_random.js matches
seeded Python runs bit for bit (node test vs tests/fixtures/ca_reference.json). Pages: home, simulator (dF/F0 /
Ca / ER / mito / IP3 images, overlays, uncaging, probe trace), model, JS = Python check, about.
No linescan (the Python has none). Nothing ran on the Mac mini.

Open questions for George: IP3R close rate 10000 s^-1 (sub-ms openings) and D_Ca 20 um^2/s with an extra
100 uM immobile buffer are the code's values, kept; with them the default cell shows single-channel blips, no waves.

## 2026-11-11 (later) -- Literature presets (George: "fix any other problems")

The default cell only showed single-channel blips. Defaults kept; cell states now "Blips (default)", "Puffs",
"High Calcium", "Low Calcium". Puffs: 30 IP3Rs/cluster (Shuai et al. 2006 doi:10.1529/biophysj.106.088880; Swillens
et al. 1999 doi:10.1073/pnas.96.24.13750), close rate 83.3 s^-1 = 12 ms open time (Shuai 2006), D_IP3 283 um^2/s
(Allbritton et al. 1992 doi:10.1126/science.1465619); Smith & Parker 2009 doi:10.1073/pnas.0810799106 for ~6 channels
per puff. DOIs checked on Crossref, numbers taken from the PubMed abstracts. Several DOIs recalled from memory were
wrong and were replaced by Crossref search results.
Tests (tests/test_presets.py + node): Blips max 1 channel open per cluster, < 0.5 um^2 above 2x rest; Puffs >= 3
open together (7), 1-20 um^2 (3.1). JS reproduces Python's event statistics exactly.
No Waves state: gating rate c^2/(c+0.3)^3 rises <= 2.6x with Ca2+; tried D_Ca 65, ER 1000 uM, 3x cluster density,
20 uM buffer with IP3 2 uM: still local events only (max ~18 um^2). Waves need a different gating model (George).
