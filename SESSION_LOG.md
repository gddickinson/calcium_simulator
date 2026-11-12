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
