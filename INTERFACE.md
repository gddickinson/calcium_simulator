# Calcium Simulator -- Interface Map

## Project Structure

```
calcium_simulator/
  __init__.py              # Package docstring + version
  main.py                  # Entry point: creates QApplication, CalciumModel, MainWindow
  calcium_model.py         # CalciumModel: reaction-diffusion solver, cell structure, IP3R gating, uncaging
  model_numerics.py        # laplacian9, diffusion_rate, substeps, equilibrium_calcium (SERCA = leak)
  param_specs.py           # PARAM_SPECS / UNCAGING_SPECS (label, unit, range, decimals), CELL_STATES (Qt-free)
  gui.py                   # MainWindow: tabbed image views, menus, simulation timer, apply/load
  gui_docks.py             # ParameterDocksMixin: parameter/uncaging/overlay docks, sync_widgets_from_model
  configs/                 # Cell states: default / high_calcium / low_calcium _state.json
  run_web.py               # Starts the web app (uvicorn) on 127.0.0.1:8205
  requirements-web.txt     # Pinned web/test dependencies
  web/
    app.py                 # FastAPI: /api/info, /api/specs, /api/states, /api/check; STATIC_ENDPOINTS
    reference.py           # Seeded Python reference runs (CASES, summarise, geometry, fixture)
    static/                # Front end (plain ES modules, no build step)
      index.html, style.css (kit theme), style_app.css, core.js (kit router/helpers), main.js, static_api.js, bg.js
      np_random.js         # NumpyRandom: MT19937 + numpy legacy rand/randint/choice/poisson
      ca_geometry.js       # createCellStructure (skimage line, binary dilation), as create_cell_structure()
      ca_engine.js         # CalciumModel JS port (bit-exact), laplacian9, equilibriumCalcium
      ca_readout.js        # dF/F0 of an equilibrium indicator, colour maps
      ca_engine.test.mjs   # node test vs tests/fixtures/ca_reference.json + physics checks
      sim.js / sim_panel.js / sim_view.js  # Simulator page: loop, controls, drawing (field, overlays, trace)
      home.js, model.js, check.js, about.js  # Other pages
  scripts/
    build_static.py        # Exports site/ for GitHub Pages
    screenshots.py         # Playwright smoke test + docs/img screenshots
    make_ca_fixture.py     # Writes tests/fixtures/ca_reference.json
  tests/
    test_calcium_model.py  # Smoke tests for CalciumModel
    test_physics.py        # Units/stability/steady-state/gating/uncaging regression tests
    test_gui.py            # Offscreen Qt: widgets drive the model, loaded states sync
    test_web.py            # API, static build, fixture current, node engine test
  .github/workflows/pages.yml  # Build + publish static site on push to main
```

## Key Classes / Functions

| Name | File | Purpose |
|------|------|---------|
| CalciumModel | calcium_model.py | 2D reaction-diffusion: IP3R gating, SERCA, leak, PMCA, MCU, buffer, IP3 |
| laplacian9, substeps, equilibrium_calcium | model_numerics.py | Stable explicit scheme, resting Ca2+ |
| MainWindow (+ ParameterDocksMixin) | gui.py, gui_docks.py | PyQt5 desktop GUI |
| PARAM_SPECS, CELL_STATES | param_specs.py | Shared by GUI docks and web API |
| CalciumModel (JS), NumpyRandom | web/static/ca_engine.js, np_random.js | Browser engine matching Python |

## Data Flow

1. `main.py` creates `CalciumModel` and passes it to `MainWindow`; the timer calls `step()` at 20 fps
2. `step()`: IP3R gating (rate x dt) -> n Ca2+ sub-steps -> uncaging events -> IP3 sub-steps -> clips
3. Web: `web/app.py` serves specs/states/reference run; the browser runs `ca_engine.js`
4. `scripts/make_ca_fixture.py` -> `tests/fixtures/ca_reference.json` -> `ca_engine.test.mjs`
