# Calcium Simulator -- Roadmap

## Current State
An interactive 2D intracellular calcium dynamics simulator with PyQt5/pyqtgraph GUI. Three source files: `calcium_model.py` (reaction-diffusion solver with stochastic IP3R gating, ER/mitochondria compartments), `gui.py` (tabbed visualization with overlay controls), and `main.py` (entry point). Includes three predefined cell state JSON files. Well-documented mathematical model with realistic biophysical parameters. Clean, focused codebase.

## Short-term Improvements
- [x] Add `requirements.txt` (PyQt5, pyqtgraph, numpy, scipy, scikit-image)
- [x] Add unit tests for `CalciumModel` -- test diffusion kernel, SERCA flux calculation, IP3R gating statistics
- [ ] Add input validation in `CalciumModel.__init__` -- reject negative rates, zero grid sizes, dt > stability limit
- [ ] Add numerical stability check -- warn if dt is too large relative to dx and diffusion coefficients (CFL condition)
- [ ] Add type hints throughout `calcium_model.py` and `gui.py`
- [ ] Add a `--headless` mode in `main.py` for running simulations without the GUI (save frames to TIFF)
- [ ] Document the mathematical model in a separate `MODEL.md` with equations and parameter references

## Feature Enhancements
- [ ] Add line profile tool -- draw a line on the image and plot Ca2+ concentration along it over time
- [ ] Implement IP3 uncaging simulation -- user clicks to release IP3 at specific locations
- [ ] Add NAADP-sensitive stores and lysosomal calcium release
- [ ] Implement SOCE (store-operated calcium entry) triggered by ER depletion
- [ ] Add ryanodine receptor (RyR) clusters for CICR modeling
- [ ] Support rectangular and non-circular cell geometries
- [ ] Add time-series recording -- plot Ca2+ at selected ROIs over simulation time
- [ ] Implement parameter sweep mode -- run multiple simulations with varying parameters and compare
- [ ] Add colormap selection and adjustable contrast in the GUI

## Long-term Vision
- [ ] Extend to 3D simulation using volumetric rendering (vispy or VTK)
- [ ] Add multi-cell simulation with gap junction coupling
- [ ] Implement real cell morphologies from segmented microscopy images
- [ ] Create a model fitting tool -- optimize parameters to match experimental Ca2+ traces
- [ ] Add SBML/CellML model export for interoperability with other simulators
- [ ] Publish as a Napari plugin for integration with the imaging community
- [ ] Support GPU-accelerated diffusion using CuPy or JAX

## Technical Debt
- [ ] `calcium_model.py` handles both simulation math and cell structure generation -- split `create_cell_structure()` into a separate `cell_geometry.py`
- [ ] `gui.py` mixes visualization setup with simulation control logic -- extract a simulation controller
- [x] `__init__.py` re-exports from a `src` subpackage that does not exist in the current flat structure (fixed to package docstring)
- [x] The three JSON state files at the root should move into a `states/` or `configs/` directory (moved to `configs/`)
- [x] No `.gitignore` for `__pycache__/` and generated outputs
- [x] No package metadata (`setup.py` or `pyproject.toml`) for installation (added `pyproject.toml`)
