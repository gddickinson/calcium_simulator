# Calcium Simulator -- Interface Map

## Project Structure

```
calcium_simulator/
  __init__.py              # Package docstring + version
  main.py                  # Entry point: creates QApplication, CalciumModel, MainWindow
  calcium_model.py         # CalciumModel: reaction-diffusion solver (222 lines)
  gui.py                   # MainWindow: PyQt5 tabbed GUI with parameter controls (537 lines)
  pyproject.toml           # Package metadata and dependencies
  requirements.txt         # Pinned dependency versions
  configs/                 # Predefined cell state JSON files
    default_state.json
    high_calcium_state.json
    low_calcium_state.json
  examples/                # Example usage
  tests/
    test_calcium_model.py  # Smoke tests for CalciumModel
```

## Key Classes

| Class | File | Purpose |
|-------|------|---------|
| CalciumModel | calcium_model.py | 2D reaction-diffusion solver: IP3R gating, SERCA, ER/mito compartments, buffering |
| MainWindow | gui.py | PyQt5 GUI: tabbed image views, parameter sliders, IP3 uncaging, state save/load |

## Data Flow

1. `main.py` creates `CalciumModel` and passes it to `MainWindow`
2. `MainWindow.timer` calls `update_simulation()` at 20fps
3. Each tick calls `CalciumModel.step()` which updates calcium, ER, mito, IP3, buffers
4. GUI reads model arrays and displays via pyqtgraph ImageView
5. JSON state files in `configs/` store/restore parameter sets
