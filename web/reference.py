"""Python reference runs of CalciumModel that the JavaScript engine (web/static/ca_engine.js) is checked against.

CASES are seeded runs (np.random.seed before building the model, as the desktop app's global RNG); run_case()
replays one and summarise() records the fields. scripts/make_ca_fixture.py writes them to
tests/fixtures/ca_reference.json for the node test, and the /api/check endpoint serves one for the in-browser check.
"""
import json
from pathlib import Path

import numpy as np

from calcium_model import CalciumModel

CONFIGS = Path(__file__).resolve().parent.parent / "configs"

# Each case: seed, model parameters, uncaging calls made before stepping, number of steps.
CASES = [
    {"slug": "local-uncaging", "seed": 7, "params": {"grid_size": 40},
     "uncage": [["local", 20, 18, 5, 1.0, 0.02], ["global", 0.3, 0.005]], "steps": 40},
    {"slug": "active-channels", "seed": 3,
     "params": {"grid_size": 32, "ip3r_open_rate": 1000.0, "buffer_total": 20, "ip3r_cluster_density": 0.05},
     "uncage": [["global", 1.0, 0.002]], "steps": 30},
]
GEOMETRY_SEED = 0  # default 200 x 200 cell

# Event-class check of the cell states: global IP3 0.5 uM over 50 ms, then EVENT_STEPS steps of 1 ms, seed 0.
EVENT_STEPS = 400
EVENT_UNCAGING = (0.5, 0.05)


def build(case):
    np.random.seed(case["seed"])
    m = CalciumModel(**case["params"])
    for call in case["uncage"]:
        if call[0] == "local":
            m.add_ip3_local(*call[1:])
        else:
            m.add_ip3_global(*call[1:])
    return m


def run_case(case):
    m = build(case)
    open_trace, substep_trace = [], []
    for _ in range(case["steps"]):
        substep_trace.append(m.calcium_substeps())
        m.step()
        open_trace.append(int(m.ip3r_open.sum()))
    return m, open_trace, substep_trace


def _floats(a):
    return [float(v) for v in np.asarray(a, dtype=float).ravel()]


def summarise(case):
    m, open_trace, substep_trace = run_case(case)
    return {
        "case": case, "eq_calcium": float(m.eq_calcium),
        "er_pixels": int(m.er.sum()), "mito_pixels": int(m.mitochondria.sum()),
        "clusters": [[int(k), int(v)] for k, v in enumerate(m.ip3r_clusters.ravel()) if v],
        "open_trace": open_trace, "substep_trace": substep_trace, "ip3_substeps": m.ip3_substeps(),
        "calcium": _floats(m.calcium), "ip3": _floats(m.ip3_conc), "open": [int(v) for v in m.ip3r_open.ravel()],
        "er_calcium_sum": float(m.er_calcium.sum()), "mito_calcium_sum": float(m.mito_calcium.sum()),
        "buffer_bound_sum": float(m.buffer_bound.sum()),
    }


def geometry():
    """The default 200 x 200 cell for np.random.seed(GEOMETRY_SEED)."""
    np.random.seed(GEOMETRY_SEED)
    m = CalciumModel()
    er = np.flatnonzero(m.er.ravel())
    return {"seed": GEOMETRY_SEED, "eq_calcium": float(m.eq_calcium), "er_pixels": int(er.size),
            "er_index_sum": int(er.sum()), "mito_pixels": int(m.mitochondria.sum()),
            "clusters": [[int(k), int(v)] for k, v in enumerate(m.ip3r_clusters.ravel()) if v],
            "calcium_substeps": m.calcium_substeps(), "ip3_substeps": m.ip3_substeps()}


def event_stats(params, seed=0, steps=EVENT_STEPS):
    """Largest number of channels open at once in one cluster, and the largest area (um^2) where cytosolic [Ca2+]
    exceeds twice its resting level, over `steps` steps after a global IP3 uncaging."""
    np.random.seed(seed)
    m = CalciumModel(**params)
    m.add_ip3_global(*EVENT_UNCAGING)
    max_open, max_area = 0, 0.0
    for _ in range(steps):
        m.step()
        max_open = max(max_open, int(m.ip3r_open.max()))
        max_area = max(max_area, float((m.calcium > 2 * m.eq_calcium).sum()) * m.dx * m.dx)
    return {"max_open": max_open, "max_area_um2": max_area}


def preset_events():
    """event_stats for the Blips and Puffs cell states (configs/)."""
    from param_specs import CELL_STATES
    out = {}
    for name in ("Blips (default)", "Puffs"):
        params = json.loads((CONFIGS / CELL_STATES[name]).read_text())
        out[name] = {"params": params, **event_stats(params)}
    return out


def fixture():
    return {"geometry": geometry(), "cases": [summarise(c) for c in CASES], "events": preset_events()}
