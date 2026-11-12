"""FastAPI web app for the Calcium Simulator: serves the front end in static/ and the JSON API under /api.

Run with `python run_web.py` (http://127.0.0.1:8205, local only). `scripts/build_static.py` exports every path in
STATIC_ENDPOINTS to JSON for the GitHub Pages site, so keep that list in step with the GET routes below. The
simulation itself runs in the browser (static/ca_engine.js, a bit-exact port of calcium_model.py); the API serves
the model's defaults, parameter specs and cell states, and a Python reference run for the in-browser check.
"""
import inspect
import json
from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from calcium_model import CalciumModel
from model_numerics import equilibrium_calcium
from param_specs import CELL_STATE_INFO, CELL_STATES, PARAM_SPECS, UNCAGING_SPECS, spec_dict
from web.reference import CASES, summarise

ROOT = Path(__file__).resolve().parent.parent
STATIC = Path(__file__).resolve().parent / "static"

app = FastAPI(title="Calcium Simulator")
app.mount("/static", StaticFiles(directory=STATIC), name="static")


def model_defaults():
    """CalciumModel's constructor defaults."""
    return {k: v.default for k, v in inspect.signature(CalciumModel.__init__).parameters.items() if k != "self"}


def cell_states():
    out = []
    for name, filename in CELL_STATES.items():
        params = json.loads((ROOT / "configs" / filename).read_text())
        out.append({"name": name, "slug": filename.removesuffix("_state.json"), "file": f"configs/{filename}",
                     "params": params, **CELL_STATE_INFO[name]})
    return out


@app.get("/", include_in_schema=False)
def index():
    return FileResponse(STATIC / "index.html")


@app.get("/api/info")
def info():
    """Project summary and counts for the home page (taken from the Python code and configs)."""
    d = model_defaults()
    return {
        "title": "Calcium Simulator",
        "description": "Interactive 2D intracellular calcium dynamics: reaction-diffusion of cytosolic Ca2+ with "
                       "stochastic IP3 receptor clusters on a reticulated ER, SERCA, PMCA and MCU transport, "
                       "buffering, and IP3 diffusion with uncaging.",
        "grid": d["grid_size"], "dx_um": d["dx"], "dt_s": d["dt"], "cell_um": d["grid_size"] * d["dx"],
        "parameters": len(PARAM_SPECS), "cell_states": len(CELL_STATES),
        "resting_calcium_uM": equilibrium_calcium(d["serca_rate"], d["serca_k"], d["leak_rate"], d["er_calcium_init"]),
        "compartments": ["cytosol", "ER", "mitochondria", "IP3"],
    }


@app.get("/api/specs")
def specs():
    """Model defaults, parameter specs (label, unit, range, step) as in the desktop GUI, and uncaging controls."""
    return {"defaults": model_defaults(), "params": spec_dict(PARAM_SPECS), "uncaging": spec_dict(UNCAGING_SPECS)}


@app.get("/api/states")
def states():
    """The predefined cell states of the Cell States menu (configs/*.json) with descriptions and DOI references."""
    return cell_states()


@app.get("/api/check")
def check():
    """A seeded Python run (web/reference.py) that the browser repeats with the JS engine and compares."""
    return summarise(CASES[0])


# GET paths (without the /api prefix) written to the static site as api/<path>.json.
STATIC_ENDPOINTS = ["/info", "/specs", "/states", "/check"]
