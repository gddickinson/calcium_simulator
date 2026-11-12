"""Web interface tests: every exported endpoint answers, the static build is complete, and the JS engine tests pass.

Run: python -m pytest -q tests/test_web.py
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from web.app import STATIC, STATIC_ENDPOINTS, app  # noqa: E402

client = TestClient(app)


def test_index_served():
    r = client.get("/")
    assert r.status_code == 200 and "static/main.js" in r.text


@pytest.mark.parametrize("path", STATIC_ENDPOINTS)
def test_endpoint(path):
    r = client.get("/api" + path)
    assert r.status_code == 200, r.text
    assert r.json() is not None


def test_static_build(tmp_path):
    sys.path.insert(0, str(ROOT / "scripts"))
    import build_static
    build_static.build(tmp_path / "site")
    assert (tmp_path / "site" / "index.html").read_text().count("APP_STATIC") == 1
    for path in STATIC_ENDPOINTS:
        assert (tmp_path / "site" / "api" / (path.strip("/") + ".json")).exists()


def test_static_assets_referenced_exist():
    html = (STATIC / "index.html").read_text()
    for name in ("style.css", "main.js"):
        assert f"static/{name}" in html and (STATIC / name).exists()


@pytest.mark.skipif(not shutil.which("node"), reason="node not installed")
def test_js_engine():
    tests = sorted(STATIC.glob("*.test.mjs"))
    for t in tests:
        r = subprocess.run(["node", str(t)], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr


# --- project-specific ---------------------------------------------------------------------------------------------
import json  # noqa: E402

import numpy as np  # noqa: E402


def test_js_fixture_is_current():
    """tests/fixtures/ca_reference.json must match the current Python model (regenerate with
    python scripts/make_ca_fixture.py). Floats compared with a tolerance: libm differs across platforms."""
    from web.reference import fixture
    saved = json.loads((ROOT / "tests" / "fixtures" / "ca_reference.json").read_text())
    now = json.loads(json.dumps(fixture()))
    assert saved["geometry"] == now["geometry"]
    for a, b in zip(saved["cases"], now["cases"]):
        for key in a:
            if isinstance(a[key], list) and a[key] and isinstance(a[key][0], float):
                np.testing.assert_allclose(a[key], b[key], rtol=1e-9, atol=1e-15, err_msg=key)
            elif isinstance(a[key], float):
                assert abs(a[key] - b[key]) <= 1e-9 * abs(b[key]), key
            else:
                assert a[key] == b[key], key


def test_states_are_the_config_files():
    states = client.get("/api/states").json()
    assert [s["name"] for s in states] == ["Default", "High Calcium", "Low Calcium"]
    for s in states:
        assert s["params"] == json.loads((ROOT / s["file"]).read_text())


def test_specs_cover_model_parameters_with_units():
    specs = client.get("/api/specs").json()
    for key, spec in specs["params"].items():
        assert key in specs["defaults"] and spec["unit"]
        assert spec["min"] <= specs["defaults"][key] <= spec["max"], key


def test_info_resting_calcium():
    info = client.get("/api/info").json()
    assert 0.11 < info["resting_calcium_uM"] < 0.12 and info["cell_um"] == 20.0
