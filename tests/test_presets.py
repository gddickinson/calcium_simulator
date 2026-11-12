"""The Blips and Puffs cell states produce their event class (measured after a global IP3 uncaging, seed 0).

Blips: single-channel events, never more than one channel open in a cluster, Ca2+ elevation confined to < 0.5 um^2.
Puffs: several channels of a cluster open together (>= 3; ~6 active IP3Rs per puff, Smith & Parker 2009,
doi:10.1073/pnas.0810799106) and the elevation spreads over >= 1 um^2 but stays local (< 20 um^2, no wave).

Run: python -m pytest -q tests/test_presets.py
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from param_specs import CELL_STATE_INFO, CELL_STATES  # noqa: E402
from web.reference import event_stats  # noqa: E402


def _params(name):
    with open(os.path.join(ROOT, "configs", CELL_STATES[name])) as f:
        return json.load(f)


def test_blips_are_single_channel_events():
    s = event_stats(_params("Blips (default)"))
    assert s["max_open"] == 1 and s["max_area_um2"] < 0.5, s


def test_puffs_are_multi_channel_local_events():
    s = event_stats(_params("Puffs"))
    assert s["max_open"] >= 3 and 1.0 <= s["max_area_um2"] < 20.0, s


def test_every_state_has_a_description_and_doi_references():
    assert set(CELL_STATE_INFO) == set(CELL_STATES)
    for info in CELL_STATE_INFO.values():
        assert info["description"]
        for ref in info["references"]:
            assert ref["doi"].startswith("10.") and ref["cite"]
    assert len(CELL_STATE_INFO["Puffs"]["references"]) >= 2
