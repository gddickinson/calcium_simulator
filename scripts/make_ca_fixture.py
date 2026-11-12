"""Write tests/fixtures/ca_reference.json: seeded runs of the Python CalciumModel (web/reference.py) that the node test
web/static/ca_engine.test.mjs compares the JavaScript engine against.

Run after changing calcium_model.py, model_numerics.py or web/reference.py:
    python scripts/make_ca_fixture.py
tests/test_web.py fails if the committed fixture is out of date.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from web.reference import fixture  # noqa: E402

FIXTURE = ROOT / "tests" / "fixtures" / "ca_reference.json"


def main():
    FIXTURE.parent.mkdir(parents=True, exist_ok=True)
    FIXTURE.write_text(json.dumps(fixture(), separators=(",", ":")))
    print(f"wrote {FIXTURE} ({FIXTURE.stat().st_size / 1e3:.0f} kB)")


if __name__ == "__main__":
    main()
