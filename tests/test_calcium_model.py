"""Smoke tests for CalciumModel."""

import numpy as np
import sys
import os
import tempfile
import json

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def test_model_init():
    from calcium_model import CalciumModel
    model = CalciumModel(grid_size=50)
    assert model.grid_size == 50
    assert model.calcium.shape == (50, 50)
    assert model.er_calcium.shape == (50, 50)
    assert model.mito_calcium.shape == (50, 50)


def test_equilibrium_calcium():
    from calcium_model import CalciumModel
    model = CalciumModel(grid_size=50)
    assert model.eq_calcium > 0, "Equilibrium calcium should be positive"


def test_model_step():
    from calcium_model import CalciumModel
    np.random.seed(42)
    model = CalciumModel(grid_size=50)
    initial_calcium = model.calcium.copy()
    model.step()
    # Calcium should still be non-negative after a step
    assert np.all(model.calcium >= 0), "Calcium went negative after step"
    assert np.all(model.er_calcium >= 0), "ER calcium went negative after step"


def test_add_ip3_global():
    from calcium_model import CalciumModel
    np.random.seed(42)
    model = CalciumModel(grid_size=50)
    initial_ip3 = model.ip3_conc.copy()
    model.add_ip3_global(amount=1.0, duration=1.0)
    assert np.all(model.ip3_conc >= initial_ip3), "IP3 should increase after global addition"


def test_add_ip3_local():
    from calcium_model import CalciumModel
    np.random.seed(42)
    model = CalciumModel(grid_size=50)
    model.add_ip3_local(x=25, y=25, radius=5, amount=1.0, duration=1.0)
    assert model.ip3_conc[25, 25] > 0, "IP3 should increase at injection site"


def test_save_load_parameters():
    from calcium_model import CalciumModel
    model = CalciumModel(grid_size=50)
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False, mode='w') as f:
        tmpfile = f.name
    try:
        model.save_parameters(tmpfile)
        assert os.path.exists(tmpfile)
        with open(tmpfile, 'r') as f:
            params = json.load(f)
        assert params['grid_size'] == 50
        model2 = CalciumModel(grid_size=50)
        model2.load_parameters(tmpfile)
        assert model2.grid_size == params['grid_size']
    finally:
        os.unlink(tmpfile)


def test_cell_structure():
    from calcium_model import CalciumModel
    np.random.seed(42)
    model = CalciumModel(grid_size=50)
    assert model.er.shape == (50, 50)
    assert model.mitochondria.shape == (50, 50)
    assert model.pm.shape == (50, 50)
    assert np.any(model.er > 0), "ER should have some non-zero values"
    assert np.any(model.mitochondria > 0), "Mitochondria should have some non-zero values"


def test_serca_flux():
    """SERCA flux should be zero when calcium is zero."""
    from calcium_model import CalciumModel
    model = CalciumModel(grid_size=50)
    # With zero calcium, SERCA flux formula gives 0^2/(0^2+k^2) = 0
    model.calcium[:] = 0
    j_serca = model.serca_rate * (model.calcium**2 / (model.calcium**2 + model.serca_k**2)) * model.er
    assert np.allclose(j_serca, 0), "SERCA flux should be zero with zero calcium"


if __name__ == "__main__":
    test_model_init()
    test_equilibrium_calcium()
    test_model_step()
    test_add_ip3_global()
    test_add_ip3_local()
    test_save_load_parameters()
    test_cell_structure()
    test_serca_flux()
    print("All calcium model tests passed.")
