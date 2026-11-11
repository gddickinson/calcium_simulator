"""Physics and numerics regression tests for CalciumModel (units, stability, steady state, gating, uncaging).

Run: python -m pytest -q tests/test_physics.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calcium_model import CalciumModel  # noqa: E402

CONFIGS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "configs")


def _diffusion_only(model):
    """Switch off every reaction so only Ca2+ diffusion acts."""
    model.serca_rate = model.leak_rate = model.pmca_rate = model.mcu_rate = 0.0
    model.buffer_total = 0.0
    model.buffer_bound = np.zeros_like(model.calcium)
    model.ip3r_clusters[:] = 0
    model.ip3r_open[:] = 0


def test_diffusion_coefficient_is_the_stated_one():
    """d c/dt = D_ca * laplacian(c): for c = 1 + x^2 (x in um) the interior rises at 2 * D_ca uM/s."""
    np.random.seed(1)
    m = CalciumModel(grid_size=41, dx=0.1, dt=1e-5, D_ca=20)
    _diffusion_only(m)
    x = np.arange(41) * m.dx
    m.calcium = 1.0 + np.tile(x ** 2, (41, 1))
    before = m.calcium.copy()
    m.step()
    rate = (m.calcium - before)[5:-5, 5:-5] / m.dt
    np.testing.assert_allclose(rate, 2 * m.D_ca, rtol=1e-6)


def test_ip3_diffusion_coefficient_is_the_stated_one():
    np.random.seed(1)
    m = CalciumModel(grid_size=41, dx=0.1, dt=1e-6, D_ip3=200, ip3_degradation_rate=0.0)
    y = np.arange(41) * m.dx
    m.ip3_conc = 1.0 + np.tile((y[:, None] - 2.0) ** 2 / 4.0, (1, 41))  # stays below the 10 uM clip
    before = m.ip3_conc.copy()
    m.step()
    rate = (m.ip3_conc - before)[5:-5, 5:-5] / m.dt
    np.testing.assert_allclose(rate, 0.5 * m.D_ip3, rtol=1e-5)


def test_default_model_stays_bounded():
    """With the default dx = 0.1 um, dt = 1 ms the explicit scheme must not blow up (it used to turn into a
    0 / 1000 uM checkerboard within 50 steps)."""
    np.random.seed(0)
    m = CalciumModel(grid_size=60)
    for _ in range(100):
        m.step()
    assert m.calcium.max() < 5.0, m.calcium.max()
    assert m.calcium.min() > 0.0
    assert np.isfinite(m.buffer_bound).all() and (m.buffer_bound <= m.buffer_total).all()


def test_local_ip3_spreads_without_overshoot():
    """Diffusion obeys the maximum principle: a local IP3 pulse never exceeds its peak while it spreads."""
    np.random.seed(0)
    m = CalciumModel(grid_size=60)
    m.ip3_conc[25:35, 25:35] = 1.0
    for _ in range(20):
        m.step()
        assert m.ip3_conc.max() <= 1.0 + 1e-12
    assert m.ip3_conc.min() >= 0.0 and m.ip3_conc[0, 0] > 0.0


def test_resting_calcium_balances_serca_and_leak():
    """The initial cytosolic [Ca2+] is the pump-leak steady state: SERCA uptake equals ER leak at ER [Ca2+]_0."""
    m = CalciumModel(grid_size=30)
    c = m.eq_calcium
    serca = m.serca_rate * c ** 2 / (c ** 2 + m.serca_k ** 2)
    leak = m.leak_rate * (m.er_calcium_init - c)
    assert abs(serca - leak) < 1e-12 * leak
    assert 0.05 < c < 0.2  # ~115 nM for the defaults: a physiological resting level


def test_loading_a_state_updates_resting_calcium():
    m = CalciumModel(grid_size=30)
    low = m.eq_calcium
    m.load_parameters(os.path.join(CONFIGS, "high_calcium_state.json"))
    m.grid_size = 30
    m.reset()
    assert m.er_calcium_init == 1000 and m.eq_calcium > low
    assert np.allclose(m.calcium, m.eq_calcium)
