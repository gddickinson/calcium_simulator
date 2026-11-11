"""GUI tests (offscreen Qt): the parameter widgets drive the model and follow loaded cell states.

Run: python -m pytest -q tests/test_gui.py
"""
import os
import sys

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

QtWidgets = pytest.importorskip("PyQt5.QtWidgets")
pytest.importorskip("pyqtgraph")


@pytest.fixture(scope="module")
def qapp():
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


_KEEP = []  # windows stay referenced until exit: pyqtgraph prints teardown noise when ImageViews die mid-session


@pytest.fixture
def window(qapp):
    from calcium_model import CalciumModel
    from gui import MainWindow
    np.random.seed(0)
    w = MainWindow(CalciumModel(grid_size=40))
    _KEEP.append(w)
    yield w
    w.hide()


def test_initial_er_widget_shows_model_value(window):
    assert window.initial_er_calcium.value() == window.calcium_model.er_calcium_init == 500


def test_apply_settings_sets_initial_conditions(window):
    m = window.calcium_model
    window.initial_er_calcium.setValue(800)
    window.initial_ip3.setValue(0.4)
    window.apply_settings()
    assert m.er_calcium_init == 800
    assert np.allclose(m.er_calcium, 800)
    assert np.allclose(m.ip3_conc, 0.4)
    assert np.allclose(m.calcium, m.eq_calcium)


def test_loading_a_cell_state_updates_the_widgets(window):
    window.load_cell_state("High Calcium")
    assert window.initial_er_calcium.value() == 1000
    window.apply_settings()  # applying afterwards must keep the loaded state, not the old widget values
    assert window.calcium_model.er_calcium_init == 1000


def test_uncaging_buttons_deliver_ip3(window):
    m = window.calcium_model
    m.ip3_degradation_rate = 0.0
    window.ip3_amount.setValue(0.5)
    window.ip3_duration.setValue(0.1)
    window.add_global_ip3()
    for _ in range(100):
        m.step()
    assert np.allclose(m.ip3_conc, 0.5, rtol=1e-6)
