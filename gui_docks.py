"""Parameter docks for the Calcium Simulator main window (mixin used by gui.MainWindow).

Each create_*_dock method builds one QDockWidget whose widgets are stored on the window; sync_widgets_from_model()
copies the model's current parameters back into them (after loading a cell state or parameter file).
"""
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QCheckBox, QDockWidget, QDoubleSpinBox, QFormLayout, QGroupBox, QLabel, QPushButton,
                             QSpinBox, QVBoxLayout, QWidget)

from param_specs import PARAM_SPECS, UNCAGING_SPECS

# Attribute name of the widget for each model parameter (the rest use the parameter name)
WIDGET_NAMES = {"er_calcium_init": "initial_er_calcium", "D_ca": "d_ca", "D_ip3": "d_ip3"}


def make_spin(spec, value):
    """Spin box with the range, step and decimals of a param_specs entry."""
    label, unit, lo, hi, step, decimals = spec[:6]
    box = QSpinBox() if decimals == 0 else QDoubleSpinBox()
    if decimals:
        box.setDecimals(decimals)
    box.setRange(lo, hi)
    box.setSingleStep(step)
    box.setValue(value)
    return box, f"{label} ({unit}):"


class ParameterDocksMixin:
    def _param_rows(self, layout, keys):
        for key in keys:
            box, text = make_spin(PARAM_SPECS[key], getattr(self.calcium_model, key))
            setattr(self, WIDGET_NAMES.get(key, key), box)
            layout.addRow(text, box)
    def create_control_dock(self):
        dock = QDockWidget("Controls", self)
        widget = QWidget()
        layout = QVBoxLayout(widget)

        self.start_button = QPushButton("Start")
        self.start_button.clicked.connect(self.toggle_simulation)
        layout.addWidget(self.start_button)

        initial_group = QGroupBox("Initial Conditions")
        initial_layout = QFormLayout(initial_group)

        # Resting cytosolic [Ca2+] is not a free input: the model sets it where SERCA uptake balances the ER leak
        self.resting_calcium = QLabel(f"{self.calcium_model.eq_calcium:.4f}")
        initial_layout.addRow("Resting [Ca2+] (μM, SERCA = leak):", self.resting_calcium)

        self._param_rows(initial_layout, ["er_calcium_init"])

        self.initial_ip3 = QDoubleSpinBox()
        self.initial_ip3.setRange(0, 10)
        self.initial_ip3.setSingleStep(0.1)
        self.initial_ip3.setValue(self.calcium_model.ip3_conc[0, 0])
        initial_layout.addRow("Initial [IP3] (μM):", self.initial_ip3)

        layout.addWidget(initial_group)

        self.apply_settings_button = QPushButton("Apply Settings")
        self.apply_settings_button.clicked.connect(self.apply_settings)
        layout.addWidget(self.apply_settings_button)

        self.reset_button = QPushButton("Reset")
        self.reset_button.clicked.connect(self.reset_simulation)
        layout.addWidget(self.reset_button)

        dock.setWidget(widget)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)

    def create_ip3r_dock(self):
        dock = QDockWidget("IP3R Parameters", self)
        widget = QWidget()
        layout = QFormLayout(widget)

        self._param_rows(layout, ["ip3r_cluster_density", "ip3r_per_cluster", "ip3r_open_rate", "ip3r_close_rate"])

        dock.setWidget(widget)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)

    def create_other_params_dock(self):
        dock = QDockWidget("Other Parameters", self)
        widget = QWidget()
        layout = QFormLayout(widget)

        self._param_rows(layout, ["D_ca", "D_ip3", "leak_rate", "serca_rate", "serca_k", "ip3_degradation_rate",
                                  "pmca_rate", "mcu_rate"])

        dock.setWidget(widget)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)

    def create_buffer_dock(self):
        dock = QDockWidget("Buffer Conditions", self)
        widget = QWidget()
        layout = QFormLayout(widget)

        self._param_rows(layout, ["buffer_total", "buffer_kd", "buffer_kon"])

        dock.setWidget(widget)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)

    def create_ip3_uncaging_dock(self):
        dock = QDockWidget("IP3 Uncaging", self)
        widget = QWidget()
        layout = QFormLayout(widget)

        for key in ("amount", "duration"):
            box, text = make_spin(UNCAGING_SPECS[key], UNCAGING_SPECS[key][6])
            setattr(self, "ip3_" + key, box)
            layout.addRow(text, box)

        self.ip3_x = QSpinBox()
        self.ip3_x.setRange(0, self.calcium_model.grid_size)
        self.ip3_x.setValue(self.calcium_model.grid_size // 2)
        layout.addRow("X position:", self.ip3_x)

        self.ip3_y = QSpinBox()
        self.ip3_y.setRange(0, self.calcium_model.grid_size)
        self.ip3_y.setValue(self.calcium_model.grid_size // 2)
        layout.addRow("Y position:", self.ip3_y)

        self.ip3_radius, text = make_spin(UNCAGING_SPECS["radius"], UNCAGING_SPECS["radius"][6])
        layout.addRow(text, self.ip3_radius)

        self.add_global_ip3_button = QPushButton("Add Global IP3")
        self.add_global_ip3_button.clicked.connect(self.add_global_ip3)
        layout.addRow(self.add_global_ip3_button)

        self.add_local_ip3_button = QPushButton("Add Local IP3")
        self.add_local_ip3_button.clicked.connect(self.add_local_ip3)
        layout.addRow(self.add_local_ip3_button)

        dock.setWidget(widget)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)

    def create_feature_visibility_dock(self):
        dock = QDockWidget("Cell Features", self)
        widget = QWidget()
        layout = QVBoxLayout(widget)

        self.show_ip3r = QCheckBox("Show IP3 Receptors")
        self.show_ip3r.stateChanged.connect(self.update_view)
        layout.addWidget(self.show_ip3r)

        self.show_er = QCheckBox("Show ER")
        self.show_er.stateChanged.connect(self.update_view)
        layout.addWidget(self.show_er)

        self.show_mito = QCheckBox("Show Mitochondria")
        self.show_mito.stateChanged.connect(self.update_view)
        layout.addWidget(self.show_mito)

        self.show_pm = QCheckBox("Show Plasma Membrane")
        self.show_pm.stateChanged.connect(self.update_view)
        layout.addWidget(self.show_pm)

        dock.setWidget(widget)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)

    def sync_widgets_from_model(self):
        """Show the model's current parameters in the dock widgets."""
        m = self.calcium_model
        self.resting_calcium.setText(f"{m.eq_calcium:.4f}")
        for key in PARAM_SPECS:
            getattr(self, WIDGET_NAMES.get(key, key)).setValue(getattr(m, key))
