"""Parameter docks for the Calcium Simulator main window (mixin used by gui.MainWindow).

Each create_*_dock method builds one QDockWidget whose widgets are stored on the window; sync_widgets_from_model()
copies the model's current parameters back into them (after loading a cell state or parameter file).
"""
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QCheckBox, QDockWidget, QDoubleSpinBox, QFormLayout, QGroupBox, QLabel, QPushButton,
                             QSpinBox, QVBoxLayout, QWidget)


class ParameterDocksMixin:
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

        self.initial_er_calcium = QDoubleSpinBox()
        self.initial_er_calcium.setRange(0, 5000)
        self.initial_er_calcium.setSingleStep(10)
        self.initial_er_calcium.setValue(self.calcium_model.er_calcium_init)
        initial_layout.addRow("Initial ER [Ca2+] (μM):", self.initial_er_calcium)

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

        self.ip3r_cluster_density = QDoubleSpinBox()
        self.ip3r_cluster_density.setRange(0, 1)
        self.ip3r_cluster_density.setSingleStep(0.01)
        self.ip3r_cluster_density.setValue(self.calcium_model.ip3r_cluster_density)
        layout.addRow("IP3R Cluster Density:", self.ip3r_cluster_density)

        self.ip3r_per_cluster = QSpinBox()
        self.ip3r_per_cluster.setRange(1, 100)
        self.ip3r_per_cluster.setValue(self.calcium_model.ip3r_per_cluster)
        layout.addRow("IP3Rs per Cluster:", self.ip3r_per_cluster)

        self.ip3r_open_rate = QDoubleSpinBox()
        self.ip3r_open_rate.setRange(0, 1000)
        self.ip3r_open_rate.setSingleStep(1)
        self.ip3r_open_rate.setValue(self.calcium_model.ip3r_open_rate)
        layout.addRow("IP3R Open Rate (s^-1):", self.ip3r_open_rate)

        self.ip3r_close_rate = QDoubleSpinBox()
        self.ip3r_close_rate.setRange(0, 100000)
        self.ip3r_close_rate.setSingleStep(100)
        self.ip3r_close_rate.setValue(self.calcium_model.ip3r_close_rate)
        layout.addRow("IP3R Close Rate (s^-1):", self.ip3r_close_rate)

        dock.setWidget(widget)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)

    def create_other_params_dock(self):
        dock = QDockWidget("Other Parameters", self)
        widget = QWidget()
        layout = QFormLayout(widget)

        self.d_ca = QDoubleSpinBox()
        self.d_ca.setRange(0, 1000)
        self.d_ca.setSingleStep(1)
        self.d_ca.setValue(self.calcium_model.D_ca)
        layout.addRow("Ca2+ Diffusion Coefficient:", self.d_ca)

        self.d_ip3 = QDoubleSpinBox()
        self.d_ip3.setRange(0, 1000)
        self.d_ip3.setSingleStep(1)
        self.d_ip3.setValue(self.calcium_model.D_ip3)
        layout.addRow("IP3 Diffusion Coefficient:", self.d_ip3)

        self.leak_rate = QDoubleSpinBox()
        self.leak_rate.setRange(0, 1)
        self.leak_rate.setSingleStep(0.0001)
        self.leak_rate.setValue(self.calcium_model.leak_rate)
        layout.addRow("ER Leak Rate:", self.leak_rate)

        self.serca_rate = QDoubleSpinBox()
        self.serca_rate.setRange(0, 10)
        self.serca_rate.setSingleStep(0.1)
        self.serca_rate.setValue(self.calcium_model.serca_rate)
        layout.addRow("SERCA Pump Rate:", self.serca_rate)

        self.serca_k = QDoubleSpinBox()
        self.serca_k.setRange(0, 10)
        self.serca_k.setSingleStep(0.1)
        self.serca_k.setValue(self.calcium_model.serca_k)
        layout.addRow("SERCA K:", self.serca_k)

        self.ip3_degradation_rate = QDoubleSpinBox()
        self.ip3_degradation_rate.setRange(0, 1)
        self.ip3_degradation_rate.setSingleStep(0.01)
        self.ip3_degradation_rate.setValue(self.calcium_model.ip3_degradation_rate)
        layout.addRow("IP3 Degradation Rate:", self.ip3_degradation_rate)

        self.pmca_rate = QDoubleSpinBox()
        self.pmca_rate.setRange(0, 1)
        self.pmca_rate.setSingleStep(0.01)
        self.pmca_rate.setValue(self.calcium_model.pmca_rate)
        layout.addRow("PMCA Rate:", self.pmca_rate)

        self.mcu_rate = QDoubleSpinBox()
        self.mcu_rate.setRange(0, 1)
        self.mcu_rate.setSingleStep(0.01)
        self.mcu_rate.setValue(self.calcium_model.mcu_rate)
        layout.addRow("MCU Rate:", self.mcu_rate)

        dock.setWidget(widget)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)

    def create_buffer_dock(self):
        dock = QDockWidget("Buffer Conditions", self)
        widget = QWidget()
        layout = QFormLayout(widget)

        self.buffer_total = QDoubleSpinBox()
        self.buffer_total.setRange(0, 1000)
        self.buffer_total.setSingleStep(10)
        self.buffer_total.setValue(self.calcium_model.buffer_total)
        layout.addRow("Total Buffer (μM):", self.buffer_total)

        self.buffer_kd = QDoubleSpinBox()
        self.buffer_kd.setRange(0.01, 10)
        self.buffer_kd.setSingleStep(0.1)
        self.buffer_kd.setValue(self.calcium_model.buffer_kd)
        layout.addRow("Buffer Kd (μM):", self.buffer_kd)

        self.buffer_kon = QDoubleSpinBox()
        self.buffer_kon.setRange(1, 1000)
        self.buffer_kon.setSingleStep(10)
        self.buffer_kon.setValue(self.calcium_model.buffer_kon)
        layout.addRow("Buffer kon (μM^-1 s^-1):", self.buffer_kon)

        dock.setWidget(widget)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)

    def create_ip3_uncaging_dock(self):
        dock = QDockWidget("IP3 Uncaging", self)
        widget = QWidget()
        layout = QFormLayout(widget)

        self.ip3_amount = QDoubleSpinBox()
        self.ip3_amount.setRange(0, 10)
        self.ip3_amount.setSingleStep(0.1)
        self.ip3_amount.setValue(0.5)
        layout.addRow("IP3 Amount (μM):", self.ip3_amount)

        self.ip3_duration = QDoubleSpinBox()
        self.ip3_duration.setRange(0.1, 10)
        self.ip3_duration.setSingleStep(0.1)
        self.ip3_duration.setValue(1)
        layout.addRow("Duration (s):", self.ip3_duration)

        self.ip3_x = QSpinBox()
        self.ip3_x.setRange(0, self.calcium_model.grid_size)
        self.ip3_x.setValue(self.calcium_model.grid_size // 2)
        layout.addRow("X position:", self.ip3_x)

        self.ip3_y = QSpinBox()
        self.ip3_y.setRange(0, self.calcium_model.grid_size)
        self.ip3_y.setValue(self.calcium_model.grid_size // 2)
        layout.addRow("Y position:", self.ip3_y)

        self.ip3_radius = QSpinBox()
        self.ip3_radius.setRange(1, 50)
        self.ip3_radius.setValue(5)
        layout.addRow("Radius:", self.ip3_radius)

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
        self.initial_er_calcium.setValue(m.er_calcium_init)
        for name in ("ip3r_cluster_density", "ip3r_per_cluster", "ip3r_open_rate", "ip3r_close_rate",
                     "leak_rate", "serca_rate", "serca_k", "ip3_degradation_rate", "pmca_rate", "mcu_rate",
                     "buffer_total", "buffer_kd", "buffer_kon"):
            getattr(self, name).setValue(getattr(m, name))
        self.d_ca.setValue(m.D_ca)
        self.d_ip3.setValue(m.D_ip3)
