import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QPushButton, QSlider, QLabel, QSizePolicy, QGroupBox, QFormLayout,
                             QSpinBox, QDoubleSpinBox, QComboBox, QScrollArea, QCheckBox,
                             QDockWidget, QTabWidget, QMenuBar, QMenu, QAction, QInputDialog, QFileDialog,
                             QMessageBox)
from PyQt5.QtCore import QTimer, Qt
import pyqtgraph as pg
import os
import json

class MainWindow(QMainWindow):
    def __init__(self, calcium_model):
        super().__init__()
        self.calcium_model = calcium_model
        self.initUI()

        self.cell_states = {
            "Default": "default_state.json",
            "High Calcium": "high_calcium_state.json",
            "Low Calcium": "low_calcium_state.json"
        }
        self.create_default_states()
        self.create_cell_state_menu()

    def initUI(self):
        self.setWindowTitle("Calcium Simulator")
        self.setGeometry(100, 100, 1400, 800)

        # Create central widget with visualization
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # Create a tab widget for different views
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)

        # Cytoplasmic calcium view
        cyto_widget = QWidget()
        cyto_layout = QVBoxLayout(cyto_widget)
        self.calcium_view = pg.ImageView()
        self.calcium_view.setLevels(min=0, max=2)
        cyto_layout.addWidget(QLabel("Cytoplasmic [Ca2+]"))
        cyto_layout.addWidget(self.calcium_view)
        self.tab_widget.addTab(cyto_widget, "Cytoplasm")

        # ER calcium view
        er_widget = QWidget()
        er_layout = QVBoxLayout(er_widget)
        self.er_calcium_view = pg.ImageView()
        self.er_calcium_view.setLevels(min=0, max=1000)
        er_layout.addWidget(QLabel("ER [Ca2+]"))
        er_layout.addWidget(self.er_calcium_view)
        self.tab_widget.addTab(er_widget, "ER")

        # Mitochondrial calcium view
        mito_widget = QWidget()
        mito_layout = QVBoxLayout(mito_widget)
        self.mito_calcium_view = pg.ImageView()
        self.mito_calcium_view.setLevels(min=0, max=10)
        mito_layout.addWidget(QLabel("Mitochondrial [Ca2+]"))
        mito_layout.addWidget(self.mito_calcium_view)
        self.tab_widget.addTab(mito_widget, "Mitochondria")

        # IP3 concentration view
        ip3_widget = QWidget()
        ip3_layout = QVBoxLayout(ip3_widget)
        self.ip3_view = pg.ImageView()
        self.ip3_view.setLevels(min=0, max=10)
        ip3_layout.addWidget(QLabel("IP3 Concentration"))
        ip3_layout.addWidget(self.ip3_view)
        self.tab_widget.addTab(ip3_widget, "IP3")

        # Cell feature overlay
        self.feature_overlay = pg.ImageItem()
        self.calcium_view.view.addItem(self.feature_overlay)

        # Create docks
        self.create_control_dock()
        self.create_ip3r_dock()
        self.create_other_params_dock()
        self.create_buffer_dock()
        self.create_ip3_uncaging_dock()
        self.create_feature_visibility_dock()

        # Create menu bar
        self.create_menu_bar()

        # Set up timer for simulation updates
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_simulation)
        self.timer.setInterval(50)  # 50 ms, 20 fps

    def create_control_dock(self):
        dock = QDockWidget("Controls", self)
        widget = QWidget()
        layout = QVBoxLayout(widget)

        self.start_button = QPushButton("Start")
        self.start_button.clicked.connect(self.toggle_simulation)
        layout.addWidget(self.start_button)

        initial_group = QGroupBox("Initial Conditions")
        initial_layout = QFormLayout(initial_group)

        self.initial_calcium = QDoubleSpinBox()
        self.initial_calcium.setRange(0, 10)
        self.initial_calcium.setSingleStep(0.01)
        self.initial_calcium.setValue(self.calcium_model.calcium[0, 0])
        initial_layout.addRow("Initial [Ca2+] (μM):", self.initial_calcium)

        self.initial_er_calcium = QDoubleSpinBox()
        self.initial_er_calcium.setRange(0, 100)
        self.initial_er_calcium.setSingleStep(0.1)
        self.initial_er_calcium.setValue(self.calcium_model.er_calcium[0, 0])
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
        self.ip3r_open_rate.setRange(0, 1)
        self.ip3r_open_rate.setSingleStep(0.01)
        self.ip3r_open_rate.setValue(self.calcium_model.ip3r_open_rate)
        layout.addRow("IP3R Open Rate:", self.ip3r_open_rate)

        self.ip3r_close_rate = QDoubleSpinBox()
        self.ip3r_close_rate.setRange(0, 100)
        self.ip3r_close_rate.setSingleStep(1)
        self.ip3r_close_rate.setValue(self.calcium_model.ip3r_close_rate)
        layout.addRow("IP3R Close Rate:", self.ip3r_close_rate)

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

    def create_menu_bar(self):
        menubar = self.menuBar()

        # File menu
        file_menu = menubar.addMenu('File')

        save_action = QAction('Save Parameters', self)
        save_action.triggered.connect(self.save_parameters)
        file_menu.addAction(save_action)

        load_action = QAction('Load Parameters', self)
        load_action.triggered.connect(self.load_parameters)
        file_menu.addAction(load_action)

        exit_action = QAction('Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # View menu
        view_menu = menubar.addMenu('View')

        toggle_calcium_action = QAction('Toggle Calcium View', self)
        toggle_calcium_action.triggered.connect(self.toggle_calcium_view)
        view_menu.addAction(toggle_calcium_action)

        toggle_ip3_action = QAction('Toggle IP3 View', self)
        toggle_ip3_action.triggered.connect(self.toggle_ip3_view)
        view_menu.addAction(toggle_ip3_action)

        # Simulation menu
        sim_menu = menubar.addMenu('Simulation')

        start_stop_action = QAction('Start/Stop', self)
        start_stop_action.setShortcut('Space')
        start_stop_action.triggered.connect(self.toggle_simulation)
        sim_menu.addAction(start_stop_action)

        reset_action = QAction('Reset', self)
        reset_action.triggered.connect(self.reset_simulation)
        sim_menu.addAction(reset_action)

    def save_simulation(self):
        # TODO: Implement save functionality
        pass

    def load_simulation(self):
        # TODO: Implement load functionality
        pass

    def toggle_calcium_view(self):
        self.calcium_view.setVisible(not self.calcium_view.isVisible())

    def toggle_ip3_view(self):
        self.ip3_view.setVisible(not self.ip3_view.isVisible())

    def update_view(self):
        # Update calcium views
        self.calcium_view.setImage(self.calcium_model.calcium.T, autoLevels=False)
        self.er_calcium_view.setImage(self.calcium_model.er_calcium.T * self.calcium_model.er.T, autoLevels=False)
        self.mito_calcium_view.setImage(self.calcium_model.mito_calcium.T * self.calcium_model.mitochondria.T, autoLevels=False)
        self.ip3_view.setImage(self.calcium_model.ip3_conc.T, autoLevels=False)

        # Create feature overlay
        overlay = np.zeros((*self.calcium_model.calcium.shape, 4), dtype=np.uint8)

        if self.show_ip3r.isChecked():
            overlay[self.calcium_model.ip3r_clusters > 0] = [255, 0, 0, 100]  # Red for closed IP3Rs
            overlay[self.calcium_model.ip3r_open > 0] = [0, 255, 0, 100]  # Green for open IP3Rs

        if self.show_er.isChecked():
            overlay[self.calcium_model.er == 1] = [0, 255, 0, 100]  # Green for ER

        if self.show_mito.isChecked():
            overlay[self.calcium_model.mitochondria == 1] = [0, 0, 255, 100]  # Blue for mitochondria

        if self.show_pm.isChecked():
            overlay[self.calcium_model.pm == 1] = [255, 255, 0, 100]  # Yellow for plasma membrane

        self.feature_overlay.setImage(overlay.transpose(1, 0, 2))

    def toggle_simulation(self):
        if self.timer.isActive():
            self.timer.stop()
            self.start_button.setText("Start")
        else:
            self.timer.start()
            self.start_button.setText("Stop")

    def update_simulation(self):
        self.calcium_model.step()
        self.update_view()

    def add_global_ip3(self):
        amount = self.ip3_amount.value()
        duration = self.ip3_duration.value()
        self.calcium_model.add_ip3_global(amount, duration)

    def add_local_ip3(self):
        x = self.ip3_x.value()
        y = self.ip3_y.value()
        radius = self.ip3_radius.value()
        amount = self.ip3_amount.value()
        duration = self.ip3_duration.value()
        self.calcium_model.add_ip3_local(x, y, radius, amount, duration)

    def apply_settings(self):
        # Update initial conditions
        self.calcium_model.initial_calcium = self.initial_calcium.value()
        self.calcium_model.initial_er_calcium = self.initial_er_calcium.value()
        self.calcium_model.initial_ip3 = self.initial_ip3.value()

        # Update IP3R parameters
        self.calcium_model.ip3r_cluster_density = self.ip3r_cluster_density.value()
        self.calcium_model.ip3r_per_cluster = self.ip3r_per_cluster.value()
        self.calcium_model.ip3r_open_rate = self.ip3r_open_rate.value()
        self.calcium_model.ip3r_close_rate = self.ip3r_close_rate.value()

        # Update other parameters
        self.calcium_model.D_ca = self.d_ca.value()
        self.calcium_model.D_ip3 = self.d_ip3.value()
        self.calcium_model.leak_rate = self.leak_rate.value()
        self.calcium_model.serca_rate = self.serca_rate.value()
        self.calcium_model.serca_k = self.serca_k.value()
        self.calcium_model.ip3_degradation_rate = self.ip3_degradation_rate.value()
        self.calcium_model.pmca_rate = self.pmca_rate.value()
        self.calcium_model.mcu_rate = self.mcu_rate.value()

        # Update buffer conditions
        self.calcium_model.set_buffer_conditions(
            self.buffer_total.value(),
            self.buffer_kd.value(),
            self.buffer_kon.value()
        )

        # Recalculate equilibrium
        self.calcium_model.eq_calcium = self.calcium_model.calculate_equilibrium_calcium()

        # Reset the simulation with new parameters
        self.calcium_model.reset()
        self.calcium_model.create_cell_structure()

        # Update the view
        self.update_view()

    def reset_simulation(self):
        self.calcium_model.reset()
        self.calcium_model.create_cell_structure()
        self.update_view()
        if self.timer.isActive():
            self.timer.stop()
            self.start_button.setText("Start")

    def create_cell_state_menu(self):
        cell_state_menu = self.menuBar().addMenu('Cell States')
        for state_name in self.cell_states:
            action = QAction(state_name, self)
            action.triggered.connect(lambda checked, name=state_name: self.load_cell_state(name))
            cell_state_menu.addAction(action)

        add_state_action = QAction('Add New State', self)
        add_state_action.triggered.connect(self.add_cell_state)
        cell_state_menu.addAction(add_state_action)

        remove_state_action = QAction('Remove State', self)
        remove_state_action.triggered.connect(self.remove_cell_state)
        cell_state_menu.addAction(remove_state_action)

    def load_cell_state(self, state_name):
        filename = self.cell_states[state_name]
        if os.path.exists(filename):
            self.calcium_model.load_parameters(filename)
            self.update_view()
        else:
            QMessageBox.warning(self, "File Not Found", f"The file {filename} does not exist.")


    def add_cell_state(self):
        name, ok = QInputDialog.getText(self, 'Add Cell State', 'Enter name for new cell state:')
        if ok and name:
            filename, _ = QFileDialog.getSaveFileName(self, 'Save Cell State', '', 'JSON Files (*.json)')
            if filename:
                self.calcium_model.save_parameters(filename)
                self.cell_states[name] = filename
                self.create_cell_state_menu()

    def remove_cell_state(self):
        name, ok = QInputDialog.getItem(self, 'Remove Cell State', 'Select state to remove:', list(self.cell_states.keys()), 0, False)
        if ok and name:
            del self.cell_states[name]
            self.create_cell_state_menu()

    def save_parameters(self):
        filename, _ = QFileDialog.getSaveFileName(self, 'Save Parameters', '', 'JSON Files (*.json)')
        if filename:
            self.calcium_model.save_parameters(filename)

    def load_parameters(self):
        filename, _ = QFileDialog.getOpenFileName(self, 'Load Parameters', '', 'JSON Files (*.json)')
        if filename:
            self.calcium_model.load_parameters(filename)
            self.update_view()

    def create_default_states(self):
        for state_name, filename in self.cell_states.items():
            if not os.path.exists(filename):
                # Create a default state file
                if state_name == "Default":
                    self.calcium_model.save_parameters(filename)
                elif state_name == "High Calcium":
                    self.calcium_model.er_calcium_init *= 2
                    self.calcium_model.save_parameters(filename)
                    self.calcium_model.er_calcium_init /= 2
                elif state_name == "Low Calcium":
                    self.calcium_model.er_calcium_init /= 2
                    self.calcium_model.save_parameters(filename)
                    self.calcium_model.er_calcium_init *= 2
