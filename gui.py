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

from gui_docks import ParameterDocksMixin
from param_specs import CELL_STATES


class MainWindow(ParameterDocksMixin, QMainWindow):
    def __init__(self, calcium_model):
        super().__init__()
        self.calcium_model = calcium_model
        self._base_dir = os.path.dirname(os.path.abspath(__file__))
        self.initUI()

        configs_dir = os.path.join(self._base_dir, "configs")
        self.cell_states = {name: os.path.join(configs_dir, f) for name, f in CELL_STATES.items()}
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
        # Update initial conditions (resting cytosolic Ca2+ follows from the pump-leak balance)
        self.calcium_model.er_calcium_init = self.initial_er_calcium.value()

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

        # Reset the simulation with new parameters (reset() recomputes the resting Ca2+)
        self.calcium_model.reset()
        self.calcium_model.create_cell_structure()
        self.calcium_model.ip3_conc[:] = self.initial_ip3.value()
        self.resting_calcium.setText(f"{self.calcium_model.eq_calcium:.4f}")

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
            self.sync_widgets_from_model()
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
            self.sync_widgets_from_model()
            self.update_view()

    def create_default_states(self):
        for state_name, filename in self.cell_states.items():
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            if not os.path.exists(filename):
                # Create a default state file
                if state_name == "Blips (default)":
                    self.calcium_model.save_parameters(filename)
                elif state_name == "High Calcium":
                    self.calcium_model.er_calcium_init *= 2
                    self.calcium_model.save_parameters(filename)
                    self.calcium_model.er_calcium_init /= 2
                elif state_name == "Low Calcium":
                    self.calcium_model.er_calcium_init /= 2
                    self.calcium_model.save_parameters(filename)
                    self.calcium_model.er_calcium_init *= 2
