# Calcium Simulator

## Introduction

The Calcium Simulator is a sophisticated tool designed to model and visualize intracellular calcium dynamics in a 2D representation of a cell. This simulator focuses on the interplay between the cytosol, endoplasmic reticulum (ER), and mitochondria, providing a detailed look at calcium signaling processes.

Key features include:
- Simulation of calcium dynamics in cytosol, ER, and mitochondria
- Modeling of IP3 receptor (IP3R) channel behavior
- Visualization of calcium concentrations and cellular structures
- Customizable parameters for tailoring simulations
- Ability to save and load different cell states

This simulator is particularly useful for researchers, students, and educators in the fields of cell biology, biophysics, and computational biology.

## Mathematical Model

The simulator is based on a reaction-diffusion model that describes the spatiotemporal evolution of calcium concentrations. The core equations are:

1. Cytosolic calcium:
d[Ca2+]_cyt/dt = D_Ca ∇²[Ca2+]_cyt + J_IP3R + J_leak - J_SERCA - J_PMCA - J_MCU - J_buffer
2. ER calcium:
d[Ca2+]_ER/dt = -J_IP3R - J_leak + J_SERCA
3. Mitochondrial calcium:
d[Ca2+]_mito/dt = J_MCU
4. IP3 concentration:
d[IP3]/dt = D_IP3 ∇²[IP3] - J_degradation

Where:
- D_Ca and D_IP3 are diffusion coefficients
- J_IP3R is the flux through IP3 receptors
- J_leak is the leak from ER to cytosol
- J_SERCA is the flux through SERCA pumps
- J_PMCA is the flux through plasma membrane Ca2+ ATPase
- J_MCU is the flux through mitochondrial Ca2+ uniporter
- J_buffer represents Ca2+ buffering

The IP3R flux is modeled stochastically, with opening and closing probabilities dependent on cytosolic calcium and IP3 concentrations.

## Installation

1. Clone this repository:
git clone https://github.com/yourusername/calcium-simulator.git

2. Navigate to the project directory:
cd calcium-simulator

3. Create a virtual environment (optional but recommended):
python -m venv venv
source venv/bin/activate  # On Windows use venv\Scripts\activate

4. Install the required packages:
pip install -r requirements.txt

## Usage

1. Run the simulator:
python main.py

2. The main window will appear with several components:
- Visualization panel showing calcium concentrations in different compartments
- Control panel for adjusting simulation parameters
- Menu bar for additional options

3. Adjust parameters in the control panel to set up your simulation.

4. Click the "Start" button to begin the simulation. You can pause it at any time by clicking the same button (now labeled "Stop").

5. Use the tabs in the visualization panel to switch between views of cytoplasmic calcium, ER calcium, mitochondrial calcium, and IP3 concentration.

6. The overlay toggles allow you to visualize different cellular structures (ER, mitochondria, plasma membrane) and IP3R states.

7. Use the "Cell States" menu to load predefined cell states or save your own.

8. The "File" menu allows you to save and load parameter sets.

## Customization

- Modify the `calcium_model.py` file to adjust the underlying mathematical model or add new features.
- Edit the `gui.py` file to change the user interface or add new controls.
- Create new cell states by adjusting parameters and saving them through the interface.

## Contributing

Contributions to improve the Calcium Simulator are welcome. Please fork the repository and submit a pull request with your changes.

## License

[MIT License](LICENSE)

## Contact

For questions or suggestions, please open an issue on the GitHub repository.

## Acknowledgements

This project was inspired by and builds upon various research in the field of calcium signaling.
