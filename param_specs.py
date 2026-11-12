"""Model parameter specs shared by the desktop GUI (gui_docks.py) and the web app (web/app.py): label, unit,
range, step and decimals for every CalciumModel parameter the user can set, plus the IP3 uncaging controls.
Qt-free so the web server can import it.
"""

# key: (label, unit, min, max, step, decimals, group). decimals = 0 means an integer spin box.
PARAM_SPECS = {
    "er_calcium_init": ("Initial ER [Ca2+]", "uM", 0, 5000, 10, 1, "Initial conditions"),
    "ip3r_cluster_density": ("IP3R cluster density", "fraction of ER pixels", 0, 1, 0.005, 3, "IP3R"),
    "ip3r_per_cluster": ("IP3Rs per cluster", "channels (Poisson mean)", 1, 100, 1, 0, "IP3R"),
    "ip3r_open_rate": ("IP3R open rate", "s^-1", 0, 1000, 1, 2, "IP3R"),
    "ip3r_close_rate": ("IP3R close rate", "s^-1", 0, 100000, 100, 1, "IP3R"),
    "D_ca": ("Ca2+ diffusion coefficient", "um^2/s", 0, 1000, 1, 2, "Diffusion"),
    "D_ip3": ("IP3 diffusion coefficient", "um^2/s", 0, 1000, 1, 2, "Diffusion"),
    "leak_rate": ("ER leak rate", "s^-1", 0, 1, 0.0001, 5, "Pumps and leaks"),
    "serca_rate": ("SERCA pump rate", "uM/s", 0, 10, 0.1, 3, "Pumps and leaks"),
    "serca_k": ("SERCA K", "uM", 0, 10, 0.05, 3, "Pumps and leaks"),
    "pmca_rate": ("PMCA rate", "s^-1", 0, 1, 0.01, 3, "Pumps and leaks"),
    "mcu_rate": ("MCU rate", "s^-1", 0, 1, 0.01, 3, "Pumps and leaks"),
    "ip3_degradation_rate": ("IP3 degradation rate", "s^-1", 0, 1, 0.01, 3, "Diffusion"),
    "buffer_total": ("Total buffer", "uM", 0, 1000, 10, 1, "Buffer"),
    "buffer_kd": ("Buffer Kd", "uM", 0.01, 10, 0.1, 3, "Buffer"),
    "buffer_kon": ("Buffer kon", "uM^-1 s^-1", 1, 1000, 10, 1, "Buffer"),
}

# IP3 uncaging controls: key: (label, unit, min, max, step, decimals, default)
UNCAGING_SPECS = {
    "amount": ("IP3 amount", "uM", 0, 10, 0.1, 2, 0.5),
    "duration": ("Duration", "s", 0.1, 10, 0.1, 2, 1.0),
    "radius": ("Radius", "grid points", 1, 50, 1, 0, 5),
}


def spec_dict(specs=PARAM_SPECS):
    """Specs as JSON-friendly dicts (for the web API)."""
    names = ("label", "unit", "min", "max", "step", "decimals")
    out = {}
    for key, values in specs.items():
        d = dict(zip(names, values[:6]))
        d["group" if specs is PARAM_SPECS else "default"] = values[6]
        out[key] = d
    return out


# Predefined cell states (Cell States menu): name -> file in configs/
CELL_STATES = {
    "Default": "default_state.json",
    "High Calcium": "high_calcium_state.json",
    "Low Calcium": "low_calcium_state.json",
}
