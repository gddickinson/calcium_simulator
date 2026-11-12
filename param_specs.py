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
    "Blips (default)": "default_state.json",
    "Puffs": "puffs_state.json",
    "High Calcium": "high_calcium_state.json",
    "Low Calcium": "low_calcium_state.json",
}

# What each cell state shows and where its numbers come from (DOIs checked against Crossref, 2026-11-11).
# There is no "Waves" state: the model's IP3R opening rate rises at most ~2.6-fold with Ca2+ (c^2/(c + 0.3)^3,
# peak at 0.6 uM), too weak for Ca2+-induced Ca2+ release to recruit neighbouring clusters with literature values.
CELL_STATE_INFO = {
    "Blips (default)": {
        "description": "The model's default parameters: clusters of ~10 IP3Rs whose openings last about one 1 ms step, "
                       "so IP3 evokes single-channel blips (one channel open per cluster at a time, sub-um^2 spread).",
        "references": [],
    },
    "Puffs": {
        "description": "Defaults except: 30 IP3Rs per cluster (puffs arise from 25-35 IP3Rs, Shuai et al. 2006; clusters "
                       "of 'two or three tens' of channels, Swillens et al. 1999); IP3R close rate 83.3 s^-1, i.e. a "
                       "12 ms mean open time at high local Ca2+, the trigger-channel opening of Shuai et al. 2006; "
                       "D_IP3 = 283 um^2/s (Allbritton et al. 1992). D_Ca stays 20 um^2/s, inside the 13-65 um^2/s "
                       "measured by Allbritton et al. 1992. Several channels of a cluster open together (puffs; about 6 "
                       "active IP3Rs per puff in Smith & Parker 2009).",
        "references": [
            {"doi": "10.1529/biophysj.106.088880", "cite": "Shuai, Rose & Parker (2006) Biophys J 91:4033"},
            {"doi": "10.1073/pnas.96.24.13750", "cite": "Swillens, Dupont, Combettes & Champeil (1999) PNAS 96:13750"},
            {"doi": "10.1126/science.1465619", "cite": "Allbritton, Meyer & Stryer (1992) Science 258:1812"},
            {"doi": "10.1073/pnas.0810799106", "cite": "Smith & Parker (2009) PNAS 106:6404"},
        ],
    },
    "High Calcium": {"description": "Defaults with the initial ER [Ca2+] doubled to 1000 uM.", "references": []},
    "Low Calcium": {"description": "Defaults with the initial ER [Ca2+] halved to 250 uM.", "references": []},
}
