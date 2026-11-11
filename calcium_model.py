import numpy as np
from scipy.ndimage import binary_dilation
from skimage.draw import line
import json

from model_numerics import LAPLACIAN_KERNEL, diffusion_rate, equilibrium_calcium, laplacian9, substeps

class CalciumModel:
    def __init__(self, grid_size=200, dx=0.1, dt=0.001,
                 ip3r_cluster_density=0.01, ip3r_per_cluster=10,
                 ip3r_open_rate=0.01, ip3r_close_rate=10,
                 D_ca=20, D_ip3=200, leak_rate=0.0002,
                 serca_rate=0.4, serca_k=0.2,
                 ip3_degradation_rate=0.1, pmca_rate=0.1, mcu_rate=0.05,
                 buffer_total=100, buffer_kd=0.5, buffer_kon=100,
                 er_calcium_init=500, mito_calcium_init=0.1):

        self.grid_size = grid_size
        self.dx = dx
        self.dt = dt

        # IP3R parameters
        self.ip3r_cluster_density = ip3r_cluster_density
        self.ip3r_per_cluster = ip3r_per_cluster
        self.ip3r_open_rate = ip3r_open_rate
        self.ip3r_close_rate = ip3r_close_rate
        self.ip3r_open = np.zeros((self.grid_size, self.grid_size), dtype=np.int32)


        # Other parameters
        self.D_ca = D_ca
        self.D_ip3 = D_ip3
        self.leak_rate = leak_rate
        self.serca_rate = serca_rate
        self.serca_k = serca_k
        self.ip3_degradation_rate = ip3_degradation_rate
        self.pmca_rate = pmca_rate
        self.mcu_rate = mcu_rate

        # Buffer parameters
        self.buffer_total = buffer_total
        self.buffer_kd = buffer_kd
        self.buffer_kon = buffer_kon

        self.er_calcium_init = er_calcium_init
        self.mito_calcium_init = mito_calcium_init

        # IP3R Ca2+ flux per open channel (s^-1, times the ER-cytosol gradient)
        self.ip3r_flux_rate = 5.0

        # Kernel for 2D diffusion: the isotropic 9-point Laplacian (times 1/dx^2), see model_numerics.py
        self.kernel = LAPLACIAN_KERNEL

        # Calculate equilibrium calcium concentration
        self.eq_calcium = self.calculate_equilibrium_calcium()

        self.reset()
        self.create_cell_structure()

    def calculate_equilibrium_calcium(self):
        """Resting cytosolic [Ca2+] (uM): SERCA uptake balances the ER leak at the initial ER [Ca2+]."""
        return equilibrium_calcium(self.serca_rate, self.serca_k, self.leak_rate, self.er_calcium_init)

    def reset(self):
        self.eq_calcium = self.calculate_equilibrium_calcium()
        self.calcium = np.full((self.grid_size, self.grid_size), self.eq_calcium, dtype=np.float64)
        self.er_calcium = np.full((self.grid_size, self.grid_size), self.er_calcium_init, dtype=np.float64)
        self.mito_calcium = np.full((self.grid_size, self.grid_size), self.mito_calcium_init, dtype=np.float64)

        self.ip3r_open = np.zeros((self.grid_size, self.grid_size), dtype=np.int32)
        self.ip3r_clusters = np.zeros((self.grid_size, self.grid_size), dtype=np.int32)
        self.ip3_conc = np.zeros((self.grid_size, self.grid_size), dtype=np.float64)
        self.buffer_bound = np.full((self.grid_size, self.grid_size),
                                    self.buffer_total * self.eq_calcium / (self.eq_calcium + self.buffer_kd))

    def create_cell_structure(self):
        # Create reticulated ER
        self.er = np.zeros((self.grid_size, self.grid_size), dtype=np.float64)
        for _ in range(50):  # Add 50 ER tubules
            x, y = np.random.randint(0, self.grid_size, 2)
            length = np.random.randint(20, 50)
            angle = np.random.rand() * 2 * np.pi
            dx, dy = int(length * np.cos(angle)), int(length * np.sin(angle))
            rr, cc = line(x, y, x + dx, y + dy)
            rr = np.clip(rr, 0, self.grid_size - 1)
            cc = np.clip(cc, 0, self.grid_size - 1)
            self.er[rr, cc] = 1
        self.er = binary_dilation(self.er, iterations=2)  # Thicken ER tubules

        # Create mitochondria
        self.mitochondria = np.zeros((self.grid_size, self.grid_size), dtype=np.float64)
        for _ in range(20):  # Add 20 mitochondria
            x, y = np.random.randint(0, self.grid_size, 2)
            self.mitochondria[max(0, x-5):min(self.grid_size, x+5),
                              max(0, y-2):min(self.grid_size, y+2)] = 1

        # Create plasma membrane
        self.pm = np.zeros((self.grid_size, self.grid_size), dtype=np.float64)
        self.pm[0, :] = self.pm[-1, :] = self.pm[:, 0] = self.pm[:, -1] = 1

        # Create IP3R clusters on ER
        self.ip3r_clusters = np.zeros((self.grid_size, self.grid_size), dtype=np.int32)
        er_sites = np.where(self.er == 1)
        num_clusters = int(len(er_sites[0]) * self.ip3r_cluster_density)
        cluster_indices = np.random.choice(len(er_sites[0]), num_clusters, replace=False)
        for idx in cluster_indices:
            self.ip3r_clusters[er_sites[0][idx], er_sites[1][idx]] = np.random.poisson(self.ip3r_per_cluster)


    def step(self):
        # IP3R dynamics
        open_prob = self.ip3r_open_rate * np.minimum(self.calcium, 1000)**2 * self.ip3_conc**2 / \
                    ((np.minimum(self.calcium, 1000) + 0.3)**3 * (self.ip3_conc + 0.2)**2)
        close_prob = self.ip3r_close_rate * np.minimum(self.calcium, 1000) / (np.minimum(self.calcium, 1000) + 0.3)

        opening = (np.random.rand(*self.ip3r_open.shape) < open_prob * (self.ip3r_clusters - self.ip3r_open)).astype(np.int32)
        closing = (np.random.rand(*self.ip3r_open.shape) < close_prob * self.ip3r_open).astype(np.int32)

        self.ip3r_open += opening - closing
        self.ip3r_open = np.clip(self.ip3r_open, 0, self.ip3r_clusters)

        # Calcium, ER, mitochondria and buffer: forward Euler in n sub-steps small enough to be stable
        n = self.calcium_substeps()
        for _ in range(n):
            self._calcium_substep(self.dt / n)

        # IP3: diffusion and degradation, sub-stepped for stability (D_ip3 is ten times D_ca)
        self._ip3_update()

        # Ensure non-negative concentrations and prevent overflow
        self.calcium = np.clip(self.calcium, 0, 1000)
        self.er_calcium = np.clip(self.er_calcium, 0, 10000)
        self.mito_calcium = np.clip(self.mito_calcium, 0, 1000)
        self.ip3_conc = np.clip(self.ip3_conc, 0, 10)
        self.buffer_bound = np.clip(self.buffer_bound, 0, self.buffer_total)

    def calcium_substeps(self):
        """Sub-steps per dt for the Ca2+ update: covers diffusion, buffer binding (k_on (B_free + Ca) + k_off),
        IP3R release and the first-order pumps/leaks, from the current state."""
        b_free = self.buffer_total - self.buffer_bound
        rate = (diffusion_rate(self.D_ca, self.dx)
                + self.buffer_kon * (float(np.max(b_free)) + float(np.max(self.calcium)) + self.buffer_kd)
                + self.ip3r_flux_rate * float(np.max(self.ip3r_open))
                + self.leak_rate + self.pmca_rate + self.mcu_rate)
        return substeps(self.dt, rate)

    def ip3_substeps(self):
        return substeps(self.dt, diffusion_rate(self.D_ip3, self.dx) + self.ip3_degradation_rate)

    def _calcium_substep(self, h):
        c = self.calcium
        er_c = self.er_calcium
        j_ip3r = self.ip3r_flux_rate * self.ip3r_open * (er_c - c) * self.er
        j_leak = self.leak_rate * (er_c - c) * self.er
        j_serca = self.serca_rate * (c * c / (c * c + self.serca_k * self.serca_k)) * self.er
        j_pmca = self.pmca_rate * c * self.pm
        j_mcu = self.mcu_rate * (c - self.mito_calcium) * self.mitochondria
        b_free = self.buffer_total - self.buffer_bound
        j_buffer = self.buffer_kon * (c * b_free - self.buffer_kd * self.buffer_bound)
        diff = self.D_ca * laplacian9(c, self.dx)

        self.calcium = c + (j_ip3r + j_leak - j_serca + diff - j_pmca - j_mcu - j_buffer) * h
        self.er_calcium = er_c - (j_ip3r + j_leak - j_serca) * h * self.er
        self.mito_calcium = self.mito_calcium + j_mcu * h
        self.buffer_bound = self.buffer_bound + j_buffer * h

    def _ip3_update(self):
        u = self.ip3_conc
        k = self.ip3_degradation_rate
        if u.max() == u.min():  # uniform field: no diffusion, degradation only
            self.ip3_conc = u - k * u * self.dt
            return
        n = self.ip3_substeps()
        h = self.dt / n
        for _ in range(n):
            u = u + (self.D_ip3 * laplacian9(u, self.dx) - k * u) * h
        self.ip3_conc = u

    def add_ip3_global(self, amount, duration):
        """Simulate global uncaging of IP3"""
        rate = amount / duration
        self.ip3_conc += rate * self.dt
        self.ip3_conc = np.clip(self.ip3_conc, 0, 10)

    def add_ip3_local(self, x, y, radius, amount, duration):
        """Simulate local uncaging of IP3"""
        xx, yy = np.meshgrid(np.arange(self.grid_size), np.arange(self.grid_size))
        mask = ((xx - x)**2 + (yy - y)**2 <= radius**2)
        rate = amount / duration
        self.ip3_conc[mask] += rate * self.dt
        self.ip3_conc = np.clip(self.ip3_conc, 0, 10)

    def set_buffer_conditions(self, total, kd, kon):
        self.buffer_total = total
        self.buffer_kd = kd
        self.buffer_kon = kon
        self.reset()  # Reset the simulation with new buffer conditions

    def save_parameters(self, filename):
        params = {
            'grid_size': self.grid_size,
            'dx': self.dx,
            'dt': self.dt,
            'ip3r_cluster_density': self.ip3r_cluster_density,
            'ip3r_per_cluster': self.ip3r_per_cluster,
            'ip3r_open_rate': self.ip3r_open_rate,
            'ip3r_close_rate': self.ip3r_close_rate,
            'D_ca': self.D_ca,
            'D_ip3': self.D_ip3,
            'leak_rate': self.leak_rate,
            'serca_rate': self.serca_rate,
            'serca_k': self.serca_k,
            'ip3_degradation_rate': self.ip3_degradation_rate,
            'pmca_rate': self.pmca_rate,
            'mcu_rate': self.mcu_rate,
            'buffer_total': self.buffer_total,
            'buffer_kd': self.buffer_kd,
            'buffer_kon': self.buffer_kon,
            'er_calcium_init': self.er_calcium_init,
            'mito_calcium_init': self.mito_calcium_init
        }
        with open(filename, 'w') as f:
            json.dump(params, f)

    def load_parameters(self, filename):
        try:
            with open(filename, 'r') as f:
                params = json.load(f)
            for key, value in params.items():
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(f"Warning: Unknown parameter '{key}' in file {filename}")
            self.reset()
            self.create_cell_structure()
        except FileNotFoundError:
            print(f"Error: File {filename} not found.")
        except json.JSONDecodeError:
            print(f"Error: File {filename} is not a valid JSON file.")
        except Exception as e:
            print(f"Error loading parameters from {filename}: {str(e)}")
