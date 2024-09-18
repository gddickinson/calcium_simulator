import numpy as np
from scipy.ndimage import convolve, binary_dilation
from skimage.draw import line
import json

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

        # Kernel for 2D diffusion
        self.kernel = np.array([[0.05, 0.2, 0.05],
                                [0.2, -1, 0.2],
                                [0.05, 0.2, 0.05]])

        # Calculate equilibrium calcium concentration
        self.eq_calcium = self.calculate_equilibrium_calcium()

        self.reset()
        self.create_cell_structure()

    def calculate_equilibrium_calcium(self):
        # Solve for equilibrium calcium concentration
        # This is a simplified calculation and might need adjustment
        a = self.serca_rate
        b = self.serca_k**2 - self.leak_rate / self.serca_rate
        c = -self.leak_rate * self.serca_k**2 / self.serca_rate

        # Quadratic formula
        eq_calcium = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)
        return eq_calcium

    def reset(self):
        self.calcium = np.full((self.grid_size, self.grid_size), self.eq_calcium, dtype=np.float64)
        self.er_calcium = np.full((self.grid_size, self.grid_size), self.er_calcium_init, dtype=np.float64)
        self.mito_calcium = np.full((self.grid_size, self.grid_size), self.mito_calcium_init, dtype=np.float64)

        self.ip3r_open = np.zeros((self.grid_size, self.grid_size), dtype=np.int32)
        self.ip3r_clusters = np.zeros((self.grid_size, self.grid_size), dtype=np.int32)
        self.ip3_conc = np.zeros((self.grid_size, self.grid_size), dtype=np.float64)
        self.buffer_bound = self.buffer_total * self.eq_calcium / (self.eq_calcium + self.buffer_kd)

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

        # Calcium dynamics
        j_ip3r = 5 * self.ip3r_open * (self.er_calcium - self.calcium) * self.er
        j_leak = self.leak_rate * (self.er_calcium - self.calcium) * self.er
        j_serca = self.serca_rate * (self.calcium**2 / (self.calcium**2 + self.serca_k**2)) * self.er
        j_pmca = self.pmca_rate * self.calcium * self.pm
        j_mcu = self.mcu_rate * (self.calcium - self.mito_calcium) * self.mitochondria

        # Buffer dynamics
        buffer_free = self.buffer_total - self.buffer_bound
        j_buffer = self.buffer_kon * (self.calcium * buffer_free - self.buffer_kd * self.buffer_bound)

        dcdt_ip3r = j_ip3r
        dcdt_leak = j_leak
        dcdt_serca = -j_serca
        dcdt_diff_ca = self.D_ca * convolve(self.calcium, self.kernel) / (self.dx**2)
        dcdt_pmca = -j_pmca
        dcdt_mcu = -j_mcu
        dcdt_buffer = -j_buffer

        self.calcium += (dcdt_ip3r + dcdt_leak + dcdt_serca + dcdt_diff_ca + dcdt_pmca + dcdt_mcu + dcdt_buffer) * self.dt
        self.er_calcium -= (dcdt_ip3r + dcdt_leak + dcdt_serca) * self.dt * self.er
        self.mito_calcium += j_mcu * self.dt

        self.buffer_bound += j_buffer * self.dt

        # IP3 dynamics
        dcdt_diff_ip3 = self.D_ip3 * convolve(self.ip3_conc, self.kernel) / (self.dx**2)
        self.ip3_conc += (dcdt_diff_ip3 - self.ip3_degradation_rate * self.ip3_conc) * self.dt

        # Ensure non-negative concentrations and prevent overflow
        self.calcium = np.clip(self.calcium, 0, 1000)
        self.er_calcium = np.clip(self.er_calcium, 0, 10000)
        self.mito_calcium = np.clip(self.mito_calcium, 0, 1000)
        self.ip3_conc = np.clip(self.ip3_conc, 0, 10)
        self.buffer_bound = np.clip(self.buffer_bound, 0, self.buffer_total)

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
