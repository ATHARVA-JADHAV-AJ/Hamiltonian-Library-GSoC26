import qutip as qt
import json
import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. CORE ARCHITECTURE (The Skeleton)
# ==========================================

class HamiltonianBase:
    """The mandatory template ensuring standardization across the library."""
    def __init__(self, name, params=None):
        self.name = name
        self.params = params or {}
        self.H = None

    def to_json(self):
        """Converts model metadata to shareable JSON format for portability."""
        data = {
            "model_name": self.name,
            "parameters": self.params,
            "matrix_shape": str(self.H.shape) if self.H else None,
            "is_hermitian": self.H.isherm if self.H else None
        }
        return json.dumps(data, indent=4)

    def plot_energy_levels(self):
        """Visualizes the eigenvalues (Energy Levels) of the Hamiltonian."""
        if self.H is None:
            print("Error: Build the Hamiltonian first.")
            return
        energies = self.H.eigenenergies()
        plt.figure(figsize=(7, 5))
        plt.plot(energies, 'ro', markersize=4)
        plt.title(f"Energy Spectrum: {self.name}")
        plt.ylabel("Energy (E)")
        plt.xlabel("Eigenvalue Index")
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.show()

# ==========================================
# 2. LIGHT-MATTER LEGENDS (Quantum Optics)
# ==========================================

class JaynesCummings(HamiltonianBase):
    """OG Model: Light-Atom interaction with Rotating Wave Approximation (RWA)."""
    def build(self, N, wc, wa, g):
        self.params = {'N': N, 'wc': wc, 'wa': wa, 'g': g}
        a = qt.tensor(qt.destroy(N), qt.qeye(2))
        sm = qt.tensor(qt.qeye(N), qt.destroy(2))
        self.H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())
        return self.H

class RabiModel(HamiltonianBase):
    """OG Model: The full light-matter interaction (no approximations)."""
    def build(self, N, wc, wa, g):
        self.params = {'N': N, 'wc': wc, 'wa': wa, 'g': g}
        a = qt.tensor(qt.destroy(N), qt.qeye(2))
        sx = qt.tensor(qt.qeye(N), qt.sigmax())
        sz = qt.tensor(qt.qeye(N), qt.sigmaz())
        self.H = wc * a.dag() * a + 0.5 * wa * sz + g * (a.dag() + a) * sx
        return self.H

class TavisCummings(HamiltonianBase):
    """OG Model: Multiple atoms interacting with a single cavity mode (RWA)."""
    def build(self, N_atoms, N_cavity, wc, wa, g):
        self.params = {'N_atoms': N_atoms, 'N_cavity': N_cavity, 'wc': wc, 'wa': wa, 'g': g}
        a = qt.tensor(qt.destroy(N_cavity), qt.qeye(2**N_atoms))
        sz_tot = 0
        sm_tot = 0
        for i in range(N_atoms):
            # Generate collective operators
            sz_i = [qt.qeye(2)] * N_atoms
            sz_i[i] = qt.sigmaz()
            sz_tot += qt.tensor([qt.qeye(N_cavity)] + sz_i)
            
            sm_i = [qt.qeye(2)] * N_atoms
            sm_i[i] = qt.destroy(2)
            sm_tot += qt.tensor([qt.qeye(N_cavity)] + sm_i)
            
        a_field = qt.tensor([qt.destroy(N_cavity)] + [qt.qeye(2)] * N_atoms)
        self.H = wc * a_field.dag() * a_field + 0.5 * wa * sz_tot + g * (a_field.dag() * sm_tot + a_field * sm_tot.dag())
        return self.H

# ==========================================
# 3. MANY-BODY & SPIN SYSTEMS
# ==========================================

class DickeModel(HamiltonianBase):
    """OG Model: N atoms interacting collectively with a single cavity mode."""
    def build(self, N_atoms, N_cavity, wc, wa, g):
        self.params = {'N_atoms': N_atoms, 'N_cavity': N_cavity, 'wc': wc, 'wa': wa, 'g': g}
        j = N_atoms / 2.0
        a = qt.tensor(qt.destroy(N_cavity), qt.qeye(int(2*j + 1)))
        jz = qt.tensor(qt.qeye(N_cavity), qt.jmat(j, 'z'))
        jx = qt.tensor(qt.qeye(N_cavity), qt.jmat(j, 'x'))
        self.H = wc * a.dag() * a + wa * jz + (g / (N_atoms**0.5)) * (a.dag() + a) * jx
        return self.H

class HeisenbergSpinChain(HamiltonianBase):
    """Many-Body Physics: A 1D chain of N spins interacting via neighbors."""
    def build(self, N_spins, Jx, Jy, Jz):
        self.params = {'N_spins': N_spins, 'Jx': Jx, 'Jy': Jy, 'Jz': Jz}
        si = qt.qeye(2); sx = qt.sigmax(); sy = qt.sigmay(); sz = qt.sigmaz()
        sx_list, sy_list, sz_list = [], [], []
        for n in range(N_spins):
            op_list = [si] * N_spins
            op_list[n] = sx; sx_list.append(qt.tensor(op_list))
            op_list[n] = sy; sy_list.append(qt.tensor(op_list))
            op_list[n] = sz; sz_list.append(qt.tensor(op_list))
        self.H = 0
        for n in range(N_spins - 1):
            self.H += Jx * sx_list[n] * sx_list[n+1] + Jy * sy_list[n] * sy_list[n+1] + Jz * sz_list[n] * sz_list[n+1]
        return self.H

# ==========================================
# 4. LATTICE & FERMIONIC MODELS
# ==========================================

class BoseHubbard(HamiltonianBase):
    """Lattice Model: Bosons hopping between sites with site-local interaction."""
    def build(self, N_sites, N_levels, t, U):
        self.params = {'N_sites': N_sites, 'N_levels': N_levels, 't': t, 'U': U}
        # Annihilation operators for each site
        ops = [qt.destroy(N_levels) for _ in range(N_sites)]
        a = [qt.tensor([ops[j] if i == j else qt.qeye(N_levels) for i in range(N_sites)]) for j in range(N_sites)]
        
        H_hop = sum([-t * (a[i].dag() * a[i+1] + a[i+1].dag() * a[i]) for i in range(N_sites-1)])
        H_int = sum([0.5 * U * a[i].dag() * a[i].dag() * a[i] * a[i] for i in range(N_sites)])
        self.H = H_hop + H_int
        return self.H

# ==========================================
# 5. INTEGRATED TEST SUITE (Portfolio Run)
# ==========================================

if __name__ == "__main__":
    print(f"--- NOISY BOY FORGE: GSoC 2026 FULL SPECTRUM TEST ---")
    print(f"System Check: {qt.about()}\n")

    portfolio = [
        JaynesCummings("Jaynes-Cummings"),
        RabiModel("Quantum-Rabi"),
        TavisCummings("Tavis-Cummings"),
        DickeModel("Dicke-Collective"),
        HeisenbergSpinChain("Heisenberg-Chain"),
        BoseHubbard("Bose-Hubbard-Lattice")
    ]

    for model in portfolio:
        # Dynamic parameter routing
        if model.name == "Dicke-Collective":
            H = model.build(N_atoms=4, N_cavity=5, wc=1.0, wa=1.0, g=0.5)
        elif model.name == "Heisenberg-Chain":
            H = model.build(N_spins=3, Jx=1.0, Jy=1.0, Jz=1.0)
        elif model.name == "Bose-Hubbard-Lattice":
            H = model.build(N_sites=3, N_levels=2, t=1.0, U=2.0)
        elif model.name == "Tavis-Cummings":
            H = model.build(N_atoms=2, N_cavity=3, wc=1.0, wa=1.0, g=0.1)
        else:
            H = model.build(N=5, wc=1.0, wa=1.0, g=0.1)

        print(f"✅ {model.name} Generated. Hilbert Space: {H.shape}")
        # print(model.to_json()) 

    print("\n--- ALL SYSTEMS OPERATIONAL ---")
    # Visualization check
    portfolio[0].plot_energy_levels()