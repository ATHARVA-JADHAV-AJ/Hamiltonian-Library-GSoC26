# ⚛️ The Noisy Boy Forge: Hamiltonian Library - GSoC 2026

### "I thought the dream was left behind. I was wrong." 

This isn't just a code repository; it is the culmination of a journey that started in the 8th standard. While others just watched the science, I decided to build the tools to simulate it. This library is my official submission for **Google Summer of Code 2026**, built and verified on a **16-CPU powerhouse**.

---

## 🛠️ The Forge Environment
To ensure research-grade accuracy, I didn't just install a package; I built the environment from source to ensure maximum performance.
* **Core Engine**: QuTiP 5.3.0.dev0+d654c48
* **Compute Power**: 16 Dedicated CPUs
* **Language**: Python 3.12.3 / Cython 3.2.4
* **OS**: Linux (x86_64) via WSL

---

## 🏛️ The 6 Physics "Legends" 
I have successfully implemented and verified the core Hamiltonians that define modern quantum research. Every model below has been tested for Hilbert space integrity on my local machine.

| Model | Physics Domain | Verified Hilbert Space |
| :--- | :--- | :--- |
| **Jaynes-Cummings** | Quantum Optics (RWA) | (10, 10) |
| **Quantum-Rabi** | Ultra-strong Coupling | (10, 10) |
| **Tavis-Cummings** | Multi-atom Cavities | (12, 12) |
| **Dicke-Collective** | Superradiance / Many-Body | (25, 25) |
| **Heisenberg-Chain** | Spin-Spin Interaction | (8, 8) |
| **Bose-Hubbard** | Lattice Quantum Gases | (8, 8) |

---

## 🚀 Project Vision
The goal of this library is to provide a **standardized interface** for these complex Hamiltonians. 
1.  **Sovereign Architecture**: Using a unified base class for all models to ensure clean, readable code.
2.  **JSON Handshake**: All models export metadata to JSON for portability across research teams.
3.  **Open Systems (Roadmap)**: Moving toward Lindblad Master Equation support to simulate real-world environmental noise.

---

## 📂 Repository Structure
```text
gsoc_submission/
└── hamiltonian_library/
    └── models.py        # The Core Forge containing all Hamiltonian Classes
