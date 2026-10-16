\# Spacetime Topological Emergent Micro-Solitons (STEMS) — Numerical Solver Suite



Official computational reproduction suite and non-linear boundary-value problem (BVP) solver for the preprint:  

\*\*"Spacetime Topological Emergent Micro-Solitons (STEMS): A Non-Linear Geometric Framework for Quantum Foundations and the Einstein-EPR Program"\*\*  

\*Dev S. Patel (2026)\* | \*\*Affiliation:\*\* Premise Narrative House



\---



\## ⚡ Terminal Execution Commands



Run these commands step-by-step in your terminal (Command Prompt, PowerShell, or Bash):



\### 0. Clone and Install Dependencies

```bash

git clone \[https://github.com/dev-s-patel/STEMS-codes-and-graphs.git](https://github.com/dev-s-patel/STEMS-codes-and-graphs.git)

cd STEMS-codes-and-graphs

pip install -r requirements.txt

```



\---



\### 1. Appendix A: Non-Linear Triad Soliton Relaxation Solver

Solves the coupled second-order non-linear differential equation for the hedgehog chiral angle $F(r)$ bypassing Derrick's theorem via quartic axial Weitzenböck torsion ($\\alpha\_4 = 1.0$), under topological boundary conditions $F(0) = 2\\pi$ and $F(\\infty) = 0$.

```bash

cd "Soliton Relaxation Solver"

python run\_stems\_soliton\_solver.py

cd ..

```

\* \*\*Algorithm:\*\* Adaptive 4th-order collocation Lobatto IIIA boundary-value solver (`scipy.integrate.solve\_bvp`).

\* \*\*Console Output:\*\* Prints core finite energy density $\\rho\_E(0)$ and solver tolerance metrics.

\* \*\*Graphic Output:\*\* Exports Figure 1 (`STEMS\_Soliton\_Convergence.png`) showing:

&#x20; 1. Hedgehog Profile Convergence ($F(r)$ descending smoothly from $2\\pi$ to $0$).

&#x20; 2. Total Static Energy Density $\\rho\_E(r)$ confirming singularity-free finite core saturation ($\\rho\_0 \\approx 0.2912$).



\---



\## 📦 System Requirements

\* Python >= 3.8

\* numpy >= 1.24.0

\* scipy >= 1.10.0

\* matplotlib >= 3.7.0

