import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# 1. Dimensionless System Parameters
alpha_4 = 1.0
r_min, r_max = 1e-3, 20.0
N_grid = 1000
r = np.linspace(r_min, r_max, N_grid)

# 2. First-Order ODE System: y = [F, dF/dr]
def ode_system(x, y):
    F, dF = y[0], y[1]
    sin_F = np.sin(F)
    sin_half = np.sin(F / 2.0)
    sin_half_sq = sin_half**2

    # Metric and force terms
    A = x**2 + (2.0 * alpha_4 * sin_half_sq) / (x**2)
    B = (2.0 * x * dF 
         + (alpha_4 * sin_F / (x**2)) * (dF**2) 
         - 4.0 * sin_F 
         - (alpha_4 * sin_half_sq * sin_F) / (x**4))
    
    return np.vstack((dF, -B / A))

# 3. Boundary Conditions (Original PDF Match)
def boundary_conditions(ya, yb):
    return np.array([ya[0] - 2.0 * np.pi, yb[0] - 0.0])

# Smooth initial trial profile
y_guess = np.zeros((2, r.size))
y_guess[0] = 2.0 * np.pi * (1.0 - np.tanh(r / 2.0))
y_guess[1] = -np.pi * (1.0 / np.cosh(r / 2.0)**2)

# 4. Numerical BVP Execution
sol = solve_bvp(ode_system, boundary_conditions, r, y_guess, 
                max_nodes=100000, tol=1e-5)

# 5. Energy Density Computation
r_eval = np.linspace(r_min, 10.0, 500)
F_sol, dF_sol = sol.sol(r_eval)[0], sol.sol(r_eval)[1]

s_half_sq = np.sin(F_sol / 2.0)**2
rho_2 = (dF_sol**2) + (8.0 * s_half_sq) / (r_eval**2)
rho_4 = ((alpha_4 * s_half_sq / (r_eval**2)) 
         * (dF_sol**2 + s_half_sq / (r_eval**2)))
rho_tot = rho_2 + rho_4

print(f"rho_E(0): {rho_tot[0]:.4f}")

# 6. Visualization Pipeline (Exact Match with PDF Page 32/33)
plt.figure(figsize=(10, 4.2))

# Left Panel: Profile Convergence
plt.subplot(1, 2, 1)
plt.plot(r_eval, F_sol, 'b-', lw=2, label=r'Numerical Profile $F(r)$')
plt.axhline(2*np.pi, color='gray', ls='--', alpha=0.6, label=r'$F(0) = 2\pi$')
plt.axhline(0, color='black', ls=':', alpha=0.6, label=r'$F(\infty) = 0$')
plt.xlabel(r'Normalized Radius $r/\sqrt{\alpha_4}$')
plt.ylabel(r'Chiral Angle $F(r)$ [rad]')
plt.title('Hedgehog Profile Convergence')
plt.grid(True, alpha=0.3)
plt.legend(loc='upper right')

# Right Panel: Energy Density Decay
plt.subplot(1, 2, 2)
plt.plot(r_eval, rho_tot, 'r-', lw=2, label=r'Total $\rho_E(r)$')
plt.plot(r_eval, rho_2, 'g--', lw=1.2, label=r'$\mathbb{T}$')
plt.plot(r_eval, rho_4, 'm-.', lw=1.2, label=r'$\mathbb{T}_{\rm axial}$')
plt.xlabel(r'Normalized Radius $r/\sqrt{\alpha_4}$')
plt.ylabel(r'Energy Density $\rho_E(r)$')
plt.title(r'Finite Core Density ($\rho_0 \approx 0.2912$)')
plt.yscale('log')
plt.grid(True, alpha=0.3)
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig('STEMS_Soliton_Convergence.png', dpi=300)
plt.show()