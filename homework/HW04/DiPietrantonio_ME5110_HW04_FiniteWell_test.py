# -*- coding: utf-8 -*-

"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 04
09/16/2025

"""
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# -------------------------
#PROBLEM 01
# -------------------------

# -------------------------
# FINITE WELL
# -------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# -------------------------
# CONSTANTS
# -------------------------
h = 6.626e-34       # Planck's constant (J·s)
h_bar = h/(2*np.pi)
m = 9.1093897e-31   # Mass of electron (kg)
eV = 1.602177e-19   # Joule to eV

# -------------------------
# WELL DEFINITION
# -------------------------
L = 1e-9            # Well width (m)
n_values = [1, 2, 3]
U0 = 50000000000 * eV         # Deep finite well
x_inside = np.linspace(0, L, 1000)
x_left = np.linspace(-L/5, 0, 200)
x_right = np.linspace(L, L+L/5, 200)

# -------------------------
# FINITE WELL ENERGY LEVELS
# -------------------------
def f(E):
    k_i = np.sqrt(2*m*E)/h_bar
    k_o = np.sqrt(2*m*(U0-E))/h_bar
    return 2*k_o*np.cos(k_i*L) + ((k_o**2 - k_i**2)/k_i)*np.sin(k_i*L)

E_guesses = [n**2 * h**2 / (8*m*L**2) for n in n_values]
E_levels = []
for guess in E_guesses:
    E_root = fsolve(f, guess)[0]
    E_levels.append(E_root)  # Joules

# -------------------------
# COEFFICIENTS AND WAVEFUNCTION
# -------------------------
def compute_coeffs(E):
    k_i = np.sqrt(2*m*E)/h_bar
    k_o = np.sqrt(2*m*(U0-E))/h_bar
    D_val = 1  # fix D scale

    # Match psi and derivative at x=0 and x=L
    M = np.array([
        [1, 0, -1, 0],
        [k_o, -k_i, 0, 0],
        [0, np.sin(k_i*L), np.cos(k_i*L), -1],
        [0, k_i*np.cos(k_i*L), -k_i*np.sin(k_i*L), k_o]
    ])
    b = np.array([D_val, 0, 0, 0])

    # Solve for A, C, G (D fixed)
    M_reduced = M[:, [0,1,3]]
    x, residuals, rank, s = np.linalg.lstsq(M_reduced, b, rcond=None)
    A, C, G = x
    return A, C, D_val, G, k_i, k_o

def wavefunction(A, C, D, G, k_i, k_o):
    psi_in = C*np.sin(k_i*x_inside) + D*np.cos(k_i*x_inside)
    psi_l = A*np.exp(k_o*x_left)
    psi_r = G*np.exp(-k_o*(x_right-L))

    x_total = np.concatenate([x_left, x_inside, x_right])
    psi_total = np.concatenate([psi_l, psi_in, psi_r])

    # Normalize
    dx = x_total[1] - x_total[0]
    psi_total /= np.sqrt(np.sum(np.abs(psi_total)**2)*dx)
    return x_total, psi_total

# -------------------------
# PLOT ENERGY LEVELS
# -------------------------
plt.figure("Finite Well Energy Levels", figsize=(8,5))
for n, E in enumerate(E_levels):
    plt.hlines(E/eV, 0, L*1e9, colors=f'C{n}', label=fr'$n={n_values[n]}$')
plt.xlabel('x (nm)')
plt.ylabel('Energy (eV)')
plt.title('Finite Well Energy Levels')
plt.legend()
plt.tight_layout()
plt.show()

# -------------------------
# PLOT WAVEFUNCTIONS
# -------------------------
plt.figure("Finite Well Wavefunctions", figsize=(8,5))
for n, E in enumerate(E_levels):
    A, C, D, G, k_i, k_o = compute_coeffs(E)
    x_plot, psi_plot = wavefunction(A, C, D, G, k_i, k_o)
    plt.plot(x_plot*1e9, psi_plot + E/eV, label=fr'$n={n_values[n]}$')  # offset by E_n
plt.xlabel('x (nm)')
plt.ylabel('Wavefunction (offset by E_n, eV)')
plt.title('Finite Well Wavefunctions')
plt.legend()
plt.tight_layout()
plt.show()

# -------------------------
# PLOT PROBABILITY DENSITIES
# -------------------------
plt.figure("Finite Well Probability Densities", figsize=(8,5))
for n, E in enumerate(E_levels):
    A, C, D, G, k_i, k_o = compute_coeffs(E)
    x_plot, psi_plot = wavefunction(A, C, D, G, k_i, k_o)
    plt.plot(x_plot*1e9, np.abs(psi_plot)**2, label=fr'$n={n_values[n]}$')
plt.xlabel('x (nm)')
plt.ylabel(r'$|\psi|^2$')
plt.title('Finite Well Probability Densities')
plt.legend()
plt.tight_layout()
plt.show()

# -------------------------
# PLOT MOMENTUM
# -------------------------
plt.figure("Momentum vs x", figsize=(8,5))
for n, E in enumerate(E_levels):
    A, C, D, G, k_i, k_o = compute_coeffs(E)
    # Momentum: inside is real, outside is imaginary
    p_inside = h_bar * k_i * np.ones_like(x_inside)
    p_left   = 1j * h_bar * k_o * np.ones_like(x_left)
    p_right  = 1j * h_bar * k_o * np.ones_like(x_right)

    p_total = np.concatenate([p_left, p_inside, p_right])
    x_total = np.concatenate([x_left, x_inside, x_right])

    plt.plot(x_total*1e9, p_total.real, label=fr'$n={n_values[n]}$')      # real part
    plt.plot(x_total*1e9, np.abs(p_total), '--', color=f'C{n}', alpha=0.5)  # magnitude

plt.xlabel('x (nm)')
plt.ylabel('Momentum (kg·m/s)')
plt.title('Finite Well Momentum')
plt.legend()
plt.tight_layout()
plt.show()
    
#BUT HOW CAN WE DEFINE MOMENTUM OUTSIDE THE WELL?










