# -*- coding: utf-8 -*-

"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 09: Problem 01
10/28/2025

"""

import numpy as np
import matplotlib.pyplot as plt

#define constants
h = 6.626*10**-34           #J-S
k = 1.381*10**-23           #J/K
R = 1.13*10**-10            #m
m_12C_amu = 12                  #amu
m_16O_amu = 15.99491461957      #amu
m_u = 1.6605402*10**-27         #kg/amu

#calculate reduced planck constant
h_bar = h/(2*np.pi)

#calculate reduced mass
m_12C = m_12C_amu*m_u
m_16O = m_16O_amu*m_u
mu = (m_12C*m_16O)/(m_12C + m_16O)

#calculate momentum of inertia
I = mu*R**2

#define rotational quantum numbers
J = np.arange(0,30)

#define temperature
T = [100, 300]

#calculate and plot distribution
plt.figure(figsize=(8,5))

for t in T:
    g_rot = 2*J + 1
    E_J = ((h_bar**2)/(2*I))*(J*(J+1))
    N_i = g_rot*np.exp(-E_J/(k*t))
    Q = N_i.sum()
    rot_dist = N_i/Q
    plt.plot(J, rot_dist, '-', label=f'T = {t} K')
    
plt.xlabel('Rotational Quantum Number, J', fontsize = 12)
plt.ylabel(r'$N_J/N$', fontsize = 12)
plt.title('Rotational Energy Level Distribution of CO')
plt.legend()
plt.show()

#add Jmax if time

