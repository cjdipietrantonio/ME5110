# -*- coding: utf-8 -*-
"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 09: Problem 03
10/28/2025

"""

import numpy as np
import matplotlib.pyplot as plt

#define constants
h = 6.626*10**-34           #J-S
k = 1.381*10**-23           #J/K
v_cm = 323                  #cm-1
c = 2.998 * 10**8           #m/s

#calculate fundamental vibrational frequency in Hz
v_hz = v_cm*100*c

#define vibrational quantum numbers
v_max = 30
v = np.arange(0,v_max)

#define temperature
T = [300, 3000]

#calculate and plot distribution
plt.figure(figsize=(8,5))

for t in T:
    E_v = (v + .5)*h*v_hz
    N_i = np.exp(-E_v/(k*t))
    Q = N_i.sum()
    vib_dist = N_i/Q
    plt.plot(v, vib_dist, '-', label=f'T = {t} K')
    
plt.xlabel('Vibrational Quantum Number, v', fontsize = 12)
plt.ylabel(r'$N_v/N$', fontsize = 12)
plt.title('Vibrational Energy Level Distribution of Diatomic Bromine')
plt.legend()
plt.show()

#add anharmonicity if time, need anharmonicity constant