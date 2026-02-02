# -*- coding: utf-8 -*-
"""

Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 10: Problem 02
11/04/2025

"""

import numpy as np
import matplotlib.pyplot as plt

#define reaction equilibrium constant
Kp = 6.1

#define pressure
P = np.linspace(0, 25, 1000)        #bar

#calculate extent of reaction over initial moles of N2O4
xi_frac = np.sqrt(Kp/((4*P) + Kp)) 

#plot
plt.figure(figsize=(8,5))
plt.plot(P, xi_frac, '-', label = r'$\xi_{eq}/n_0 = \sqrt{\frac{K_p(T)}{K_p(T) + 4P}}$')
plt.xlabel('Pressure, P (bar)', fontsize = 12)
plt.ylabel(r'Fraction Dissociated, $\xi_{eq}/n_0$', fontsize = 12)
plt.title(r'$\mathrm{N_2O_4 \rightleftharpoons 2NO_2}$ at 373 K ($K_p = 6.1$)', fontsize = 12)
plt.legend(fontsize = 12)
plt.show()
