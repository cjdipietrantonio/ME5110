# -*- coding: utf-8 -*-
"""

Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 10: Problem 01
11/04/2025

"""

import numpy as np
from scipy.optimize import fsolve

#define reaction equilibrium constant
Kp = 10 **(-0.97)

#define pressure
P = [1, 10]             #bar

#define equation for extent of reaction

def xi_calc(xi, Kp, P):
    return(((4*(xi**3))/(((1-2*xi)**2)*(1+xi))) - (Kp/P))

for p in P:
    xi_guess = (Kp/(4*p))**(1/3)
    xi_sol = fsolve(xi_calc, xi_guess, args=(Kp, p))[0]
    
    #find number of moles of each species
    N_CO2 = 1 - (2*xi_sol)
    N_CO = 2*xi_sol
    N_O2 = xi_sol
    N_tot = 1 + xi_sol
    
    #find mole fraction for each species
    x_CO2 = N_CO2/N_tot
    x_CO = N_CO/N_tot
    x_O2 = N_O2/N_tot
    
    #calculate sum of mole fractions
    x_sum = x_CO2 + x_CO + x_O2
    
    print(f'\nAt p = {p} bar:')
    print(f'xi = {xi_sol:.5f} kmol')
    print(f'N_CO2 = {N_CO2:.5f} kmol')
    print(f'N_CO = {N_CO:.5f} kmol')
    print(f'N_O2 = {N_O2:.5f} kmol')
    print(f'x_sum = {x_sum:.5f}')