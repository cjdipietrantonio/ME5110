# -*- coding: utf-8 -*-
"""

Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 10: Problem 03
11/04/2025

"""

import numpy as np
from scipy.optimize import fsolve

#define reaction equilibrium constant
Kp = 6.23*10**(-3)

#define pressure
P = 30              #bar

#define initial guess for fraction dissociated
x_guess = 0.25

#calculate fraction dissociated
def x_calc(x):
    return (((4*x)*(1-x))/((1-(2*x))**2)) - (Kp*P**2)

x = fsolve(x_calc, x_guess)[0]

print(f'Equilibirum Fraction Dissociated, x = {x:.5f}')

#calculate equilibrium mole fractions
x_CO = 1/2
print(f'x_CO = {x_CO:.5f}')

x_H2 = (1-(2*x))/(2*(1-x))
print(f'x_H2 = {x_H2:.5f}')

x_CH3OH = x/(2*(1-x))
print(f'x_CH3OH = {x_CH3OH:.5f}')
