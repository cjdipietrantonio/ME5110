# -*- coding: utf-8 -*-
"""

Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 10: Problem 06
11/04/2025

"""

import numpy as np
from scipy.optimize import fsolve

#define reaction equilibrium constant
Kp = 10 **(0.158)

#define pressure
P = 100             #bar

#define equation for extent of reaction
def xi_calc(xi):
    num = (xi*((1.5 - (.5*xi))**.5))
    den = ((P**.5)*(1-xi)*((.5 - (.5*xi))**.5))
    return (num/den) - Kp

#define initial guess for xi
xi_guess = .8293799

#solve for xi
xi_sol = fsolve(xi_calc, xi_guess)[0]

print(f'Xi = {xi_sol}')

