# -*- coding: utf-8 -*-
"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 07: Problem 03
10/14/2025

"""

import numpy as np
from scipy.optimize import fsolve

#define constants
T1 = 301.3          #K
v1 = 0.05           #m3/mol
v2 = 0.001          #m3/mol
b = 6.4754*10**-5   #m3/mol
n = 2               #mol
a = 0.5547          #Pa-m6/mol2
R = 8.3143          #J/mol-K

#define the function
def f(T2):
    return 1.131*np.log(T2/T1) + 0.019225*(T2 - T1) - (2.7805*10**-6)*((T2**2) - (T1**2)) + np.log((v2 - b)/(v1 - b))

#define initial guess
T2_guess = 631      #K

#numerically solve for T2
T2_sol = fsolve(f, T2_guess)


print('T2 =', T2_sol[0], 'K')


#calculating change in internal energy
du_1 = R*(0.131*(T2_sol[0] - T1) + ((1.9225*10**-2)/2)*((T2_sol[0]**2) - (T1**2)) - ((5.561*10**-6)/3)*((T2_sol[0]**3) - (T1**3)))

du_2 = a*((1/v1)-(1/v2))

du = du_1 + du_2

U = n*du
print('Delta U = ', U, 'J')

