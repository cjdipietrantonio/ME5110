# -*- coding: utf-8 -*-
"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 07: Problem 05
10/14/2025

"""

import numpy as np
from scipy.optimize import fsolve

#define constants
T = 373.15      #K
v_mix = 0.000111    #m3/mol
R = 8.3143      #J/mol-K

#define mole fractions
x_CO2 = 0.695
x_C2H4 = 0.305

#van der Waals parameters
a_CO2 = 0.3653          #Pa-m6/mol2
b_CO2 = 4.2798*10**-5   #m3/mol
a_C2H4 = 0.4566         #Pa-m6/mol2
b_C2H4 = 5.7491*10**-5  #m3/mol

#solve vdw for molar volume of gas
def v_vdw(P, a, b):
    f = lambda v: P - ((R*T)/(v-b) - (a/v**2))
    return fsolve(f, R*T/P)[0]

#find difference between given molar volume and molar volume found from additive volume method
def volume_err(P):
    v_CO2 = v_vdw(P, a_CO2, b_CO2)
    v_C2H4 = v_vdw(P, a_C2H4, b_C2H4)
    vmix_calc = x_CO2*v_CO2 + x_C2H4*v_C2H4
    return vmix_calc - v_mix

#initial guess for P
P_guess = 1.739*10**7

#iteratively solve for mixture pressure
P_mix = fsolve(volume_err, P_guess)[0]

#convert to bar
P_mix_bar = P_mix/(10**5)

#find final component molar volumes
v_C02_final = v_vdw(P_mix, a_CO2, b_CO2)
v_C2H4_final = v_vdw(P_mix, a_C2H4, b_C2H4)


print('Mixture Pressure = ', P_mix_bar, 'bar')
print('CO2 Molar Volume', v_C02_final, 'm3/mol')
print('C2H4 Molar Volume', v_C2H4_final, 'm3/mol')

