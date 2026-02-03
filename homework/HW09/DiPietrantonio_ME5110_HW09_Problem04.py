# -*- coding: utf-8 -*-
"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 09: Problem 04
10/28/2025

"""

import numpy as np
from scipy.optimize import fsolve

#define constants
h = 6.626*10**-34                   #J-S
k = 1.381*10**-23                   #J/K
v_cm = 2330                         #cm-1
c = 2.998 * 10**8                   #m/s
T = 2000                            #K

#calculate fundamental vibrational frequency in Hz
v_hz = v_cm*100*c

#calculate reduced planck constant
h_bar = h/(2*np.pi)

#define max vibrational quantum number
v_max = 5

#define toelrance
tolerance = 1*10**-6

#define variale for old temp values
p_old = 0

#define number of iterations
it = 0

#add vibrational energy levels until probability converges
while True:
    
    #define vibrational quantum numbers
    v = np.arange(0,v_max)

    #define function for vibrational distribution
    def p_calc(v):
        E_v = (v + .5)*h*v_hz
        N_i = np.exp(-E_v/(k*T))
        Q = np.sum(N_i)
        vib_dist = (np.sum(N_i[v<=2])/Q)
        return vib_dist

    #calculate probability
    p_sol = p_calc(v)
    print(f'v_max = {v_max}, p = {p_sol}')
    
    #check convergance
    if abs(p_sol - p_old) < tolerance:
        break
    
    #update p_old and v_max for next iteration
    p_old = p_sol
    v_max += 1 
    it += 1
    
print(f'Final probability converged to {p_sol} or {p_sol*100:.2f}%')
print(f'{v_max} Vibrational Energy Levels Required')
print(f'{it} iterations required')
print(f'final tolerance of: {abs(p_sol - p_old)}')