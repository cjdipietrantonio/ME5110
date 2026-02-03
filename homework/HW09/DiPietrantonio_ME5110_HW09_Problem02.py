# -*- coding: utf-8 -*-

"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 09: Problem 02
10/28/2025

"""

import numpy as np
from scipy.optimize import fsolve

#define constants
h = 6.626*10**-34                   #J-S
k = 1.381*10**-23                   #J/K
R = 1.275*10**-10                   #m
m_1H_amu = 1.00782503223            #amu
m_35Cl_amu = 34.968852682           #amu
m_u = 1.6605402*10**-27             #kg/amu

#calculate reduced planck constant
h_bar = h/(2*np.pi)

#calculate reduced mass
m_1H = m_1H_amu*m_u
m_35Cl = m_35Cl_amu*m_u
mu = (m_1H*m_35Cl)/(m_1H + m_35Cl)

#calculate momentum of inertia
I = mu*R**2

#define characteristic temperature for rotation
theta_rot = (h_bar**2)/(2*I*k)

#define initial guess for T
T_guess = theta_rot*4.5


#define max rotational quantum number
J_max = 10


#define toelrance
tolerance = 1*10**-6

#define variale for old temp values
T_old = 0

#define number of iterations
it = 0

#add rotational energy levels until T converges
while True:
    
    #define rotational quantum numbers
    J = np.arange(0,J_max)

    #define function for rotational distribution
    def T_calc(T):
        g_rot = 2*J + 1
        E_J = ((h_bar**2)/(2*I))*(J*(J+1))
        N_i = g_rot*np.exp(-E_J/(k*T))
        Q = np.sum(N_i)
        return (np.sum(N_i[J<=3])/Q) - 0.5

    #iteratively solve for temperature
    T_sol = fsolve(T_calc, T_guess)[0]
    print(f'J_max = {J_max},T = {T_sol} K')
    
    #check convergance
    if abs(T_sol - T_old) < tolerance:
        break
    
    #update T_old, T_guess, and J_max for next iteration
    T_old = T_sol
    T_guess = T_sol
    J_max += 1 
    it += 1
    
print(f'Final Temperature converged to {T_sol} K')
print(f'{J_max} Rotational Energy Levels Required')
print(f'{it} iterations required')
print(f'final tolerance of: {abs(T_sol - T_old)}')


    
    
    