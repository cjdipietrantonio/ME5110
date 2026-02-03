# -*- coding: utf-8 -*-
"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 07: Problem 04
10/14/2025

"""
#define constants
T1 = 600            #K
P1 = 1*10**7        #Pa
T2 = 445            #K
P2 = 2*10**6        #Pa
a = 1*10**-8        #m3-K/Pa-mol
b = 8*10**-5        #m3/mol

#calculate change in molar enthalpy
h12 = 30*(T2-T1) + 0.01*((T2**2) - (T1**2)) + b*(P2 - P1) + ((a/T2)*((P2**2) - (P1**2)))

print('h12 = ', h12, 'J/mol')
