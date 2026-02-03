# -*- coding: utf-8 -*-

"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 07: Problem 01
10/14/2025

"""

import numpy as np
import sympy as sp

#SOLVING NUMERICALLY
#define cubic polynomial coefficeints
coefficients = [1, -3, -3, -1]

#solve numerically for for roots of cubic
num_roots = np.roots(coefficients)

#define empty list of real numerical roots
real_roots_num = []

#cince x = vc/B, x must be real and positive, filter list of roots
for r in num_roots:
    if np.isreal(r) and r.real > 0:
        real_roots_num.append(r.real)

print("All roots (numerical):", num_roots)
print("Real positive root (physical x):", real_roots_num)

#SOLVING ALGEBRAICALLY
#define variable x
x = sp.symbols('x', real=True)

#define cubic polynomial
cubic = (x**3) - (3*x**2) - (3*x) - 1

#solve for algebraic roots of cubic
alg_roots = sp.solve(cubic, x)

#define empty list of real algebraic roots
real_roots_alg = []

#filter for real algebraic roots
for r in alg_roots:
        if r.is_real and r.evalf() > 0:
            real_roots_alg.append(r)

print("All roots (algebraic):")
for r in alg_roots:
    print(r)
print("Real positive root (algebraic):", real_roots_alg)

#SOLVING Zc
#define constants
A = 1
B = 1
R = 1

#define critical properties
Pc = (0.1213*(A**(2/3))*(B**(-5/3))*(R**(1/3))) - (0.09129*(A**(1/2))*(B**(-3/2))*(R**(1/2)))
Tc = 0.345*(A/(B*R))**(2/3) 
Vc = 3.8473*B

Zc = (Pc*Vc)/(R*Tc)

print("Critical compressibility factor Zc =", Zc)



