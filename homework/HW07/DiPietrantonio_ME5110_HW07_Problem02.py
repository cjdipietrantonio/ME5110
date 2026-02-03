# -*- coding: utf-8 -*-

"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 07: Problem 02
10/14/2025

"""

import numpy as np
import matplotlib.pyplot as plt


#define constants and critical properties
R = 8.3143           #J/mol-K
Tc = 126.2           #K
Pc = 4.4 * 10**6     #Pa

a = (27*(R**2)*(Tc**2))/(64*Pc)
b = (R*Tc)/(8*Pc)

#define range of molar volumes
v = np.linspace(1.001*b, 100*b, 5000)         #m3/mol

#inversion tempertature
T_inv = ((2*a)/(v**2))/  (((v*R)/((v-b)**2))-(R/(v-b)))

#pressure from van der Waals
P = ((R*T_inv)/(v-b)) - (a/(v**2)) #Pa
P_bar = P/(10**5)                  #Bar

#plot inversion line
plt.figure(figsize=(8,5))
plt.plot(P_bar, T_inv, lw=2)
plt.xlabel('Pressure (Bar)')
plt.ylabel('Inversion Temperature (K)')
plt.title('Inversion Line for N2 (van der Waals EOS)')
plt.xlim(left=0) #start x-axis at 0 to avoid plotting points at negative pressure
plt.show()

