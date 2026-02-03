# -*- coding: utf-8 -*-

"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 04
09/16/2025

"""

import matplotlib.pyplot as plt
import numpy as np

# -------------------------
#PROBLEM 03b
# -------------------------

#Define Variables
De = 7.31*10**-19                   #Units: J
a = 1.81*10**10                     #Units: 1/m
Re = 1.275*10**-10                  #Equilibrium bond length (Units: m)
k = 2*De*a**2                       #Force constant (Units: N/m)


#Generate Actual Bond Lengths
R = np.linspace(0.5*Re, 1.5*Re, 10000)                  #Units: m
x = R - Re

#Define Morse Potential
U_morse = De*(np.exp(-2*a*x) - 2*np.exp(-a*x) +1)       #Units: J


#Define Corresponding Harmonic Oscillator Potential
U_harmonic = 0.5*k*x**2                                 #Units: J


#Conversions
R_pm = R * 10**12                                       #Actual Bond Lengths (Units: pm)


#Create Plot
plt.figure("Morse Potential vs Harmonic Oscillator Potential for HCl", figsize=(8,5))
plt.plot(R_pm, U_morse, label='Morse Potential', color ='blue')
plt.plot(R_pm, U_harmonic, linestyle=':', label='Harmonic Oscillator Potential', color='Red')
plt.xlabel('Bond Length, R (pm)')
plt.ylabel('Potential Energy, U (J)')
plt.title('Morse Potential vs Harmonic Oscillator Potential for HCl')
plt.legend(loc='upper right', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.show()
