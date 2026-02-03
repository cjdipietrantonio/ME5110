# -*- coding: utf-8 -*-
"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 09: Problem 07
10/28/2025

"""

import numpy as np

#define constants
T = 273                             #K
theta_vib = 308.8                   #K
theta_rot = 0.054                   #K
R = 8.3143                          #J/mol-K

#define e_tr
e_tr = (3/2)*(R*T)
print('e_tr = ', e_tr)

#define Cp_tr
Cp_tr = (5/2)*R
print('Cp_tr = ', Cp_tr)

#define e_rot
e_rot = R*T
print('e_rot = ', e_rot)

#define Cv_rot
Cv_rot = R
print('Cv_rot = ', Cv_rot)

#define Cp_rot
Cp_rot = Cv_rot
print('Cp_rot = ', Cp_rot)

#define Q_vib
Q_vib = np.exp(-theta_vib/(2*T))/(1 - np.exp(-theta_vib/T))
print('Q_vib = ', Q_vib)

#define derivative of Q_vib
dQ_vib = np.exp(-theta_vib/(2*T))*(theta_vib/(T**2))*((1/(2*(1 - np.exp(-theta_vib/T))))
                                                      + (np.exp(-theta_vib/T)/((1 - np.exp(-theta_vib/T))**2)))
print('dQ_vib = ', dQ_vib)

#define e_vib
e_vib = (R*(T**2))*(1/Q_vib)*dQ_vib
print('e_vib = ', e_vib)

#define Cv_vib
Cv_vib = ((R*(theta_vib**2))*np.exp(-theta_vib/T))/((T**2)*((1 - np.exp(-theta_vib/T))**2))
print('Cv_vib = ', Cv_vib)

#define Cp_vib
Cp_vib = Cv_vib
print('Cp_vib = ', Cp_vib)

#calculate Cp
Cp = Cp_tr + Cp_rot + Cp_vib
print('Cp = ', Cp)