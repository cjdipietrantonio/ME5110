# -*- coding: utf-8 -*-

"""
Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Homework 09: Problem 06
10/28/2025

"""

import numpy as np

#define constants
h = 6.626*10**-34                   #J-S
P = 1.01325*10**5                   #Pa
T = 30000                           #K
k = 1.381*10**-23                   #J/K
m_1H_amu = 1.00782503223            #amu
m_u = 1.6605402*10**-27             #kg/amu
m_1H = m_1H_amu*m_u                 #kg
R = 8.3143                          #J/mol-K

#define energy levels
e_1 = 1.6338477*10**(-18)              #J
e_2 = 1.6338484*10**(-18)              #J
e_3 = 1.6338550*10**(-18)              #J

#define degeneracies
g_0 = 2
g_1 = 2
g_2 = 2
g_3 = 4

#define e_tr
e_tr = (3/2)*R*T

#define Cv_tr
Cv_tr = (3/2)*R

#calculate log term for translational molar specific entropy equation
ln_qtr_n = np.log((((2*np.pi*m_1H)/(h**2))**(3/2))*(((k*T)**(5/2))/P))
print('ln(Qtr/N) = ', ln_qtr_n)

#calculate translational molar specific entropy
s_tr = (R*(ln_qtr_n + 1)) + e_tr/T
print('s_tr = ', s_tr)

#calculate electronic partition function
Q_el = g_0 + (g_1)*np.exp((-e_1)/(k*T)) + (g_2)*np.exp((-e_2)/(k*T)) + (g_3)*np.exp((-e_3)/(k*T))
print('Q_el = ', Q_el)

#calculate first derivative of electronic partition function
dQ_el = (g_1)*np.exp((-e_1)/(k*T))*(e_1/(k*(T**2))) + (g_2)*np.exp((-e_2)/
        (k*T))*(e_2/(k*(T**2))) + (g_3)*np.exp((-e_3)/(k*T))*(e_3/(k*(T**2)))
print('dQ_el =' , dQ_el)

#calculate electronic molar specific energy
e_el = ((R*(T**2))/Q_el)*dQ_el
print('e_el = ', e_el)

#calculate second derivative of electronic partition function
d2Q_el = (g_1)*np.exp((-e_1)/(k*T))*((e_1**2)/((k**2)*(T**4)) + ((-2*e_1)/(k*(T**3)))) + (g_2)*np.exp((-e_2)/
         (k*T))*((e_2**2)/((k**2)*(T**4)) + ((-2*e_2)/(k*(T**3)))) + (g_3)*np.exp((-e_3)/(k*T))*((e_3**2)/
         ((k**2)*(T**4)) + ((-2*e_3)/(k*(T**3))))
print('d2Q_el = ', d2Q_el)

#calculate first derivative of the natural log of the electronic partition function
dlnQ_el = (1/Q_el)*dQ_el
print('dlnQ_el = ', dlnQ_el)

#calculate second derivative of the natural log of the electronic partition function
d2lnQ_el = ((Q_el*d2Q_el) - (dQ_el)**2)/(Q_el**2)
print('d2lnQ_el = ', d2lnQ_el)

#calculate the electronic molar specific volume constant specific heat
Cv_el = (R*(T**2))*d2lnQ_el + (2*R*T)*dlnQ_el
print('Cv_el = ', Cv_el)

#calculate the electronic molar specific entropy
s_el = (R*np.log(Q_el)) + (e_el/T)
print('s_el = ', s_el)

#calculate total molar specific internal energy
e = e_tr + e_el
print('e = ', e)

#calculate total molar specific volume constant specfic heat
Cv = Cv_tr + Cv_el
print('Cv = ', Cv)

#calculate total molar specific entropy
s = s_tr + s_el
print('s = ', s)




