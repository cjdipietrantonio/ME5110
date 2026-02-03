# -*- coding: utf-8 -*-
"""

Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Final Exam, Project
11/30/2025

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#define JANNAF table temperatures and equilibrium constants
T = [500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]       #K
log_kf = [-0.501, -1.382, -2.029, -2.524, -2.916, -3.234, -3.496, -3.716, -3.903, -4.064, -4.203]


#define pressures
P = [1, 15, 30, 75, 150, 300]       #bar


#create dataframe for temperature and equilibrium constants
df = pd.DataFrame({
    'T': T,
    'logKf': log_kf
})


#add Kf and Kp to dataframe
df['Kf'] = 10**df['logKf']
df['Kp'] = 1/((df['Kf'])**2)


#define equation for the equilibium ammonia conversion factor
def xi(Kp, P):
    A = (4/P)*np.sqrt(Kp/27)
    return np.sqrt(A/(1+A))


#calculate xi as a function of temperature for each pressure
for p in P:
    df[f'xi_{p}bar'] = xi(df['Kp'], p)
    
    
#plot xi against temperature at each pressure
markers = ['o', 's', '^', 'd', 'v', 'x']

plt.figure(figsize=(12,7))

for p, m in zip(P, markers):
    plt.plot(df['T'], df[f'xi_{p}bar'], marker=m, markersize=8, label=f'{p} bar', markevery=[1, 2, 3, 4, 5])
    
plt.xlabel('Temperature (K)', fontsize=14)
plt.ylabel(r'Equilibrium Ammonia Conversion Factor, $\xi$', fontsize=14)
plt.title('NH3 Equilibrium Conversion vs. Temperature', fontsize=16, fontweight='bold')
plt.xlim(500, 1500)
plt.xticks([500, 750, 1000, 1250, 1500], fontsize = 14)
plt.yticks(fontsize=14)
plt.legend(title= r'$\mathbf{Pressure}$', fontsize=14, title_fontsize=14)
plt.show()






