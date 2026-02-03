# -*- coding: utf-8 -*-
"""

Christian DiPietrantonio
ME 5110: Advanced Thermodynamics
Final Project: Problem 04
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
    

#define constraint points
xi_min = 0.8
Tmax = 800          #K
#Tmin_15 = 624.36   #K
Tmin_15 = np.interp(xi_min, df['xi_15bar'], df['T'])    #K
#Tmin_30 = 670.88   #K
Tmin_30 = np.interp(xi_min, df['xi_30bar'], df['T'])    #K

 
#plot xi against temperature at each pressure
markers = ['o', 's', '^', 'd', 'v', 'x']

plt.figure(figsize=(12,7))

for p, m in zip(P, markers):
    plt.plot(df['T'], df[f'xi_{p}bar'], marker=m, markersize=8, label=f'{p} bar', markevery=[1, 2, 3, 4, 5], zorder=1)


#initialize lists for shading coordinates
T_fill = []
xi_fill = []


#define coordinates for minimum temperatures on the p = 15 bar and 30 bar curves
T_fill.extend([Tmin_15, Tmin_30])
xi_fill.extend([xi_min, xi_min])


#defne coordinates for the p = 30 bar curve within the operating range 
T_fill.extend(df['T'][(df['T'] >= Tmin_30) & (df['T'] <= Tmax)])
xi_fill.extend(df['xi_30bar'][(df['T'] >= Tmin_30) & (df['T'] <= Tmax)])


#define coordinates for Tmax on the p = 15 curve 
T_fill.append(Tmax)
xi_fill.append(df.loc[df['T'] == Tmax, 'xi_15bar'].values[0])


#defne coordinates for the p = 15 bar curve within the operating range and reverse their order
T_fill.extend(df['T'][(df['T'] >= Tmin_15) & (df['T'] <= Tmax)][::-1])
xi_fill.extend(df['xi_15bar'][(df['T'] >= Tmin_15) & (df['T'] <= Tmax)][::-1])


#fill the shape traced by the coordinates defined above
plt.fill(T_fill, xi_fill,
         facecolor='green', alpha=0.9,
         edgecolor='black', linestyle='--', linewidth=2,
         label='Operating Window', zorder=2)

plt.xlabel('Temperature (K)', fontsize=14)
plt.ylabel(r'Equilibrium Ammonia Conversion Factor, $\xi$', fontsize=14)
plt.title('NH3 Equilibrium Conversion vs. Temperature', fontsize=16, fontweight='bold')
plt.xlim(500, 1500)
plt.xticks([500, 750, 1000, 1250, 1500], fontsize = 14)
plt.yticks(fontsize=14)
plt.legend(title= r'$\mathbf{Pressure}$', fontsize=14, title_fontsize=14)
plt.show()

#make a close up plot of operating window
plt.figure(figsize=(12,7))


#plot p = 15 bar and p = 30 bar curves
for p, m in zip([15,30], ['s','^' ]):
    if p == 15:
        color = 'orange'
    elif p ==30:
        color = 'green'
    plt.plot(df['T'], df[f'xi_{p}bar'], marker=m, markersize=8, label=f'{p} bar', markevery=[1, 2, 3, 4, 5], color=color, zorder=1)


#fill the operating window
plt.fill(T_fill, xi_fill,
         facecolor='green', alpha=0.9,
         edgecolor='black', linestyle='--', linewidth=2,
         label='Operating Window', zorder=2)

plt.xlabel('Temperature (K)', fontsize=14)
plt.ylabel(r'Equilibrium Ammonia Conversion Factor, $\xi$', fontsize=14)
plt.title('NH3 Cracking Operating Pressure/Temperature Window', fontsize=16, fontweight='bold')
plt.xlim(Tmin_15 - 25, Tmax + 25)
plt.ylim(0.75, 1)
plt.xticks(fontsize = 14)
plt.yticks(fontsize=14)
plt.legend(title= r'$\mathbf{Pressure}$', fontsize=14, title_fontsize=14, loc='lower right')
plt.show()





