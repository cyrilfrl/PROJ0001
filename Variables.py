
import numpy as np

u_g = 1             # (m/s) vitesse superficielle d'entrée du gaz
R = 8.314           # PV = NRT
T = 700 + 273.15    # (kelvin)
eta = 0.3           # efficacité du catalyseur
epsilon = 0.5       # porosité du réacteur
rho_cat = 1100      # (kg/m^3) masse volumique du catalyseur
rho_CaO = 1620      # (kg/m^3) masse volumique du CaO

M_CaO = 56
M_CH4 = 16
M_H20 = 18
M_H2 = 2
M_CO = 28
M_CO2 = 44

# constantes d'état du réacteur
h = 0.29                     # hauteur du réacteur
D_r = 0.024                    # diamètre du réacteur
V = ((np.pi*(D_r**2))/4)*h     # volume du réacteur
V_input = ((np.pi*(D_r**2))/4) # volume d'entrée des gaz dans le réacteur

# constantes des pelets
d_p = 3*10e-3   # (m) diamètre pellets

P = 3                           # (bar)
F_tot = 1/22.4                  # (mol/L)
P_H2 = 10e-15                   # () pression partiell H2

K_1 = 4.707*10e12*np.exp(-224000/R*T) # (bar^2)
K_3 = 1.142*10e-2*np.exp(37300/(R*T))   # (dimensionless)
K_2 = K_1 * K_3                         # (bar^2)

k_1 = (1.842/3600)*10e-4*np.exp((-240100/R)*(1/T - 1/648)) # (kmol.(bar)^1/2 /(kg.s))  
k_2 = (2.193/3600)*10e-5*np.exp((-243900/R)*(1/T - 1/648)) # (kmol.(bar)^1/2 /(kg.s))
k_3 = (7.558/3600)*np.exp((-67130/R)*(1/T - 1/648))        # (kmol/(kg.s.bar))

K_CH4 = 0.179*np.exp((38280/R)*(1/T - 1/823))   # (/bar)
K_H20 = 0.4152*np.exp((-88680/R)*(1/T - 1/823)) # (dimensionless)
K_H2 = 0.0296*np.exp((82900/R)*(1/T - 1/648))   # (/bar)
K_CO = 40.91*np.exp((70650/R)*(1/T - 1/648))    # (/bar)
