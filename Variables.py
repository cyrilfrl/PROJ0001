
import numpy as np

# constantes recherche solution équation diff
n = 1500

# constantes d'état du réacteur
h_r = 0.29                     # hauteur du réacteur
D_r = 0.024                    # diamètre du réacteur
V = ((np.pi*(D_r**2))/4)*h_r   # volume du réacteur
V_input = ((np.pi*(D_r**2))/4) # volume d'entrée des gaz dans le réacteur
u_s = 10**(-3)      # (m/s) vitesse d'entrée du CaO
u_g = 1             # (m/s) vitesse d'entrée des gaz
R = 8.314           # PV = NRT
T = 700 + 273.15    # (kelvin)
eta = 0.3           # efficacité du catalyseur
epsilon = 0.5       # porosité du réacteur
rho_cat = 1100      # (kg/m^3) masse volumique du catalyseur
rho_CaO = 1620      # (kg/m^3) masse volumique du CaO
W_CaO = 83.6*10**(-3) # (kg/h)
W_cat = 16.4*10**(-3) # (kg/h)
rho_s = (W_cat + W_CaO)/(W_cat/rho_cat + W_CaO/rho_CaO)
mu = 2.8*10**(-3)   # (Bar.s) (sans le 10^-5 donné en Pa) viscosité

# constantes des substances
M_CaO = 56
M_CH4 = 16
M_H2O = 18
M_H2 = 2
M_CO = 28
M_CO2 = 44
M_k = 303
N_k = -13146
M_b = 1.6
N_b = 5649
C_ps = 0.98 # (kJ/(kmol.K)) capacité thermique solide
C_pg = 8.45 # (kJ/(kmol.K)) capacité thermique gaz
k_g = 2.59*10**(-4) # conductivité thermique du gaz
k_s = 10**(-3) # conductivité thermique du solide
k_z0 = k_g*(epsilon + (1-epsilon)/(0.139*epsilon-0.0339 + (2/3)*(k_g/k_s)))
T_W = 700 + 273.15

H_cbn = -178.8*10**(3) #(kJ/kmol) soit (J/mol) CaO + CO2 ⇔ CaCO3, concerne la réaction de carbonatation (capture)
H_R1 = 206*10**(3)     # (kJ/kmol) CH4 + H2O ⇔ CO + 3H2
H_R2 = 164.9*10**(3)   # (kJ/kmol) CH4 + 2H2O ⇔ CO2 + 4H2
H_R3 = -41.1*10**(3)   # (kJ/kmol) CO + H2O ⇔ CO2 + H2

# constantes des pelets
d_p = 3*10**(-3)   # (m) diamètre pellets

#//# useless?
P = 3                           # (bar)
F_tot = 1/22.4                  # (mol/L)    
P_H2 = 10**(-15)                # () pression partielle H2
