
# modules import
import numpy as np 
import RechercheRacine
import SimReacteur
import test

###########################################
#  définition des constantes et formules  #
###########################################
u_g = 1     # vitesse d'entrée du gaz
R = 8.314
T = 700 + 273.15
eta = 0.3
epsilon = 0.5
rho_cat = 1100
rho_CaO = 1620

M_CaO = 56
M_CH4 = 16
M_H20 = 18
M_H2 = 2
M_CO = 28
M_CO2 = 44

# constantes de pression
h = 0.29 # hauteur du réacteur
D = 0.024 # diamètre du réacteur
V = ((np.pi*(D**2))/4)*h # volume du réacteur
V_input = ((np.pi*(D**2))/4) # volume d'entrée des gaz dans le réacteur

P = 3 # (bar)
F_tot = 1/22.4 # (mol/L)
P_H2 = 10e-11

K_1 = 4.707 * 10e12*np.exp(-224000/R*T)
K_3 = 1.142 * 10e-2*np.exp(37300/R*T)
K_2 = K_1 * K_3

k_1 = (1.842/3600)*10e-4*np.exp((-240100/R)*(1/T - 1/648))
k_2 = (2.193/3600)*10e-5*np.exp((-243900/R)*(1/T - 1/648))
k_3 = (7.558/3600)*np.exp((-67130/R)*(1/T - 1/648))

K_CH4 = 0.179*np.exp(38280/R)*(1/T - 1/823)
K_H20 = 0.4152*np.exp(-88680/R)*(1/T - 1/823)
K_H2 = 0.0296*np.exp(82900/R)*(1/T - 1/648)
K_CO = 40.91*np.exp(70650/R)*(1/T - 1/648)

# bar et kelvin

# odefunction

def dcidz(r_i, r_cbn):
    return (eta(1-epsilon)*rho_cat*r_i - (1-epsilon)*rho_CaO*r_cbn)/u_g


def odefunction(z,C):
    # définition variables dépendantes des arguments de la fonction
    C_CH4 = C[0]
    P_CH4 = C_CH4*V # (mol)
    C_H20 = C[1]
    P_H20 = C_H20*V # (mol)
    C_CO = C[3]
    P_CO = C_CO*V
    C_CO2 = C[4]
    P_CO2 = C_CO2*V
    X = C[5]
    T = C[6]
    P = C[7]

    DEN = 1 + K_CO*P_CO + K_H2*P_H2 + K_CH4*P_CH4 + K_H20*P_H20/P_H2
    R_1 = (k_1/(P_H2**2.5))*(P_CH4*P_H20 - (P_H2**3*P_CO)/K_1)/DEN**2
    R_2 = (k_2/(P_H2**3.5))*(P_CH4*(P_H20**2) - (P_H2**4*P_CO2)/K_2)/DEN**2
    R_3 = (k_3/P_H2)*(P_CO*P_H20 - (P_H2*P_CO2)/K_3)/DEN**2

    r_CH4 = -R_1 -R_2
    r_H20 = -R_1 - 2*R_2 -R_3
    r_H2 = 3*R_1 + 4*R_2 + R_3
    r_CO = R_1 - R_3
    r_CO2 = R_2 + R_3

    M_k = 303
    N_k = -13146
    k_c = M_k*np.exp(N_k/T)
    M_b = 1.6
    N_b = 5649
    b = M_b*np.exp(N_b/T)
    X_u = k_c*b
    r_cbn = (k_c/M_CaO)*(1-X/X_u)**2

    r = [r_CH4, r_H20, r_H2, r_CO, r_CO2]
    dC_dz = [dcidz(r_i, r_cbn) for r_i in r]

    return dC_dz

z, C = calculConcentrationsEuler([z0, z1], C0)

