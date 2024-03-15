
# modules import
import numpy as np 
#import RechercheRacine
import Variables as var
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# odefunction
def dcidz(r_i, r_cbn):
    return (var.eta*(1-var.epsilon)*var.rho_cat*r_i - (1-var.epsilon)*var.rho_CaO*r_cbn)/var.u_g

def odefunction(z,C):
    # définition variables dépendantes des arguments de la fonction
    X = C[5]              # (dimensionless) conversion fractionnaire du CaO en CaCO3
    T = C[6]
    P = C[7]

    # réevaluation des pression partielles (P * C/C_tot avec PV=NRT et P_part = P_tot*n/n_tot, sachant C = n/V où V s'annule)
    C_tot = np.sum(C[0:5])
    C_CH4 = C[0]
    P_CH4 = P*(C_CH4/C_tot)
    C_H2O = C[1]
    P_H2O = P*(C_H2O/C_tot)
    C_H2 = C[2]
    P_H2 = P*(C_H2/C_tot)
    C_CO = C[3]
    P_CO = P*(C_CO/C_tot)
    C_CO2 = C[4]
    P_CO2 = P*(C_CO2/C_tot)
    
    K_1 = 4.707*(10**12)*np.exp(-224000/(var.R*T))
    K_3 = 1.142*(10**(-2))*np.exp(37300/(var.R*T))   # (dimensionless)
    K_2 = K_1 * K_3                         # (bar^2)

    k_1 = (1.842/3600)*10**(-4)*np.exp((-240100/var.R)*(1/T - 1/648)) # (kmol.(bar)^1/2 /(kg.s))  
    k_2 = (2.193/3600)*10**(-5)*np.exp((-243900/var.R)*(1/T - 1/648)) # (kmol.(bar)^1/2 /(kg.s))
    k_3 = (7.558/3600)*np.exp((-67130/var.R)*(1/T - 1/648))        # (kmol/(kg.s.bar))

    K_CH4 = 0.179*np.exp((38280/var.R)*(1/T - 1/823))   # (/bar)
    K_H2O = 0.4152*np.exp((-88680/var.R)*(1/T - 1/823)) # (dimensionless)
    K_H2 = 0.0296*np.exp((82900/var.R)*(1/T - 1/648))   # (/bar)
    K_CO = 40.91*np.exp((70650/var.R)*(1/T - 1/648))    # (/bar)

    DEN = 1 + K_CO*P_CO + K_H2*P_H2 + K_CH4*P_CH4 + K_H2O*P_H2O/P_H2

    R_1 = (k_1/((P_H2)**2.5))*(P_CH4*P_H2O - (P_H2**3*P_CO)/K_1)/DEN**2       # vitesse de réaction reformage (méthane + eau -> monoxyde de carbone + dihydrogène)    
    R_2 = (k_2/((P_H2)**3.5))*(P_CH4*(P_H2O**2) - (P_H2**4*P_CO2)/K_2)/DEN**2 # méthane + eau -> dioxyde de carbone + dihydrogène
    R_3 = (k_3/(P_H2))*(P_CO*P_H2O - (P_H2*P_CO2)/K_3)/DEN**2                 # monoxyde de carbone + eau -> dioxyde de carbone pour dihydrogène

    r_CH4 = -R_1 - R_2          # (kmol/(kg.s))
    r_H2O = -R_1 - 2*R_2 - R_3  # (kmol/(kg.s))
    r_H2 = 3*R_1 + 4*R_2 + R_3  # (kmol/(kg.s))
    r_CO = R_1 - R_3            # (kmol/(kg.s))
    r_CO2 = R_2 + R_3           # (kmol/(kg.s))

    k_c = var.M_k*np.exp(var.N_k/T)             # (/s) vitesse apparente de carbonatation
    b = var.M_b*np.exp(var.N_b/T)               # (/s) vitesse apparente de carbonatation
    X_u = k_c*b                             # (dimensionless) conversion ultime
    r_cbn = (k_c/var.M_CaO)*(1-X/X_u)**2    # (kmol/(kg.s)) taux de consommation par carbonatation de CO2 (conversion fractionnaire CaO en CaCO3)   
    
    # calcul du bilan énergétique
    rho_g = (1/(var.R*T))*(var.M_CH4*P_CH4 + var.M_CO*P_CO + var.M_CO2*P_CO2 + var.M_H2*P_H2 + var.M_H2O*P_H2O)*100 # (kg/m^3) masse volumique de la phase gazeuse, *100 pour la conversion d'unités    
    Re_p = (var.u_g*var.epsilon*rho_g*var.d_p)/var.mu
    
    if Re_p > 20 and var.d_p/var.D_r < 0.3:
        h_w = 2.03*(k_g/var.D_r)*(Re_p**(0.8))*np.exp((-6*var.d_p)/(var.D_r))
        
    elif Re_p < 20:
        h_w = 6.15*(var.k_z0/var.D_r)
        
    else:
        h_w = "error"
    
    # calcul des dérivées
    dC_dz = np.zeros(8)
    r = [r_CH4, r_H2O, r_H2, r_CO, r_CO2]   # uniques paramètres variant d'une équation à l'autre
    for i in range(5):
        dC_dz[i] = dcidz(r[i], r_cbn)
        
    # X
    dX_dz = (var.M_CaO*r_cbn)/var.u_s
    dC_dz[5] = dX_dz
    
    # T
    dT_dz = (-(1 - var.epsilon) * var.rho_cat * var.eta * (R_1*var.H_R1 + R_2*var.H_R2 + R_3*var.H_R3) - (1 - var.epsilon) * var.rho_CaO*r_cbn*var.H_cbn + h_w * (var.T_W - T) * 4 / var.D_r) / ((1 - var.epsilon)*var.rho_s*var.u_s*var.C_ps + rho_g*var.u_g*var.C_pg)
    dC_dz[6] = dT_dz
    
    # P
    dP_dz = (-rho_g*(var.u_g**2)/var.d_p)*((1- var.epsilon)/var.epsilon)*((150*(1 - var.epsilon)*var.mu)/(var.d_p*rho_g*var.u_g) + 1.75)*(10**(-5))
    dC_dz[7] = dP_dz
    
    return dC_dz

def calculConcentrationsEuler(Z, C_0):
    z_0, z_f = Z
    z_t = z_0 # z at time t
    c_t = C_0
    n_euler = 500000 # nombre de points
    
    # Hyperparameters variables
    z = np.linspace(z_0, z_f, n_euler)    # coordinates between z_0 and z_f with n values (datapoints)
    #h = 0.001                          # step size (pas)
    h = (z_f-z_0)/n_euler                 # cannot define interval/step size/number of steps simultaneously, one var must be the consequence of the others
    print(h)
    #print("Step size:", h)
    C = np.zeros([8, n_euler])            # tableau de taille 8 x n, n étant le nombre d'éléments recherchés
    
    # Single loop
    for i in range(n_euler):
        z_t += h
        derivee = odefunction(z_t, c_t)
        c_tt = c_t + h*derivee  # alternative to "c++"" (next c_t)
        C[:,i] = c_tt
        c_t = c_tt
        
    return [z, C]

def calculConcentrationsIVP(z, C_0): # solve Initial Problem Value (IVP)
    z_0, z_f = z
    z = np.linspace(z_0, z_f, var.n)

    sol = solve_ivp(lambda t, C: odefunction(t, C), [z_0, z_f], C_0, t_eval=z, rtol=10e-5) # time span and initial concentrations (tried 0.0002 but not sufficient because had "0.00023864609591485773 not less than 0.0002")
    z = sol.t
    C = sol.y
    
    return z, C

##################################################
#                   EVALUATION                   #
##################################################

p_e = 1/22.4*var.R*973.15/100000 # pression au début du réacteur
P_CH4 = (1/22.4)/(4)
P_H20 = (1/22.4)*(3/4)
C_0 = [P_CH4, P_H20, 10e-3, 0.00, 0.00, p_e, 973.15, 3]
z = [0, var.h_r] # h_r étant longueur du réacteur, le calcul se fait de 0 à 0.29 (m)

# Euler
z, C = calculConcentrationsEuler(z, C_0)

# IVP
z = [0, var.h_r] # longueur du réacteur
z, C = calculConcentrationsIVP(z, C_0)

fig, (ax1, ax2) = plt.subplots(2)

ax1.set_title('Concentrations')
ax1.grid(True)
ax1.set_xlabel('Distance (z)')
ax1.set_ylabel('Concentration')

labels = ['CH4', 'H2O', 'CH2', 'CO', 'CO2']
for i in range(5):
    ax1.plot(z, C[i], label=f'{labels[i]}')

ax1.legend(loc="upper right")

# X, T, P plot
ax2.set_title('X, T, P')

labels = ['X', 'T', 'P']
for i in range(5, 8):
    ax2.plot(z, C[i], label=f'{labels[i-5]}')

ax2.legend(loc="upper right")

plt.tight_layout()
plt.grid(True)
plt.legend()
plt.show()



