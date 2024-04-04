
# modules import
import numpy as np 
#import RechercheRacine
import Variables as var
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# odefunction
def dcidz(r_i, r_cbn, ug=var.u_g):
    return (var.eta*(1-var.epsilon)*var.rho_cat*r_i - (1-var.epsilon)*var.rho_CaO*r_cbn)/ug

def odefunction(z,C, CCUS=True, ug=var.u_g, us=var.u_s):
    
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
    K_2 = K_1 * K_3                                  # (bar^2)

    k_1 = (1.842/3600)*10**(-4)*np.exp((-240100/var.R)*(1/T - 1/648)) # (kmol.(bar)^1/2 /(kg.s))  
    k_2 = (2.193/3600)*10**(-5)*np.exp((-243900/var.R)*(1/T - 1/648)) # (kmol.(bar)^1/2 /(kg.s))
    k_3 = (7.558/3600)*np.exp((-67130/var.R)*(1/T - 1/648))           # (kmol/(kg.s.bar))

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


    if CCUS == True:
        u_s_adapted = us
        k_c = var.M_k*np.exp(var.N_k/T)             # (/s) vitesse apparente de carbonatation
        b = var.M_b*np.exp(var.N_b/T)               # (/s) vitesse apparente de carbonatation
        X_u = k_c*b                                 # (dimensionless) conversion ultime
        r_cbn = (k_c/var.M_CaO)*(1-X/X_u)**2        # (kmol/(kg.s)) taux de consommation par carbonatation de CO2 (conversion fractionnaire CaO en CaCO3)   
        
    elif CCUS == False:
        u_s_adapted = 0 # on ne fait pas avancer le pellet
        k_c = 0     # vitesse de carbonatation nulle puisque pas de carbonatation
        b = 0       # temps de conversion ultime nul puisque n'avance pas
        X_u = 0
        r_cbn = 0

    # calcul du bilan énergétique
    rho_g = (1/(var.R*T))*(var.M_CH4*P_CH4 + var.M_CO*P_CO + var.M_CO2*P_CO2 + var.M_H2*P_H2 + var.M_H2O*P_H2O)*100 # (kg/m^3) masse volumique de la phase gazeuse, *100 pour la conversion d'unités    
    Re_p = (ug*var.epsilon*rho_g*var.d_p)/var.mu
    
    if Re_p > 20 and var.d_p/var.D_r < 0.3:
        h_w = 2.03*(var.k_g/var.D_r)*(Re_p**(0.8))*np.exp((-6*var.d_p)/(var.D_r))
        
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
    if u_s_adapted != 0:
        dX_dz = (var.M_CaO*r_cbn)/u_s_adapted
    elif u_s_adapted == 0:
        dX_dz = 0
    dC_dz[5] = dX_dz
    
    # T
    if CCUS == True:
        dT_dz = (-(1 - var.epsilon) * var.rho_cat * var.eta * (R_1*var.H_R1 + R_2*var.H_R2 + R_3*var.H_R3) - (1 - var.epsilon) * var.rho_CaO*r_cbn*var.H_cbn + h_w * (var.T_W - T) * 4 / var.D_r) / ((1 - var.epsilon)*var.rho_s*u_s_adapted*var.C_ps + rho_g*ug*var.C_pg)
    
    elif CCUS == False:
        dT_dz = (-(1 - var.epsilon) * var.rho_cat * var.eta * (R_1*var.H_R1 + R_2*var.H_R2 + R_3*var.H_R3) + h_w * (var.T_W - T) * 4 / var.D_r) / ((1 - var.epsilon)*var.rho_s*u_s_adapted*var.C_ps + rho_g*ug*var.C_pg)
    
    dC_dz[6] = dT_dz
    
    # P
    dP_dz = (-rho_g*(ug**2)/var.d_p)*((1- var.epsilon)/var.epsilon)*((150*(1 - var.epsilon)*var.mu)/(var.d_p*rho_g*ug) + 1.75)*(10**(-5))
    dC_dz[7] = dP_dz
    
    return dC_dz

def calculConcentrationsEuler(Z, C_0, CCUS, ug=var.u_g, us=var.u_s):
    z_0, z_f = Z
    z_t = z_0 # z at time t
    c_t = C_0
    n_euler = 500000 # nombre de points
    
    # Hyperparameters variables
    z = np.linspace(z_0, z_f, n_euler)    # coordinates between z_0 and z_f with n values (datapoints)
    #h = 0.001                          # step size (pas)
    h = (z_f-z_0)/n_euler                 # cannot define interval/step size/number of steps simultaneously, one var must be the consequence of the others
    #print("Step size:", h)
    C = np.zeros([8, n_euler])            # tableau de taille 8 x n, n étant le nombre d'éléments recherchés, ici n_euler est défini
    
    # Single loop
    for i in range(n_euler):
        z_t += h
        derivee = odefunction(z_t, c_t, CCUS, ug, us)
        c_tt = c_t + h*derivee  # alternative to "c++"" (next c_t)
        C[:,i] = c_tt
        c_t = c_tt
        
    return [z, C]

def calculConcentrationsIVP(z, C_0, CCUS=True, ug=var.u_g, us=var.u_s): # solve Initial Problem Value (IVP)
    z_0, z_f = z
    z = np.linspace(z_0, z_f, var.n)

    sol = solve_ivp(lambda t, C: odefunction(t, C, CCUS, ug, us), [z_0, z_f], C_0, t_eval=z, rtol=10e-5) # time span and initial concentrations (tried 0.0002 but not sufficient because had "0.00023864609591485773 not less than 0.0002")
    z = sol.t
    C = sol.y
    
    return z, C

def completePlot(z1, z2, C1, C2, specific_label="_variante"):
    fig, (ax1, ax2) = plt.subplots(2)
    ax1.set_title('Concentrations')
    ax1.grid(True)
    ax1.set_xlabel('Distance (z)')
    ax1.set_ylabel('Concentration')

    labels = ['CH4', 'H2', 'CH2', 'CO', 'CO2']

    for i in range(5):
        if i != 1: # pour exclure C[1] avec H2O
            ax1.plot(z1, C1[i], label=f'{labels[i]}')

    for i in range(5):
        if i != 1: # pour exclure C[1] avec H2O
            ax1.plot(z2, C2[i], label=f'{labels[i]}_{specific_label}')

    ax1.legend(loc="upper right")

    ax2.set_title('T (Kelvins)')
    ax2.plot(z1, C1[6], label='T')
    ax2.plot(z2, C2[6], label=f'T_{specific_label}')
    ax2.legend(loc="upper right")

    plt.tight_layout()
    plt.grid(True)
    plt.legend()
    plt.show()

    

##################################################
#                   EVALUATION                   #
##################################################

p_e = 1/22.4*var.R*973.15/100000 # pression au début du réacteur
P_CH4 = (1/22.4)/(4)
P_H2O = (1/22.4)*(3/4)
C_0 = [P_CH4, P_H2O, 10e-3, 0.00, 0.00, p_e, 973.15, 3]
z = [0, var.h_r] # h_r étant longueur du réacteur, le calcul se fait de 0 à 0.29 (m)


##################################################
#                VARIANTES EVAL                  #
##################################################

# clés pour print les graphiques de concentrations spécifiques
labels_tokenized = {'CH4': 0, 'H2O': 1, 'H2': 2, 'CO': 3, 'CO2': 4, 'TEMP': 6, 'GENERAL':7}
labels = ['CH4', 'H2O', 'H2', 'CO', 'CO2','TEMP']

""" Q 3.1"""
"""
# Carbon Capture Enabled vs Disabled Comparison
z = [0, var.h_r] # longueur du réacteur
z1, C1 = calculConcentrationsIVP(z, C_0, CCUS=False)
z2, C2 = calculConcentrationsIVP(z, C_0, CCUS=True)

# plot des variations observables pour les différentes variations du modèle (procédé similaire pour toutes les variations étudiées)
while True:
    key = input("Quel variable désirez-vous plotter? (exit pour sortir): ")
    if key == "exit":
        break
    
    else:
        i = labels_tokenized[key]
        fig, (ax1, ax2) = plt.subplots(2)
        ax1.set_title(f'Concentrations {labels[i]}')
        ax1.grid(True)
        ax1.set_xlabel('Distance (z)')
        ax1.set_ylabel('Concentration')

        ax1.plot(z1, C1[i], label=f'{labels[i]}', color='red')
        ax1.plot(z2, C2[i], label=f'{labels[i]}_CCEN', color='green') # Carbon capture enabled
        ax1.legend(loc="upper right")
        
        ax2.set_title('T (Kelvins)')
        ax2.plot(z1, C1[6], label='T', color='red')
        ax2.plot(z2, C2[6], label='T', color='green')
        ax2.legend(loc="upper right")

        plt.tight_layout()
        plt.grid(True)
        plt.legend()
        plt.show()
"""
"""Q 3.2"""
if (choice := input("BCOMP/VCOMP?: ")) == 'BCOMP': # making use of the walrus operator to choose between 2 values comparison or vector value comparison (e.g. 1 to 10 with 0.1 as step size)

    # Carbon Capture Enabled vs Disabled Comparison
    z = [0, var.h_r] # longueur du réacteur


    us_valeurs = []
    us_valeurs.append(float(input("Quelle doit être la première valeur de u_s?: ")))
    us_valeurs.append(float(input("Quelle doit être la seconde valeur de u_s?: ")))

    #//# us_valeurs = [10**(-3), 10]

    z1, C1 = calculConcentrationsIVP(z, C_0, CCUS=True, us=us_valeurs[0])
    z2, C2 = calculConcentrationsIVP(z, C_0, CCUS=True, us=us_valeurs[1])

    # plot des variations observables pour les différentes variations du modèle (procédé similaire pour toutes les variations étudiées)
    while True:
        key = input("Quel variable désirez-vous plotter? (exit pour sortir): ")
        if key == "exit":
            break
        
        elif key == "GENERAL":
            completePlot(z1, z2, C1, C2, specific_label="Variations u_s")
        
        else:
            i = labels_tokenized[key]
            fig, (ax1, ax2) = plt.subplots(2)
            ax1.set_title('Concentrations')
            ax1.grid(True)
            ax1.set_xlabel('Distance (z)')
            ax1.set_ylabel('Concentration')

            ax1.plot(z1, C1[i], label=f'{labels[i]}', color='red')
            ax1.plot(z2, C2[i], label=f'{labels[i]}_CCEN', color='green') # Carbon capture enabled
            ax1.legend(loc="upper right")
            
            ax2.set_title('T (Kelvins)')
            ax2.plot(z1, C1[6], label='T', color='red')
            ax2.plot(z2, C2[6], label='T', color='green')
            ax2.legend(loc="upper right")

            plt.tight_layout()
            plt.grid(True)
            plt.legend()
            plt.show()

elif choice == 'VCOMP':
    initial_value = float(input("Initial value?: ")) # 0.06 minimal
    final_value = float(input("Final value?: "))
    nb_steps = int(input("How many steps in between?: "))
        
    z = [0, var.h_r] # longueur du réacteur
    UsSteps = np.linspace(int(initial_value), int(final_value), nb_steps)
    UsConcentrations = []

    for temp_u_s in UsSteps:
        zr, Cr = calculConcentrationsIVP(z, C_0, CCUS=True, us=temp_u_s) # zr and Cr stand respectively for z result and C result
        UsConcentrations.append(Cr)
    
    # Plot des différents graphiques
    fig, ax = plt.subplots(3, 2, figsize=(12, 8), constrained_layout=True)
    
    for i in range(6):
        concentrations = [k[i][-1] for k in UsConcentrations]

        ax[i%3, i%2].plot(UsSteps, concentrations, label=f'{labels[i]}')            
        ax[i%3, i%2].set_title(f'Concentrations de sortie en {labels[i]} selon us')
        ax[i%3, i%2].set_ylabel(f'Concentrations {labels[i]}')
        ax[i%3, i%2].set_xlabel('Valeur u_s')
        ax[i%3, i%2].grid(True)
        ax[i%3, i%2].legend(loc="upper right")

    plt.show()



"""Q 3.3"""
"""
if (choice := input("BCOMP/VCOMP?: ")) == 'BCOMP':
    z = [0, var.h_r] # longueur du réacteur

    ug_valeurs = []
    ug_valeurs.append(float(input("Quelle doit être la première valeur de u_g?: ")))
    ug_valeurs.append(float(input("Quelle doit être la seconde valeur de u_g?: ")))

    z1, C1 = calculConcentrationsIVP(z, C_0, CCUS=True, ug=ug_valeurs[0])
    z2, C2 = calculConcentrationsIVP(z, C_0, CCUS=True, ug=ug_valeurs[1])

    # plot des variations observables pour les différentes variations du modèle (procédé similaire pour toutes les variations étudiées)
    while True:
        key = input("Quel variable désirez-vous plotter? (exit pour sortir): ")
        if key == "exit":
            break
        
        elif key == "GENERAL":
            completePlot(z1, z2, C1, C2, specific_label="Variations u_g")
        
        else:
            i = labels_tokenized[key]
            fig, (ax1, ax2) = plt.subplots(2)
            ax1.set_title('Concentrations')
            ax1.grid(True)
            ax1.set_xlabel('Distance (z)')
            ax1.set_ylabel('Concentration')

            ax1.plot(z1, C1[i], label=f'{labels[i]}_{ug_valeurs[0]} (m/s)', color='red')
            ax1.plot(z2, C2[i], label=f'{labels[i]}_{ug_valeurs[1]} (m/s)', color='green') # Carbon capture enabled
            ax1.legend(loc="upper right")
            
            ax2.set_title('T (Kelvins)')
            ax2.plot(z1, C1[6], label='T', color='red')
            ax2.plot(z2, C2[6], label='T', color='green')
            ax2.legend(loc="upper right")

            plt.tight_layout()
            plt.grid(True)
            plt.legend()
            plt.show()

elif choice == 'VCOMP':
    initial_value = float(input("Initial value?: ")) # 0.06 minimal
    final_value = float(input("Final value?: "))
    nb_steps = int(input("How many steps in between?: "))
        
    z = [0, var.h_r] # longueur du réacteur
    UgSteps = np.linspace(int(initial_value), int(final_value), nb_steps)
    #print("UsSteps: ", UsSteps)
    UgConcentrations = []
    #print("UsConcentrations: ", UsConcentrations)

    for temp_u_g in UgSteps:
        zr, Cr = calculConcentrationsIVP(z, C_0, CCUS=True, ug=temp_u_g) # zr and Cr stand respectively for z result and C result
        #print(Cr)
        UgConcentrations.append(Cr)
        
    #print("UsConcentrations: ", UsConcentrations)
    
    while True:
        key = input("Quel variable désirez-vous plotter? (exit pour sortir): ")
        if key == "exit":
            break
                
        else:
            if (graph_choice := input("MESH/LAST: ")) == "LAST": # print a 3D graph of all values or the last value only (output)                
                i = labels_tokenized[key]
                concentrations = [] #np.zeros(nb_steps)
                
                for k in UgConcentrations:
                    concentrations.append(k[i][-1]) # on néglige les points avant le dernier (car ce qui se passe dedans ne nous intéresse pas)

                fig, ax = plt.subplots()
                ax.plot(UgSteps, concentrations, label=f'{labels[i]}')
                ax.set_title('Comparaison des concentrations en fonction des valeurs de u_s')
                ax.set_ylabel(f'Value {labels[i]}')
                ax.set_xlabel('Value u_s')
                ax.grid(True)
                ax.legend(loc="upper right")
                plt.show()

            elif graph_choice == "MESH":
                from matplotlib import cbook, cm
                from matplotlib.colors import LightSource
                
                # Load and format data
                z = dem['elevation']
                nrows, ncols = z.shape
                x = np.linspace(dem['xmin'], dem['xmax'], ncols)
                y = np.linspace(dem['ymin'], dem['ymax'], nrows)
                x, y = np.meshgrid(x, y)

                region = np.s_[5:50, 5:50]
                x, y, z = x[region], y[region], z[region]

                # Set up plot
                fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

                ls = LightSource(270, 45)
                # To use a custom hillshading mode, override the built-in shading and pass
                # in the rgb colors of the shaded surface calculated from "shade".
                rgb = ls.shade(z, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
                surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=rgb,
                                    linewidth=0, antialiased=False, shade=False)

                plt.show()
"""

"""Q 3.4"""
"""
p_e = 1/22.4*var.R*973.15/100000 # pression au début du réacteur
ratios = input("Combien de valeurs du ratio CH4/H2O essayer?: ") # nombre de ratios à essayer
tol = 10**(-3) # tolerance (mentioned in report)
ratios = np.linspace(tol, 1-50*tol, int(ratios))

concentrationsSortie = [] # une "case" pour chaque ratio
index = 1

import time
time_per_calculation = [] # délai écoule (s)

for ratio in ratios:
    if index % 50 == 0: # print une valeur sur x pour voir où en est le calcul sans retarder
        percentage_achieved = (100/len(ratios))*index
        print(f'Calculated {percentage_achieved:.2f}%')
    index+=1
    
    P_CH4 = (1/22.4)/(ratio)
    P_H2O = (1/22.4)*(1-ratio)

    C_0 = [P_CH4, P_H2O, 10e-3, 0.00, 0.00, p_e, 973.15, 3]
    z = [0, var.h_r] # h_r étant longueur du réacteur, le calcul se fait de 0 à 0.29 (m)

    #ti = time.time()
    
    zr, Cr = calculConcentrationsIVP(z, C_0)
    #tf = time.time()
    #timespan = tf - ti
    #time_per_calculation.append(timespan)
    #print(f'{timespan:.2f} secondes écoulées')
    concentrationsSortie.append(Cr)

#fig, ax = plt.subplots()
#ax.plot(time_per_calculation)
#plt.show()

while True:
    key = input("Quel variable désirez-vous plotter? (exit pour sortir): ")
    if key == "exit":
        break

    else:
        i = labels_tokenized[key]
        
        concentrations_plot = []
        for k in concentrationsSortie:
            concentrations_plot.append(k[i][-1])
        
        plt.plot(concentrations_plot)
        plt.show()
"""

"""Q 4.1"""
"""
p_e = 1/22.4*var.R*973.15/100000
P_CH4 = (1/22.4)/(4)
P_H2O = (1/22.4)*(3/4)
C_0 = [P_CH4, P_H2O, 10e-3, 0.00, 0.00, p_e, 973.15, 3]

def optimal_us(reforming_rate):
    initial_value = float(input("Initial value?: ")) # 0.06 minimal
    final_value = float(input("Final value?: "))
    nb_steps = int(input("How many steps in between?: "))
        
    z = [0, var.h_r] # longueur du réacteur
    UsSteps = np.linspace(initial_value, final_value, nb_steps)
    #print("UsSteps: ", UsSteps)
    UsConcentrations = []
    for temp_u_s in UsSteps:
        zr, Cr = calculConcentrationsIVP(z, C_0, CCUS=True, us=temp_u_s)
        UsConcentrations.append(Cr)
    

    concentrations_CO2 = [k[4][-1] for k in UsConcentrations] # où on s'intéresse seulement aux points de sortie
    concentrations_CO = [k[3][-1] for k in UsConcentrations]
    nb_comp = len(concentrations_CO) # nombre de points de comparaisons CO/CO2

    concentrations_comparison = [(concentrations_CO[i]/concentrations_CO2[i])*100 for i in range(nb_comp)] # sur 100% de CO, combien de % ont été reformés en CO2?
     
    
    fig, ax = plt.subplots(2, 1)
    ax[0].plot(UsSteps, concentrations_CO2, label=f'CO2')
    ax[0].plot(UsSteps, concentrations_CO, label=f'CO')
    ax[0].set_title('Comparaison des concentrations en fonction des valeurs de u_s')
    ax[0].set_ylabel(f'Value CO2')
    ax[0].set_xlabel('Value u_s')
    ax[0].grid(True)
    ax[0].legend(loc="upper right")
    
    ax[1].plot(UsSteps, concentrations_comparison, label=f'Reformage CO2')
    ax[1].grid(True)
    ax[1].set_ylabel('Y%')
    ax[1].set_xlabel('Conversion rate')
    plt.show()

u_s = optimal_us(reforming_rate=0.075) # decimal for 7.5%
print(u_s)
"""


"""
(Y%)
u_g fixé
C_0 = [P_CH4, P_H2O, 10e-3, 0.00, 0.00, p_e, 973.15, 3]
u_s donné

z = [0, var.h_r]
z2, C2 = calculConcentrationsIVP(z, C_0, CCUS=True)
"""


