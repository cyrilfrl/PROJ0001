
# modules import
import matplotlib.pyplot as plt
import numpy as np
import RechercheRacine as root_search
from scipy.integrate import solve_ivp
import Variables as var

# plotting
def easy_plotting(absi, ordo, title, x_label, y_label, color='blue'):
    fig, ax = plt.subplots()
    ax.plot(absi, ordo, c=color)
    ax.set_title(title)
    ax.grid(True)
    plt.show()

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
    
    if P_H2 <= 0: # éviter epsilon machine
        P_H2 = 1e-5
    
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
        k_c = var.M_k*np.exp(var.N_k/T)             # (/s) vitesse apparente de carbonatation
        b = var.M_b*np.exp(var.N_b/T)               # (/s) vitesse apparente de carbonatation
        X_u = k_c*b                                 # (dimensionless) conversion ultime
        r_cbn = (k_c/var.M_CaO)*(1-X/X_u)**2        # (kmol/(kg.s)) taux de consommation par carbonatation de CO2 (conversion fractionnaire CaO en CaCO3)   
        
    elif CCUS == False:
        us = 0 # on ne fait pas avancer le pellet
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
    if us != 0:
        dX_dz = (var.M_CaO*r_cbn)/us
    elif us == 0:
        dX_dz = 0
    dC_dz[5] = dX_dz
    
    # T
    rho_g = (1/(var.R*T))*100*(var.M_CH4*P_CH4 + var.M_H2O*P_H2O + var.M_H2*P_H2 + var.M_CO*P_CO + var.M_CO2*P_CO2)
    rho_s = (var.W_cat + var.W_CaO)/((var.W_cat/var.rho_cat)+(var.W_CaO/var.rho_CaO))
    
    divider = (1 - var.epsilon)*rho_s*us*var.C_ps + rho_g*ug*var.C_pg
    if divider <= 1e-6: # éviter epsilon machine, résolution bug proposée par professeur tp
        divider = 1e-6

    dT_dz = (-(1 - var.epsilon) * var.rho_cat*var.eta * (R_1 * var.H_R1 + R_2 * var.H_R2 + R_3 * var.H_R3) - (1 - var.epsilon)*var.rho_CaO*r_cbn*var.H_cbn + h_w*(var.T_W-T)*(4/var.D_r))/divider
    dC_dz[6] = dT_dz
        
    # P
    dP_dz = (-rho_g*(ug**2)/var.d_p)*((1- var.epsilon)/var.epsilon)*((150*(1 - var.epsilon)*var.mu)/(var.d_p*rho_g*ug) + 1.75)*(10**(-5))
    rho_gas = rho_g
    dC_dz[7] = dP_dz
    
    return dC_dz

def calculConcentrationsEuler(Z, C_0, CCUS=True, ug=var.u_g, us=var.u_s):
    z_0, z_f = Z
    z_t = z_0 # z at time t
    c_t = C_0
    n_euler = 500000 # nombre de points
    
    # Hyperparameters variables
    z = np.linspace(z_0, z_f, n_euler)    # coordinates between z_0 and z_f with n values (datapoints)
    #h = 0.001                            # step size (pas)
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

    sol = solve_ivp(lambda t, C: odefunction(t, C, CCUS, ug, us), [z_0, z_f], C_0, t_eval=z, rtol=10e-5, method='RK45') # time span and initial concentrations (tried 0.0002 but not sufficient because had "0.00023864609591485773 not less than 0.0002")
    z = sol.t
    C = sol.y
    
    return z, C

##################################################
#                   EVALUATION                   #
##################################################

p_e = 1/22.4*var.R*973.15/100000 # pression au début du réacteur
P_CH4 = (1/22.4)/(4)
P_H2O = (1/22.4)*(3/4)
C_0 = [P_CH4, P_H2O, 10e-3, 0.00, 0.00, p_e, 973.15, 3]
C_0 = [P_CH4, P_H2O, 1e-5, 0.00, 0.00, 0.00, 973.15, 3]
z = [0, var.h_r] # h_r étant longueur du réacteur, le calcul se fait de 0 à 0.29 (m)
# clés pour print les graphiques de concentrations spécifiques
labels_tokenized = {'CH4': 0, 'H2O': 1, 'H2': 2, 'CO': 3, 'CO2': 4, 'TEMP': 6, 'GENERAL':7}
labels = ['CH4', 'H2O', 'H2', 'CO', 'CO2','TEMP']


##################################################
#                VARIANTES EVAL                  #
##################################################

while True:
    print('')
    choice = input("NOCC/US_VARI/UG_VARI/ENTRY_RATIO/OPTI_US/OPTI_TEMP: ")
    
    if choice == "NOCC":
        """ Q 3.1"""
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
                ax1.set_title(f'Concentrations {labels[i]} avec ou sans carbon capture')
                ax1.grid(True)
                ax1.set_xlabel('Distance z (m)')
                ax1.set_ylabel('Concentration (mol/L)')

                ax1.plot(z1, C1[i], label=f'{labels[i]}', color='red')
                ax1.plot(z2, C2[i], label=f'{labels[i]}_CCEN', color='green') # Carbon capture enabled
                ax1.legend(loc="upper right")
                
                ax2.set_title('Température avec ou sans carbon capture')
                ax2.plot(z1, C1[6], label='T', color='red')
                ax2.plot(z2, C2[6], label='T_CCEN', color='green')
                ax2.legend(loc="upper right")
                ax2.set_xlabel('Distance z (m)')
                ax2.set_ylabel('Température (K)')
                plt.tight_layout()
                plt.grid(True)
                plt.legend()
                plt.show()        
    
    elif choice == "US_VARI":
        """Q 3.2"""
        initial_value = 0.01
        final_value = 3
        nb_steps = 200
        z = [0, var.h_r] # longueur du réacteur
        Steps = np.linspace(float(initial_value), float(final_value), nb_steps)
        concentrations_unfiltered = []
        C_0 = [P_CH4, P_H2O, 10**(-3), 0.00, 0.00, 0.00, 973.15, 3]

        for temp_u_s in Steps:
            zr, Cr = calculConcentrationsIVP(z, C_0, CCUS=True, us=temp_u_s) # zr and Cr stand respectively for z result and C result
            concentrations_unfiltered.append(Cr)

        # Plot des différents graphiques
        fig, ax = plt.subplots(3, 2, figsize=(12, 8), constrained_layout=True)

        for i in [j for j in range(6) if j != 1]:
            concentrations = [k[i][-1] for k in concentrations_unfiltered]

            label = labels[i]
            ax[i%3, i%2].plot(Steps, concentrations, label=label)            
            if label != "TEMP":
                ax[i%3, i%2].set_title(f'Concentrations de sortie en {label} selon us')
                ax[i%3, i%2].set_ylabel(f'Concentrations {label} (mol/L)')

            elif label == "TEMP":
                ax[i%3, i%2].set_title(f'Température de sortie selon us')
                ax[i%3, i%2].set_ylabel(f'Température (K)')

            ax[i%3, i%2].set_xlabel('Valeur u_s (m/s)')
            ax[i%3, i%2].grid(True)
            ax[i%3, i%2].legend(loc="upper right")


        for i in [j for j in range(6) if j != 1]:
            concentrations = [k[i][-1] for k in concentrations_unfiltered]
            ax[1%3, 1%2].plot(Steps, concentrations, label=labels[i])
            ax[1%3, 1%2].set_title(f'Concentrations de toutes les substances selon us')
            ax[1%3, 1%2].set_ylabel('Concentrations (mol/L)')
            ax[1%3, 1%2].set_xlabel('Valeur u_s (m/s)')
            ax[1%3, 1%2].grid(True)
            ax[1%3, 1%2].legend(loc="upper right")

        plt.show()
    
    elif choice == "UG_VARI":
        """Q 3.3"""
        initial_value = 0.01
        final_value = 20
        nb_steps = 200
        z = [0, var.h_r] # longueur du réacteur
        Steps = np.linspace(float(initial_value), float(final_value), nb_steps)
        concentrations_unfiltered = []
        C_0 = [P_CH4, P_H2O, 10**(-3), 0.00, 0.00, 0.00, 973.15, 3]

        for temp_u_g in Steps:
            zr, Cr = calculConcentrationsIVP(z, C_0, CCUS=True, ug=temp_u_g) # zr and Cr stand respectively for z result and C result
            concentrations_unfiltered.append(Cr)

        # Plot des différents graphiques
        fig, ax = plt.subplots(3, 2, figsize=(12, 8), constrained_layout=True)

        for i in [j for j in range(6) if j != 1]:
            concentrations = [k[i][-1] for k in concentrations_unfiltered]

            label = labels[i]
            ax[i%3, i%2].plot(Steps, concentrations, label=label)            
            if label != "TEMP":
                ax[i%3, i%2].set_title(f'Concentrations de sortie en {label} selon ug')
                ax[i%3, i%2].set_ylabel(f'Concentrations {label} (mol/L)')

            elif label == "TEMP":
                ax[i%3, i%2].set_title(f'Température de sortie selon ug')
                ax[i%3, i%2].set_ylabel(f'Température (K)')

            ax[i%3, i%2].set_xlabel('Valeur u_g (m/s)')
            ax[i%3, i%2].grid(True)
            ax[i%3, i%2].legend(loc="upper right")


        for i in [j for j in range(6) if j != 1]:
            concentrations = [k[i][-1] for k in concentrations_unfiltered]
            ax[1%3, 1%2].plot(Steps, concentrations, label=labels[i])
            ax[1%3, 1%2].set_title(f'Concentrations de toutes les substances selon ug')
            ax[1%3, 1%2].set_ylabel('Concentrations (mol/L)')
            ax[1%3, 1%2].set_xlabel('Valeur u_g (m/s)')
            ax[1%3, 1%2].grid(True)
            ax[1%3, 1%2].legend(loc="upper right")

        plt.show()
        
    elif choice == "ENTRY_RATIO":
        """Q 3.4"""
        ratios = input("Combien de valeurs du ratio CH4/H2O essayer?: ") # nombre de ratios à essayer
        tol = 10**(-3) # tolerance (mentioned in report)
        ratios = np.linspace(tol, 1-50*tol, int(ratios))

        concentrationsSortie = [] # une "case" pour chaque ratio
        index = 1

        for ratio in ratios:
            if index % 10 == 0: # print une valeur sur x pour voir où en est le calcul sans retarder
                percentage_achieved = (100/len(ratios))*index
                print(f'Calculated {percentage_achieved:.2f}%')
            index+=1
            
            P_CH4 = (1/22.4)/(ratio)
            P_H2O = (1/22.4)*(1-ratio)

            C_0 = [P_CH4, P_H2O, 10e-3, 0.00, 0.00, p_e, 973.15, 3]

            zr, Cr = calculConcentrationsIVP(var.z, C_0)
            concentrationsSortie.append(Cr)

        # Plot des différents graphiques
        fig, ax = plt.subplots(3, 2, figsize=(12, 8), constrained_layout=True)

        for i in range(6):
            concentrations = [k[i][-1] for k in concentrationsSortie]

            label = labels[i]
            ax[i%3, i%2].plot(ratios, concentrations, label=label)            
            if label != "TEMP":
                ax[i%3, i%2].set_title(f'Concentrations de sortie en {label} selon ratio CH4/H20')
                ax[i%3, i%2].set_ylabel(f'Concentrations {label} (mol/L)')

            elif label == "TEMP":
                ax[i%3, i%2].set_title(f'Température de sortie selon ratio CH4/H20')
                ax[i%3, i%2].set_ylabel(f'Température (K)')

            ax[i%3, i%2].set_xlabel('Valeur ratio')
            ax[i%3, i%2].grid(True)
            ax[i%3, i%2].legend(loc="upper right")

        plt.show()

    
    elif choice == 'OPTI_US':
        """Q 4.1"""
        ug = float(input("Valeur de ug à utiliser (m/s)?: "))
        Y = float(input("Rendement souhaité? (décimal e.g. 7.5% = 0.075): "))
        print(f'Optimisation pour {Y*100}% ...')
        
        def optimization_us(C, ug, Y):
            def rendement(us):
                partial_ode = lambda z, C: odefunction(z=var.z, C=C, ug=ug, us=us)
                sol = solve_ivp(partial_ode, var.z, C, rtol = 1e-6)
                Cr = sol.y # Cr stands for C"result" because C was already taken
                
                P_CO2 = Cr[4][-1] # pression partielle en sortie de CO2
                P_TOT = sum([Cr[i][-1] for i in range(5) if i != 1]) # somme pressions partielles hors eau
                R_CO2 = P_CO2/P_TOT - Y # rendement en CO2
                
                return R_CO2
            
            specific_yield_us = root_search.bissection(rendement, 10e-6, 0.3, tol=10e-6)

            graph_yn = input("Graph Y/N: ")
            if graph_yn == "Y":
                absi = np.linspace(10e-2, 2, 50)
                ordo = [rendement(i) for i in absi]
                                
                easy_plotting(absi, ordo,
                            'Rendement en CO2 en fonction de la vitesse du CaO',
                            'Vitesse du CaO us (m/s)',
                            'Rendement CO2 (%)',
                            color='darkviolet')

            return specific_yield_us
        
        C_0 = [P_CH4, P_H2O, 1e-5, 0.00, 0.00, 0.00, 973.15, 3]
        optimal_value = optimization_us(C_0, ug, Y)[0]
        print(f'Valeur optimale: {optimal_value} m/s')

    elif choice == "OPTI_TEMP":
        """Q 4.3"""
        def optimization_us(C, ug, Y):
            def rendement(us):
                partial_ode = lambda z, C: odefunction(z=var.z, C=C, ug=ug, us=us)
                sol = solve_ivp(partial_ode, var.z, C, rtol = 1e-6)
                Cr = sol.y # Cr stands for C"result" because C was already taken
                
                P_CO2 = Cr[4][-1] # pression partielle en sortie de CO2
                P_TOT = sum([Cr[i][-1] for i in range(5) if i != 1]) # somme pressions partielles hors eau
                R_CO2 = P_CO2/P_TOT - Y # rendement en CO2
                
                return R_CO2
            
            specific_yield_us = root_search.bissection(rendement, 10e-6, 0.3, tol=10e-6)

            return specific_yield_us

        Y = float(input("Rendement souhaité? (décimal e.g. 7.5% = 0.075): "))
        absi = np.linspace(965, 1000, 50)
        ordo = []

        for temperature in absi:
            try:
                us = optimization_us([P_CH4, P_H2O, 1e-5, 0.00, 0.00, 0.00, temperature, 3], var.u_g, Y)[0]
                
                if not isinstance(us, str):
                    ordo.append(us)
                    
                else:
                    ordo.append(None)
                    
            except:
                ordo.append(None)

        fig, ax = plt.subplots()
        ax.plot(absi, ordo, c='darkviolet')
        Y *= 100
        ax.set_title(f'Us optimal pour {Y} % de rendement souhaité')
        ax.grid(True)
        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel('Us optimal (m/s)')
        plt.show()    
        
    elif choice == "DOCU":
        options = {
            "NOCC":"NO Carbon Capture",
            "US_VAR":"Variations des concentrations selon valeur de Us",
            "UG_VARI": "Variation des concentrations selon valeur de Ug",
            "ENTRY_RATIO":"Between 0 and 1, variate the proportions of CH4/H20 but keep the same amount of matter",
            "OPTI_US":"For a given value of ug (speed of gas) return optimal value of us to set to get a given yield of CO2 reconversion, une fonction complémentaire pour avoir \
                un rendu visuel du graphique de la fonction dont il faut trouver les racines est proposée pendant le déroulement de la recherche",
            "OPTI_TEMP":"For a given desired CO2 yield compare different optimal us (speed of CaO) depending on initial temperature"
        }
        for key in options:
            print(f'{key}: {options[key]}')

