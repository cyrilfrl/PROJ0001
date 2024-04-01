

# modules import
import numpy as np 
#import RechercheRacine
import Variables as var
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#On definit cette fonction a utiliser dans odefunction.
#Les seuls parametres qui changent sont r_i en fct de l'element et r_cbn qui n'est pas constant
    
def dcidz(r_i, r_cbn):
    return (var.eta*(1-var.epsilon)*var.rho_cat*r_i - (1-var.epsilon)*var.rho_CaO*r_cbn)/var.u_g
    
def odefunction(z,C):
    
    # définition variables dépendantes des arguments de la fonction
    
    X = C[5]              # (dimensionless) conversion fractionnaire du CaO en CaCO3
    T = C[6]
    P = C[7]

    #les pressions partielles dependent directement des arguments de C = arguments de odefunction!
    #Calcul des pression partielles (P * C/C_tot avec PV=NRT et P_part = P_tot*n/n_tot, sachant C = n/V où V s'annule)
    
    C_tot = np.sum(C[0:5]) #Concentration totale = somme des concentations partielles = somme des 5 premiers elements de C[]
    
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
    
    #definition des variables qui dependent de T la temperature
    
    #constantes d'equilibre K
    
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
    
    #vitesses de reaction

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
    #r_cbn = (k_c/var.M_CaO)*(1-X/X_u)**2    # (kmol/(kg.s)) taux de consommation par carbonatation de CO2 (conversion fractionnaire CaO en CaCO3)   
    r_cbn = 0       #valeur test
    
    ################################################
    
    # calculs du bilan énergétique
    
    rho_g = (1/(var.R*T))*(var.M_CH4*P_CH4 + var.M_CO*P_CO + var.M_CO2*P_CO2 + var.M_H2*P_H2 + var.M_H2O*P_H2O)*100 # (kg/m^3) masse volumique de la phase gazeuse, *100 pour la conversion d'unités    
    Re_p = (var.u_g*var.epsilon*rho_g*var.d_p)/var.mu
    
    #page 10- h_w est une variable qui depend de Re_p qui varie en foncition de T(?) et on la definit dans divers cas de figure
    
    if Re_p > 20 and var.d_p/var.D_r < 0.3:
        h_w = 2.03*(var.k_g/var.D_r)*(Re_p**(0.8))*np.exp((-6*var.d_p)/(var.D_r))
        
    elif Re_p < 20:
        h_w = 6.15*(var.k_z0/var.D_r)
        
    else:
        h_w = "error"
    
    # calcul des dérivées de Ci, P , X et T grace aux formules donnees
    
    dC_dz = np.zeros(8) #on initialise le tableau
    
    r = [r_CH4, r_H2O, r_H2, r_CO, r_CO2]   # uniques paramètres variant d'une équation à l'autre, en f(T)
    for i in range(5): #i va de 0 a 4 pour parcourir les 5 elements de r
        dC_dz[i] = dcidz(r[i], r_cbn) 
        #les 5 premieres derivees sont mises dans dC_dz a retourner.
        
    # Derivee de X
    dX_dz = (var.M_CaO*r_cbn)/var.u_s
    dC_dz[5] = dX_dz
    
    # Derivee de T
    dT_dz = (-(1 - var.epsilon) * var.rho_cat * var.eta * (R_1*var.H_R1 + R_2*var.H_R2 + R_3*var.H_R3) - (1 - var.epsilon) * var.rho_CaO*r_cbn*var.H_cbn + h_w * (var.T_W - T) * 4 / var.D_r) / ((1 - var.epsilon)*var.rho_s*var.u_s*var.C_ps + rho_g*var.u_g*var.C_pg)
    dC_dz[6] = dT_dz
    
    # Derivee de P
    dP_dz = (-rho_g*(var.u_g**2)/var.d_p)*((1- var.epsilon)/var.epsilon)*((150*(1 - var.epsilon)*var.mu)/(var.d_p*rho_g*var.u_g) + 1.75)*(10**(-5))
    dC_dz[7] = dP_dz
    
    return dC_dz

#On a maintenant la derivee des concentrations par rapport a la distance z parcourue dans le gaz dans le reacteur
#On definit Euler pour calculer

def calculConcentrationsEuler(Z, C_0):
    
    z_0, z_f = Z  #z0 le debut du reacteur, zf la fin du reacteur. Z est un tuple qui englobe les 2.
    
    z_t = z_0  # z_t la position au temps t. au debut, c'est z0 (pour les CI)
    c_t = C_0 # c_t la concentration au temps t. c'est un array de taille 8
    
    # nombre de points n = nombre d'iterations qui definiront le pas
    n_euler = 500000 
    
    # Hyperparameters variables
    
    # on cree z, un array de n=n_euler points qui commence de z_0 et termine a z_f 
    z = np.linspace(z_0, z_f, n_euler)   
    
    h = (z_f-z_0)/n_euler  
    print(f'Le pas h = {h}') # h = 5.8e-5
    
    C = np.zeros([8, n_euler])            # initialisation tableau de taille 8 x n, n étant le nombre d'éléments recherchés
    
    ############CALCUL###########
    #methode d'euler = calcul des concentrations a partir des derivees et concentrations initiales
    
    for i in range(n_euler): # i va de 0 a n_euler 
        z_t += h  #la position actuelle = abscisse est mise a jour avec le pas
        derivee = odefunction(z_t, c_t) #pour utiliser la formule euler. retourne un array
        c_tt = c_t + h*derivee  #formule d'euler : calcul de la prochaine concentration. c_tt = array + h* array
        C[:,i] = c_tt #chaque colonne de C est remplacee par l'array c_tt
        c_t = c_tt #on remplace c_t par c_actuel pour continuer la boucle 
        
    return [z, C]

def calculConcentrationsIVP(z, C_0): # solve Initial Problem Value (IVP)

#cette fonction de numpy permet de resoudre directement l'EDO. On s'en sert pour comparer nos resultats

    z_0, z_f = z  
    z = np.linspace(z_0, z_f, var.n) #avec n = 300, z = un array de n points qui commence de z_0 et termine a z_f

    sol = solve_ivp(lambda t, C: odefunction(t, C), [z_0, z_f], C_0, t_eval=z, rtol=10e-5) 
    #avec lambda on peut coder une fonction en une ligne : t,c = arguments, return = odefunction(t,C). 
    #t_eval = intervalle de resolution et rtol = tolerace qu'on impose au solver
    z = sol.t
    C = sol.y
    
    return z, C


##################################################
#                   EVALUATION                   #
##################################################

#test de la fonction avec des C.I. imposees par la consigne

p_e = 1/22.4*var.R*973.15/100000 # pression au début du réacteur
P_CH4 = (1/22.4)/(4)
P_H20 = (1/22.4)*(3/4)  ######## H20 et pas H2O? Cyril t'es trop con je te jure
C_0 = [P_CH4, P_H20, 10e-3, 0.00, 0.00, p_e, 973.15, 3] #ici aussi Cyril
z = [0, var.h_r] # h_r étant longueur du réacteur, le calcul se fait de 0 à 0.29 (m)

# Test Euler
z, C = calculConcentrationsEuler(z, C_0)

# Test IVP
z = [0, var.h_r] # longueur du réacteur
z, C = calculConcentrationsIVP(z, C_0)

#plots

fig, (ax1, ax2) = plt.subplots(2)  #on cree une figure avec 2 subplots : ax1 et ax2

ax1.set_title('Concentrations SANS capture C02')
ax1.grid(True)
ax1.set_xlabel('Distance (z)')
ax1.set_ylabel('Concentration')

#ici on cree une boucle qui va attribuer le nom a la courbe des 5 premiers elements de C.

labels = ['CH4', 'H2O', 'CH2', 'CO', 'CO2']
for i in range(5):
    ax1.plot(z, C[i], label=f'{labels[i]}')
    
    """
    À chaque itération de la boucle, cette instruction trace une courbe sur 
    la figure ax1 avec les données z en abscisse et les données correspondantes de C[i] en ordonnée. 
    L'étiquette de cette courbe est définie comme le nom de l'espèce chimique correspondante, 
    récupéré à partir de la liste labels à l'indice i.
    
    """

ax1.legend(loc="upper right")

# X, T, P plot
ax2.set_title('X, T, P SANS capture C02')

#pareil pour le graphe 2 sauf qu'ici on doit prendre les 3 derniers elements de C

labels = ['X', 'T', 'P']
for i in range(5, 8):  #5 a 8 car ce sont les elements de C
    ax2.plot(z, C[i], label=f'{labels[i-5]}')
    
    """
    Lorsque i est égal à 5, i - 5 vaut 0. Donc, vous utilisez labels[0], ce qui correspond au premier élément de la liste labels, soit 'X'.
Lorsque i est égal à 6, i - 5 vaut 1. Ainsi, vous accédez à labels[1], qui est le deuxième élément de la liste labels, soit 'T'.
Lorsque i est égal à 7, i - 5 vaut 2. Vous accédez donc à labels[2], qui est le troisième élément de la liste labels, soit 'P'.
    """

ax2.legend(loc="upper right")

plt.tight_layout()
plt.grid(True)
plt.legend()
plt.show()

########### plot du CO2 et X

fig, (ax1,ax2) = plt.subplots(2)
    
ax1.plot(z,C[4],label='CO2')
ax1.set_title('Concentration CO2 SANS capture C02')
ax1.set_xlabel('z distance parcourue dans reacteur')
ax1.set_ylabel('Concentration de CO2')
ax1.grid(True)
ax1.legend(loc="lower right")

plt.tight_layout()


ax2.plot(z,C[5],label='X')
ax2.set_title('Conversion fractionnaire de CaO en CaCO3')
ax2.set_xlabel('z distance parcourue dans reacteur')
ax2.set_ylabel('X')
ax2.grid(True)
ax2.legend(loc="lower right")

plt.show()
plt.tight_layout()









