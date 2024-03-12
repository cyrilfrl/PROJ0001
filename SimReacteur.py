
# modules import
import numpy as np 
import RechercheRacine
import SimReacteur
import Variables


# odefunction

def dcidz(r_i, r_cbn):
    return (eta(1-epsilon)*rho_cat*r_i - (1-epsilon)*rho_CaO*r_cbn)/u_g

def odefunction(z,C):
    # définition variables dépendantes des arguments de la fonction
    C_CH4 = C[0]
    P_CH4 = C_CH4*V # (mol)
    C_H20 = C[1]
    P_H20 = C_H20*V # (mol)
    C_H2 = C[2]
    P_H2 = C_H2*V
    C_CO = C[3]
    P_CO = C_CO*V
    C_CO2 = C[4]
    P_CO2 = C_CO2*V
    X = C[5] # (dimensionless) conversion fractionnaire du CaO en CaCO3
    T = C[6]
    P = C[7]

    DEN = 1 + K_CO*P_CO + K_H2*P_H2 + K_CH4*P_CH4 + K_H20*P_H20/P_H2
    R_1 = (k_1/(P_H2**2.5))*(P_CH4*P_H20 - (P_H2**3*P_CO)/K_1)/DEN**2       # vitesse de réaction reformage (méthane + eau -> monoxyde de carbone + dihydrogène)
    R_2 = (k_2/(P_H2**3.5))*(P_CH4*(P_H20**2) - (P_H2**4*P_CO2)/K_2)/DEN**2 # méthane + eau -> dioxyde de carbone + dihydrogène
    R_3 = (k_3/P_H2)*(P_CO*P_H20 - (P_H2*P_CO2)/K_3)/DEN**2                 # monoxyde de carbone + eau -> dioxyde de carbone pour dihydrogène

    r_CH4 = -R_1 - R_2          # (kmol/(kg.s))
    r_H20 = -R_1 - 2*R_2 - R_3  # (kmol/(kg.s))
    r_H2 = 3*R_1 + 4*R_2 + R_3  # (kmol/(kg.s))
    r_CO = R_1 - R_3            # (kmol/(kg.s))
    r_CO2 = R_2 + R_3           # (kmol/(kg.s))

    M_k = 303
    N_k = -13146
    k_c = M_k*np.exp(N_k/T)             # (/s) vitesse apparente de carbonatation
    M_b = 1.6                           # (s)  temps pour arriver à la moitié de la conversion ultime X_u
    N_b = 5649
    b = M_b*np.exp(N_b/T)
    X_u = k_c*b                         # conversion ultime
    r_cbn = (k_c/M_CaO)*(1-X/X_u)**2    # (kmol/(kg.s)) taux de consommation par carbonatation de CO2 (conversion fractionnaire CaO en CaCO3)

    r = [r_CH4, r_H20, r_H2, r_CO, r_CO2]
    dC_dz = [dcidz(r_i, r_cbn) for r_i in r]

    return dC_dz

def calculConcentrationsEuler(Z, C0):
    z0, zf = Z

    #(kmol/m^3)

    return z, C 

#z, C = calculConcentrationsEuler([z0, z1], C0)


def euler(sirmodel, y0, h, Ti, Tf):
    pas = int(round(float(Tf - Ti)/h))# arrondir si nombre a virgules
    f_ = lambda y, t: sirmodel(t, y, beta, gamma)#ca permet de definir les variables d'une fonction
    y = zeros((pas+1, len(y0)))#initialise le tableau a deux dimensions avec des zeros (shape)
    t = linspace(0, pas*h, len(y))#len compte la taille du tableau, linspace crée l'intervalle de temps
    y[0] = y0#initialisation du tableau y avec les conditions initiales
    for i in range(pas):#utiliser la formule du cours en fonction du pas
        y[i+1] = y[i] + h*f_(y[i], t[i])
    return y, t

def EulerSEIR(R0cible, gamma, sigma, eta, y0):
    ti= 0
    tf = 400
    h=0.01
    
    y1=np.zeros(5)
    y2= np.zeros(5)

    y1[0] = y0[0]
    y1[1] = y0[1]
    y1[2] = y0[2]
    y1[3] = y0[3]
    y1[4] = y0[4]
    i=1


    intervalle = int((tf-ti)/h + 1)

    temps = np.zeros(intervalle)
    s = np.zeros_like(temps)
    x = np.zeros_like(temps)
    r = np.zeros_like(temps)
    e = np.zeros_like(temps)
    b  = np.zeros_like(temps)
    
    s[0] = y0[0]
    x[0] = y0[1]
    r[0] = y0[2]
    e[0] = y0[3]
    b[0] = y0[4]
    temps[0] = ti
    t = ti+h


    while i<intervalle-1:
       
       dy= seirmodel(t,y1, gamma,sigma,eta, R0cible)
       y2[0]= y1[0] +dy[0]*h
       y2[1]= y1[1] +dy[1]*h
       y2[2]= y1[2] +dy[2]*h
       y2[3]= y1[3] +dy[3]*h
       y2[4]= y1[4] +dy[4]*h
      
       s[i] = y2[0]
       e[i]= y2[1]
       x[i] = y2[2]
       r[i] = y2[3]
       b[i]= y2[4]
       temps[i] = t
      
       
       y1[0]=y2[0]
       y1[1]=y2[1]
       y1[2]=y2[2]
       y1[3]=y2[3]
       y1[4]=y2[4]
       
       i=i+1 
       
       t= t+h
       
       

    #Calcul pour le temps t= 400
    dy = seirmodel(tf,y1,gamma, sigma, eta, R0cible)
    y2[0]= y1[0] +dy[0]*h
    y2[1]= y1[1] +dy[1]*h
    y2[2]= y1[2] +dy[2]*h
    y2[3]= y1[3] +dy[3]*h
    y2[4]= y1[4] +dy[4]*h

    s[-1] = y2[0]
    e[-1]= y2[1]
    x[-1] = y2[2]
    r[-1]= y2[3]
    b[-1]= y2[4]
    temps[-1] = tf
    
    # A choisir d'afficher ou non
    """plt.figure(figsize=(10,10))
    plt.plot(temps, s, 'green', label='Susceptibles')
    plt.plot(temps, e, 'grey', label='Exposés')
    plt.plot(temps, x, 'black', label='Infectés')
    plt.plot(temps, r, 'red', label='Guéris')
    plt.title(" Evolution de la propagation du coronavirus")
    plt.legend(loc='best')
    plt.xlabel('Temps (jour)')
    plt.ylabel('population')"""
    
    return max(x)