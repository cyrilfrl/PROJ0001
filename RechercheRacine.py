
import math
import numpy as np

def bissection(f, x0, x1, tol):
    # statut
    statut = 0 # on part du principe que l'algo trouvera une racine, donc on définit le statut comme étant correct d'entrée de jeu

    # hypothèses
    y0 = f(x0)
    y1 = f(x1)
    
    if y0*y1 > 0:
        erreur = "ErrorMessage: no root in the function"
        statut = 1
        return [erreur, statut]

    # cas x0 après x1 (mauvaise utilisation)
    if x1 < x0:
        print('Attention, merci de fournir les bornes dans leur ordre croissant')
        x0, x1 = x1, x0
    
    # max iterations
    N = math.ceil(math.log((math.fabs(x1 - x0))/(2*tol), math.e)/math.log(2, math.e))
        
    # convergence loop
    if f(x0) > f(x1):
        k = 0 # nombre d'itérations actuelles
    
        while k < N:
            # calcul de l'abscisse centrale
            x_mid = x0 + math.fabs(x1-x0)/2
            f_mid = f(x_mid)
            
            if f_mid < 0:
                x1 = x_mid
                
            elif f_mid > 0:
                x0 = x_mid
            
            # on augmente le nombre d'itérations de 1
            k += 1
            
            # on vérifie si convergence atteinte
            if math.fabs(f_mid) <= tol:
                return [x_mid, statut]
        
    elif f(x1) > f(x0):
        k = 0 # nombre d'itérations actuelles
    
        while k < N:
            x_mid = x0 + math.fabs(x1-x0)/2
            f_mid = f(x_mid)
            debug.append(f_mid)
            
            if f_mid < 0:
                x0 = x_mid
                
            elif f_mid > 0:
                x1 = x_mid
                
            k += 1
            
            # on vérifie si convergence atteinte
            if math.fabs(f_mid) <= tol:
                return [x_mid, statut]

    statut = 1
    return[None,statut]

def secante(f, x0, x1, tol) :
    statut = 0

    # Hypothesis verification
    y0 = f(x0)
    y1 = f(x1)
    k = 0

    if tol <= 10**(-10):
        erreur = "ErrorMessage: invalid tolerance"
        statut = 1
        return [4, statut]
    
    if abs(y0 - y1) <= 10**(-10):
        erreur = "ErrorMessage: point unique"
        statut = 1
        return [4, statut]
    
    N = int(np.log(abs(x1-x0)/(2*tol)))/(np.log(2)) + 1
    while abs(y1) > tol:
        x2 = x1
        x1 = (x1 - (y1 * (x1-x0))/(y1 - y0))
        x0 = x2
        y0 = y1
        y1 = f(x1)
        
        if(abs(y0-y1) <= 10**(-10)):
            erreur = "ErrorMessage: doesn't converge"
            statut = -1
            return [x1, statut]
        
        k += 1
        if k > N:
            statut = -1
            return [x1, statut]
    
    return [x1, statut]
