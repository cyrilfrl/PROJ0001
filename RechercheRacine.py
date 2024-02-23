
import math
import numpy as np

def bissection(f, x0, x1, tol):
    # cas x0 après x1
    if x1 < x0:
        x0, x1 = x1, x0
    
    # status verification (hypotheses)
    statut = 0
    half_ecart = (x1-x0)/2
    y0 = f(x0)
    y1 = f(x1)
    if y0*y1 > 0:
        erreur = "ErrorMessage: no root in the function"
        statut = 1
        return [erreur, statut]

    elif y0*y1 == 0:
        if y0 == 0:
            answer = x0
        elif y1 == 0:
            answser = x1

        return [0, answer]
    
    # calcul du nombre max d'itérations
    N = math.ceil(math.log((x1 - x0)/(2*tol), math.e)/math.log(2, math.e))
    k = 0
    while True:
        k += 1
        y = f(x0+half_ecart)
        
        if np.sign(y0) == np.sign(y):
            x0 += half_ecart
            y0 = y
        
        elif np.sign(y1) == np.sign(y):
            x1 -= half_ecart
            y1 = y

        # checking if reached tolerance threshold
        half_ecart = (x1-x0)/2 # coordinates between x0 and x1
        if abs(f(x0 + half_ecart)) <= tol:
            return [x0 + half_ecart, statut]

    if k >= N:
        erreur = "ErrorMessage: la fonction ne converge pas"
        statut = -1
        return [erreur, statut]

def secante(f ,x0, x1, tol):
    statut = 0 # initialisation
    y0 = f(x0)
    y1 = f(x1)

    #cas ou x1>x0
    if x1<x0:
        x0, x1 = x1, x0
         
    #cas ou x0 est une racine
    if abs(y0) < tol : 
        return [x0 , statut]
    
    #cas ou x1 est une racine
    elif abs(y1) < tol :
        return [x1 , statut]
    
    # calcul du nombre max d'itérations
    N = math.ceil(math.log((x1 - x0)/(2*tol), math.e)/math.log(2, math.e))
    k = 0
    while k <= N:
        k += 1
        #si x0 est inferieur a x1! donc il faut prendre le nouvel intervalle a droite en changeant x0
        x_suiv = x1 - ( (y1*(x1-x0)) / (y1 - y0) )
        x0 = x1
        y0 = f(x0)
        x1 = x_suiv
        y1 = f(x1)
        
        if abs(y0) < tol:
            return [x0, statut]
        
    if k >= N:
        erreur = "ErrorMessage: la fonction ne converge pas"
        statut = -1
        return [erreur, statut]