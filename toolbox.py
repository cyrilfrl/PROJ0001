
from math import*
import matplotlib.pyplot as plt

def sign(x):
    return x >= 0

def bissection(f, x0, x1, tol):
    # cas x0 après x1
    if x0 > x1:
        x0, x1 = x1, x0
    
    valeurs = []
    #//# if x0 > x1 x1 = x0 et x0 = x1
    # status verification (hypotheses)
    statut = 0
    half_ecart = (x1-x0)/2
    y0 = f(x0)
    y1 = f(x1)
    if y0*y1 > 0:
        erreur = "laidron avait raison"
        statut = 1
        return [erreur, statut]
    elif y0*y1 == 0:
        if y0 == 0:
            answer = x0
        elif y1 == 0:
            answser = x1

        return [0, answer]
    
    # max iterations
    #//#k = ln((b-a)/2epsilon)/ln(2)
    
    while True:
        y = f(x0+half_ecart)
        
        if sign(y0) == sign(y):
            x0 += half_ecart
            y0 = y
        
        elif sign(y1) == sign(y):
            x1 -= half_ecart
            y1 = y

        # checking if reached tolerance threshold
        half_ecart = (x1-x0)/2 # coordinates between x0 and x1
        valeurs.append(f(x0 + half_ecart))
        if abs(f(x0 + half_ecart)) <= tol:
            break

    plt.plot(valeurs)
    plt.show()
    return [x0 + half_ecart, statut]

def secante(f, x0, x1, tol):
    pass
