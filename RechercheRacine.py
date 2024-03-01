
import math
import numpy as np

def bissection(f, x0, x1, tol):
    # status verification (hypotheses)
    statut = 0
    half_ecart = (x1-x0)/2
    y0 = f(x0)
    y1 = f(x1)

    # cas x0 après x1
    if x1 < x0:
        x0, x1 = x1, x0
    
    # cas fonction décroissante
    if(y0>y1):
        x0, x1 = x1, x0
    
    if y0*y1 > 0:
        erreur = "ErrorMessage: no root in the function"
        statut = 1
        return [erreur, statut]

    elif tol <= 10**(-10):
        statut = 1
        return[None,statut]

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

def secante(f, x0, x1, tol) :
    statut = 0

    # Hypothesis verification
    y0 = f(x0)
    y1 = f(x1)
    k = 0

    if tol <= 10**(-10):
        erreur = "ErrorMessage: invalid tolerance"
        statut = 1
        return [erreur, statut]
    
    if abs(y0 - y1) <= 10**(-10):
        erreur = "ErrorMessage: point unique"
        statut = 1
        return [erreur, statut]
    
    N = int(np.log(abs(x1-x0)/(2*tol)))/(np.log(2)) + 1
    while abs(y1) > tol :
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
            
    statut = 0
    
    return [x1, statut]

import math
def f(x):
    return x**3 - math.cos(x)

x0 = 0.4
x1 = 1
print(bissection(f, x0, x1, 0.01))
print(secante(f, x0, x1, 0.01))

import matplotlib.pyplot as plt
absi = np.linspace(0, 1.5, 100)
plt.plot(absi, [f(x) for x in absi])
#plt.show()