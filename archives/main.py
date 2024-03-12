
from RechercheRacine import secante, bissection
import numpy as np

def f(x):
    #return x**3 - x**2 -5*x + 2
    return x**4 - x**3 - 3*x**2

# Calculating roots of f
x0 = 2
x1 = 3
tol = 0.01

x, statut = secante(f, x0, x1, tol)
print(f'{x}')

x, statut = bissection(f, x0, x1, tol)
print(f'{x}')
