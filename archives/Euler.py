# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 11:22:03 2024

@author: Leon Carmona
"""

import numpy as np
from main import odefunction as f

#la fonction consideree est dC/dz
    
def calculConcentrationsEuler(z0,zf,C0):
    
    #definition du pas (à fixer plus tard)
    
    h = 0.5 
    n = (zf - z0) / h #le nombre de points 
    
    #initialisation - condition initiale
    
    C0 = np.zeros(8)
    C = np.array #deja defini dans odefunction
    
    y0 = f(z0,C0)
    
    #calcul des approximations des points suivants
    
    for i in C:
        
        z_i = z0 + i*h
        y_i = f(z_i,C)
 
        y_1 = y0+h*f(z_i,y_i)
    
    
    
    
    

