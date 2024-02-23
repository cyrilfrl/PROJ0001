
from math import*
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
        
        if np.sign(y0) == np.sign(y):
            x0 += half_ecart
            y0 = y
        
        elif np.sign(y1) == np.sign(y):
            x1 -= half_ecart
            y1 = y

        # checking if reached tolerance threshold
        half_ecart = (x1-x0)/2 # coordinates between x0 and x1
        if abs(f(x0 + half_ecart)) <= tol:
            break

    return [x0 + half_ecart, statut]

            # -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 13:56:19 2024

@author: Leon Carmona
"""

def secante(f ,x0, x1, tol):
    
    statut = 0      #initialisation par defaut, aucune erreur
    y0 = f(x0)
    y1 = f(x1)
    values = []
    #il reste le cas ou y1*y2 = 0

    #cas ou x1>x0
    if x1<x0:
        x0, x1 = x1, x0
         
    #cas ou x0 est une racine
    if abs(y0) < tol : 
        return [x0 , statut]
    
    #cas ou x1 est une racine
    elif abs(y1) < tol :
        return [x1 , statut]
    
    #cas ou il n'y aucune racine (deux ordonnees sont du meme np.signe)
    elif y0*y1 > 0 :    
        statut = -1
        erreur = "Il n'y a aucune racine dans cet intervalle"
        return [erreur, statut]
         

    #cas ou il y a une racine (au moins une des deux ordonnees est negative)
    elif y0*y1 < 0 :    
        
        while abs(y0) > tol: #verifier ceci!!
            
            #si x0 est inferieur a x1! donc il faut prendre le nouvel intervalle a droite en changeant x0
            x_suiv = x1 - ( (y1*(x1-x0)) / (y1 - y0) )
            x0 = x1
            y0 = f(x0)
            x1 = x_suiv
            y1 = f(x1)
            values.append(x0)
        
            if abs(y0) < tol:
                return [x0, statut]

