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
        print("size issue")
        x0, x1 = x1, x0
         
    #cas ou x0 est une racine
    if abs(y0) < tol : 
        return [x0 , statut]
    
    #cas ou x1 est une racine
    elif abs(y1) < tol :
        return [x1 , statut]
    
    #cas ou il n'y aucune racine (deux ordonnees sont du meme signe)
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
                return [x0, statut], values #statut = 0 par defaut
            
            
                    
            
                
            
            
    
        
        
        
        
                
        
       
        
        
    
    
        
   
        
        
    
    
    
    
    
   






        #x, statut = secante(f, x0, x1, tol)







