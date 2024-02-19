# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 13:56:19 2024

@author: Leon Carmona
"""

def secante(f ,x0, x1, tol):
    
    statut = 0      #initialisation par defaut, aucune erreur
    y0 = f(x0)
    y1 = f(x1)
    resultat = [None, None]

  #il reste le cas ou x1>x0, et si y1*y2 = 0

     #cas ou x1>x0
    if x1>x0:
         a = x0
         x0 = x1
         x1 = a
    
    #cas ou x0 est une racine
    elif abs(y0) < tol : 
        racines = [x0 , statut]
        return racines
    
    #cas ou x1 est une racine
    elif abs(y1) < tol : 
        racines = [x1 , statut]
        return racines
    
    #cas ou il n'y aucune racine (deux ordonnees sont du meme signe)
    elif y0*y1 > 0 :
        statut = -1
        erreur = "Il n'y a aucune racine dans cet intervalle"
        return [print(erreur), statut]
         

    #cas ou il y a une racine (au moins une des deux ordonnees est negative)
    elif y0*y1 < 0 :    
        
        while abs(y0) > tol: #verifier ceci!!
            
            #si x0 est inferieur a x1! donc il faut prendre le nouvel intervalle a droite en changeant x0
            x_suiv = x1 - ( (y1*(x1-x0)) / (y1 - y0) )
            x0 = x1
            y0 = f(x0)
            x1 = x_suiv
            y1 = f(x1)
        
            if abs(y0) < tol:
                resultat = [x0, statut]
                return print(resultat) #statut = 0 par defaut
                break 
            
            
                    
            
                
            
            
    
        
        
        
        
                
        
       
        
        
    
    
        
   
        
        
    
    
    
    
    
   






        #x, statut = secante(f, x0, x1, tol)







