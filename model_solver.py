# -*- coding: utf-8 -*-
"""
Created on Wed May 10 09:34:17 2017

@author: jaumep
"""
from __future__ import division

import numpy as np
import os
from scipy import optimize
from model_functions import BD_epsilon, BD_sigma, StreamPlot, QuiverPlot, TimeEvol


k = 4
inf_rate = 4/k


if (os.path.isdir("Rho_plots") == False):
    os.makedirs("Rho_plots")

os.chdir("Rho_plots")

print ("")
print ("")

print ("       *********************************")
print ("       *---------MODEL PLOTS-----------*")
print ("       *********************************")

print ("")
print ("OPTION 1: Computing diagrams for a range of (epsilon,sigma)")
print ("OPTION 2: Computing diagrams for individual values")
print ("")
option = input("  - What option do you choose? : ")
print ("")
while ((option != 1) and (option != 2)):
    print ("")
    print ("ERROR!!! Please enter the value again")
    print ("")
    option = input("  - What option do you choose? : ")

if (option == 1):

    min_epsilon = input()
    max_epsilon = input()
    min_sigma = input()
    max_sigma = input()
    
    Epsilon_Values = np.linspace(min_epsilon,max_epsilon,5)
    Sigma_Values = np.linspace(min_sigma,max_sigma,4)

else:
    print ("  - Enter the parameters you want to work with:")
    epsilon = input("      - epsilon (fitness) = ")
    sigma = input("      - sigma (coupling) = ")

print ("")
print ("    - PLOT 1: Bifurcation Diagram")
print ("    - PLOT 2: Streamplot")
print ("    - PLOT 3: Quiverplot")
print ("    - PLOT 4: Time evolution")
print ("")

plot = input("  - Enter plot type: ")



###############################################################################
#                           BIFURCATION DIAGRAMS                              #
###############################################################################

if (plot == 1):
    if (option == 1):
        
        for sigma in Sigma_Values:
            BD_epsilon(k, inf_rate,sigma)
    
        for epsilon in Epsilon_Values:
            BD_sigma(k, inf_rate,epsilon)
    
    else:
        BD_epsilon(k, inf_rate,sigma)
        BD_sigma(k, inf_rate,epsilon)
    
###############################################################################
###############################################################################
#                                                                             #
#                           PHASE PORTRAIT PLOTS                              #
#                                                                             #
###############################################################################
###############################################################################
elif (plot != 1):
    
    def dRho_dt(Rho):

        den_w = (epsilon*(Rho[0])**sigma + (1-epsilon)*(Rho[1])**sigma)
        w1 = (epsilon*Rho[0]**sigma)/den_w
        w2 = ((1-epsilon)*Rho[1]**sigma)/den_w
    
        f = np.zeros(2)
    
        f[0] = Rho[0]*(inf_rate*k*w1*(1-Rho[0]) - 1)
        f[1] = Rho[1]*(inf_rate*k*w2*(1-Rho[1]) - 1)

        return f
    
#    RESULTS = {}
#    POINT = {}
#    z = 0
    
    
    
#    if (option == 1):
#        
#        for sigma in Sigma_Values:        
#            for epsilon in Epsilon_Values:
#                
#                POINT[z] = [epsilon,sigma]               
#                RESULTS[z] = []
#                
#                for r in range(100):
#                    
#                    guess = [np.random.rand(),np.random.rand()]  
#        #            print('GUESS IS %s' %(guess))
#                    sol,info,ier,mesg = optimize.fsolve(dRho_dt, [guess[0],guess[1]], xtol = 1e-6,full_output = True)
#                    if (ier == 1): 
#                        RESULTS[z].append(list(sol))
#                        z += 1
#                    else:
#                        continue
#        #                print (mesg)
#        
#    else:
#         
#        for r in range(100):
#            
#            POINT[z] = [epsilon,sigma]                 
#            RESULTS[z] = []
#            
#            guess = [np.random.rand(),np.random.rand()]  
#
#            sol,info,ier,mesg = optimize.fsolve(dRho_dt, [guess[0],guess[1]], xtol = 1e-6,full_output = True)
#            if (ier == 1): 
#                RESULTS[z].append(list(sol))
#                z += 1
#            else:
#                continue
            
###############################################################################
#                           STREAMPLOT                                        #
############################################################################### 

    if (plot == 2):
        if (option == 1):
            
            for sigma in Sigma_Values:        
                for epsilon in Epsilon_Values:
                    StreamPlot(k,inf_rate,epsilon,sigma)
        else:
            
            StreamPlot(k,inf_rate,epsilon,sigma)
            
###############################################################################
#                           QUIVER PLOT                                       #
############################################################################### 
        
    elif (plot == 3):
        if (option == 1):
            
            for sigma in Sigma_Values:        
                for epsilon in Epsilon_Values:
                    QuiverPlot(k,inf_rate,epsilon,sigma)
        else:
            
            QuiverPlot(k,inf_rate,epsilon,sigma)
            
###############################################################################
#                           TIME EVOLUTION                                    #
############################################################################### 
    
    elif (plot == 4):
        if (option == 1):
        
            for sigma in Sigma_Values:        
                for epsilon in Epsilon_Values:
                    TimeEvol(k,inf_rate,epsilon,sigma)
        else:
        
            TimeEvol(k,inf_rate,epsilon,sigma)
        