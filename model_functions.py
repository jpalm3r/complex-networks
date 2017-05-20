# -*- coding: utf-8 -*-
"""
Created on Wed May 10 18:18:23 2017

@author: jaumep
"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import os
from scipy import integrate, optimize
import plotly.plotly as py
import plotly.offline as pyoff
import plotly.figure_factory as ff
import plotly.graph_objs as go
import sys

def BD_epsilon(k,inf_rate,coupling):

    RhoA_Values = np.linspace(0.01,0.99,100)
    Rho2A_Values = [(((rho1-1)/(inf_rate*k*(1-rho1)-1)) + 1) for rho1 in RhoA_Values]  
    
    
    if (os.path.isdir("Epsilon_plots") == False):
        os.makedirs("Epsilon_plots")
    os.chdir("Epsilon_plots")
        
    Epsilon_Values = []
    
    for rho_a in RhoA_Values:
        
        Om = (((rho_a-1)*(1-inf_rate*k)-1)**coupling)/(((rho_a)**coupling)*(inf_rate*k*(1-rho_a)-1)**(coupling+1))
        f = Om/(1+Om)
        Epsilon_Values.append(f)
     
    plt.ylabel(r"$\rho_{i}^{a *}$")
    plt.xlabel(r"$\epsilon$")
    plt.xlim(0.3,0.7)
    plt.ylim(0.0,1.0)
    title = r"$\sigma$ = " + str(coupling)[:3]
    plt.title(title)
    
    plt.plot(Epsilon_Values,RhoA_Values, 'b', label = r"$\rho_{1}^{a *}$", linewidth = 1.5)
    plt.plot(Epsilon_Values,Rho2A_Values, 'r', label = r"$\rho_{2}^{a *}$", linewidth = 1.5)
    
    
    plt.legend(bbox_to_anchor=(0.98, 0.5), loc=5, borderaxespad=0.)
    
    fname = "rho_vs_epsilon_cc-" + str(coupling)[:4] +".png"
    plt.savefig(fname, dpi = 500)
    plt.close()
    
    
    os.chdir("../")
 
###############################################################################
#                          ------------------                                 #
###############################################################################  
   
def BD_sigma(k,inf_rate,epsilon):
    
    RhoA_Values = np.linspace(0.01,0.99,100)
    RhoA_List = list(RhoA_Values)
    
    if (os.path.isdir("Sigma_plots") == False):
        os.makedirs("Sigma_plots")
    os.chdir("Sigma_plots")
           
    Coupling_Values = []
    
    for rho_a in RhoA_Values:
        
        comp_log_num = (inf_rate*k*(1-rho_a) - 1)*(epsilon/(1-epsilon))        
        comp_log_den = ((inf_rate*k - 1)*(1-rho_a) - 1)/(rho_a*(inf_rate*k*(1-rho_a) - 1))
        num = np.log(comp_log_num)
        den = np.log(comp_log_den)
        f = num/den
        if(((comp_log_num > 0) and (comp_log_num > 0)) and (f > 0)):
            
            Coupling_Values.append(f)
        else:
            RhoA_List.remove(rho_a)
        
    Rho2A_Values = [(((rho1-1)/(inf_rate*k*(1-rho1)-1)) + 1) for rho1 in RhoA_List]

    
        
    
    plt.ylabel(r"$\rho_{i}^{a *}$")
    plt.xlabel(r"$\sigma$")
    
    plt.xlim(0.5,2)
    plt.ylim(0.0,1)
    title = r"$\epsilon$ = " + str(epsilon)[:4]
    plt.title(title)
    
    plt.plot(Coupling_Values,RhoA_List, 'c', label = r"$\rho_{1}^{a *}$", linewidth = 1.5)
    plt.plot([0 for x in range(len(RhoA_List))],RhoA_List, 'c', linewidth = 1.5)
    plt.plot(Coupling_Values,Rho2A_Values, 'm', label = r"$\rho_{2}^{a *}$", linewidth = 1.2)
    
    plt.legend(bbox_to_anchor=(0.98, 0.15), loc=5, borderaxespad=0.)
    
    fname = "rho_vs_sigma_cc-" + str(epsilon)[:4] +".png"
    plt.savefig(fname, dpi = 500)
    plt.close()
    
    os.chdir("../")

###############################################################################
#                          ------------------                                 #
###############################################################################  
    
def StreamPlot(k,inf_rate,epsilon,sigma): #NOT WORKING !!!!!
    
    # 1) Find system solutions    
    
    def dRho_dt(Rho, t=0):
    
        den_w = (epsilon*(Rho[0])**sigma + (1-epsilon)*(Rho[1])**sigma)
        w1 = (epsilon*Rho[0]**sigma)/den_w
        w2 = ((1-epsilon)*Rho[1]**sigma)/den_w
    #        Rho = [rho1,rho2]
        return np.array([Rho[0]*(inf_rate*k*w1*(1-Rho[0]) - 1),
                      Rho[1]*(inf_rate*k*w2*(1-Rho[1]) - 1)])

    
    RESULTS = []    
    z = 0
    one_solution = False
    
    x_total = 10
    y_total = 10
    x0 = np.linspace(0,1,x_total)
    y0 = np.linspace(0,1,y_total)
    
    for r1 in range(10):
        for r2 in range(10):
                              
            guess = [x0[r1],y0[r2]]  
    
            sol,info,ier,mesg = optimize.fsolve(dRho_dt, guess, xtol = 1e-6,full_output = True)
            
            sol_l = list(sol)
            
            if (ier == 1):                
                if (not one_solution):
                    one_solution = True
                    RESULTS.append(sol_l)
                else:
                    solution_found = False
                        
                    for j in RESULTS:
                        solution_found = (abs(j[0]-sol_l[0]) <= 0.001) and (abs(j[1]-sol_l[1]) <= 0.001)
                        
                    if(not solution_found):
                        RESULTS.append(sol_l)
                        z += 1
                    else:
                        continue
    # 2) Draw Streamplot
    
    N = 50
    rho1_start, rho1_end = 0.01, 1.0
    rho2_start, rho2_end = 0.01, 1.0
    rho1 = np.linspace(rho1_start, rho1_end, N)
    rho2 = np.linspace(rho2_start, rho2_end, N)
    RHO1, RHO2 = np.meshgrid(rho1, rho2)
    
    den_w = (epsilon*(RHO1)**sigma + (1-epsilon)*(RHO2)**sigma)
    w1 = (epsilon*RHO1**sigma)/den_w
    w2 = ((1-epsilon)*RHO2**sigma)/den_w
    #Rho = [rho1,rho2]
    U = RHO1*(inf_rate*k*w1*(1-RHO1) - 1)
    V = RHO2*(inf_rate*k*w2*(1-RHO2) - 1)
    
    rho_on_axis = (inf_rate*k - 1)/(inf_rate*k)
    
    # Add source point
    rho1_axis = go.Scatter(x=[rho_on_axis], y=[0],
                              mode='markers',
                              marker=go.Marker(size=14,symbol='diamond',color='red'))
                              
    rho2_axis = go.Scatter(x=[0], y=[rho_on_axis],
                              mode='markers',
                              marker=go.Marker(size=14,symbol='diamond',color='red')) 
    
    fig = ff.create_streamline(rho1, rho2, U, V, density=1.2, arrow_scale=.01)
    
    # Add source point to figure
    fig['data'].append(rho1_axis)
    fig['data'].append(rho2_axis) 
    
    for res in RESULTS:
        fixed_point = go.Scatter(x=[res[0]], y=[res[1]],
                              mode='markers',
                              marker=go.Marker(size=14,symbol='circle',color='red'))
        fig['data'].append(fixed_point) 
    
    # Image returned in html and open in browser
    pyoff.plot(fig, filename='Streamline.html')

###############################################################################
#                          ------------------                                 #
###############################################################################  

def FindFP(k,inf_rate,epsilon,sigma, guess0):
    
    def dRho_dt(Rho):
        
        den_w = (epsilon*(Rho[0])**sigma + (1-epsilon)*(Rho[1])**sigma)
        w1 = (epsilon*Rho[0]**sigma)/den_w
        w2 = ((1-epsilon)*Rho[1]**sigma)/den_w

        f = np.zeros(2)
        
        f[0] = Rho[0]*(inf_rate*k*w1*(1-Rho[0]) - 1)
        f[1] = Rho[1]*(inf_rate*k*w2*(1-Rho[1]) - 1)

        return f
    
#    def Jacob(Rho):
#        
#        den_w = (epsilon*(Rho[0])**sigma + (1-epsilon)*(Rho[1])**sigma)
#        w1 = (epsilon*Rho[0]**sigma)/den_w
#        w2 = ((1-epsilon)*Rho[1]**sigma)/den_w
#        
#        dw1_drho1 = (sigma*epsilon*(1-epsilon)*(Rho[0]**(sigma-1))*Rho[1]**sigma)/(epsilon*Rho[0]**sigma + (1-epsilon)*Rho[1]**sigma)**2 
#        dw1_drho2 = (-epsilon*sigma*(Rho[0]**sigma)*(1-epsilon)*sigma*Rho[1]**(sigma-1))/(epsilon*Rho[0]**sigma + (1-epsilon)*Rho[1]**sigma)**2 
#        dw2_drho1 = (-epsilon*sigma*(Rho[1]**sigma)*(1-epsilon)*sigma*Rho[0]**(sigma-1))/(epsilon*Rho[0]**sigma + (1-epsilon)*Rho[1]**sigma)**2 
#        dw2_drho2 = (sigma*epsilon*(1-epsilon)*(Rho[1]**(sigma-1))*Rho[0]**sigma)/(epsilon*Rho[0]**sigma + (1-epsilon)*Rho[1]**sigma)**2 
#        
#        a11 = inf_rate*k*Rho[0]*(1-Rho[0])*dw1_drho1 + inf_rate*k*w1*(1-2*Rho[0])-1
#        a12 = inf_rate*k*Rho[0]*(1-Rho[0])*dw1_drho2
#        a21 = inf_rate*k*Rho[1]*(1-Rho[1])*dw2_drho1        
#        a22 = inf_rate*k*Rho[1]*(1-Rho[1])*dw2_drho2 + inf_rate*k*w2*(1-2*Rho[1])-1
#        
#        J = [[],[]]
#        J[0] = [a11, a12]
#        J[1] = [a21, a22]         
#        
#        return J
    
    sol, ier, msg = optimize.fsolve(dRho_dt,[guess0[0],guess0[1]])
       
    
    return (sol, ier, msg)

###############################################################################
#                          ------------------                                 #
###############################################################################  

def TimeEvol(k,inf_rate,epsilon,sigma):
    
    def dRho_dt(Rho, t=0):
    
        den_w = (epsilon*(Rho[0])**sigma + (1-epsilon)*(Rho[1])**sigma)
        w1 = (epsilon*Rho[0]**sigma)/den_w
        w2 = ((1-epsilon)*Rho[1]**sigma)/den_w
    #        Rho = [rho1,rho2]
        return np.array([Rho[0]*(inf_rate*k*w1*(1-Rho[0]) - 1),
                      Rho[1]*(inf_rate*k*w2*(1-Rho[1]) - 1)])
                        
    
    t = np.linspace(0,30,1000)
    Rho02 = [0.9,0.4] # Initial conditions
    Rho2 = integrate.odeint(dRho_dt, Rho02, t)    
    
    rho1, rho2 = Rho2.T # .T means the transposed
    
    f1 = plt.figure()
    plt.plot(t, rho1, 'r-', label=r'$\rho_{1}^{a}$')
    plt.plot(t, rho2, 'b-', label=r'$\rho_{2}^{a}$')
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.xlim(0,15)
    plt.ylim(0,1)
    plt.ylabel(r'$\rho_{i}^{a}$')
    plt.title('Evolution of network density. ($\sigma = %s; \epsilon = %s $)' %(str(sigma)[:4],str(epsilon)[:4]))
    
    f1.savefig('evolution1_e-%s_s-%s.png' %(str(epsilon)[:4],str(sigma)[:4]), dpi = 500)

###############################################################################
#                          ------------------                                 #
###############################################################################    

def FilterSolutions(RESULTS,z):
    
    SOLUTIONS = []
    if (z != 0):
        SOLUTIONS.append(RESULTS[0])
        for i in range(1,z-1):
            sol = RESULTS[i]
            for j in SOLUTIONS:            
                if ((abs(j[0]-sol[0]) <= 0.001) and (abs(j[1]-sol[1]) <= 0.001)):
                    SOLUTIONS.append(sol)
    
    return SOLUTIONS
    
###############################################################################
#                          ------------------                                 #
###############################################################################     

def QuiverPlot(k,inf_rate,epsilon,sigma):

    def dRho_dt(Rho, t=0):
    
        den_w = (epsilon*(Rho[0])**sigma + (1-epsilon)*(Rho[1])**sigma)
        w1 = (epsilon*Rho[0]**sigma)/den_w
        w2 = ((1-epsilon)*Rho[1]**sigma)/den_w
    #        Rho = [rho1,rho2]
        return np.array([Rho[0]*(inf_rate*k*w1*(1-Rho[0]) - 1),
                      Rho[1]*(inf_rate*k*w2*(1-Rho[1]) - 1)])

    
    RESULTS = []    
    z = 0
    one_solution = False
    
    x_total = 10
    y_total = 10
    x0 = np.linspace(0,1,x_total)
    y0 = np.linspace(0,1,y_total)
    
    for r1 in range(10):
        for r2 in range(10):
                              
            guess = [x0[r1],y0[r2]]  
    
            sol,info,ier,mesg = optimize.fsolve(dRho_dt, guess, xtol = 1e-6,full_output = True)
            
            sol_l = list(sol)
            
            if (ier == 1):                
                if (not one_solution):
                    one_solution = True
                    RESULTS.append(sol_l)
                else:
                    solution_found = False
                        
                    for j in RESULTS:
                        solution_found = (abs(j[0]-sol_l[0]) <= 0.001) and (abs(j[1]-sol_l[1]) <= 0.001)
                        
                    if(not solution_found):
                        RESULTS.append(sol_l)
                        z += 1
                    else:
                        continue
                                
#    SOLUTIONS = FilterSolutions(RESULTS,z)    

    t = np.linspace(0,30,1000)
    
    Rho01 = [0.1,0.1]
    Rho02 = [0.9,0.4]
    Rho03 = [0.4,0.7]
    Rho04 = [0.6,0.2]

    Rho1 = integrate.odeint(dRho_dt, Rho01, t)
    Rho2 = integrate.odeint(dRho_dt, Rho02, t)
    Rho3 = integrate.odeint(dRho_dt, Rho03, t)
    Rho4 = integrate.odeint(dRho_dt, Rho04, t)

    
    f2 = plt.figure()
    
    # plot trajectorY       
    plt.plot(Rho1[:,0], Rho1[:,1], lw=1.3,color='r')
    plt.plot(Rho01[0],Rho01[1],marker='x',color='r')
    plt.plot(Rho2[:,0], Rho2[:,1], lw=1.3,color='r')
    plt.plot(Rho02[0],Rho02[1],marker='x',color='r')    
    plt.plot(Rho3[:,0], Rho3[:,1], lw=1.3,color='r')
    plt.plot(Rho03[0],Rho03[1],marker='x',color='r')
    plt.plot(Rho4[:,0], Rho4[:,1], lw=1.3,color='r')    
    plt.plot(Rho04[0],Rho04[1],marker='x',color='r')
    
    rho_on_axis = (inf_rate*k - 1)/(inf_rate*k)
    
    plt.plot([0],[rho_on_axis],marker='D',color='r')
    plt.plot([rho_on_axis],[0],marker='D',color='r')
    for fp in RESULTS:
        plt.plot(fp[0],fp[1],marker='D',color='r')
    
    
    x = np.linspace(0,1,20)
    y = np.linspace(0,1,20)

    X1, Y1 = np.meshgrid(x,y) # create a grid
    DX1, DY1 = dRho_dt([X1,Y1]) # compute growth rate on the gridt
    M = (np.hypot(DX1, DY1)) # Norm of the growth rate 
    M[M == 0] = 1. # Avoid zero division errors 
    DX1 /= M # Normalize each arrows
    DY1 /= M
    
    plt.title(r'Trajectories and direction fields. ($\sigma = %s; \epsilon = %s$)' %(str(sigma)[:4],str(epsilon)[:4]))
    plt.quiver(X1, Y1, DX1, DY1, M, pivot='mid', cmap=plt.cm.ocean)
    plt.xlabel(r'$\rho_{1}^{a}$')
    plt.ylabel(r'$\rho_{2}^{a}$')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    f2.savefig('evolution2_e-%s_s-%s.png' %(str(epsilon)[:4],str(sigma)[:4]), dpi = 500)

###############################################################################
#                          ------------------                                 #
###############################################################################  
    
def W(k,inf_rate,epsilon,sigma):
    
    def dRho_dt(Rho, t=0):
    
        den_w = (epsilon*(Rho[0])**sigma + (1-epsilon)*(Rho[1])**sigma)
        w1 = (epsilon*Rho[0]**sigma)/den_w
        w2 = ((1-epsilon)*Rho[1]**sigma)/den_w
    #        Rho = [rho1,rho2]
        return np.array([Rho[0]*(inf_rate*k*w1*(1-Rho[0]) - 1),
                      Rho[1]*(inf_rate*k*w2*(1-Rho[1]) - 1)])
    
    t = np.linspace(0,30,1000)
    Rho01 = [0.1,0.3]
    Rho = integrate.odeint(dRho_dt, Rho01, t)
    rho1, rho2 = Rho.T # .T means the transposed
    
    w1 = (epsilon*rho1**sigma)/(epsilon*rho1 + (1-epsilon)*rho2**sigma)    

    plt.plot(rho1,w1)
    plt.plot(rho2,w1)
    plt.show()
