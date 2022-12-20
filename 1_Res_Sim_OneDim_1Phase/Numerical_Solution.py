# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% import packages and libraries
import numpy as np
import math
from scipy.sparse.linalg import spsolve



#%% Functions

def P_1D_1PH_Hom_Numerical(tf, Po, dt, Nx, Act, T, J, Q, theta):
    """
    This function returns the numerical solution to the diffusivity equation using finite differences
    for a 1D, single phase, homogeneous resrvoir
    
    Pressure is a function of distance x and time t --> P(x,t)

    Parameters
    ----------
    tf : Simulation end time --> last time step
    Po : Initial pressure --> value from initial condition
    dt : delta time.
    Nx : Discretization in x-direction --> # of gridblocks in x-dir
    Act : Matrix of Accumulation * total compressibility
    T : Transmissibility matrix.
    J : Boundary condition matrix - Dirichlet.
    Q : Boundary condition vector - Neumann.
    theta : 0=Implicit, 0.5=C-N, 1=Explicit method.

    Returns
    -------
    P_matrix : Pressure in each gridblock in matrix form.

    """
    total_time_steps =  int(tf/dt)    
    time = np.zeros((total_time_steps + 1))             #initializing time vector
    P_matrix = np.zeros((Nx, total_time_steps + 1))     #matrix to save pressure and time 
    P = Po                                              #initializing iterable current pressure
    
    for i in range(total_time_steps + 1):
        Pn = P 
        ##print(np.round(Pn,4))
        P_matrix[:,i] = P[:,0]  
        #Get P^n+1 using a solver from scipy to solve an eq of the form Ax=b
        P = np.transpose([spsolve(Act + (1 - theta) * (T+J), (Act - theta * (T + J)) @ Pn + Q)])              
        if i!=total_time_steps:                     #Save the results at time n                               
            time[i+1] = time[i] + dt                #time counter
   
    return P_matrix
