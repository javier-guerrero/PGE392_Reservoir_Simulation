# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% Import libraries
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix, csr_matrix
import matplotlib.pyplot as plt
from time import perf_counter, localtime, strftime
#Import functions and depedences from other files
from Get_Arrays import get_arrays                              #Get matrix T, A, Act, Q, J
from Preprocess import Num, Res, PVT, PP, Well, IC, BC, CF     #Import input data
from Initialization import initialization                      #Calculate Pw, Po, Sw and So
from Well_Arrays import Well_Arrays                            #Update arrays when we have wells
from Productivity_Index import Jindex
from Postprocess import plots

#%% Initial conditions
#Po=IC.P+(PVT.rhoosc/144)*(Res.D-IC.Dref)

#%%Initialization
Sw, So, Pw, Po = initialization(Res.D, PP, PVT, Res)

#%% Numerical Method to solve the eq. Ax=b
# theta=Num.one_phase_method
theta=0

#%% Call the function to get arrays        
T, A, Act, Q, J, G=get_arrays(Num, PVT, Res,PP, BC)
J,Q=Well_Arrays(J,Q, Well,Res,PVT,Num)

#%% Solve P from a equation of type Ax=b
def P_SinglePhase_Numerical(Po, Num, Act, T, J, Q, theta, G=0):
    '''
    Function returns Numerical Solution for  Pressure
    input:          Po, Num, Act, T, J, Q, theta, G
    output:         Plots of P a different times, Acumulative production

    Parameters
    ----------
    Po : Initial conitions --> initial pressure
    Num : Numerical inputs
    Act : Accumulation matrix
    T : Transmissibility matrix
    J : Dirichlet BC  and Well conditions matrix
    Q : Rates vector
    theta : Value for solution method (explicit/implicit/CN)
    G : Gravity vector, optional -> The default is 0.

    Returns
    -------
    P_matrix : Solution for Pressure
    time : time

    '''
    start = perf_counter()
    print(strftime("%H:%M:%S", localtime()), ' Simulation started')
    total_time_steps =  int(Num.t_final/Num.dt)    
    time = np.zeros((total_time_steps+1))               # initializing time vector
    P_matrix= np.zeros((Num.N, total_time_steps+1))     # matrix to save P(t) 
    P_matrix[:,0]=Po
    for i in range(1,total_time_steps+1):
        time[i]= time[i-1]+ Num.dt 
        P = np.transpose([spsolve(Act+(1-theta)*(T+J) ,(Act-theta*(T+J))@P_matrix[:,i-1]+Q+G)])   
        P_matrix[:,i] = P[:,0]
        print(strftime("%H:%M:%S", localtime()),' Time', time[i], 'days', 'Simulation percentage', round((time[i])*100/Num.t_final,2), '%', 'Paverage= ', round(np.mean(P_matrix[:,i]),1))
    end = perf_counter() 
    print("Time elapsed during the calculation:", round(end - start,2), 'sec')  
    return P_matrix, time

#%% Call function to calculate P
Num.t_final = 2000
Num.dt = 1
Wellindex = np.where(Q>0)
P, t = P_SinglePhase_Numerical(Po, Num, Act, T, J, Q, theta, G)    

#%%

#%%
plots(P, t, Res, Num, PVT, Well)
    
