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
from scipy.sparse.linalg import inv
import matplotlib.pyplot as plt
from time import perf_counter, localtime, strftime

#Import functions and  from other files
from BCM_kr import BCM_RelativePermeability, CapillaryPressure
from Preprocess import Num, Res, PVT, petro, Well, IC, BC, CF  # Import input data
from Initialization import initialization                      # Calculate Pw, Po, Sw and So
from T_interblock import T_interblock                          # Interblock transmissibility
from Get_Arrays import grid_arrays                              # Get matrix T, A, Act, Q, J
from Get_Arrays_SS import get_arrays_SS                        # Get matrix T, A, Act, Q, J
from Productivity_Index import Jindex                          # Productivity index
from Well_Arrays import Well_Arrays, Well_Arrays_SS            # Update arrays when we have wells
 

#%% Initialization
Sw, So, Pw, Po, rhoo, rhow = initialization(Res.D, petro, PVT, Res)
PVT.rhoo = rhoo
PVT.rhow = rhow

#Calculate relative permeability
krw, kro = BCM_RelativePermeability(Sw, petro)
Swo = np.reshape(Sw, (Num.N,1))
Po = np.reshape(Po, (Num.N,1))
#
Num.t_initial = 0
Num.t_final = 5000


#%% Solve P from a equation of type Ax=b
start = perf_counter()    
print(strftime("%H:%M:%S", localtime()), ' Simulation starts')
total_time_steps =  int((Num.t_final - Num.t_initial)/Num.dt)    
time = np.zeros((total_time_steps+1))                # initializing time vector
P_matrix = np.zeros((Num.N, total_time_steps+1))     # matrix to save P(t) 
Sw_matrix = np.zeros((Num.N, total_time_steps+1))    # matrix to save Sw(t)
Qw_matrix = np.zeros((Num.N, total_time_steps+1))    # matrix to save Sw(t)
Sw_hyst = np.empty((Num.N,2))
Sw_hyst[:,0] = Swo[:,0]


#Initial Pressure and initial Saturation
P_matrix[:,0] = Po[:,0] 
Pn = Po
Sw_matrix[:,0] = Swo[:,0]
Swn = Swo
 
for i in range(1, total_time_steps+1):
    time[i] = time[i-1]+ Num.dt 
    if Num.MF_method == 'IMPES':
        T, Tw, To, A, Act, Q, J, Jw, Jo, G, BetaW, Qw, Qo, G1, Pc, d11, d12 = grid_arrays(Num, PVT, Res, petro, BC, Pn, Swn, Sw_hyst)
        J, Jw, Jo, Q, Qw, Qo = Well_Arrays(Swn, J, Jw, Jo, Q, Qw, Qo, Well, Res, PVT, Num, petro)
        Q = Qw + Qo + 800.0 * J @ np.ones((Num.N, 1))
        P = np.transpose([spsolve(Act+(T+J) , Act@Pn+Q+G)]) 
        Pc = np.reshape(np.asarray(Pc),(Num.N,1))
        Sw = (Swn + inv(A) * (-Tw@(P - Pc - (PVT.rhow / 144) * Res.D.reshape(Num.N,1)) + Qw + Jw @(800-P)))-(BetaW * Swn) * (PVT.cw + PVT.cf)@(P - Pn)
        P_matrix[:,i] = P[:,0]
        Qw_matrix[:,i] = np.reshape(Qw.toarray(), (Num.N))
        Sw_matrix[:,i] = Sw.T
        #print(Pc)
        Pn = P
        Swn = np.asarray(Sw)
        #To calculate hysteresis
        for j in range(0, Num.N):
            if Sw_matrix[j,i] > Sw_matrix[j,i-1] and Sw_hyst[j,1] == 0:  # [i,1] is a flag
                Sw_hyst[j,0] = Sw_matrix[j,i]
                Sw_hyst[j,1] = 1.0

            elif Sw_matrix[j,i] < Sw_matrix[j,i-1]:
                Sw_hyst[j,0] = Sw_matrix[j,i]    
        print(strftime("%H:%M:%S", localtime()),' Time', time[i], 'days', 'Simulation Progress:', round((time[i])*100/(Num.t_final-Num.t_initial),2), '%')
end = perf_counter() 
print("Time elapsed during the calculation:", round(end - start,2), 'sec')  
J=J.toarray()
# Q=Q.toarray()
Qw=Qw.toarray()
Qo=Qo.toarray()
G=G.toarray()
To=To.toarray()
Tw = Tw.toarray()
T = T.toarray()