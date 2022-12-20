# -*- coding: utf-8 -*-
"""
PGE-392K-  Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

"""
ASSIGNMENT:
Develop a reservoir simulator to solve for flow in a 1D, homogeneous reservoir 
using a finite difference method that is second order accurate in space. 
The program should allow the user to define properties in an input ‘.yml’ file.
The program should also be flexible enough to allow for different boundary conditions 
and numerical methods (explicit, implicit, Crank-Nicholson)
"""

#%% import packages and libraries
import numpy as np
from Preprocess import Num, Res, PVT, petro, Well, IC, BC, CF #Import input data
from scipy.sparse.linalg import spsolve
from Grid_Arrays import Create_Grid_Arrays  #Create arrays to solve equations
from Postprocess import Plot_Analytical_Numerical   #To plot the solutions numerical and analytical

#%% Function to calculate P  ==> PROBLEM 2
def P_1D_SinglePhase_Homogeneous_Numerical(tf, Po, dt, Nx, Act, T, J, Q, theta):
    total_time_steps =  int(tf/dt)    
    time = np.zeros((total_time_steps+1))           #initializing time vector
    P_matrix= np.zeros((Nx, total_time_steps+1))    #matrix to save pressure and time 
    P= Po                                           #initializing iterable current pressure
    for i in range(total_time_steps+1):
        Pn =P 
        ##print(np.round(Pn,4))
        P_matrix[:,i] = P[:,0]  
        #Get P^n+1 using a solver from scipy to solve an eq of the form Ax=b
        P = np.transpose([spsolve(Act+(1-theta)*(T+J) ,(Act-theta*(T+J))@Pn+Q)])              
        if i!=total_time_steps:                     #Save the results at time n                               
            time[i+1] = time[i] + dt                #time counter
    return P_matrix

#%% CHANGE INPUT parameters for PROBLEMS 3 TO 5
'''
PROBLEM 2 - Edit inputs to match Example 3.4 and validate Explicit, Implicit and C-N solutions w/ Ex 3.4 & 3.5
PROBLEM 3 - Increase number of block N=25 --> num['N']=25 and choose a time step stable for all methods.
            Compare to ANALYTICAL solution. Plot PD vs xD for tD=0.25
PROBLEM 4 - Repeat Problem 3 but use a timestep that is 5% larger  than value required for stability of Explicit method
PROBLEM 5 - Repeat Problem 3, but change the Boundary conditions to Dirichlet-Neumann and Dirichlet-Dirichlet where the
            Dirichlet BC on the right equials initial condition
UNCOMMENT/COMMENT THE NEXT CODE LINES DEPENDING ON THE PRBLEM TO SOLVE
'''

xD=np.linspace(Num.dx/2, Res.L-Num.dx/2, Num.N)/Res.L

tD=0.25
Num.tf=tD/((Res.k*petro.kroo/(Res.phi*PVT.muo*Res.ct))*(1/Res.L**2)*CF.mDft_cp_to_ft3_psiday)
print(Num.tf)
Num.dt=0.2

#Get Parameters A and T for matrixes
Accumulation=(Res.w*Res.h)*Res.phi*Num.dx/Num.dt
Transmissibility=(Res.k*petro.kroo*(Res.w*Res.h)/(PVT.muo*Num.dx))*CF.mDft_cp_to_ft3_psiday

#Call Grid_Arrays to get the matrix and arrays for the solution
T, J, Q, A, ct, Act, Po = Create_Grid_Arrays(Accumulation, Transmissibility, Res.ct, Num.N, BC.type, BC.value, IC.P)

#%% Warning for stability criteria when explicit method
etha=(Res.k*petro.kroo/(Res.phi*PVT.muo*Res.ct))*(Num.dt/Num.dx**2)*CF.mDft_cp_to_ft3_psiday
if etha > 0.5:
    print('**Warning!   For Explicit Method** \n Stability cruiteria etha<=0.5 \n Current value etha= ', etha)
    dt=0.5/((Res.k*petro.kroo/(Res.phi*PVT.muo*Res.ct))*(1/Num.dx**2)*CF.mDft_cp_to_ft3_psiday)
    print('dt should be <=', np.round(dt,4)) #Suggested dt value at defined  dx
    

#%% SOLUTION TO PROBLEMS 3 TO 5

# SOLUTION TO PROBLEM 3 --> Increase Nx to 25 blocks
Analyticaltype='Dirichlet-Neumann'; BC.type= [['Neumann'],['Dirichlet']]; BC.value = [0,1000];  
fig1 = Plot_Analytical_Numerical(1,xD, tD, BC.type, BC.value, IC.P, Analyticaltype, Accumulation, Transmissibility, Res.ct, Num.N, Num.dt, Num.tf, etha)

# SOLUTION TO PROBLEM 5 --> Modify Boundary Conditions Configuration
Analyticaltype='Dirichlet-Neumann'; BC.type= [['Dirichlet'],['Neumann']]; BC.value = [1000,0];  
fig2 = Plot_Analytical_Numerical(2,xD, tD, BC.type, BC.value, IC.P, Analyticaltype, Accumulation, Transmissibility, Res.ct, Num.N, Num.dt, Num.tf, etha)

Analyticaltype='Dirichlet-Dirichlet'; BC.type= [['Dirichlet'],['Dirichlet']]; BC.value = [1000,3000];  
fig3 = Plot_Analytical_Numerical(3,xD, tD, BC.type, BC.value, IC.P, Analyticaltype, Accumulation, Transmissibility, Res.ct, Num.N, Num.dt, Num.tf, etha)

# SOLUTION TO PROBLEM 4 --> Increase tD to make Eplicit method unstable
Num.dt=0.25
etha=(Res.k*petro.kroo/(Res.phi*PVT.muo*Res.ct))*(Num.dt/Num.dx**2)*CF.mDft_cp_to_ft3_psiday
#Get Parameters A and T for matrixes
Accumulation=(Res.w*Res.h)*Res.phi*Num.dx/Num.dt
Transmissibility=(Res.k*petro.kroo*(Res.w*Res.h)/(PVT.muo*Num.dx))*CF.mDft_cp_to_ft3_psiday

#Call Grid_Arrays to get the matrix and arrays for the solution
T, J, Q, A, ct, Act, Po = Create_Grid_Arrays(Accumulation, Transmissibility, Res.ct, Num.N, BC.type, BC.value, IC.P)

Analyticaltype='Dirichlet-Neumann'; BC.type= [['Neumann'],['Dirichlet']]; BC.value = [0,1000];  
fig4 = Plot_Analytical_Numerical(4,xD, tD, BC.type, BC.value, IC.P, Analyticaltype, Accumulation, Transmissibility, Res.ct, Num.N, Num.dt, Num.tf, etha)



