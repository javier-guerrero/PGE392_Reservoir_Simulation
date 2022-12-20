# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""
#%% import packages and libraries
import numpy as np

#%% Functions
def Create_Grid_Arrays(Accumulation, Transmissibility, ct, Nx, BC_type, BC_value, IC):
    """
    This function takes the output of the PREPROCESS  function and creates the required arrays (matrix / vectors)
    to pass them as inputs for the solution methods (analytical or numerical)

    Parameters
    ----------
    Accumulation : Accumulation term 
    Transmissibility : Transmissibility 
    ct : Total compressibility 
    Nx : # of gridblocs in x-dir
    BC_type : boundary configuration
    BC_value : values for the boundary conditions
    IC : Initial COndition

    Returns
    -------
    T : Transmissibility matrix.
    J : Boundary condition matrix - Dirichlet.
    Q : Boundary condition vector - Neumann.
    A : Accumulation matrix
    ct : Total compressibility matrix
    Act : A * ct --> matrix
    Po : Initial pressure --> value from initial condition

    """
    # Create one-Diagonal matrix A and ct
    A = np.diag(np.full(Nx, Accumulation), k=0)
    ct = np.diag(np.full(Nx, ct), k=0)
    Act = A*ct
    # Initializate matrix T, J and Q
    T = np.zeros((Nx, Nx))
    J = np.zeros((Nx, Nx))
    Q = np.zeros((Nx, 1))
    
    # Populate the tree diagonals of T
    T = np.diag(np.full(Nx, 2*Transmissibility), k=0)\
        +np.diag(np.full(Nx-1, -Transmissibility), k=-1)\
            +np.diag(np.full(Nx-1, -Transmissibility), k=1)
    T[0,0] = Transmissibility
    T[Nx-1, Nx-1] = Transmissibility
    
    # INITIAL CONDITIONS AND BOUNDARY CONDITIONS
    
    #Initial conditions P(x, t=0)
    Po = np.full((Nx, 1), IC)
    
    #Boundary Conditions P(x=0, t) 
    if 'Dirichlet' in BC_type[0]:
        J[0,0] = 2 * Transmissibility
        Q[0] = 2 * Transmissibility * BC_value[0]
   
    if 'Neumann' in BC_type[0]:
        Q[0] = BC_value[0]
    
    #Boundary Conditions P(x=L, t)
    if 'Dirichlet' in BC_type[1]:
        J[Nx-1,Nx-1] = 2 * Transmissibility
        Q[Nx-1] = 2 * Transmissibility * BC_value[1]
    
    if 'Neumann' in BC_type[1]:
        Q[Nx-1] = BC_value[1]
        
    return T, J, Q, A, ct, Act, Po
