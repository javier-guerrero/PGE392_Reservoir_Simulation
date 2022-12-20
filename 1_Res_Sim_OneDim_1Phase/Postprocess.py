# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% import packages and libraries
import numpy as np
import matplotlib.pyplot as plt
from Grid_Arrays import Create_Grid_Arrays
from Analytical_Solution import PD_1D_1PH_Hom_Analytical
from Numerical_Solution import P_1D_1PH_Hom_Numerical

#%% Functions 
def Plot_Analytical_Numerical(i, xD, tD, BCtype, BCvalue, ICP, Analyticaltype, Accumulation, Transmissibility, Resct, NumN, Numdt, Numtf, etha):
    """
    This function creates a comparison plot between the analytical solution and
    numerical solution using the explicit, implicit and Crank-Nicholson methods
    
    As this function calls the analytical and numerical funtions to create the solutions
    then we need to create the arrays to solve the diffusivity eqn in matrix form --> calls Grid_Arrays

    Parameters
    ----------
    As this fn calls 3 other funs (Grid_Arrays, Analytical_Solution, Numerical_+Solution)
    then the parameters are the same used in these fns

    Returns
    -------
    Figure with analytical solution (continous line) and numerical solutions (points)

    """
    T, J, Q, A, ct, Act, Po = Create_Grid_Arrays(Accumulation, Transmissibility, Resct, NumN, BCtype, BCvalue, ICP) 
   
    if BCtype == [['Neumann'],['Dirichlet']]:
        PDi = 0
    elif BCtype == [['Dirichlet'],['Neumann']]:
        PDi = 0
    else:
        PDi = (ICP - BCvalue[1]) / (BCvalue[0] - BCvalue[1])
    if BCtype == [['Neumann'],['Dirichlet']]:
        # the Analytical solution is coded assuming Dirichlet in the left and Nuemann in the right 
        # Therefore, specify (1-xD) to invert the configutration ==> Nuemann in the left and Dirichlet in the Right boundary
        PD_Analytical = PD_1D_1PH_Hom_Analytical(1-xD, tD, PDi, Analyticaltype) 
    else:
        PD_Analytical = PD_1D_1PH_Hom_Analytical(xD, tD, PDi, Analyticaltype)
    
    plt.figure(i)
    plt.style.use('seaborn')
    plt.plot(xD, PD_Analytical, 'k-', label='Analytical')
    
    methods = ['Implicit', 'Crank-Nicholson', 'Explicit']
    theta = np.linspace(0, 1, 3)
    for i in range(len(theta)):
        P_Num = P_1D_1PH_Hom_Numerical(Numtf, Po, Numdt, NumN, Act, T, J, Q, theta[i])
        if BCtype == [['Neumann'],['Dirichlet']]:
            PD_Num = (P_Num-ICP)/(BCvalue[1]-ICP)
        elif BCtype == [['Dirichlet'],['Neumann']]:
            PD_Num = (P_Num-ICP)/(BCvalue[0]-ICP)
        else:
            PD_Num = (P_Num-BCvalue[1])/(BCvalue[0]-BCvalue[1])
        plt.plot(xD, PD_Num[:,int(Numtf/ Numdt)], 'o--',label=methods[i])
    plt.title(f'ANALYTCAL AND NUMERICAL SOLUTIONS (tD={np.round(tD,2)} and t={round(Numtf,1)} days) \n Boundary Conditions {BCtype} \n Nx={NumN} dt={Numdt} days  $\eta$={round(etha,3)}')
    plt.ylabel('$P_D$'); plt.xlabel('$x_{D}$');
    plt.xlim([0,1])
    plt.grid(True)
    plt.legend()
    
    # return plt.figure()

    