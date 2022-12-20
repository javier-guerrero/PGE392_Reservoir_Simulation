# -*- coding: utf-8 -*-
"""
PGE-392K-  NUmerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""
#%% import packages and libraries
import numpy as np
import math
from scipy.sparse.linalg import spsolve

#%% Parameters
pi=np.pi


#%% Functions
def PD_1D_1PH_Hom_Analytical(xD, tD, PDi, method, XiD=0):
    """
    This function returns the analytical solution to the difussivity eqn in 1D, 
    single pahse, homogeneous reservoir 
    
    The analytical solution is in dimensionelss variables PD(xD, tD)

    Parameters
    ----------
    xD : Dimensionlees distance
    tD : Dimensionless time
    PDi : Initial condition --> value od PD at initial conditions
    method : boundary configuration
    XiD : Initial value od xD. The default is 0.

    Returns
    -------
    PD : Dimensionless Pressure

    """
    error=1e6
    Tolerance=1e-6
    
    if method=='Dirichlet-Dirichlet':
        n=1
        sum1=0
        sum2=0;
        
        while error>Tolerance:
            new1 = (((-1)**n)/n) * np.sin(n*pi*xD) * np.exp(-n**2*pi**2*tD)
            new2 = (((-1)**n)/n) * np.sin(n*pi*(xD-1)) * np.exp(-n**2*pi**2*tD)
            sum1 += new1
            sum2 += new2
            error = max(abs(new2)) + max(abs(new1))
            n += 1
            if n > 100:
                break
        PD = 1-xD - (2 * PDi / pi) * sum1 - (2 * (1-PDi) / pi) * sum2
    
    if method=='Dirichlet-Neumann':
        n=0
        sum1=0
        
        while error>Tolerance:
            new = (4/(pi*(2*n+1))) * np.exp(-(2*n+1)**2*pi**2*tD/4) * np.sin((2*n+1)*pi*xD/2)
            sum1 += new
            error = max(abs(new))
            n += 1
            if n > 100:
                break
        if tD == 0:
            PD = np.zeros(len(xD))
        else:
            PD = 1-sum1
    
    if method == 'Dirchlet-Semi-Infinite':
        #XiD=np.linspace(0,3,11)
        PD = np.ones(len(XiD))
        for i, var in enumerate(XiD):
            PD[i] = math.erfc(var)
    
    return PD

