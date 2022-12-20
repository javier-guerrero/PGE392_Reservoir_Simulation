# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""
from Productivity_Index import Jindex
#from scipy.sparse import csr_matrix

def Well_Arrays(J,Q, Well,Res,PVT,Num):
    '''
    Update the vectors Q and J according to the wells given in the input file

    Parameters
    ----------
    J : Dirichlet BC  and Well conditions (diagonal matrix)
    Q : Neumann BC and Well conditions (vector)
    Well :      Well inputs (yml)
    Res :       Reservoir Inputs (yml)
    PVT :       PVT Inputs (yml)
    Num :       Numerical inputs (yml)

    Returns
    -------
    J : Updated J 
    Q : Updated Q

    '''
    for i in Well.well_id:
        for j in Well.blocks[i][0]:
            if Well.type[i]==1:
                if Well.Jindex!=None:
                    J[j,j]=J[j,j]+Well.Jindex[i]
                    Q[j]=Q[j]+J[j,j]*Well.rates[i]
                else:
                    J[j,j]=J[j,j]+Jindex(j,i, Well,Res,PVT,Num)
                    Q[j]=Q[j]+J[j,j]*Well.rates[i]      
            else:
                if Well.kind[i]==0:
                    Q[j]=Q[j]-Well.rates[i]/(len(Well.blocks[i][0]))
                else:
                    Q[j]=Q[j]+Well.rates[i]/(len(Well.blocks[i][0]))
    
    return J,Q
