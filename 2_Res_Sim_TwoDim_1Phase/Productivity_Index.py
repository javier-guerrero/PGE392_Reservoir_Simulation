# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""
#%% Import Libraries and packages
import numpy as np

#%% Define functions

def Jindex(WellGrid, WellID, Well, Res, PVT, Num):  
    '''
    Calculate the productivity index of a well located in a i block
    input:          WellGrid, WellID, Num, Res, PVT, Well
    output:         Productivity Index

    Parameters
    ----------
    WellGrid :  Gridblocks where well is located
    WellID :    Well ID
    Well :      Well inputs (yml)
    Res :       Reservoir Inputs (yml)
    PVT :       PVT Inputs (yml)
    Num :       Numerical inputs (yml)

    Returns
    -------
    Productivity index

    '''
    i = WellGrid
    j = WellID
   
    if Well.direction[j]=='x-direction':
        k1=Res.kz[i]; k2=Res.ky[i]; d1=Num.dy[i]; d2=Num.dz[i]; d3=Num.dx[i]
    
    elif Well.direction[j]=='y-direction':
        k1=Res.kz[i]; k2=Res.kx[i]; d1=Num.dz[i]; d2=Num.dx[i]; d3=Num.dy[i]
    
    elif Well.direction[j]=='z-direction':
        k1=Res.ky[i]; k2=Res.kx[i]; d1=Num.dy[i]; d2=Num.dx[i]; d3=Num.dz[i]
        
    #Peaceman model for anysotropy Eq 5.17
    req   = 0.28*np.sqrt(np.sqrt(k1/k2)*d1**2.0 + np.sqrt(k2/k1)*d2**2.0)/ ((k1/k2)**0.25 + (k2/k1)**0.25)
    Jwell = 6.33E-3*(2*np.pi*np.sqrt(k1*k2)*d3)/(PVT.muo*(np.log(req/Well.rw[j]) + Well.skin[j]))
    #print(k1,k2,d1,d2,d3,req,Jwell)
    
    return Jwell;