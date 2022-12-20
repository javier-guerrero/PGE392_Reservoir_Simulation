# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% Import Libraries and Packages
import numpy as np
from BCM_kr import BCM_RelativePermeability

#%% Define Functions
def Jindex(WellGrid,WellID, Well,Res,PVT,Num, Sw, petro):  
    '''
    Function:       Calculate the productivity index of a well located in a i block
    input:          WellGrid, WellID, Num, Res, PVT, Well
    output:         Productivity Index

    Parameters
    ----------
    WellGrid : grid in which the well is located
    WellID : Well ID
    Well : Well inputs
    Res : Reservoir inputs
    PVT : pvt inputs
    Num : Numerical Inputs
    Sw : Water Saturation
    PP : petrophysical inputs

    Returns
    -------
    Jwell_oil : Prod Index of the well for oil
    Jwell_water : Prod Index of the well for water
    '''
    
    i = WellGrid 
    j = WellID
   
    if Well.direction[j] == 'x-direction':
        k1=Res.kz[i] 
        k2=Res.ky[i] 
        d1=Num.dy[i] 
        d2=Num.dz[i] 
        d3=Num.dx[i]
    
    elif Well.direction[j] == 'y-direction':
        k1=Res.kz[i] 
        k2=Res.kx[i] 
        d1=Num.dz[i] 
        d2=Num.dx[i]
        d3=Num.dy[i]
    
    elif Well.direction[j] == 'z-direction':
        k1=Res.ky[i] 
        k2=Res.kx[i] 
        d1=Num.dy[i] 
        d2=Num.dx[i]
        d3=Num.dz[i]
    
    krw, kro=BCM_RelativePermeability(Sw, petro)
    
    #Peaceman model for anysotropy Eq 5.17
    req   = 0.28*np.sqrt(np.sqrt(k1/k2)*d2**2.0 + np.sqrt(k2/k1)*d1**2.0)/ ((k1/k2)**0.25 + (k2/k1)**0.25) #Eqn 5.17
    
    #Productivity index for each phase
    Jwell_oil = 6.33E-3*(2*np.pi*kro*np.sqrt(k1*k2)*d3)/(PVT.muo*(np.log(req/Well.rw[j]) + Well.skin[j]))   #Eqn 5.18
    Jwell_water = 6.33E-3*(2*np.pi*krw*np.sqrt(k1*k2)*d3)/(PVT.muw*(np.log(req/Well.rw[j]) + Well.skin[j])) #Eqn 5.18
    #print('grid',i,'well',j)
    #print(req,(2*(kro/5)*np.pi*np.sqrt(k1*k2)*d3)/((np.log(req/Well.rw[j]) + Well.skin[j])), Jwell_water)
    
    return Jwell_oil, Jwell_water;