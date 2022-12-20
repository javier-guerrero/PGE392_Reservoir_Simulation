"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

# Import Librariees and packages
import numpy as np
from BCM_kr import BCM_RelativePermeability

# Define Functions
def T_interblock(i, j, dz, dx, dy, kx, ky,  direction, P, Sw, petro, PVT, Res):
    '''
    Parameters
    ----------
    i : block i
    j : block j
    dz : block depth - Thickness
    dx : delta x
    dy : delta y
    kx : Perm in x-dir
    ky : perm in y-dir
    direction : Transmissibility direction: x-dir or y-dir
    P : pressure --> for direction of flow
    Sw : Water Saturation
    PP : petrophysical properties
    PVT : PVT properties
    
    Dependencies
    ------------
    BCM_RelativePermeability --> fn to compute relative permeability

    Returns
    -------
    T_interblock_w : Interblock Transmissibility for Water 
    T_interblock_o : Interblock Transmissibility for Oil
    
    TO-DO: add upwinding for gas when we use 3 phase flow
    
    '''
    
    CF = 6.33e-3 # Conversion factor
    
    #Harmonic mean geometric properties Eq. 9.17a
    if direction == 'x':
        T_interblock = 2 * dz[i] * kx[i] * (dy[i] * kx[j] * dy[j] / (kx[i] * dy[i] * dx[j] + kx[j] * dy[j] * dx[i]))
    elif direction == 'y':
        T_interblock = 2 * dz[i] * ky[i] * (dx[i] * ky[j] * dx[j] / (ky[i] * dx[i] * dy[j] + ky[j] * dx[j] * dy[i]))
    
    # The direction of the potentails has been specified 
    INITPOT = np.asarray([4,8,3,5,9,2,6,7,1,])*1e-5
    POT_i  = P[i,0]  -  PVT.rhoo / 144.0 * (Res.D[i] - Res.woc) + INITPOT[i]
    POT_j  = P[j,0]  - PVT.rhoo / 144.0 * (Res.D[j] - Res.woc) + INITPOT[j]
    #print(POT_i, POT_j)
    
    # UPWINDING fluid properties Eq. 9.17b
    if POT_i >= POT_j:
        krw, kro = BCM_RelativePermeability(Sw[i,0], petro)
        # For this assignment --> properties are known and NOT time depedent
        Bw  = PVT.Bw 
        Bo  = PVT.Bo
        # Bg = PVT.Bg  
        muw = PVT.muw 
        muo = PVT.muo 
        #mug]PVT.mug

    elif POT_i < POT_j:
        krw, kro = BCM_RelativePermeability(Sw[j,0], petro)
        #For this assignment --> properties are known and NOT time depedent
        Bw  = PVT.Bw 
        Bo  = PVT.Bo
        # Bg=PVT.Bg  
        muw = PVT.muw
        muo = PVT.muo 
        #mug]PVT.mug
       
    #Interblock phase transmisibility Eq.9.17c
    T_interblock_w = CF * (krw / (muw * Bw)) * T_interblock
    T_interblock_o = CF * (kro / (muo * Bo)) * T_interblock
    # T_interblock_o = CF * (kro / (muo)) * T_interblock # to match the same results as Dr. Balhoff
    # T_interblock_g = CF * (kro / (muo * Bo)) * T_interblock
    
    
    return T_interblock_w, T_interblock_o # T_interblock_g

