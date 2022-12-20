# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

def initialization(depth, PP, PVT, Res):
    '''
    Function:       Initialization function to calculate phase pressure and Saturations @Given Depth
    inputs:         depth(as an array), PP, PVT and Res (as classes)
    Outputs:        Sw, So, Pw, Po @given depth

    Parameters
    ----------
    depth : depth 
    PP : Petrophysical Inputs
    PVT : PVT Inputs
    Res : Reservoir Inputs
    

    Returns
    -------
    Sw : Water Saturation
    So : Oil Saturation
    Pw : Pressure of water phase
    Po : Pressure of oleic phase

    '''
    #Constants (Rs was not given in the input)
    R=      10.73               #psi-ft^3/lbmol-R
    Tsc=    60+459.67           #R
    Psc=    14.7                #psi
    rs=     500                 #scf/STB
    
    #Density Calculations
    rhow=   PVT.rhowsc/ PVT.Bw                                      #Eq.1.6
    rhogsc= PVT.Mg*Psc/(R*Tsc)                                      #Eq.1.6
    rhoo=   (PVT.rhoosc + rhogsc*rs/5.61)/PVT.Bob                   #Eq. 1.13  5.61 ft^3 in a barrel (bbl)
    
    #Phase Pressure calculations
    Pw=     Res.pwoc+ rhow* (depth - Res.woc)/144.0                 #Eq.1.38   144 in^2 in 1 ft^2 
    Po =    (Res.pwoc + PP.Pce) + rhoo/ 144.0 * (depth - Res.woc)   #Eq.1.39 
    
    #Water saturation calculation
    Sw =    PP.Swr + (1-PP.Swr)*((Po-Pw)/PP.Pce) **(-PP.lam)        #From Corey-Brooks draining curve
    
    #Conditional for depths below the woc where Sw=1
    #Sw[depth>=Res.woc] = 1.0
    So=1-Sw
    
    return Sw, So, Pw, Po
