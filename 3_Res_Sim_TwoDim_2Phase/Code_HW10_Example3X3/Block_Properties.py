# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

def Block_Properties(i, Sw, Res, PVT, petro, Num, dPcow):
    '''
    Computes Block Properties
    
    For this assignment --> Block fluid properties are considered as homogenous and NO time dependent
    '''
    
    #Water properties
    Bw = PVT.Bw
    cw = PVT.cw
    muw = PVT.muw
    Sw = Sw[i,0]
    
    #Oil properties
    Bo = PVT.Bo
    co = PVT.co
    muo = PVT.muo
    So = 1 - Sw
    
    #Gas properties 
    Bg = PVT.Bg
    cg = PVT.cg
    Sg = 0
    
    #Capillary effects are not considered in this assignment
    # dPcow = 0
    
    #Using Eq. 9.8
    Beta_o = Bo
    Beta_w = Bw / (1 - Sw * cw * dPcow)
    Beta_g = Bg / (1 - Sg * cg * dPcow)
    A = Num.dx[i] * Num.dy[i] * Num.dz[i] * Res.phi[i] / Num.dt
    
    #Using equation Eq. 9.10
    #print(Beta_w,Sw,cw,PVT.cf,Bw,So,co,PVT.cf,Beta_g,Sg,cg,PVT.cf,Bg)
    ct = Beta_w * Sw * (cw + PVT.cf) / Bw + So * (co + PVT.cf) + Beta_g * Sg * (cg + PVT.cf) / Bg
   
    #Coefficients for each block on the accumulation term Eq. 9.7
    C1 = Beta_w * A * Sw * (cw + PVT.cf)
    C2 = A * (1 - Sw) * (co + PVT.cf)
    C3 = Beta_g * A * (Sg * (cg + PVT.cf) + So * PVT.Rs * (co +PVT.cf))
    
    return Beta_o, Beta_w, A, ct, C1, C2, C3


# def Block_Properties(i, Sw, Res, PVT, petro, Num,dPcow):
#     #Block fluid properties (for this assignments these properties are considered as homogenous an no time depedant)
    
#     #Oil properties
#     Bo = PVT.Bo
#     co=PVT.co
#     muo= PVT.muo
#     Sw=Sw[i,0]
#     #Water properties
#     Bw = PVT.Bw
#     cw=PVT.cw
#     muw=PVT.muw
#     So=1-Sw
#     #Gas properties 
#     Bg= PVT.Bg
#     cg=PVT.cg
#     Sg=0
#     #Capillary effects are not considered in this assignment
#     #Using Eq. 9.8
#     Beta_o=Bo
#     Beta_w=Bw/(1-Sw*cw*dPcow)
#     Beta_g=Bg/(1-Sg*cg*dPcow)
#     A=Num.dx[i]*Num.dy[i]*Num.dz[i]*Res.phi[i]/Num.dt
#     #Using equation Eq. 9.10
#     #print(Beta_w,Sw,cw,PVT.cf,Bw,So,co,PVT.cf,Beta_g,Sg,cg,PVT.cf,Bg)
#     ct=Beta_w*Sw*(cw+PVT.cf)/Bw+So*(co+PVT.cf)+Beta_g*Sg*(cg+PVT.cf)/Bg
#     #Coefficients for each block on the accumulation term Eq. 9.7
#     C1=Beta_w*A*Sw*(cw+PVT.cf)
#     C2=A*(1-Sw)*(co+PVT.cf)
#     C3=Beta_g*A*(Sg*(cg+PVT.cf)+So*PVT.Rs*(co+PVT.cf))
#     return Beta_o, Beta_w, A, ct, C1, C2, C3