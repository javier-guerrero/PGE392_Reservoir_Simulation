# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""
import numpy as np
def BCM_RelativePermeability(Sw, petro):
    '''
    Computes the relative permeability using Brooks-Corey model

    Parameters
    ----------
    Sw : Water Saturation
    petro : Petrophysical properties

    Returns
    -------
    krw : Water Relative permeability
    kro : Oil Relative Permeability

    '''
    S=(Sw-petro.Swr)/(1-petro.Swr-petro.Sor)
    krw=petro.krwo*S**petro.n_w
    kro=petro.kroo*(1-S)**petro.n_o
    return krw, kro

def CapillaryPressure1(petro, Sw):
    '''
    Computes the Capillary pressure

    Parameters
    ----------
    Sw : Water Saturation
    petro : Petrophysical properties

    Returns
    -------
    Pc : Capillary pressure
    dPc : Derivative of Capillary pressure
    eps : 0.0001 small value added to Sw to avoid division by 0 in Pc and dPc and avoid inf and -inf values

    '''
    #Normalized Sw*=Swx and Se
    Se=(Sw-petro.Swr)/(1-petro.Swr-petro.Sor) #Eq. 1.33
    #Capillary Pressure(Pc)
    Pc=petro.Pe*(np.power(Se,(-1/petro.lam))-1) #Imbibition Eq. 1.32b
    dPc=-(petro.Pe/petro.lam)*(np.power(Se,(-(1+petro.lam)/petro.lam)))
    return Pc, dPc

def CapillaryPressure2(petro, Sw, Sw_hyst):
    '''
    Computes the Capillary pressure

    Parameters
    ----------
    petro : Petrophysical properties
    Sw : Water Saturation
    Sw_hyst: Water saturation at hysteresys

    Returns
    -------
    Pc_S : Capillary pressure for Scanning curve
    dPc_S : Derivative of Capillary pressure for Scanning curve
    eps : 0.0001 small value added to Sw to avoid division by 0 in Pc and dPc and avoid inf and -inf values

    '''
    
    #Normalized Sw*=Swx and Se
    Swx=(Sw-petro.Swr)/(1-petro.Swr) #Eq. 1.33
    Se=(Sw-petro.Swr)/(1-petro.Swr-petro.Sor) #Eq. 1.33
    #Capillary Pressure(Pc)
    Pc_D=petro.Pe*Swx**(-1/petro.lam) #Drainage Eq.1.32a
    Pc_I=petro.Pe*(Se**(-1/petro.lam)-1) #Imbibition Eq. 1.32b
    #rint(Sw_hyst, type(Sw_hyst), petro.epspc, type(petro.epspc))
    f=(((1-petro.Sor)-Sw_hyst+petro.epspc)/((1-petro.Sor)-Sw_hyst))*(Sw-Sw_hyst)/(Sw-Sw_hyst+petro.epspc) #Eq.1.37
    if Sw<=Sw_hyst:
        Pc_S=Pc_D
    else:
        Pc_S=f*Pc_I+(1-f)*Pc_D
        
    #Pc_S=[Pc_D[i] if Sw[i]<=Sw_hyst else f[i]*Pc_I[i]+(1-f[i])*Pc_D[i] for i in range(len(Sw))] #Scanning Eq.1.36
    #Derivative of Capillary Pressure (dPc)
    dPc_D=(1/(1-petro.Swr))*petro.Pe*(-1/petro.lam)*Swx**((-1/petro.lam)-1)# Derivative of Eq.1.32a
    dPc_I=(1/(1-petro.Swr-petro.Sor))*petro.Pe*((-1/petro.lam)*Se**(-1/petro.lam-1))# Derivative of Eq.1.32b
    if Sw<=Sw_hyst:
        dPc_S=dPc_D
    else:
        dPc_S=f*dPc_I+(1-f)*dPc_D
    #dPc_S=[dPc_D[i] if Sw[i]<=Sw_hyst else f[i]*dPc_I[i]+(1-f[i])*dPc_D[i] for i in range(len(Sw))]   
    return Pc_S, dPc_S 

def CapillaryPressure(petro, Sw, Sw_hyst):
    Pc_S = 0 
    dPc_S = 0
    return Pc_S, dPc_S 
    