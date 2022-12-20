# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% Import Libraries and packages
import numpy as np


#%% Define functions
def RelativePerm_3Phase_StoneI(petro, Sw, So, Sg, phase):
    '''
    Computes Relative Permeability for 2-phase or 3-phase systems

    Parameters
    ----------
    petro : petrophysical properties
    Sw : Water saturation
    So : Oil Saturation
    Sg : Gas Saturation
    phase : Number of phases present in the reservoir

    Returns
    -------
    krw : Water rel perm
    kro : oil rel perm
    krg : gas rel perm
    krow : Oil rel perm in O-W system
    krog : Oil rel perm in O-G  system
    '''
    
    if phase == 2: # 2-phase system O-W (Sg=0)
        S_norm = (Sw-petro['Swr'])/(1-petro['Swr']-petro['Sor']) # normalized saturation
        krw, kro = petro['krwo']*S_norm**petro['n_w'], petro['kroo']*(1-S_norm)**petro['n_o'] # water/oil rel perm
        krg, krow, krog = np.zeros(krw.shape), np.zeros(krw.shape), np.zeros(krw.shape) # others are equal to zero
    
    elif phase == 3:
        SwD, SgD = (Sw-petro['Swr'])/(1-petro['Swr']-petro['Sorw']), (Sg-petro['Sgr'])/(1-petro['Sgr']-petro['Sorg']-petro['Swr']) # normalized saturations
        SoD = (So-petro['Sorw'])/(1-petro['Swr']-petro['Sorw']-petro['Sgr']) # normalized saturations
        a = Sg/(1-petro['Swr']-petro['Sorg']) # "a" factor --> see CH. 1 notes
        Som = (1-a)*petro['Sorw'] + a*petro['Sorg'] # weighting --> see CH. 1 notes
        Sos, Sws, Sgs = (So-Som)/(1-petro['Swr']-Som), (Sw-petro['Swr'])/(1-petro['Swr']-Som), Sg/(1-petro['Swr']-Som) # other values needed
        krow, krog = petro['krowo']*(1-SwD)**petro['Now'], petro['krogo']*(1-SgD)**petro['Nog'] # oil/water and oil/gas rel perm
        krw, krg = petro['krwo']*SwD**petro['Nw'], petro['krgo']*SgD**petro['Ng'] # water and gas rel perm
        kro = (Sos*krow*krog)/(petro['krogo']*(1-Sws)*(1-Sgs)) # oil rel perm
        kro[SoD>1] = petro['krowo'] # check for other values
        kro[SoD<0] = 0 # check for other values
        krg[Sg<petro['Sgr']] = 0
    
    return krw, kro, krg, krow, krog # return all rel perms


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