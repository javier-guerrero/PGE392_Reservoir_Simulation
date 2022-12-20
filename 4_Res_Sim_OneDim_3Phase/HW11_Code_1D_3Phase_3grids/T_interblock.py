# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% Import Libraries and packages
import yaml
from yaml.loader import SafeLoader
import numpy as np
from numpy.linalg import matrix_power
import matplotlib.pyplot as plt
from scipy import special, sparse, integrate
import scipy.sparse.linalg as lg
import sympy as sym
import scipy as sp
import scipy.optimize as opt


#%% Define functions
def T_Interblock(i1,i2,direction,res,petro,pvt,numeric,krw,kro,krg,potential,Bw,Bo,Bg):
    '''
    Parameters
    ----------
    i1 : block i
    i2 : block j
    direction : Transmissibility direction: x-dir or y-dir
    res : reservoir props
    petro : petrophysical properties
    PVT : PVT properties
    krw, kro, krg : Relative perms
    Potential : fn to determine flow direction 
    Bw, Bo, Bg : Formation volume factors
    

    Returns
    -------
    Tw : Interblock Transmissibility for Water 
    To : Interblock Transmissibility for Oil
    Tg : Interblock Transmissibility for Gas

    '''
    conv =  0.00633 # Convert md/cp to ft^2/psi-day
    
    # Harmonc Average for Geometric (reservoir) properties
    if direction == 'x':
        num = 2*res['kx'][i1]*numeric['D'][i1,1]*numeric['D'][i1,2]*res['kx'][i2]*numeric['D'][i2,1]*numeric['D'][i2,2]
        den1 = res['kx'][i1]*numeric['D'][i1,1]*numeric['D'][i1,2]*numeric['D'][i2,0]
        den2 = res['kx'][i2]*numeric['D'][i2,1]*numeric['D'][i2,2]*numeric['D'][i1,0]
    elif direction == 'y':
        num = 2*res['ky'][i1]*numeric['D'][i1,0]*numeric['D'][i1,2]*res['ky'][i2]*numeric['D'][i2,0]*numeric['D'][i2,2]
        den1 = res['ky'][i1]*numeric['D'][i1,0]*numeric['D'][i1,2]*numeric['D'][i2,1]
        den2 = res['ky'][i2]*numeric['D'][i2,0]*numeric['D'][i2,2]*numeric['D'][i1,1]
    elif direction == 'z':
        num = 2*res['kz'][i1]*numeric['D'][i1,0]*numeric['D'][i1,1]*res['kz'][i2]*numeric['D'][i2,0]*numeric['D'][i2,1]
        den1 = res['kz'][i1]*numeric['D'][i1,0]*numeric['D'][i1,1]*numeric['D'][i2,2]
        den2 = res['kz'][i2]*numeric['D'][i2,0]*numeric['D'][i2,1]*numeric['D'][i1,2]
    else:
        print('Error in Interblock Function')
    harmonic = num/(den1+den2)
    
    # UPWINDNG --> flow direction
    
    idx_for_iw = np.asarray(np.where(potential[i1,0] >= potential[i2,0]))[1,:]
    idx_for_iw2 = np.asarray(np.where(potential[i1,0] < potential[i2,0]))[1,:]
    idx_for_kw, idx_for_kw2 = i1[0,idx_for_iw], i2[0,idx_for_iw2]
    upwind_w = np.zeros((len(i1[0,:]),1))
    upwind_w[idx_for_iw] = krw[idx_for_kw]/pvt['muw']/Bw[idx_for_kw]
    upwind_w[idx_for_iw2] = krw[idx_for_kw2]/pvt['muw']/Bw[idx_for_kw2]
    
    idx_for_io = np.asarray(np.where(potential[i1,1] >= potential[i2,1]))[1,:]
    idx_for_io2 = np.asarray(np.where(potential[i1,1] < potential[i2,1]))[1,:]
    idx_for_ko, idx_for_ko2 = i1[0,idx_for_io], i2[0,idx_for_io2]
    upwind_o = np.zeros((len(i1[0,:]),1))
    upwind_o[idx_for_io] = kro[idx_for_ko]/pvt['muo'][idx_for_ko]/Bo[idx_for_ko]
    upwind_o[idx_for_io2] = kro[idx_for_ko2]/pvt['muo'][idx_for_ko2]/Bo[idx_for_ko2]
    
    idx_for_ig = np.asarray(np.where(potential[i1,2] >= potential[i2,2]))[1,:]
    idx_for_ig2 = np.asarray(np.where(potential[i1,2] < potential[i2,2]))[1,:]
    idx_for_kg, idx_for_kg2 = i1[0,idx_for_ig], i2[0,idx_for_ig2]
    upwind_g = np.zeros((len(i1[0,:]),1))
    upwind_g[idx_for_ig] = krg[idx_for_kg]/pvt['mug'][idx_for_kg]/Bg[idx_for_kg]
    upwind_g[idx_for_ig2] = krg[idx_for_kg2]/pvt['mug'][idx_for_kg2]/Bg[idx_for_kg2]
    
    Tw = np.transpose(upwind_w)*harmonic*conv
    To = np.transpose(upwind_o)*harmonic*conv
    Tg = np.transpose(upwind_g)*harmonic*conv
    
    return Tw, To, Tg


def potentials(pvt,res,numeric,po,pc,Bw,Bo,Bg,z):
     R = 10.7316 # gas constant
     rho_g = pvt['Mg']*po/z/R/(res['T']+460) # gas density [lbm/ft^3]
     rhogsc = pvt['Mg']/379.4 # gas density at standard cond. [lbm/ft3]
     rho_w, rho_o = pvt['rhowsc']/Bw, (pvt['rhoosc']+pvt['Rs']*rhogsc/5.615)/Bo # water and oil density [lbm/ft3]
     D_mid = res['D']
     phi_o = po-rho_o*(D_mid-res['woc']).reshape((numeric['N'],1))/144
     phi_w = po-rho_w*(D_mid-res['woc']).reshape((numeric['N'],1))/144-pc
     phi_g = po-rho_g*(D_mid-res['woc']).reshape((numeric['N'],1))/144-pc
     phi_o, phi_w, phi_g = phi_o.reshape(numeric['N']), phi_w.reshape(numeric['N']), phi_g.reshape(numeric['N'])
     potential = np.array([phi_w,phi_o,phi_g]).T
     
     return potential