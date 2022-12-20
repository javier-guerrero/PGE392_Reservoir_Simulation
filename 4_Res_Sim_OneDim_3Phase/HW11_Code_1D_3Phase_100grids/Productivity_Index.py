# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% Import Libraries and packages
import numpy as np


#%% Define functions
def Jindex(well_no,block,well,res,pvt,petro,numeric,Sw,krw,kro,krg):
    conv =  0.00633 # Convert [md/cp] to [ft^2/psi-day]
    dx,dy,dz = numeric['D'][block,0], numeric['D'][block,1], numeric['D'][block,2]
    
    if well['direction'][well_no] == 'z': # vertical well!
        dw1, dw2, dw3, k1, k2 = dx, dy, dz, res['kx'][block], res['ky'][block]
    elif well['direction'][well_no] == 'x': # horizontal well!
        dw1, dw2, dw3, k1, k2 = dz, dy, dx, res['kz'][block], res['ky'][block]
    elif well['direction'][well_no] == 'y': # horizontal well!
        dw1, dw2, dw3, k1, k2 = dx, dz, dy, res['kx'][block], res['kz'][block]
        
    num = np.sqrt(np.sqrt(k2/k1)*dw1**2 + np.sqrt(k1/k2)*dw2**2)
    den = (k2/k1)**(1/4) + (k1/k2)**(1/4)
    req = 0.28*num/den
    
    if (well['rates'][well_no] > 0 and well['type'][well_no] == 0):
        krw[block] = 1
    else:
        pass
    
    Jw = conv*2*np.pi*dw3*np.sqrt(k1*k2)*krw[block]/(pvt['muw']*(np.log(req/well['rw'][well_no])+well['skin'][well_no]))
    Jo = conv*2*np.pi*dw3*np.sqrt(k1*k2)*kro[block]/(pvt['muo'][block]*(np.log(req/well['rw'][well_no])+well['skin'][well_no]))
    Jg = conv*2*np.pi*dw3*np.sqrt(k1*k2)*krg[block]/(pvt['mug'][block]*(np.log(req/well['rw'][well_no])+well['skin'][well_no]))
    return Jw, Jo, Jg