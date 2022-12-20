# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% Import Libraries and packages
import numpy as np
from PVT_props_fns_pressure import P_dependent_properties


#%% Define functions
def block_properties(pvt, res, petro, numeric, p, Sw, So, Sg):
    '''
    Computes Block Properties
    
    For this assignment --> Block fluid properties are considered as homogenous and NO time dependent
    '''
    
    pc = np.zeros(Sw.shape)
    Bw = pvt['Bw'] * np.ones((numeric['N'],1))
    Bo, Bg, co_star, z = P_dependent_properties(p,res,pvt)
    betaw, betao, betag = Bw, Bo, Bg
    volume = numeric['D'][:,0] * numeric['D'][:,1] * numeric['D'][:,2]
    A = volume * res['phi'] / numeric['dt'] # [ft^3/day]
    A = A.reshape((numeric['N'],1))
    ct = betaw / Bw * Sw * (pvt['cw'] + pvt['cf']) + So*(pvt['co']+pvt['cf']) + betag/Bg*Sg*(pvt['cg']+pvt['cf'])
    C1 = betaw / Bw * A * Sw * (pvt['cw']+pvt['cf'])
    C2 = A * So * (co_star + pvt['cf'])
    C3 = betag / Bg * A * (Sg * (pvt['cg'] + pvt['cf']) + So*pvt['Rs'] * (pvt['co'] + pvt['cf']))
    
    return A, ct, betaw, betao, betag, pc, C1, C2, C3, Bw, Bo, Bg, z, co_star