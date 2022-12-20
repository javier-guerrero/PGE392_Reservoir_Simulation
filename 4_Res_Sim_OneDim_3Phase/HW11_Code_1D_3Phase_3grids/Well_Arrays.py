# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""


#%% Import Libraries and packages

import numpy as np
from Productivity_Index import Jindex

#%% Define functions
def Well_Arrays(well,res,pvt,petro,numeric,Q,J,Qw,Qo,Qg,Jw,Jo,Jg,Sw,So,betaw,betao,betag,Bw,Bo,Bg,po,time,krw,kro,krg):
    '''
    Updates vectors Q and J according to the wells given in the input file
    
    Inputs:         J,Q, Well
    Outputs:        Updated J, Q

    '''
    ro = np.zeros(well['array'].shape)
    if (well['pwf'][0] < well['min_BHP'][0] and well['type'][0] == 0 and time > numeric['dt']):
        well['type'][0] = 1
        well['rates'][0] =  well['min_BHP'][0]
        print('\n Constant Rate Producer converted to BHP \n')
    
    if (time > 730 and well['rates'][1] == 0 and well['rates'][2] == 0):
        well['rates'][1], well['rates'][2] = well['inj_rate'][1], well['inj_rate'][2]
        print('\n Injectors Online \n')

    for k in range(len(well['well_id'])):
        for j in range(well['perf_blocks'][k]):
            block = int(well['array'][k,j])
            well['Jw'][k,j], well['Jo'][k,j], well['Jg'][k,j] = Jindex(k,block,well,res,pvt,petro,numeric,Sw,krw,kro,krg)
            
            if well['type'][k] == 0: # constant rate well
            
                if well['rates'][k] > 0: # Injector well
                    well['wc'][k,j] = 1.0
                elif well['rates'][k] < 0: # Producer well
                    well['wc'][k,j] = (krw[block]/pvt['muw']/Bw[block])/(krw[block]/pvt['muw']/Bw[block] + kro[block]/pvt['muo'][block]/Bo[block])
                
                Qo_sc = well['rates'][k]*(1-well['wc'][k,j])
                Qw_sc = well['rates'][k]*well['wc'][k,j]
                Qg_sc = Qo_sc*(well['Jg'][k,j]/well['Jo'][k,j])*(Bo[block]/Bg[block])
                Qw[block], Qo[block], Qg[block] = Qw_sc*betaw[block], Qo_sc*betao[block], Qg_sc*betag[block]
                Q[block] = Qw[block] + Qo[block] + Qg[block]
                well['Q_field'][k,j] = Qo.toarray()[block]/betao[block] + Qw.toarray()[block]/betaw[block]
                ro[k] = -Qo.toarray()[block]/Bo[block]
                
                
            elif well['type'][k] == 1: # constant BHP well
                Qw[block], Qo[block], Qg[block] = well['Jw'][k,j]*well['rates'][k], well['Jo'][k,j]*well['rates'][k], well['Jg'][k,j]*well['rates'][k]
                Q[block] = Q[block] + Qw[block] + Qo[block] + Qg[block]
                Jw[block,block], Jo[block,block], Jg[block,block] = well['Jw'][k,j], well['Jo'][k,j], well['Jg'][k,j] 
                J[block,block] = J[block,block] + Jw[block,block] + Jo[block,block] + Jg[block,block]
                well['Q_field'][k,j] = -well['Jo'][k,j]*(po[block]-well['rates'][k])/Bo[block] - well['Jw'][k,j]*(po[block]-well['rates'][k])/Bw[block]
                ro[k] = well['Jo'][k,j]*(po[block]-well['rates'][k])/Bo[block]
                
            if well['type'][k] == 0:
                well['pwf'][k,j] = po[block] + Qo.toarray()[block]/well['Jo'][k,j]
                well['GOR'][k] = (Qg[block]+pvt['Rs'][block]*Qo[block])/Bg[block]/Qo_sc
            else:
                well['pwf'][k,j] = well['rates'][k]
                qgrc = -well['Jg'][k,j]*(po[block]-well['rates'][k])
                qorc = -well['Jo'][k,j]*(po[block]-well['rates'][k])
                qosc = qorc/Bo[block]
                well['GOR'][k] = (qgrc + qorc*pvt['Rs'][block])/Bg[block]/qosc
                
    return Q, J, Qw, Qo, Qg, Jw, Jo, Jg, ro