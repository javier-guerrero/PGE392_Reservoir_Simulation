# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""
#%% Import Libraries and Packages
from Productivity_Index import Jindex
from scipy.sparse import csr_matrix
import numpy as np
from BCM_kr import BCM_RelativePermeability

#%% Define Functions
def Well_Arrays(Sw, J, Jw, Jo, Q, Qw, Qo, Well, Res, PVT, Num, petro):
    '''
    Updates vectors Q and J according to the wells given in the input file
    
    Inputs:         J,Q, Well
    Outputs:        Updated J, Q

    '''
    J = J.toarray() 
    Q = Q.toarray()
    Qw = Qw.toarray() 
    Qo = Qo.toarray()
    Jw = Jw.toarray() 
    Jo = Jo.toarray()
    
    for i in Well.well_id:
        # NEXT 2 LINES ARE FOR THE NEW VERSION
        qw = 0
        qo = 0
        
        for j in Well.blocks[i][0]:
            krw,kro=BCM_RelativePermeability(Sw[j], petro)
            # BHP Well
            if Well.type[i] == 1: 
                if Well.Jindex != None:
                    J[j,j] = J[j,j] + Well.Jindex[i]
                    Q[j] = Q[j] + J[j,j] * Well.rates[i]
                else:
                    Jwell_oil, Jwell_water = Jindex(j, i, Well, Res, PVT, Num, Sw[j,0], petro)
                    Jw[j,j] = Jw[j,j] + Jwell_water
                    Jo[j,j] = Jo[j,j] + Jwell_oil
                    # Qw[j] = Qw[j] + Jwell_water * Well.rates[i]
                    # Qo[j] = Qo[j] + Jwell_oil * Well.rates[i]  
                    # Qw[j] = Jwell_water * Well.rates[i]
                    # Qo[j] = Jwell_oil * Well.rates[i] 
            #Constant Rate Well
            else: 
                # OLD VERSION --> with fractional flow
                # if Well.kind[i] == 0: #0-Producer 1-Injector
                #     krw, kro = BCM_RelativePermeability(Sw[j], petro)
                #     M= (kro * PVT.Bw) / (krw * PVT.Bo)
                #     fw= 1.0 / (1.0 + M)
                #     #print('fw',fw)
                #     Qw[j] = Qw[j] - fw * Well.rates[i] / (len(Well.blocks[i][0])) * PVT.Bw
                #     Qo[j] = Qo[j] - (1 - fw) * Well.rates[i] / (len(Well.blocks[i][0])) * PVT.Bo
                
                # NEW VERSION
                Jwell_oil, Jwell_water = Jindex(j, i, Well, Res, PVT, Num, Sw[j,0], petro)
                qw += ((Jwell_water / PVT.Bw) / (Jwell_water / PVT.Bw + Jwell_oil / PVT.Bo)) * Well.rates[i]
                qo += ((Jwell_oil / PVT.Bo) / (Jwell_water / PVT.Bw + Jwell_oil / PVT.Bo)) * Well.rates[i]
                if Well.kind[i] == 0: # 0-producer 1-Injector
                    Qw[j] = Qw[j] - qw
                    Qo[j] = Qo[j] - qo * 1.5
                
                
                else:# Injector well (It considers that the fluid injected is 100% water)
                    # OLD VERSION
                    # Qw[j] = Qw[j] + Well.rates[i] / (len(Well.blocks[i][0]))
                    # Qo[j] = Qo[j] + 0
                    
                    #NEW VERSION
                    Qw[j] = Qw[j] + Well.rates[i] 
                    Qo[j] = Qo[j] + 0
    Q = Qw + Qo
    J = Jw + Jo
            
    J = csr_matrix(J)  
    Q=csr_matrix(Q)
    Qw=csr_matrix(Qw) 
    Qo=csr_matrix(Qo) 
    Jo=csr_matrix(Jo)  
    Jw=csr_matrix(Jw)
    
    return J, Jw, Jo, Q, Qw, Qo

def Well_Arrays_SS(Sw,J,Q, Well,Res,PVT,Num, petro):
    '''
    Updates vectors Q and J according to the wells given in the input file
    for SS Method
    
    Inputs:         J,Q, Well
    Outputs:        Updated J, Q

    '''
    Qw = np.zeros((Num.N,1))
    Q1 = np.zeros((2*Num.N,1))
    J = J.toarray()
    Q = Q.toarray();
    
    for i in Well.well_id:
        for j in Well.blocks[i][0]:
            krw, kro = BCM_RelativePermeability(Sw[j], petro)
            if Well.type[i] == 1:
                if Well.Jindex != None:
                    J[j,j] = J[j,j] + Well.Jindex[i]
                    Q[j] = Q[j] + J[j,j] * Well.rates[i]
                else:
                    J[j,j] = J[j,j] + Jindex(j, i, Well, Res, PVT, Num)  # Calls Jindex fn
                    Q[j] = Q[j] + J[j,j] * Well.rates[i]      
            else: 
                if Well.kind[i] == 0: #0=Producer 1=Injector
                    Q[j] = Q[j] - Well.rates[i] / (len(Well.blocks[i][0]))
                    Qw[j] = Q[j] * (krw / (PVT.muw * PVT.Bw) / (krw/(PVT.muw * PVT.Bw) + kro / (PVT.muo * PVT.Bo)))
                else:
                    Q[j] = Q[j] + Well.rates[i] / (len(Well.blocks[i][0]))
                    Qw[j] = Q[j]
                Q1[2*j] = Qw[j]
                Q1[2*j+1] = Q[j] - Qw[j]
    
    
    Qo = Q - Qw
    J = csr_matrix(J) 
    Q = csr_matrix(Q) 
    Qw = csr_matrix(Qw) 
    Qo = csr_matrix(Qo) 
    Q1 = csr_matrix(Q1)
    
    return J, Q1, Qw, Qo
