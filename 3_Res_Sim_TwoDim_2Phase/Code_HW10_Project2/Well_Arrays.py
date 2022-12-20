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
def Well_Arrays(Sw,J,Jw,Jo,Q,Qw,Qo,Well,Res,PVT,Num, PP,P,Qom,Qwm):
    '''
    Updates vectors Q and J according to the wells given in the input file
    
    Inputs:         J,Q, Well
    Outputs:        Updated J, Q

    '''
    J = J.toarray() 
    Jw = Jw.toarray() 
    Jo = Jo.toarray()
    Q = Q.toarray()
    Qw = Qw.toarray() 
    Qo = Qo.toarray()
    
    J_test=np.zeros((Num.N,Num.N))
    Jw_test=np.zeros((Num.N,Num.N))
    Jo_test=np.zeros((Num.N,Num.N))
    BHP=np.zeros((Num.N,1))
    for i in Well.well_id:
        for j in Well.blocks[i][0]:
            krw,kro=BCM_RelativePermeability(Sw[j], PP)
            if Well.type[i]==1: # BHP Well
                if Well.Jindex!=None: #In case that J is given
                    J[j,j]=J[j,j]+Well.Jindex[i]
                    Q[j]=Q[j]+J[j,j]*Well.rates[i]
                else:
                    Jwell_oil, Jwell_water=Jindex(j,i, Well,Res,PVT,Num, Sw[j,0], PP);
                    Jw[j,j]=Jw[j,j]+Jwell_water
                    Jo[j,j]=Jo[j,j]+Jwell_oil
                    Jw_test[j,j]=Jw_test[j,j]+Jwell_water
                    Jo_test[j,j]=Jo_test[j,j]+Jwell_oil
                    BHP[j,0]=Well.rates[i]
            else: #Constant Rate Well
                Jwell_oil, Jwell_water=Jindex(j,i, Well,Res,PVT,Num, Sw[j,0], PP);
                qw=((Jwell_water/PVT.Bw)/(Jwell_water/PVT.Bw+Jwell_oil/PVT.Bo))*(Well.rates[i]/len(Well.blocks[i][0]))
                qo=((Jwell_oil/PVT.Bo)/(Jwell_water/PVT.Bw+Jwell_oil/PVT.Bo))*(Well.rates[i]/len(Well.blocks[i][0]))
                Jw_test[j,j]=Jw_test[j,j]+Jwell_water
                Jo_test[j,j]=Jo_test[j,j]+Jwell_oil
                #print(Jw[j,j],Jo[j,j])
                if Well.kind[i]==0:#0-Producer 1-Injector
                    Qw[j]=Qw[j]-qw*PVT.Bw
                    Qo[j]=Qo[j]-qo*PVT.Bo
                else:# Injector well (It considers that the fluid injected is 100% water)
                    Qw[j]=Qw[j]+(Well.rates[i]/len(Well.blocks[i][0]))*PVT.Bw
                    Qo[j]=Qo[j]+0
    Q=Qw+Qo
    J=Jw+Jo
    J_test=Jw_test+Jo_test
            
    J=csr_matrix(J);  Q=csr_matrix(Q); Qw=csr_matrix(Qw); Qo=csr_matrix(Qo); 
    Jo=csr_matrix(Jo);  Jw=csr_matrix(Jw); J_test=csr_matrix(J_test);
    return J,Jw,Jo,Q, Qw,Qo,J_test,BHP, Jw_test, Jo_test

def Well_Arrays_SS(Sw,J,Q, Well,Res,PVT,Num, PP):
    Qw=np.zeros((Num.N,1))
    Q1=np.zeros((2*Num.N,1))
    J=J.toarray(); Q=Q.toarray();
    for i in Well.well_id:
        for j in Well.blocks[i][0]:
            krw,kro=BCM_RelativePermeability(Sw[j], PP)
            if Well.type[i]==1:
                if Well.Jindex!=None:
                    J[j,j]=J[j,j]#+Well.Jindex[i]
                    Q[j]=Q[j]+J[j,j]*Well.rates[i]
                    Qw[j]=Q[j]*(krw/(PVT.muw*PVT.Bw)/(krw/(PVT.muw*PVT.Bw)+kro/(PVT.muo*PVT.Bo)))
                else:
                    J[j,j]=J[j,j]#+Jindex(j,i, Well,Res,PVT,Num)
                    Q[j]=Q[j]+J[j,j]*Well.rates[i] 
                    Qw[j]=Q[j]*(krw/(PVT.muw*PVT.Bw)/(krw/(PVT.muw*PVT.Bw)+kro/(PVT.muo*PVT.Bo)))
            else: 
                if Well.kind[i]==0:#0-Producer 1-Injector
                    Q[j]=Q[j]-Well.rates[i]/(len(Well.blocks[i][0]))
                    Qw[j]=Q[j]*(krw/(PVT.muw*PVT.Bw)/(krw/(PVT.muw*PVT.Bw)+kro/(PVT.muo*PVT.Bo)))
                else:
                    Q[j]=Q[j]+Well.rates[i]/(len(Well.blocks[i][0]))
                    Qw[j]=Q[j]
                Q1[2*j]=Qw[j]
                Q1[2*j+1]=Q[j]-Qw[j]
    
    
    Qo=Q-Qw
    J=csr_matrix(J);  Q=csr_matrix(Q); Qw=csr_matrix(Qw); Qo=csr_matrix(Qo); 
    Q1=csr_matrix(Q1)
    return J,Q1, Qw,Qo
