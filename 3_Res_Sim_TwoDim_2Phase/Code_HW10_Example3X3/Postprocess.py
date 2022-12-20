# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""
import matplotlib.pyplot as plt
import numpy as np
from Productivity_Index import Jindex
#%%
def plots(P, t, Res, Num, PVT, Well):
    
    OOIP=np.sum(Res.phi*Num.dx*Num.dy*Num.dz)/5.615/PVT.Bob #STB
    Pav=np.mean(P,axis=0)
    
    Qstdwell=np.zeros((len(Well.well_id),len(t)))
    for k in range(1,len(t)):
        for i in Well.well_id:
            Qwell=0
            for j in Well.blocks[i][0]:
                if P[j,k]>2000:
                    Qwell+=Jindex(j,i, Well,Res,PVT,Num)*(P[j,k]-2000)
                #print(Jindex(j,i, Well,Res,PVT,Num)*(2000-P[j,k]))
            Qstdwell[i,k]=Qwell/5.615/PVT.Bob
    
    plt.figure(0)
    x=np.linspace(0,Res.L, Num.Nx)
    y=np.linspace(0,Res.w, Num.Ny)
    X, Y = np.meshgrid(x, y)
    Z=np.reshape(P[:,0],(Num.Ny,Num.Nx))
    plt.pcolormesh(X, Y, Z, cmap='viridis', shading='auto')
    plt.colorbar(label='P [psi]');
    plt.clim(np.max(P[:,0]),np.min(P[:,-1]))
    plt.title(f'Pressure at t={t[0]} days')
    plt.xlabel('x [ft]'); plt.ylabel('y [ft]');
        
    i=np.where(Pav<=2750)[0][0]
    plt.figure(1)
    Z=np.reshape(P[:,i],(Num.Ny,Num.Nx))
    plt.pcolormesh(X, Y, Z, cmap='viridis', shading='auto')
    plt.colorbar(label='P [psi]');
    #plt.clim(np.max(P[:,0]),np.min(P[:,-1]))
    plt.title(f'Pressure Average = {round(Pav[i],0)} psi at t={t[i]} days')
    plt.xlabel('x [ft]'); plt.ylabel('y [ft]');
    
    i=np.where(Pav<=2500)[0][0]
    plt.figure(2)
    Z=np.reshape(P[:,i],(Num.Ny,Num.Nx))
    plt.pcolormesh(X, Y, Z, cmap='viridis', shading='auto')
    plt.colorbar(label='P [psi]');
    #plt.clim(np.max(P[:,0]),np.min(P[:,-1]))
    plt.title(f'Pressure Average = {round(Pav[i],0)} psi at t={t[i]} days')
    plt.xlabel('x [ft]'); plt.ylabel('y [ft]');
    
    i=np.where(Pav<=2250)[0][0]
    plt.figure(3)
    Z=np.reshape(P[:,i],(Num.Ny,Num.Nx))
    plt.pcolormesh(X, Y, Z, cmap='viridis', shading='auto')
    plt.colorbar(label='P [psi]');
    #plt.clim(np.max(P[:,0]),np.min(P[:,-1]))
    plt.title(f'Pressure Average = {round(Pav[i],0)} psi at t={t[i]} days')
    plt.xlabel('x [ft]'); plt.ylabel('y [ft]');
    
    i=2000
    plt.figure(4)
    Z=np.reshape(P[:,i],(Num.Ny,Num.Nx))
    plt.pcolormesh(X, Y, Z, cmap='viridis', shading='auto')
    plt.colorbar(label='P [psi]');
    #plt.clim(np.max(P[:,0]),np.min(P[:,-1]))
    plt.title(f'Pressure Average = {round(Pav[i],0)} psi at t={t[i]} days')
    plt.xlabel('x [ft]'); plt.ylabel('y [ft]');
    
    plt.figure(5)
    plt.plot(t,Pav)  
    plt.title('Pressure Average vs Time')
    plt.xlabel('time [days]'); plt.ylabel('P average [psi]');
    
    plt.figure(6)
    for i in Well.well_id:
        if i!=4:
            plt.plot(t[1:],Qstdwell[i,1:], label=f'Well{Well.well_id[i]}') 
    plt.title('Qstd per well with respect to time')
    plt.xlabel('time [days]'); plt.ylabel('Qstd [STB/day]');
    plt.legend()   
    
    plt.figure(7)
    plt.plot(t[1:],Qstdwell[4,1:], label=f'Well{Well.well_id[i]}') 
    plt.title('Qstd per well with respect to time')
    plt.xlabel('time [days]'); plt.ylabel('Qstd [STB/day]');
    plt.legend()  
    
    plt.figure(8)
    for i in Well.well_id:
        plt.plot(t[1:],np.full(len(t[1:]),2000), label=f'Well{Well.well_id[i]}') 
    plt.title('BHP pressure per well with respect to time')
    plt.xlabel('time [days]'); plt.ylabel('BHP [psi]');
    plt.legend()   
    
    plt.figure(9)
    TotalQstd=np.sum(Qstdwell,axis=0)
    CumulativeProduction=np.cumsum(TotalQstd)*100/OOIP
    plt.plot(t,CumulativeProduction)
    plt.title('Cumulative Oil Production (as % OOIP)')
    plt.xlabel('time [days]'); plt.ylabel('% production recovered');
    