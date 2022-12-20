
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% Import libraries
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import inv
import matplotlib.pyplot as plt
from time import perf_counter, localtime, strftime

#Import functions and depedences from other files
from Get_Arrays import get_arrays                              #Get matrix T, A, Act, Q, J
from Get_Arrays_SS import get_arrays_SS                        #Get matrix T, A, Act, Q, J
from Preprocess import Num, Res, PVT, petro, Well, IC, BC, CF     #Import input data
from Initialization import initialization                      #Calculate Pw, Po, Sw and So
from Well_Arrays import Well_Arrays, Well_Arrays_SS          #Update arrays when we have wells
from Productivity_Index import Jindex
from BCM_kr import BCM_RelativePermeability, CapillaryPressure
from T_interblock import T_interblock 
from Initialization import initialization 

#%%
for i in Well.well_id:
        for j in Well.blocks[i][0]:
            print('Well #:',i, 'Perforated Blocks:', j)
#%% Initialization
Sw, So, Pw, Po,rhoo, rhow = initialization(Res.D, petro, PVT, Res)
PVT.rhoo = rhoo
PVT.rhow = rhow

#Plots of initial pressure and initial saturation
x = np.linspace(0,Res.L, Num.Nx)
y = np.linspace(0,Res.w, Num.Ny)
X, Y = np.meshgrid(x, y)
Z = np.reshape(Po,(Num.Ny,Num.Nx))

plt.figure(0)
plt.style.use('seaborn')
plt.pcolormesh(X, Y, Z, cmap='jet')
plt.colorbar(label='P [psi]')
plt.title('Initial Pressure', fontsize=16)
plt.xlabel('x [ft]', fontsize=14)
plt.ylabel('y [ft]', fontsize=14)
plt.grid()
#%%
plt.figure(1)
plt.style.use('seaborn')
Z1 = np.reshape(Sw,(Num.Ny,Num.Nx))
plt.pcolormesh(X, Y, Z1, cmap='winter_r')
plt.colorbar(label='Sw');
plt.title('Initial Saturation', fontsize=16)
plt.xlabel('x [ft]', fontsize=14)
plt.ylabel('y [ft]', fontsize=14)
plt.grid()

#%%
#Calculate relative permeability#krw,kro=BCM_RelativePermeability(Sw, PP)
Swo=np.reshape(Sw,(Num.N,1))
Po=np.reshape(Po,(Num.N,1))
#
Num.t_initial=0
Num.t_final=9116
Num.dt=1

#%% Solve P from a equation of type Ax=b
start = perf_counter()    
print(strftime("%H:%M:%S", localtime()), ' ::::::::Simulation Starts:::::::')
total_time_steps =  int((365+365*24/10))    
time = np.zeros((total_time_steps+1))               # initializing time vector
P_matrix= np.zeros((Num.N, total_time_steps+1))     # matrix to save P(t) 
Sw_matrix= np.zeros((Num.N, total_time_steps+1))     # matrix to save Sw(t)
Qw_matrix= np.zeros((Num.N, total_time_steps+1))     # matrix to save Sw(t)
Qo_matrix= np.zeros((Num.N, total_time_steps+1))     # matrix to save Sw(t)
Q_matrix= np.zeros((Num.N, total_time_steps+1))     # matrix to save Sw(t)
Q_perwell=np.zeros((5, total_time_steps+1))
Qw_perwell=np.zeros((5, total_time_steps+1))
Qo_perwell=np.zeros((5, total_time_steps+1))
Pwf_perwell=np.zeros((5,total_time_steps+1))
Pwf_matrix= np.zeros((Num.N, total_time_steps+1)) 
q=np.zeros(Num.N)
qw=np.zeros(Num.N)
qo=np.zeros(Num.N)
Pwf=np.zeros(Num.N)
Sw_hyst=np.empty((Num.N,2))
Sw_hyst[:,0]=Swo[:,0]
#Initial Pressure and initial Saturation
P_matrix[:,0]=Po[:,0]; Pn=Po; Sw_matrix[:,0]=Swo[:,0]; Swn=Swo; 
for i in range(1,total_time_steps+1):
    time[i]= time[i-1]+ Num.dt 
    if Num.MF_method == 'IMPES':
        T,Tw,To, A, Act, Q, J,Jw,Jo, G, BetaW, Jw, Jo, Qw, Qo, G1, Pc=get_arrays(Num, PVT, Res,petro, BC, Pn, Swn, Sw_hyst)
        J,Jw,Jo,Q, Qw, Qo, J_test,BHP,Jw_test, Jo_test=Well_Arrays(Swn,J,Jw,Jo,Q, Qw, Qo, Well,Res,PVT,Num, petro,Pn, Qo_matrix[:,i],Qw_matrix[:,i])
        Q=Qw+Qo+J@(BHP) #Create a vector BHP
        J_test=J_test.toarray()
        Q=csr_matrix(Q)
        J=csr_matrix(J)
        Pn=csr_matrix(Pn)
        P = np.transpose([spsolve(Act+(T+J) ,Act@Pn+Q+G)]) 
       
        #This is to calculate the flowrate of oil and water
        for a in Well.well_id:
            for b in Well.blocks[a][0]:
                if Well.type[a]==1: #BHP
                    qw[b]=Jw_test[b,b]*(P[b,0]-BHP[b])
                    qo[b]=Jo_test[b,b]*(P[b,0]-BHP[b])
                    q[b]=-(qw[b]+qo[b])
                elif Well.type[a]==0: #Constant Rate
                    qw[b]=Qw[b,0]
                    qo[b]=Qo[b,0]
                    q[b]=qw[b]+qo[b]
        

                    
        for j in range(len(P)):
            if J_test[j,j]!=0:
                Pwf[j]=P[j]-abs(q[j])/J_test[j,j]
            else:
                Pwf[j]=P[j]
        #print(Well.type)
        for well in Well.well_id:
            Pwf_perwell[well,i-1]=np.average(Pwf[Well.blocks[well][0]])
            Q_perwell[well,i-1]=np.sum(qo[Well.blocks[well][0]]/PVT.Bo+qw[Well.blocks[well][0]]/PVT.Bo)
            Qo_perwell[well,i-1]=np.sum(qo[Well.blocks[well][0]]/PVT.Bo)
            Qw_perwell[well,i-1]=np.sum(qw[Well.blocks[well][0]]/PVT.Bw)
            
            if well==0:
                if np.average(Pwf[Well.blocks[well][0]])<800:
                    Well.type[well]=1
                    Well.rates[well]=800
                if Q_perwell[well,i-1]>5615:
                    Well.type[well]=0
                    Well.rates[well]=5615
            elif well==1 and time[i]<=365:
                if np.average(Pwf[Well.blocks[well][0]])<1000:
                    Well.type[well]=1
                    Well.rates[well]=1000
            elif well==1:
                if np.average(Pwf[Well.blocks[well][0]])>7000:
                     Well.type[well]=1
                     Well.rates[well]=7000
                if np.average(Pwf[Well.blocks[well][0]])>7000:
                    Well.type[well]=1
                    Well.rates[well]=7000
            elif well==2 and time[i]<5116:
                if np.average(Pwf[Well.blocks[well][0]])<900:
                    Well.type[well]=1
                    Well.rates[well]=900
                if Q_perwell[well,i-1]>5615:
                    Well.type[well]=0
                    Well.rates[well]=5615
            elif well==2:
                if np.average(Pwf[Well.blocks[well][0]])>7000:
                     Well.type[well]=1
                     Well.rates[well]=7000
                    
            elif well==3:
                if np.average(Pwf[Well.blocks[well][0]])<650:
                    Well.type[well]=1
                    Well.rates[well]=650
                if Q_perwell[well,i-1]>5615:
                    Well.type[well]=0
                    Well.rates[well]=5615
               
            elif well==4 and time[i]<=365:
                if np.average(Pwf[Well.blocks[well][0]])<1000:
                    Well.type[well]=1
                    Well.rates[well]=1000
            elif well==4:
                if np.average(Pwf[Well.blocks[well][0]])>7000:
                     Well.type[well]=1
                     Well.rates[well]=7000
      
        if i==366:
            Well.kind[1]=1
            Well.rates[1]=1500*5.62
            Well.type[1]=0
            Well.kind[4]=1
            Well.rates[4]=1500*5.62
            Well.type[4]=0
            Num.dt=10
            print(Well.type)
            print(Well.kind)
        if time[i]==5116:
            Well.kind[2]=1
            Well.rates[2]=1000*5.62
            Well.type[2]=0
            
            
            
        
        Pc=np.reshape(Pc.toarray(),(Num.N,1))
        Sw=(Swn+inv(A)*(-Tw@(P-Pc-(PVT.rhow/144)*Res.D.reshape(Num.N,1))+Qw+Jw@(BHP-P)))-(BetaW*Swn)*(PVT.cw+PVT.cf)@(P-Pn)
        Sw[Sw > 1.0] = 1.0
        P_matrix[:,i] = P[:,0]; 
        Qw_matrix[:,i]=np.reshape(Qw.toarray(), (Num.N))
        Qo_matrix[:,i]=np.reshape(Qo.toarray(), (Num.N))
        Q_matrix[:,i]=np.reshape(Q.toarray(), (Num.N))
        Sw_matrix[:,i] = Sw.T;
        Pwf_matrix[:,i] = Pwf; 
        Pn=P
        Swn=np.asarray(Sw)
        #To calculate hysteresis
        for j in range(0, Num.N):
            if Sw_matrix[j,i]> Sw_matrix[j,i-1] and Sw_hyst[j,1] == 0:  # [i,1] is a flag
                Sw_hyst[j,0] = Sw_matrix[j,i]
                Sw_hyst[j,1] = 1.0

            elif Sw_matrix[j,i] < Sw_matrix[j,i-1]:
                Sw_hyst[j,0] = Sw_matrix[j,i]    
        print(strftime("%H:%M:%S", localtime()),i,' Time', time[i], 'days', 'Simulation Progress:', round((time[i])*100/(Num.t_final-Num.t_initial),2), '%')
end = perf_counter() 
print("Time elapsed during the calculation:", round(end - start,2), 'sec')

#%% Plot Pressure at 1 year
plt.figure(2)
plt.style.use('seaborn')
Z3=np.reshape(P_matrix[:,365],(Num.Ny,Num.Nx))
plt.pcolormesh(X, Y, Z3, cmap='jet')
plt.colorbar(label='P [psi]');
plt.title('Pressure at 1 year', fontsize=16)
plt.xlabel('x [ft]', fontsize=14) 
plt.ylabel('y [ft]', fontsize=14)
plt.grid()

#%% Plot Pressure at 14 years
plt.figure(3)
plt.style.use('seaborn')
Z3=np.reshape(P_matrix[:,841],(Num.Ny,Num.Nx))
plt.pcolormesh(X, Y, Z3, cmap='jet')
plt.colorbar(label='P [psi]');
plt.title('Pressure at 14 years', fontsize=16)
plt.xlabel('x [ft]', fontsize=14) 
plt.ylabel('y [ft]', fontsize=14)
plt.grid()

#%% Plot Pressure at 25 years
plt.figure(4)
plt.style.use('seaborn')
Z3=np.reshape(P_matrix[:,i-1],(Num.Ny,Num.Nx))
plt.pcolormesh(X, Y, Z3, cmap='jet')
plt.colorbar(label='P [psi]');
plt.title('Pressure at 25 years', fontsize=16)
plt.xlabel('x [ft]', fontsize=14)
plt.ylabel('y [ft]', fontsize=14)
plt.grid()
#%% Plot Sw at 25 years (final)
plt.figure(5)
plt.style.use('seaborn')
Z3=np.reshape(Sw_matrix[:,i-1],(Num.Ny,Num.Nx))
plt.pcolormesh(X, Y, Z3, cmap='winter_r')
plt.colorbar(label='P [psi]');
plt.title('Water Saturation at 25 years', fontsize=16)
plt.xlabel('x [ft]', fontsize=14)
plt.ylabel('y [ft]', fontsize=14)
plt.grid()

#%% Plot Sw at 1 years 
plt.figure(6)
plt.style.use('seaborn')
Z3=np.reshape(Sw_matrix[:,365],(Num.Ny,Num.Nx))
plt.pcolormesh(X, Y, Z3, cmap='winter_r')
plt.colorbar(label='Sw');
plt.title('Water Saturation at 1 year', fontsize=16)
plt.xlabel('x [ft]', fontsize=14)
plt.ylabel('y [ft]', fontsize=14)
plt.grid()

#%% Plot Sw at 14 years 
plt.figure(7)
plt.style.use('seaborn')
Z3=np.reshape(Sw_matrix[:,841],(Num.Ny,Num.Nx))
plt.pcolormesh(X, Y, Z3, cmap='winter_r')
plt.colorbar(label='Sw');
plt.title('Water Saturation at 14 years', fontsize=16)
plt.xlabel('x [ft]', fontsize=14)
plt.ylabel('y [ft]', fontsize=14)
plt.grid()

#%% Plot BHP vs time
plt.figure(8)
plt.style.use('seaborn')
#Pwf_matrix[:,0]=P_matrix[:,1]
plt.plot(time,Pwf_perwell[0,:], label='Well-0')
plt.plot(time,Pwf_perwell[1,:], label='Well-1')
plt.plot(time,Pwf_perwell[2,:], label='Well-2')
plt.plot(time,Pwf_perwell[3,:], label='Well-3')
plt.plot(time,Pwf_perwell[4,:], label='Well-4')
plt.title('Well BHP', fontsize=16)
plt.xlabel('Time [days]', fontsize=14)
plt.ylabel('Well BHP (psi)', fontsize=14)
plt.xlim(0,364)
plt.ylim(0,3000)
plt.legend()
#%% Plot Liquid Rate vs time
plt.figure(9)
plt.style.use('seaborn')
plt.plot(time,np.abs(Q_perwell[0,:])/5.61, label='Well-0')
plt.plot(time,np.abs(Q_perwell[1,:])/5.61, label='Well-1')
plt.plot(time,np.abs(Q_perwell[2,:])/5.61, label='Well-2')
plt.plot(time,np.abs(Q_perwell[3,:])/5.61, label='Well-3')
plt.plot(time,np.abs(Q_perwell[4,:])/5.61, label='Well-4')
plt.title('Liquid Rate SC (Oil + Water) [bbl/day]', fontsize=16)
plt.xlabel('Time [days]', fontsize=14)
plt.ylabel('Liquid Rate SC (BBL/day)', fontsize=14)
# plt.xlim(0,364)
plt.legend()
#%% PLot Oil rate vs time
plt.figure(10)
plt.style.use('seaborn')
plt.plot(time,np.abs(Qo_perwell[0,:])/5.61, label='Well-0')
plt.plot(time,np.abs(Qo_perwell[2,:])/5.61, label='Well-2')
plt.plot(time,np.abs(Qo_perwell[3,:])/5.61, label='Well-3')
plt.title('Oil Rate SC [bbl/day]', fontsize=16)
plt.xlabel('Time [days]', fontsize=14)
plt.ylabel('Oil rate', fontsize=14)
#plt.xlim(0,364)
plt.legend()

#%% Plot Water cut vs time
wc=np.divide(Qw_perwell,Q_perwell)*100

plt.figure(11)
plt.style.use('seaborn')
plt.plot(time/365,np.abs(wc[0,:]), label='Well-0')
plt.plot(time/365,np.abs(wc[1,:]), label='Well-1')
plt.plot(time/365,np.abs(wc[2,:]), label='Well-2')
plt.plot(time/365,np.abs(wc[3,:]), label='Well-3')
plt.plot(time/365,np.abs(wc[4,:]), label='Well-4')
plt.title('Water Cut', fontsize=14)
plt.xlabel('Time [days]', fontsize=14)
plt.ylabel('% water cut', fontsize=14)
plt.ylim(0,100)
#plt.xlim(0,364)
plt.legend()
#%%
J=J.toarray()