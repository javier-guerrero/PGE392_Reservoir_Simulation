# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""


#%% Import Libraries and Packages
import numpy as np
import time
import matplotlib.pyplot as plt
plt.style.use('seaborn')

from preprocess import preprocess
from BCM_kr import RelativePerm_3Phase_StoneI
from Block_properties import block_properties
from PVT_props_fns_pressure import P_dependent_properties
from T_interblock import T_Interblock, potentials
from Well_Arrays import Well_Arrays
from Grid_Arrays import Grid_Arrays
from Solver import Solve, So_eqn, Sw_eqn

#%% Preprocessing
file = "Thomas.yml" # Import YML file
res, pvt, petro, numeric, well, bcs = preprocess(file) # Call preprocess function


#%%  Initiall conditions
poi = bcs['P_init'] * np.ones(numeric['N'])
Swi = petro['Swi'] * np.ones(numeric['N'])
Soi =  1 - Swi
Sgi = np.zeros(numeric['N'])


#%%  Instantiate matrices and vectors
n_step = int(numeric['t_final'] / numeric['dt']) + 1 # total number of time steps
tt = np.zeros(n_step) # initialize time
t, k = numeric['dt'], 1 # initialize time and increment
po = np.zeros((numeric['N'],n_step)) # initial condition
Sw = np.zeros((numeric['N'],n_step)) # initial condition
So = np.zeros((numeric['N'],n_step)) # initial condition
po[:,0] = poi # initial condition
Sw[:,0] = Swi # initial condition
So[:,0] = Soi # initial condition
qqo = np.zeros(n_step)
pwf = np.zeros((len(well['well_id']), n_step))
gor = np.zeros((len(well['well_id']), n_step))


#%% Main = time loop
print('\n --------Simulation Starts -----------')
start = time.time()
while t <= numeric['t_final']:
    po_old, tt[k] = np.reshape(po[:,k-1],(numeric['N'],1)), t # set old pressures and time value
    Sw_old = np.reshape(Sw[:,k-1],(numeric['N'],1)) # set old water saturation
    So_old = np.reshape(So[:,k-1],(numeric['N'],1)) # set old oil saturation
    Sg_old = 1 - Sw_old - So_old # set old gas saturation
    print('Timestep: ', tt[k], 'Simulation Progress: ', round(k / n_step *100, 2), '[%]')
    
    krw, kro, krg, _, _ = RelativePerm_3Phase_StoneI(petro,Sw_old,So_old,Sg_old,int(pvt['phase'])) # call relative permeabilities
    A, ct, betaw, betao, betag, pc, C1, C2, C3, Bw, Bo, Bg, z, co_star = block_properties(pvt,res,petro,numeric,po_old,Sw_old,So_old,Sg_old) # call block properties
    potential = potentials(pvt,res,numeric,po_old,pc,Bw,Bo,Bg,z) # call potentail function
    T, A, J, Q, G, Tw, To, Tg, Jw, Jo, Jg, Qw, Qo, Qg, ct = Grid_Arrays(A,ct,Sw_old,numeric,res,petro,pvt,bcs,krw,kro,krg,potential,betaw,betao,betag,Bw,Bo,Bg,po_old,pc,z) # call grid arrays
    Q, J, Qw, Qo, Qg, Jw, Jo, Jg, ro = Well_Arrays(well,res,pvt,petro,numeric,Q,J,Qw,Qo,Qg,Jw,Jo,Jg,Sw_old,So_old,betaw,betao,betag,Bw,Bo,Bg,po_old,tt[k],krw,kro,krg) # call well arrays
    po[:,k] = np.reshape(Solve(numeric, T, A, J, Q, G, po_old), (numeric['N'],)) # solve pressure equation
    Sw[:,k] = Sw_eqn(Sw_old,A,ct,betaw,Tw,Qw,Jw,po[:,k],po_old,numeric,pvt)
    So[:,k] = So_eqn(So_old,A,ct,betao,To,Qo,Jo,po[:,k],po_old,numeric,co_star,pvt)
    
    qqo[k] = ro[0]
    pwf[0,k], pwf[1,k], pwf[2,k] = well['pwf'][0], well['pwf'][1], well['pwf'][2]
    gor[0,k], gor[1,k], gor[2,k] = well['GOR'][0], well['GOR'][1], well['GOR'][2]
    
    t += numeric['dt'] # increment time
    k += 1 # increment index

end = time.time()
delta_t = (end - start) / 60
print('\n Simulation Complete - Elapsed time: ', round(delta_t,2), ' minutes')


#%% BHP  Plot
fig = plt.figure()
fig.set_dpi(300)
ax = plt.subplot(111)
ax.set_xlim([0,n_step])
ax.set_ylim(0,1100)
ax.set_xlabel('Time [days]',fontsize=14)
ax.set_ylabel('BHP [psi]',fontsize=14)
ax.set_title('Well Bottomhole Pressure vs. Time', fontsize=16, fontweight='bold')
plt.plot(tt[1:n_step],pwf[0,1:n_step], color='green', label='Prod-1')
plt.plot(tt[730:n_step],pwf[1,730:n_step], 'r-', label='Inj-1')
plt.plot(tt[730:n_step],pwf[2,730:n_step], 'b--', label='Inj-2')
plt.plot(tt[1:730],pwf[1,1:730], 'r-')
plt.plot(tt[1:730],pwf[2,1:730], 'b--')
plt.plot((730,730), (-100,1200), 'k--')
plt.text(800,800,'$Injection$ $Begins$', fontsize=12, fontweight='bold')
ax.legend(frameon=True,fontsize=14)
plt.show()

#%% Oil Rate Plot
fig = plt.figure()
fig.set_dpi(300)
ax = plt.subplot(111)
ax.set_xlabel('Time [day]',fontsize=14)
ax.set_ylabel('Surface Oil Rate [STB]',fontsize=14)
ax.set_title('Oil Rate vs. Time', fontsize=16, fontweight='bold')
plt.plot(tt[1:n_step],qqo[1:n_step]/5.615, color='green')#, label='$Q_o$')
#plt.plot(tt,qqg, label='$Q_g^{sc}$')
# ax.legend(frameon=False, fontsize=14)
plt.show()

#%% GOR  Plot
fig = plt.figure()
fig.set_dpi(300)
ax = plt.subplot(111)
ax.set_xlim([0,n_step])
ax.set_xlabel('Time [days]',fontsize=14)
ax.set_ylabel('Gas Oil Ratio [-]', fontsize=14)
ax.set_title('Gas Oil Ratio [GOR]', fontsize=16, fontweight='bold')
ax.set_ylim([0,300])
plt.plot(tt[1:n_step],gor[0,1:n_step]/5.615, color='red')#, label='GOR')
# ax.legend(frameon=False, fontsize=14)
plt.show()
