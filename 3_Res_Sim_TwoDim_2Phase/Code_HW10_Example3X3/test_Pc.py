#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 12:03:24 2022

@author: javier
"""

import numpy as np
from Preprocess import Num, Res, PVT, PP, Well, IC, BC, CF     #Import input data
from BCM_kr import BCM_RelativePermeability, CapillaryPressure
import matplotlib.pyplot as plt

#%%
PP.Swr = 0.2
PP.Sor = 0.2

Sw = np.linspace(PP.Swr, (1 - PP.Sor), 61)
# print(Sw)

# Pcow,dPcow = CapillaryPressure(PP, Sw)
                               
# Pcow2,dPcow2 = CapillaryPressure(PP, (Sw + 0.001))

# Sw2 = Sw + 0.001

# Pcow3,dPcow3 = CapillaryPressure(PP, Sw2)

# delta_pc = Pcow2 - Pcow3
# delta_dpc = dPcow2 - dPcow3

#%%

plt.figure(1)
plt.plot(Sw, Pcow)
# plt.plot(Sw, Pcow2, 'r--')
# plt.xlim([0,1])
# plt.ylim([0,1])
plt.xlabel('Sw')
plt.ylabel('Pc')
plt.grid()


plt.figure(2)
plt.plot(Sw, dPcow)
# plt.plot(Sw, dPcow2, 'r--')
# plt.xlim([0,1])
# plt.ylim([0,1])
plt.xlabel('Sw')
plt.ylabel('dPc')
plt.grid()


#%%
for i in range(1, 6):
    print(i)
    
#%%
Sw8 = BCM_RelativePermeability(0.570605, petro)
print(Sw8)

                                 