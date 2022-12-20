# -*- coding: utf-8 -*-
"""
PGE-392K-  NUmerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% import packages and libraries
import yaml
import numpy as np


#%% Function to read the file .yml  
filename="HW4_Variables.yml"

def read_well_parameters(filename):
    '''
    Reads input parameters related to reservoir properties.
    
    argument:
        filename   - well activity file
    returns:
        parameters - Python dictonary containing well parameters
    '''
    with open(filename) as f:
        parameters = yaml.load(f, Loader=yaml.Loader)
    return parameters


#%% Convert dict to class  
class Dict2Class(object):
    def __init__(self, my_dict):
        for key in my_dict:
            setattr(self, key, my_dict[key])

#%% Read parameters and separate it in classes            
parameters = read_well_parameters(filename)
Res = Dict2Class(parameters['reservoir'])
PVT = Dict2Class(parameters['pvt'])
petro = Dict2Class(parameters['petrophysical'])
Num = Dict2Class(parameters['numerical'])
Well = Dict2Class(parameters['well'])


#%% Read the array of Depth, Permeability, and Porosity
Depth = np.loadtxt("HW4_Depth.txt", comments="#", unpack=False) + Res.Dtop
Permeability = np.loadtxt("HW4_Perm.txt", comments="#", unpack=False)
Porosity = np.loadtxt("HW4_porosity.txt", comments="#", unpack=False)

#%% Add the array values of Depth, Permeability and Porosity to the Dict
Res.k=Permeability
Res.D=Depth
Res.phi=Porosity

#%% New class for Initial and Boundary Conditions
class IC:
    def __init__(self):
        self.P=[]
class BC:
    def __init__(self):
        self.P=[]
class CF:
    def __init__(self):
        self.P=[]
        
CF.mDft_cp_to_ft3_psiday = 6.3285e-3
IC.P = 3000
BC.type= [['Neumann'],['Dirichlet']] 
BC.value = [0,1000]   

#%%  Add parameters of example 3.4 of the notes
Res.L=4000
Res.w=1000
Res.h=20
Res.k=100
Res.phi=0.2
Res.ct=1e-5
PVT.muo=5
petro.kroo=1
petro.swr=0.1
Num.N=25
Num.dx=Res.L/Num.N
Num.dt=5 #0.5/((Res.k*PP.kroo/(Res.phi*PVT.muo*Res.ct))*(1/Num.dx**2)*CF.mDft_cp_to_ft3_psiday)
Num.tf=15
Num.theta=0



