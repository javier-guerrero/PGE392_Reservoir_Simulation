# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""
#%% Import Libraries and Packages
import yaml
import numpy as np

#%% Function to read the file .yml
filename="Thomas.yml"

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
#%%Convert dict to class
class Dict2Class(object):
    def __init__(self, my_dict):
        for key in my_dict:
            setattr(self, key, my_dict[key])
#%% Read parameters and separate it in classes            
parameters=read_well_parameters(filename)
Res = Dict2Class(parameters['reservoir'])
PVT = Dict2Class(parameters['pvt'])
PP = Dict2Class(parameters['petrophysical'])
Num = Dict2Class(parameters['numerical'])
Well = Dict2Class(parameters['well'])
BC = Dict2Class(parameters['BC'])
IC = Dict2Class(parameters['IC'])

#%% Use this if you want to read this values from a txt file
if isinstance(Res.D, str):
    Depth = np.loadtxt("Thomas_Depth_shifted.txt", comments="#", unpack=False)
    Res.D=Depth
if isinstance(Res.kx, str):
    Permeability= np.loadtxt("Thomas_Perm.txt", comments="#", unpack=False)
    Res.kx=Permeability
if isinstance(Res.phi, str):
    Porosity = np.loadtxt("Thomas_porosity.txt", comments="#", unpack=False)
    Res.phi=Porosity/100
#%% New class for Initial and Boundary Conditions

#For convertion factors and constants
class CF:
    def __init__(self):
        self.mDft_cp_to_ft3_psiday=6.33e-3


if Num.dx==None:
    Num.dx=[np.asarray(Res.L/Num.Nx)]
else:
    Num.dx=[np.asarray(Num.dx)]
if Num.dy==None:
    Num.dy=[np.asarray(Res.w/Num.Ny)]
else:
    Num.dy=np.asarray(Num.dy)
if Num.dz==None:
    Num.dz=[np.asarray(Res.h/Num.Nz)]
else:
    Num.dz=np.asarray(Num.dz)


Num.N=Num.Nx*Num.Ny*Num.Nz
Num.dx=np.asarray(Num.dx)
Num.dy=np.asarray(Num.dy)
Num.dz=np.asarray(Num.dz)

#If equally spaced grid
if len(Num.dx)==1:
    Num.dx=np.full(Num.N,Num.dx)
if len(Num.dy)==1:
    Num.dy=np.full(Num.N,Num.dy)
if len(Num.dz)==1:
    Num.dz=np.full(Num.N,Num.dz)

Res.kx=np.asarray(Res.kx)
Res.phi=np.asarray(Res.phi)

#If reservoir is homogeneous
if len(Res.kx)==1:
    Res.kx=np.full(Num.N,Res.kx)
if len(Res.phi==1):
    Res.phi=np.full(Num.N,Res.phi)
   
Res.ky=1*Res.kx
Res.kz=1*Res.kx

Res.D=np.asarray(Res.D)+Res.Dtop

if len(Res.D)==1:
    Res.D=np.full(Num.N,Res.D)

BC.value=np.asarray(BC.value)


Well.direction=[]
for i in Well.well_id:
    if Well.x_start[i]!=Well.x_end[i]:
        direction='x-direction'
    elif Well.y_start[i]!=Well.y_end[i]:
        direction='y-direction'
    else:
        direction='z-direction'
    Well.direction.append(direction)

#%% Create a matrix with the locations of the blocks in which is well reside only considers vertical and horizontal
Well.x_start = np.asarray(Well.x_start)*Res.L
Well.x_end = np.asarray(Well.x_end)*Res.L
Well.y_start = np.asarray(Well.y_start)*Res.w
Well.y_end = np.asarray(Well.y_end)*Res.w

#%% Determine which grid(s) every well resides based on the x,y,z location (start and end) of the well
Well.blocks=[]
for i in Well.well_id:
    Well.blocks.append([])
    Well.x_start[i]
    xgrid_start=np.where(np.cumsum(Num.dx)>=Well.x_start[i])[0][0]
    ygrid_start=np.where(np.cumsum(Num.dy.reshape(Num.Nx, Num.Ny),axis=0)>=Well.y_start[i])[0][0]
    zgrid=0
    Well.x_end[i]
    xgrid_end=np.where(np.cumsum(Num.dx)>=Well.x_end[i])[0][0]
    ygrid_end=np.where(np.cumsum(Num.dy.reshape(Num.Nx, Num.Ny),axis=0)>=Well.y_end[i])[0][0]
    grid_start=int(xgrid_start+ygrid_start*Num.Nx+zgrid*Num.Nx*Num.Ny)
    grid_end=int(xgrid_end+ygrid_end*Num.Nx+zgrid*Num.Nx*Num.Ny)
    if Well.direction[i]=='x-direction':
        Block=np.arange(grid_start, grid_end+1,1)
    elif Well.direction[i]=='z-direction':
        Block=np.arange(grid_start, grid_end+1,1)
    elif Well.direction[i]=='y-direction':
        Block=np.arange(grid_start, grid_end+1,1) 
    Well.blocks[i].append(Block)

for i in Well.well_id:
    print('Well',i, 'is drilled in Block', Well.blocks[i][0] )



