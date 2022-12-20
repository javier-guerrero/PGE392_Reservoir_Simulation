# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""
import yaml
import numpy as np

#%% Function to read the file .yml
filename="Thomasv2.yml"
def read_reservoir_parameters(filename):
    '''
    Reads reservoir inputs from yml file and 
    stores data into python dictionaries
    
    argument:
        filename   - well activity file
    returns:
        clean_data - Python list containing well activity information
    '''
    with open(filename) as f:
        parameters = yaml.load(f, Loader=yaml.Loader)
    
    return parameters

#%%Convert dict to class
class Dict2Class(object):
    '''
    Converts a python dictionary into a class
    '''
    def __init__(self, my_dict):
        for key in my_dict:
            setattr(self, key, my_dict[key])


#%% Read parameters and separate it in classes            
parameters = read_reservoir_parameters(filename) # Reads yml file

Res = Dict2Class(parameters['reservoir'])       # Reservoir inputs
PVT = Dict2Class(parameters['pvt'])             # PVT and fluids inputs
petro = Dict2Class(parameters['petrophysical']) # Petrophysical inputs
Num = Dict2Class(parameters['numerical'])       # Numerical inputs
Well = Dict2Class(parameters['well'])           # Well inouts
BC = Dict2Class(parameters['BC'])               # Boundary conditions 
IC = Dict2Class(parameters['IC'])               # Initial Conditions


#%% Use this if you want to read this values from a txt file
if isinstance(Res.D, str):
    Depth = np.loadtxt("Thomas_Depth.txt", comments="#", unpack=False)
    Res.D = Depth
if isinstance(Res.kx, str):
    Permeability = np.loadtxt("Thomas_Permeability.txt", comments="#", unpack=False)
    Res.kx = Permeability
if isinstance(Res.phi, str):
    Porosity = np.loadtxt("Thomas_porosity.txt", comments="#", unpack=False)
    Res.phi = Porosity/100


#%% New class for Initial and Boundary Conditions

#For convertion factors and constants
class CF:
    def __init__(self):
        self.mDft_cp_to_ft3_psiday=6.33e-3

# Computes delta_x
if Num.dx == None:
    Num.dx = [np.asarray(Res.L/Num.Nx)]
else:
    Num.dx = [np.asarray(Num.dx)]

# Computes delta_x
if Num.dy == None:
    Num.dy = [np.asarray(Res.w/Num.Ny)]
else:
    Num.dy = np.asarray(Num.dy)
# Computes delta_x
if Num.dz == None:
    Num.dz = [np.asarray(Res.h/Num.Nz)]
else:
    Num.dz = np.asarray(Num.dz)


Num.N = Num.Nx * Num.Ny * Num.Nz
Num.dx = np.asarray(Num.dx)
Num.dy = np.asarray(Num.dy)
Num.dz = np.asarray(Num.dz)
# Num.t_initial = 0

#If equally spaced grid --> checks if dx, dy, dz are single constant values
if len(Num.dx) == 1:
    Num.dx = np.full(Num.N, Num.dx)
if len(Num.dy) == 1:
    Num.dy=np.full(Num.N, Num.dy)
if len(Num.dz) == 1:
    Num.dz=np.full(Num.N, Num.dz)

Res.kx = np.asarray(Res.kx)
Res.phi = np.asarray(Res.phi)

#If reservoir is Homogeneous and Isotropic
if len(Res.kx) == 1:  # Checks if perm_x is a single value 
    Res.kx = np.full(Num.N, Res.kx)
if len(Res.phi==1):     # Checks if phi is a single value 
    Res.phi = np.full(Num.N, Res.phi)
   
# To match the reservoir anisotropy in the example for HW10
Res.ky = Res.kx
Res.kz = Res.kx

Res.D = np.asarray(Res.D) + Res.Dtop  # Computes real depth of each gridblock

if len(Res.D) == 1:
    Res.D = np.full(Num.N,Res.D)


BC.value = np.asarray(BC.value)

#%% Determine direction of the well 
# This accounts only for vertical or horizontal wells (not deviated or slanted)

Well.direction=[]
for i in Well.well_id:
    if Well.x_start[i]!= Well.x_end[i]:
        direction = 'x-direction'
    elif Well.y_start[i]!= Well.y_end[i]:
        direction = 'y-direction'
    else:
        direction = 'z-direction'
    
    Well.direction.append(direction)

#%%Create a matrix with the locations of the blocks in which the weel is located reside only considers vertical and horizontal
Well.x_start = np.asarray(Well.x_start) * Res.L 
Well.x_end = np.asarray(Well.x_end) * Res.L
Well.y_start = np.asarray(Well.y_start) * Res.w 
Well.y_end = np.asarray(Well.y_end) * Res.w


#%% Determine grid location for every well resides based on the x,y,z location (start and end) of the well
Well.blocks=[]
for i in Well.well_id:
    Well.blocks.append([])
    Well.x_start[i]
    xgrid_start=np.where(np.cumsum(Num.dx)>Well.x_start[i])[0][0]
    ygrid_start=np.where(np.cumsum(Num.dy.reshape(Num.Nx, Num.Ny),axis=0)>Well.y_start[i]--Num.dy[0]/2)[0][0]
    zgrid=0
    Well.x_end[i]
    xgrid_end=np.where(np.cumsum(Num.dx)>Well.x_end[i])[0][0]
    ygrid_end=np.where(np.cumsum(Num.dy.reshape(Num.Nx, Num.Ny),axis=0)>Well.y_end[i]--Num.dy[0]/2)[0][0]
    grid_start=int(xgrid_start+ygrid_start*Num.Nx+zgrid*Num.Nx*Num.Ny)
    grid_end=int(xgrid_end+ygrid_end*Num.Nx+zgrid*Num.Nx*Num.Ny)
    if Well.direction[i]=='x-direction':
        Block=np.arange(grid_start, grid_end+1,1)
    elif Well.direction[i]=='z-direction':
        Block=np.arange(grid_start, grid_end+1,1)
    elif Well.direction[i]=='y-direction':
        Block=np.arange(grid_start, grid_end+1,Num.Nx) 
    Well.blocks[i].append(Block)
#%%
for i in Well.well_id:
    print('Well',i, 'Block', Well.blocks[i][0] )

#%%
petro.epspc = float(petro.epspc)

# import yaml
# import numpy as np
# #%% Function to read the file .yml
# filename="Thomasv2.yml"
# def read_well_parameters(filename):
#     with open(filename) as f:
#         parameters = yaml.load(f, Loader=yaml.Loader)
#     return parameters
# #%%Convert dict to class
# class Dict2Class(object):
#     def __init__(self, my_dict):
#         for key in my_dict:
#             setattr(self, key, my_dict[key])
# #%% Read parameters and separate it in classes            
# parameters=read_well_parameters(filename)
# Res=Dict2Class(parameters['reservoir']); PVT=Dict2Class(parameters['pvt']); PP=Dict2Class(parameters['petrophysical']);
# Num=Dict2Class(parameters['numerical']); Well=Dict2Class(parameters['well']);
# BC=Dict2Class(parameters['BC']); IC=Dict2Class(parameters['IC']);
# #%% Use this if you want to read this values from a txt file
# if isinstance(Res.D, str):
#     Depth = np.loadtxt("Thomas_Depth.txt", comments="#", unpack=False)
#     Res.D=Depth
# if isinstance(Res.kx, str):
#     Permeability= np.loadtxt("Thomas_Permeability.txt", comments="#", unpack=False)
#     Res.kx=Permeability
# if isinstance(Res.phi, str):
#     Porosity = np.loadtxt("Thomas_Porosity.txt", comments="#", unpack=False)
#     Res.phi=Porosity/100
# #%% New class for Initial and Boundary Conditions

# #For convertion factors and constants
# class CF:
#     def __init__(self):
#         self.mDft_cp_to_ft3_psiday=6.33e-3


# if Num.dx==None:
#     Num.dx=[np.asarray(Res.L/Num.Nx)]
# else:
#     Num.dx=[np.asarray(Num.dx)]
# if Num.dy==None:
#     Num.dy=[np.asarray(Res.w/Num.Ny)]
# else:
#     Num.dy=np.asarray(Num.dy)
# if Num.dz==None:
#     Num.dz=[np.asarray(Res.h/Num.Nz)]
# else:
#     Num.dz=np.asarray(Num.dz)


# Num.N=Num.Nx*Num.Ny*Num.Nz
# Num.dx=np.asarray(Num.dx)
# Num.dy=np.asarray(Num.dy)
# Num.dz=np.asarray(Num.dz)

# #If equally spaced grid
# if len(Num.dx)==1:
#     Num.dx=np.full(Num.N,Num.dx)
# if len(Num.dy)==1:
#     Num.dy=np.full(Num.N,Num.dy)
# if len(Num.dz)==1:
#     Num.dz=np.full(Num.N,Num.dz)

# Res.kx=np.asarray(Res.kx)
# Res.phi=np.asarray(Res.phi)

# #If reservoir is homogeneous
# if len(Res.kx)==1:
#     Res.kx=np.full(Num.N,Res.kx)
# if len(Res.phi==1):
#     Res.phi=np.full(Num.N,Res.phi)
   
# Res.ky=Res.kx
# Res.kz=Res.kx

# Res.D=np.asarray(Res.D)+Res.Dtop

# if len(Res.D)==1:
#     Res.D=np.full(Num.N,Res.D)

# BC.value=np.asarray(BC.value)


# Well.direction=[]
# for i in Well.well_id:
#     if Well.x_start[i]!=Well.x_end[i]:
#         direction='x-direction'
#     elif Well.y_start[i]!=Well.y_end[i]:
#         direction='y-direction'
#     else:
#         direction='z-direction'
#     Well.direction.append(direction)
# #%%
# #Create a matrix with the locations of the blocks in which is well reside only considers vertical and horizontal
# Well.x_start=np.asarray(Well.x_start)*Res.L; Well.x_end=np.asarray(Well.x_end)*Res.L; Well.y_start=np.asarray(Well.y_start)*Res.w; Well.y_end=np.asarray(Well.y_end)*Res.w;
# #%% Determine which grid(s) every well resides based on the x,y,z location (start and end) of the well
# Well.blocks=[]
# for i in Well.well_id:
#     Well.blocks.append([])
#     Well.x_start[i]
#     xgrid_start=np.where(np.cumsum(Num.dx)>Well.x_start[i])[0][0]
#     ygrid_start=np.where(np.cumsum(Num.dy.reshape(Num.Nx, Num.Ny),axis=0)>Well.y_start[i]--Num.dy[0]/2)[0][0]
#     zgrid=0
#     Well.x_end[i]
#     xgrid_end=np.where(np.cumsum(Num.dx)>Well.x_end[i])[0][0]
#     ygrid_end=np.where(np.cumsum(Num.dy.reshape(Num.Nx, Num.Ny),axis=0)>Well.y_end[i]--Num.dy[0]/2)[0][0]
#     grid_start=int(xgrid_start+ygrid_start*Num.Nx+zgrid*Num.Nx*Num.Ny)
#     grid_end=int(xgrid_end+ygrid_end*Num.Nx+zgrid*Num.Nx*Num.Ny)
#     if Well.direction[i]=='x-direction':
#         Block=np.arange(grid_start, grid_end+1,1)
#     elif Well.direction[i]=='z-direction':
#         Block=np.arange(grid_start, grid_end+1,1)
#     elif Well.direction[i]=='y-direction':
#         Block=np.arange(grid_start, grid_end+1,Num.Nx) 
#     Well.blocks[i].append(Block)
# #%%
# PP.epspc=float(PP.epspc)



