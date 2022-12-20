# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% Import Libraries and packages
import yaml
from yaml.loader import SafeLoader
import numpy as np


#%% Define functions
def preprocess(file):
    '''
    Reads reservoir inputs from yml file and 
    stores data into python dictionaries
    
    argument:
        filename   - well activity file
    returns:
        clean_data - Python list containing well activity information
    '''
    with open(file, 'r') as thomas:
        d = yaml.load(thomas, Loader=SafeLoader)   # get .yml file
        Res, PVT, petro = d.get('reservoir'), d.get('pvt'), d.get('petrophysical') # seperate dictionaries
        Num, Well, BC = d.get('numerical'), d.get('well'), d.get('boundary_conditions')  # seperate dictionaries
        N = Num['Nx']*Num['Ny']*Num['Nz'] # Total number of unknowns
        Num['N'] = N
        Num['block'] = np.linspace(0,N-1,N, dtype=int)
        D = Grid_increment(Res['L']/Num['Nx'],Res['w']/Num['Ny'],Res['h']/Num['Nz'],Res,Num) # Call grid_increment
        Num['D'] = D # Add this to numerical dictionary
        
        for key in ['D', 'k', 'phi']:  # check for additional .txt files
            if type(Res[key]) == str:
                Res[key] = np.genfromtxt(Res[key]) # generate variables from .txt
            elif (type(Res[key]) == int or type(Res[key]) == float):
                Res[key] = np.ones(N)*Res[key]
                
        Res['kx'],Res['ky'],Res['kz'] = np.copy(Res['k'])*Res['kx_val'],np.copy(Res['k'])*Res['ky_val'],np.copy(Res['k'])*Res['kz_val']
        
        Well['array'], Well['direction'], Well['perf_blocks'] = well_blocks(Well, Num)
        Well['pwf'], Well['GOR'] = np.zeros(Well['array'].shape), np.zeros(len(Well['well_id']))
        Well['Q_o'] = np.copy(Well['pwf'])
        Well['Q_g'] = np.copy(Well['pwf'])
        Well['Q_field'] = np.copy(Well['pwf'])
        Well['J'] = np.zeros((len(Well['well_id']), max(Well['perf_blocks'])))
        Well['Jw'] = np.zeros((len(Well['well_id']), max(Well['perf_blocks'])))
        Well['Jo'] = np.zeros((len(Well['well_id']), max(Well['perf_blocks'])))
        Well['Jg'] = np.zeros((len(Well['well_id']), max(Well['perf_blocks'])))
        Well['wc'] = np.zeros((len(Well['well_id']), max(Well['perf_blocks'])))
        return Res, PVT, petro, Num, Well, BC # return dictionaries and cell centers
    
def well_blocks(well, numeric):
    '''
    Determine grid location for every well resides 
    based on the x,y,z location (start and end) of the well

    '''
    i_s, i_e = np.zeros(len(well['well_id'])), np.zeros(len(well['well_id']))
    for pos in zip(well['x_start'], well['y_start']):
        xp, yp = pos[0], pos[1] # get x,y coordinates of start position
        ind = well['x_start'].index(xp) # get vector index
        j = int(np.round_(xp*numeric['Nx']))+1 # x-pos in j,k,l space
        k = int(np.round_(yp*numeric['Ny']))+1 # y-pos in j,k,l space
        i_s[ind] = j + (k-1)*numeric['Nx'] - 1 # "i" block for start location
        
    for pos in zip(well['x_end'], well['y_end']):
        xp, yp = pos[0], pos[1] # get x,y coordinates of start position
        ind = well['x_end'].index(xp) # get vector index
        j = int(np.round_(xp*numeric['Nx']))+1 # x-pos in j,k,l space
        k = int(np.round_(yp*numeric['Ny']))+1 # y-pos in j,k,l space
        i_e[ind] = j + (k-1)*numeric['Nx'] - 1 # "i" block for end location

    di = i_e - i_s # compute difference of start/end block
    dir_ = np.empty(len(well['well_id']), dtype=str) # Initialize direction array
    dir_[di<max(numeric['Nx'],numeric['Ny'])] = 'x'
    dir_[di>max(numeric['Nx'],numeric['Ny'])] = 'y'
    dir_[di==0] = 'z'
    
    len_ = np.empty(len(well['well_id']))
    len_[dir_=='z'] = int(1)
    
    array = np.zeros((len(well['well_id']),int(max(len_)))) # Initialize well array
    array[dir_=='z',0] = i_s[dir_=='z']
    perf = np.count_nonzero(array, axis=1)
    perf[perf==0] = 1
    return array, dir_, perf


def Grid_increment(dx,dy,dz,res,numeric):
    N = numeric['N']
    dx = np.ones(numeric['Nx'])*dx
    dy = np.ones(numeric['Ny'])*dy
    dx_ = np.tile(dx, (1, len(dy))).T
    dy_ = np.reshape(np.repeat(dy,len(dx)),dx_.shape)
    dz_ = np.ones((N,1))*res['h']
    D = np.concatenate((dx_,dy_, dz_), axis=1)
    
    return D