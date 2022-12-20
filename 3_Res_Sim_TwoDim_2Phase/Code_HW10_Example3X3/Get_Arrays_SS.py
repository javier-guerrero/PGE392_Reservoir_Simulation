"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""


#%% Import Libraries and packages
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from T_interblock import T_interblock    #for calculating transmissibility
from Block_Properties import Block_Properties


#%% Define functions
def get_arrays_SS(Num, PVT, Res, petro, BC, P, Sw):
    #Setting up matrix T, To, Tw, J, A, Act
    T = lil_matrix((2*Num.N, 2*Num.N))
    To = lil_matrix((2*Num.N, 2*Num.N))
    Tw = lil_matrix((2*Num.N, 2*Num.N))
    J = lil_matrix((Num.N, Num.N))
    A = lil_matrix((2*Num.N, 2*Num.N))
    Act = lil_matrix((2*Num.N, 2*Num.N))
    C = lil_matrix((2*Num.N, 2*Num.N))
    
    #Setting up the arrays Q, G
    Q = lil_matrix((Num.N,1))
    G = lil_matrix((2*Num.N,1))
    
    for i in range(Num.N):
        Beta_o, Beta_w, A1, ct, C1, C2, C3 = Block_Properties(i, Sw, Res, PVT, petro, Num)
        A[i,i] = A1
        Act[i,i] = A[i,i] * ct
        C[2*i,2*i] = C1
        C[2*i+1,2*i+1] = -A1 
        C[2*i,2*i+1] = A1
        C[2*i+1,2*i] = C2
        
        #Interblock Trasmisibility in x direction
        if (i+1) % Num.Nx != 0: #Not in the right boundary
            T_interblock_w, T_interblock_o = T_interblock(i,i+1, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'x', P, Sw, petro, PVT)
            Tw[2*i,2*i+2] = -T_interblock_w
            Tw[2*(i+1),2*i] = Tw[2*i,2*i+2]
            To[2*i+1,2*i+2] = -T_interblock_o
            To[2*i+3,2*i] = To[2*i+1,2*i+2]
        
        #Interblock Trasmisibility in y direction
        if i // Num.Nx > 0:
            T_interblock_w, T_interblock_o = T_interblock(i,i-Num.Nx, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'y', P, Sw, petro, PVT)
            To[i,i-Num.Nx]=-T_interblock_o
            Tw[i,i-Num.Nx]=-T_interblock_w
        
        if i // Num.Nx < Num.Ny-1:
            T_interblock_w, T_interblock_o = T_interblock(i,i+Num.Nx, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'y', P, Sw, petro, PVT)
            To[i,i+Num.Nx] = -T_interblock_o
            Tw[i,i+Num.Nx] = -T_interblock_w
        To[2*i+1,2*i] = abs(np.sum(To[2*i+1,:]))
        Tw[2*i,2*i] = abs(np.sum(Tw[2*i,:]))
        
        #Left boundary conditions P(x=0, t)
        if (i+1) % Num.Nx == 1:
            if 'Neumann' in BC.type[0]:
                Q[i]==BC.value[0]
            elif 'Dirichlet' in BC.type[0]:
                J[i,i] = 2 * T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'x', P, Sw, petro, PVT)
                Q[i] = 2 * T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky,  'x', P, Sw, petro, PVT) * BC.value[0]
       
        #Right boundary conditions P(x=L, t)
        if (i+1) % Num.Nx == 0:
            if 'Neumann' in BC.type[1]:
                Q[i] == BC.value[1]
            elif 'Dirichlet' in BC.type[1]:
                J[i,i] = 2 * T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'x', P, Sw, petro, PVT)
                Q[i] = 2 * T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'x', P, Sw, petro, PVT) * BC.value[1]
        
        #Bottom boundary conditions P(y=0, t)
        if i // Num.Nx == 0:
            if len(BC.type) > 2:
                if 'Neumann' in BC.type[2]:
                    Q[i] == BC.value[2]
                elif 'Dirichlet' in BC.type[2]:
                    J[i,i] = 2 * T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'y', P, Sw, petro, PVT.muo)
                    Q[i] = 2 * T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'y', P, Sw, petro, PVT) * BC.value[2]
            else:
                None
        
        #Top boundary conditions P(y=w, t)
        if i // Num.Nx == Num.Ny-1:
            if len(BC.type) > 2:
                if 'Neumann' in BC.type[3]:
                    Q[i] == BC.value[3]
                elif 'Dirichlet' in BC.type[3]:
                    J[i,i] = 2 * T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'y', P, Sw, petro, PVT)
                    Q[i] = 2 * T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'y', P, Sw, petro, PVT) * BC.value[3]
            else:
                None
       
    #G=(PVT.rhoosc/144)*(T@Res.D) 
    T = To + Tw
           
    #Convert matrices and arrays to sparse   
    G = csr_matrix(G) 
    Act = csr_matrix(Act) 
    A = csr_matrix(A)
    T = csr_matrix(T)
    Tw = csr_matrix(Tw) 
    To = csr_matrix(To) 
    C = csr_matrix(C)
   
    return T, Tw, To, A, Act, Q, J, G, C;
