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
from BCM_kr import BCM_RelativePermeability, CapillaryPressure

#%% Define Functions
def grid_arrays(Num, PVT, Res, petro, BC, P, Sw, Sw_hyst):
    #Setting up matrix T, To, Tw, J, A, Act
    T = lil_matrix((Num.N, Num.N))  # Total transmissibility matrix
    To = lil_matrix((Num.N, Num.N)) # Transmissibility matrix for oil phase
    Tw = lil_matrix((Num.N, Num.N)) # Transmissibility matrix for water phase
    Jw = lil_matrix((Num.N, Num.N))
    Jo = lil_matrix((Num.N, Num.N))
    J = lil_matrix((Num.N, Num.N))  # J matrix
    A = lil_matrix((Num.N, Num.N))  # Accumulation matrix
    Act = lil_matrix((Num.N, Num.N))# Accumulation * compressibility matrix
    d11 = lil_matrix((Num.N, Num.N))
    d12 = lil_matrix((Num.N, Num.N)) 
    
    #Setting up the arrays Q, G
    Q = lil_matrix((Num.N, 1))   # Rates vector
    Qo = lil_matrix((Num.N,1))
    Qw = lil_matrix((Num.N,1))
    G = lil_matrix((Num.N, 1))   # Gravity vextor
    
    for i in range(Num.N):
        Pcow, dPcow = CapillaryPressure(petro, Sw, Sw_hyst[i,0])
        Beta_o, Beta_w, A1, ct, C1, C2, C3 = Block_Properties(i, Sw, Res, PVT, petro, Num, dPcow[i])
        A[i,i] = A1
        Act[i,i] = A[i,i] * ct
        
        #Interblock Trasmisibility in x direction
        if (i + 1) % Num.Nx != 0: #Not in the right boundary
            T_interblock_w, T_interblock_o = T_interblock(i, i+1, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'x', P, Sw, petro, PVT, Res)
            To[i,i+1] = -T_interblock_o
            To[i+1,i] = To[i,i+1]
            Tw[i,i+1] = -T_interblock_w
            Tw[i+1,i] = Tw[i,i+1]
        
        #Interblock Trasmisibility in y direction
        if i // Num.Nx > 0:  
            T_interblock_w, T_interblock_o = T_interblock(i,i-Num.Nx, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'y', P, Sw, petro, PVT, Res)
            To[i,i-Num.Nx] = -T_interblock_o
            Tw[i,i-Num.Nx] = -T_interblock_w
        
        if i // Num.Nx < Num.Ny-1:
            T_interblock_w, T_interblock_o = T_interblock(i,i+Num.Nx, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'y', P, Sw, petro, PVT, Res)
            To[i,i+Num.Nx] = -T_interblock_o
            Tw[i,i+Num.Nx] = -T_interblock_w
        
        # Diagonal To and Tw matrix
        To[i,i] = abs(np.sum(To[i,:]))
        Tw[i,i] = abs(np.sum(Tw[i,:]))
        
        
        #Left boundary conditions P(x=0, t)
        if (i+1) % Num.Nx == 1:
            if 'Neumann' in BC.type[0]:
                Q[i] == BC.value[0]
            elif 'Dirichlet' in BC.type[0]:
                T_interblock_w, T_interblock_o = T_interblock(i, i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'x', P, Sw, petro, PVT, Res)
                Jw[i,i] = 2 * T_interblock_w
                Jo[i,i] = 2 * T_interblock_o
                Qo[i] = 2 * T_interblock_o * BC.value[0]
                Qw[i,i] = 2 * T_interblock_w * BC.value[0]
        
        #Right boundary conditions P(x=L, t)
        if (i+1)%Num.Nx == 0:
            if 'Neumann' in BC.type[1]:
                Q[i] == BC.value[1]
            elif 'Dirichlet' in BC.type[1]:
                T_interblock_w, T_interblock_o = T_interblock(i, i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'x', P, Sw, petro, PVT, Res)
                Jw[i,i] = 2 * T_interblock_w
                Jo[i,i] = 2 * T_interblock_o
                Qo[i] = 2 * T_interblock_o * BC.value[1]
                Qw[i,i] = 2 * T_interblock_w * BC.value[1]
        
        #Bottom boundary conditions P(y=0, t)
        if i // Num.Nx == 0:
            if len(BC.type) > 2:
                if 'Neumann' in BC.type[2]:
                    Q[i] == BC.value[2]
                elif 'Dirichlet' in BC.type[2]:
                    T_interblock_w, T_interblock_o = T_interblock(i, i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'x', P, Sw, petro, PVT, Res)
                    Jw[i,i] = 2 * T_interblock_w
                    Jo[i,i] = 2 * T_interblock_o
                    Qo[i] = 2 * T_interblock_o * BC.value[2]
                    Qw[i,i] = 2 * T_interblock_w * BC.value[2]
            else:
                None
        
        #Top boundary conditions P(y=w, t)
        if i // Num.Nx == Num.Ny-1:
            if len(BC.type) > 2:
                if 'Neumann' in BC.type[3]:
                    Q[i] == BC.value[3]
                elif 'Dirichlet' in BC.type[3]:
                    T_interblock_w, T_interblock_o = T_interblock(i, i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky, 'x', P, Sw, petro, PVT, Res)
                    Jw[i,i] = 2 * T_interblock_w
                    Jo[i,i] = 2 * T_interblock_o
                    Qo[i] = 2 * T_interblock_o * BC.value[3]
                    Qw[i,i] = 2 * T_interblock_w * BC.value[3]
            else:
                None
        Vp = Num.dx[i] * Num.dy[i] * Num.dz[i] * Res.phi[i]  # Pore Volume
        d11[i,i] = Vp * Sw[i,0] * (PVT.cw + PVT.cf) / (PVT.Bw * Num.dt)
        d12[i,i] = Vp / (PVT.Bw * Num.dt)*(1.0 - Sw[i,0] * Res.phi[i] * PVT.cw * dPcow[i] )
       
    BetaW = np.diag(np.full((Num.N), Beta_w)) 
    BetaO = np.diag(np.full((Num.N), Beta_o))  
    J = Jw + Jo
    G2 = (BetaW * Tw) @ Pcow
    G1 = np.reshape((((PVT.rhoo / 144) * (BetaO * To) + (PVT.rhow / 144) * (BetaW * Tw)) @ Res.D), (Num.N,1))
    G = G1 + G2
    T = Beta_o * To + Beta_w[0] * Tw  
   
    
    #Convert matrices and arrays to sparse   
    G = csr_matrix(G)
    Act = csr_matrix(Act) 
    A = csr_matrix(A)
    T = csr_matrix(T)
    Tw = csr_matrix(Tw) 
    To = csr_matrix(To)
    J = csr_matrix(J)
    Jo = csr_matrix(Jo) 
    Jw = csr_matrix(Jw)
    Q = csr_matrix(Q)
    Qo = csr_matrix(Qo) 
    Qw = csr_matrix(Qw)
    d12 = d12.tocsr() 
    d11 = d11.tocsr() 
    
    return T, Tw, To, A, Act, Q, J, Jw, Jo, G, BetaW, Qw, Qo, G1, Pcow, d11, d12
