"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu

# JOG: Pending to add front and back boundary for 3d reservoirs
# JOG: Pending change matrix to sparse
"""


#%% import packages and libraries
import numpy as np
from T_interblock import T_interblock    #for calculating transmissibility
# from Preprocess import Num, PVT, Res, PP, BC
  
#%% Define Functions 
# fluid, reservoir and simulation parameters   
def get_arrays(Num, PVT, Res, PP, BC):
    '''
    Function:       Calculate arrays for calculation
    input:          Num, PVT, Res,PP, BC
    output:         T, A, Act, Q, J, G
    """

    Parameters
    ----------
    Num : Numerical inputs
    PVT : PVT Inputs
    Res : Reservoir Inputs
    PP : Petrophysical Inputs
    BC : Boundary conditions 

    Returns
    -------
    T : Transmissibility matrix.
    A : Accumulation matrix
    Act : A * ct --> matrix
    Q : Boundary condition vector - Neumann.
    J : Boundary condition matrix - Dirichlet.
    G : Gravity vector

    '''
    #Setting up matrix T, B, and Q
    ct = PVT.ct*np.ones(Num.N)
    T = np.zeros((Num.N, Num.N))
    J = np.zeros((Num.N, Num.N))
    A = np.zeros((Num.N, Num.N))
    Act= np.zeros((Num.N, Num.N))
    Q = np.zeros(Num.N)
    G= np.zeros(Num.N)
    
    for i in range(Num.N):
        #Interblock Trasmisibility in x direction
        if (i+1)%Num.Nx!=0: #Not in the right boundary
            T[i,i+1]=-T_interblock(i,i+1, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky,  'x', PP.kroo,PVT.muo)
            T[i+1,i]=T[i,i+1]
        #Interblock Trasmisibility in y direction
        if i//Num.Nx>0: 
            T[i,i-Num.Nx]=-T_interblock(i,i-Num.Nx, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky,  'y', PP.kroo,PVT.muo)
        if i//Num.Nx<Num.Ny-1:
            T[i,i+Num.Nx]=-T_interblock(i,i+Num.Nx, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky,  'y', PP.kroo,PVT.muo)
        T[i,i]=abs(np.sum(T[i,:]))
        
        #Left boundary conditions P(x=0, t)
        if (i+1)%Num.Nx==1:
            if 'Neumann' in BC.type[0]:
                Q[i]==BC.value[0]
            elif 'Dirichlet' in BC.type[0]:
                J[i,i]=2*T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky,  'x', PP.kroo,PVT.muo)
                Q[i]=2*T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky,  'x', PP.kroo,PVT.muo)*BC.value[0]
        #Right boundary conditions P(x=L, t)
        if (i+1)%Num.Nx==0:
            if 'Neumann' in BC.type[1]:
                Q[i]==BC.value[1]
            elif 'Dirichlet' in BC.type[1]:
                J[i,i]=2*T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky,  'x', PP.kroo,PVT.muo)
                Q[i]=2*T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky,  'x', PP.kroo,PVT.muo)*BC.value[1]
        #Bottom boundary conditions P(y=0, t)
        if i//Num.Nx==0:
            if len(BC.type)>2:
                if 'Neumann' in BC.type[2]:
                    Q[i]==BC.value[2]
                elif 'Dirichlet' in BC.type[2]:
                    J[i,i]=2*T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky,  'y', PP.kroo,PVT.muo)
                    Q[i]=2*T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky,  'y', PP.kroo,PVT.muo)*BC.value[2]
            else:
                None
        #Top boundary conditions P(y=w, t)
        if i//Num.Nx==Num.Ny-1:
            if len(BC.type)>2:
                if 'Neumann' in BC.type[3]:
                    Q[i]==BC.value[3]
                elif 'Dirichlet' in BC.type[3]:
                    J[i,i]=2*T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky,  'y', PP.kroo,PVT.muo)
                    Q[i]=2*T_interblock(i,i, Num.dz, Num.dx, Num.dy, Res.kx, Res.ky,  'y', PP.kroo,PVT.muo)*BC.value[3]
            else:
                None
        G=(PVT.rhoosc/144)*(T@Res.D)
           
        A[i,i]=Num.dx[i]*Num.dy[i]*Num.dz[i]*Res.phi[i]/Num.dt
        Act[i,i]=A[i,i]*ct[i]
    
    return T, A, Act, Q, J, G;
