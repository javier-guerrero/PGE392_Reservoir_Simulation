"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

def T_interblock(i,j, dz, dx, dy, kx, ky,  direction, krao,mu):
    '''
    This function computes the interblock transmissibility between 2 blocks 

    Parameters
    ----------
    i : block index --> first block or from where the flow comes (e.g. if flow goes from left to right --> this index id for the block on the left)
    j : block index --> second block or where the flow goes
    dx : Delta x --> discretization in x-dir
    dy : Delta y --> discretization in y-dir
    dz : Delta z --> discretization in z-dir
    kx : Permeability in the x-dir
    ky : Permeability in the y-dir
    direction : Direction of the connection x-dir, y-dir, z-dir
    krao : Phase Endpoint relative permeability --> for single phase only 1 phase flows then we can use the endpoint  
    mu : phase viscosity --> in this case oil viscosity

    Returns
    -------
    T_interblock : Transmissibility

    '''
    CF = 6.33e-3  # Conversion factor from mD-ft / cp to ft3/psi-day
    
    if direction=='x':
        T_interblock=CF*(krao/mu)*2*dz[i]*kx[i]*dy[i]*kx[j]*dy[j]/(kx[i]*dy[i]*dx[j]+kx[j]*dy[j]*dx[i]) # Units: ft3/psi-day
   
    elif direction=='y':
        T_interblock=CF*(krao/mu)*2*dz[i]*ky[i]*dx[i]*ky[j]*dx[j]/(ky[i]*dx[i]*dy[j]+ky[j]*dx[j]*dy[i]) # Units: ft3/psi-day
    
    return T_interblock

