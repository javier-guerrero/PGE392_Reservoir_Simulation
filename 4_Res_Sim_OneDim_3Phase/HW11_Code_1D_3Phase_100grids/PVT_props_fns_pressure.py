# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""



#%% Define functions
def P_dependent_properties(p,res,pvt):
    '''
    Compute Pressure dependent properties using 
    fitted equations 

    '''
    # Constants
    R =   10.73               #psi-ft^3/lbmol-R
    Tsc = 60+459.67           #R
    Psc = 14.7                #psi
    
    pvt['Rs'] = (3.43e-5 * p ** 2 + 9.21e-2 * p) / 5.615 # scf/ft^3
    Bo = 6.119e-5 * p + 1.083 # rb/STB
    z = 1.336e-8 * p ** 2 - 8.804e-5 * p +0.9952 # unitless
    pvt['muo'] = 4.846e-7 * p ** 2 - 1.533e-3 * p + 2.375
    pvt['mug'] = 5.929e-10 * p ** 2 - 6.764e-7 * p + 0.01309
    T = res['T'] + 460
    Bg = 0.0282 * z * T / p
    co_star = -6.119e-5 / Bo
    pvt['co'] = co_star + 0.0282 * z * T / 5.616 / p / Bo * (2 * 3.43e-5 * p + 9.21e-2)
    pvt['cg'] = 1/p - (2 * 1.336e-8 * p - 8.804e-5) / z
    
    return Bo, Bg, co_star, z