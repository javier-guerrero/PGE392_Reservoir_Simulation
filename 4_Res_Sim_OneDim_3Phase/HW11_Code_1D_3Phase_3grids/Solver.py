# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% Import Libraries and packages
import numpy as np
from numpy.linalg import matrix_power
from scipy import special, sparse, integrate
import scipy.sparse.linalg as lg
import sympy as sym
import scipy as sp
import scipy.optimize as opt

#%% Define functions
def Solve(numeric, T, A, J, Q, G, p_old):
    if numeric['1-phase_method'] == "EX":
        theta = 1
    elif numeric['1-phase_method'] == "IM":
        theta = 0
    elif numeric['1-phase_method'] == "CN":
        theta = 0.5
    else:
        print("Error: check the solution method!")
        
    lhs = (1-theta)*(T+J) + A
    rhs = (A-theta*(T+J))*p_old + Q + sparse.lil_matrix(G)
    p = lg.cg(lhs,rhs)[0]
    p = np.reshape(p,(p.shape[0],1))
    return p

def So_eqn(So_old,A,ct,betao,To,Qo,Jo,po,po_old,numeric,co_star,pvt):
    a1 = So_old.reshape((numeric['N'],)) + ct.multiply(sparse.linalg.inv(A))*(-sparse.diags(betao.reshape(numeric['N'],),0)*To*po+Qo.toarray().reshape((numeric['N'],))-Jo*po)
    a2 = -(So_old*(co_star+pvt['cf'])*(po.reshape((numeric['N'],1))-po_old)).reshape((numeric['N'],))
    return a1 + a2
    
def Sw_eqn(Sw_old,A,ct,betaw,Tw,Qw,Jw,po,po_old,numeric,pvt):
    a1 = Sw_old.reshape((numeric['N'],)) + ct.multiply(sparse.linalg.inv(A))*(-sparse.diags(betaw.reshape(numeric['N'],),0)*Tw*po+Qw.toarray().reshape((numeric['N'],))-Jw*po)
    a2 = -(Sw_old*(pvt['cw']+pvt['cf'])*(po.reshape((numeric['N'],1))-po_old)).reshape((numeric['N'],))
    return a1 + a2