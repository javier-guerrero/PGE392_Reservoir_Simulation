# -*- coding: utf-8 -*-
"""
PGE-392K-  Numerical Simulation of Reservoirs 
@author: Javier Guerrero - JOG496
email: javier.guerrero@utexas.edu
"""

#%% Import Libraries and packages

import numpy as np
from scipy import special, sparse, integrate
from T_interblock import T_Interblock

#%% Define functions
def Grid_Arrays(A, ct, Sw, numeric, res, petro, pvt, bcs, krw, kro, krg, potential, betaw, betao, betag, Bw, Bo, Bg, po, pc, z):
    Tw, To, Tg = sparse.csr_matrix((numeric['N'],numeric['N'])), sparse.csr_matrix((numeric['N'],numeric['N'])), sparse.csr_matrix((numeric['N'],numeric['N']))
    T = sparse.csr_matrix((numeric['N'],numeric['N']))
    A = np.reshape(A,(numeric['N'],))
    ct = np.reshape(ct,(numeric['N'],))
    A = sparse.diags(A,0)
    ct = sparse.diags(ct,0)
    A = A.multiply(ct)
    Jw, Jo, Jg = sparse.csr_matrix((numeric['N'],numeric['N'])), sparse.csr_matrix((numeric['N'],numeric['N'])), sparse.csr_matrix((numeric['N'],numeric['N']))
    J = sparse.csr_matrix((numeric['N'],numeric['N']))
    Qw, Qo, Qg = sparse.csr_matrix((numeric['N'],1)), sparse.csr_matrix((numeric['N'],1)), sparse.csr_matrix((numeric['N'],1))
    Q = sparse.csr_matrix((numeric['N'],1))
    G = sparse.csr_matrix((numeric['N'],1))
    
    lb = (np.remainder(numeric['block'],numeric['Nx']) == 0)
    rb = (np.remainder(numeric['block'],numeric['Nx']) == numeric['Nx']-1)
    bb = (np.isin(numeric['block'],np.linspace(0,numeric['Nx']-1,numeric['Nx'])))
    tb = (np.isin(numeric['block'],np.linspace(numeric['N']-numeric['Nx'],numeric['N']-1,numeric['Nx'])))
    ins = np.logical_not(np.logical_or(rb,np.logical_or(lb,np.logical_or(tb,bb))))
    
    true_lb = np.asarray(np.nonzero(np.array(lb==True,dtype=int)))
    true_rb = np.asarray(np.nonzero(np.array(rb==True,dtype=int)))
    true_bb = np.asarray(np.nonzero(np.array(bb==True,dtype=int)))
    true_tb = np.asarray(np.nonzero(np.array(tb==True,dtype=int)))
    
    not_lb = np.asarray(np.nonzero(np.array(lb==False,dtype=int)))
    not_rb = np.asarray(np.nonzero(np.array(rb==False,dtype=int)))
    not_bb = np.asarray(np.nonzero(np.array(bb==False,dtype=int)))
    not_tb = np.asarray(np.nonzero(np.array(tb==False,dtype=int)))
    ins = np.asarray(np.nonzero(np.array(ins==True,dtype=int)))
    
    Tw[not_lb,not_lb-1], To[not_lb,not_lb-1], Tg[not_lb,not_lb-1] = tuple([-1*x for x in T_Interblock(not_lb,not_lb-1,'x',res,petro,pvt,numeric,krw,kro,krg,potential,Bw,Bo,Bg)])
    Tw[not_rb,not_rb+1], To[not_rb,not_rb+1], Tg[not_rb,not_rb+1] = tuple([-1*x for x in T_Interblock(not_rb,not_rb+1,'x',res,petro,pvt,numeric,krw,kro,krg,potential,Bw,Bo,Bg)])
    Tw[not_tb,not_tb+numeric['Nx']], To[not_tb,not_tb+numeric['Nx']], Tg[not_tb,not_tb+numeric['Nx']] = tuple([-1*x for x in T_Interblock(not_tb,not_tb+numeric['Nx'],'y',res,petro,pvt,numeric,krw,kro,krg,potential,Bw,Bo,Bg)])
    Tw[not_bb,not_bb-numeric['Nx']], To[not_bb,not_bb-numeric['Nx']], Tg[not_bb,not_bb-numeric['Nx']] = tuple([-1*x for x in T_Interblock(not_bb,not_bb-numeric['Nx'],'y',res,petro,pvt,numeric,krw,kro,krg,potential,Bw,Bo,Bg)])
    diagw = np.asarray(-Tw.sum(axis=1)) # sum the rows
    diagw = sparse.diags(diagw[:,0],0) # make sparse diagonal matrix
    Tw += diagw # append on main diagonal
    
    diago = np.asarray(-To.sum(axis=1)) # sum the rows
    diago = sparse.diags(diago[:,0],0) # make sparse diagonal matrix
    To += diago # append on main diagonal
    
    diagg = np.asarray(-Tg.sum(axis=1)) # sum the rows
    diagg = sparse.diags(diagg[:,0],0) # make sparse diagonal matrix
    Tg += diagg # append on main diagonal
    
    betaw, betao, betag = betaw.reshape((numeric['N'],1)), betao.reshape((numeric['N'],1)), betag.reshape((numeric['N'],1))
    Bw, Bo, Bg = Bw.reshape((numeric['N'],1)), Bo.reshape((numeric['N'],1)), Bg.reshape((numeric['N'],1))
    
    if bcs['BC_xi'] == 'Dirichlet':
        Txw, Txo, Txg = T_Interblock(true_lb,true_lb,'x',res,petro,pvt,numeric,krw,kro,krg,potential,Bw,Bo,Bg)
        Jw[true_lb,true_lb], Jo[true_lb,true_lb], Jg[true_lb,true_lb] = 2*Txw, 2*Txo, 2*Txg
        Qw[true_lb], Qo[true_lb], Qg[true_lb] = Jw[true_lb,true_lb]*bcs['BC_xi_val'], Jo[true_lb,true_lb]*bcs['BC_xi_val'], Jg[true_lb,true_lb]*bcs['BC_xi_val']
        Q[true_lb] = Qw[true_lb] + Qo[true_lb] + Qg[true_lb]
    else:
        Qw[true_lb[0,:]], Qo[true_lb[0,:]], Qg[true_lb[0,:]] = bcs['BC_xi_val'], bcs['BC_xi_val'], bcs['BC_xi_val']
        Q[true_lb[0,:]] = Qw[true_lb[0,:]].multiply(betaw[true_lb[0,:]]) + Qo[true_lb[0,:]].multiply(betao[true_lb[0,:]]) + Qg[true_lb[0,:]].multiply(betag[true_lb[0,:]])
    if bcs['BC_xe'] == 'Dirichlet':
        Txw, Txo, Txg = T_Interblock(true_rb,true_rb,'x',res,petro,pvt,numeric,krw,kro,krg,potential,Bw,Bo,Bg)
        Jw[true_rb,true_rb], Jo[true_rb,true_rb], Jg[true_rb,true_rb] = 2*Txw, 2*Two, 2*Txg
        Qw[true_rb], Qo[true_rb], Qg[true_rb] = Jw[true_rb,true_rb]*bcs['BC_xe_val'], Jo[true_rb,true_rb]*bcs['BC_xe_val'], Jg[true_rb,true_rb]*bcs['BC_xe_val']
        Q[true_rb] = Qw[true_rb] + Qo[true_rb] + Qg[true_rb]
    else:
        Qw[true_rb], Qo[true_rb], Qg[true_rb] = bcs['BC_xe_val'], bcs['BC_xe_val'], bcs['BC_xe_val']
        Q[true_rb[0,:]] = Qw[true_rb[0,:]].multiply(betaw[true_rb[0,:]]) + Qo[true_rb[0,:]].multiply(betao[true_rb[0,:]]) + Qg[true_rb[0,:]].multiply(betag[true_rb[0,:]])
    if (numeric['Ny']==1 and numeric['Nz']==1):
        pass
    else:
        if bcs['BC_yi'] == 'Dirichlet':
            Tyw, Tyo, Tyg = T_Interblock(true_bb,true_bb,'y',res,petro,pvt,numeric,krw,kro,krg,potential,Bw,Bo,Bg)
            Jw[true_bb,true_bb], Jo[true_bb,true_bb], Jg[true_bb,true_bb] = 2*Tyw, 2*Tyo, 2*Tyg
            Qw[true_bb], Qo[true_bb], Qg[true_bb] = Jw[true_bb,true_bb]*bcs['BC_yi_val'], Jo[true_bb,true_bb]*bcs['BC_yi_val'], Jg[true_bb,true_bb]*bcs['BC_yi_val']
            Q[true_bb] = Qw[true_bb] + Qo[true_bb] + Qg[true_bb]
        else:
            Qw[true_bb], Qo[true_bb], Qg[true_bb] = bcs['BC_yi_val'], bcs['BC_yi_val'], bcs['BC_yi_val']
            Q[true_bb[0,:]] = Qw[true_bb[0,:]].multiply(betaw[true_bb[0,:]]) + Qo[true_bb[0,:]].multiply(betao[true_bb[0,:]]) + Qg[true_bb[0,:]].multiply(betag[true_bb[0,:]])
        if bcs['BC_ye'] == 'Dirichlet':
            Tyw, Tyo, Tyg = T_Interblock(true_tb,true_tb,'y',res,petro,pvt,numeric,krw,kro,krg,potential,Bw,Bo,Bg)
            Jw[true_tb,true_tb], Jo[true_tb,true_tb], Jg[true_tb,true_tb] = 2*Tyw, 2*Tyo, 2*Tyg
            Qw[true_tb], Qo[true_tb], Qg[true_tb] = Jw[true_tb,true_tb]*bcs['BC_ye_val'], Jo[true_tb,true_tb]*bcs['BC_ye_val'], Jg[true_tb,true_tb]*bcs['BC_ye_val']
            Q[true_tb] = Qw[true_tb] + Qo[true_tb] + Qg[true_tb]
        else:
            Qw[true_tb], Qo[true_tb], Qg[true_tb] = bcs['BC_ye_val'], bcs['BC_ye_val'], bcs['BC_ye_val']
            Q[true_tb[0,:]] = Qw[true_tb[0,:]].multiply(betaw[true_tb[0,:]]) + Qo[true_tb[0,:]].multiply(betao[true_tb[0,:]]) + Qg[true_tb[0,:]].multiply(betag[true_tb[0,:]])
    
    T = Tw.multiply(betaw) + To.multiply(betao) + Tg.multiply(betag)
    R = 10.7316 # gas constant
    rho_g = pvt['Mg']*po/z/R/(res['T']+460) # gas density [lbm/ft^3]
    rhogsc = pvt['Mg']/379.4 # gas density at standard cond. [lbm/ft3]
    rho_w, rho_o = pvt['rhowsc']/Bw, (pvt['rhoosc']+pvt['Rs']*rhogsc/5.615)/Bo # water and oil density [lbm/ft3]
    gravw = rho_w*betaw*res['D'].reshape((numeric['N']),1)/144
    gravo = rho_o*betao*res['D'].reshape((numeric['N']),1)/144
    gravg = rho_g*betag*res['D'].reshape((numeric['N']),1)/144
    betaw_diag = sparse.diags(betaw.reshape(numeric['N'],),0)
    betag_diag = sparse.diags(betag.reshape(numeric['N'],),0)
    
    G = Tw*gravw + To*gravo+ Tg*gravg + betaw_diag*Tw*pc+ betag_diag*Tg*pc
    J = Jw + Jo + Jg
    
    return T, A, J, Q, G, Tw, To, Tg, Jw, Jo, Jg, Qw, Qo, Qg, ct