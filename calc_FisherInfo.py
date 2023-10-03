#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 14:14:58 2022

@author: tobias
"""

import numpy as np
from attractor_with_input import deriv_one_point_source

def FisherInfo_from_corr(corr_mat, x_list, xi, A_one_point, w_one_point, T):
    """ Computes the Fisher information given the correlation matrix and the 
    parameters for the one-point source term."""
    
    U_prime = deriv_one_point_source(xi, x_list, A_one_point, w_one_point)
    return np.dot(U_prime/T, np.dot(corr_mat, U_prime/T))

def FisherInfo_from_corr_compare_noise_corr(corr_mat, var_mat, x_list, xi, A_one_point, w_one_point, T):
    """ Computes the Fisher information given the correlation matrix and the 
    parameters for the one-point source term."""
    
    U_prime = deriv_one_point_source(xi, x_list, A_one_point, w_one_point)
    f_prime = np.dot(corr_mat, U_prime/T)


    # corr_mat_inv = np.linalg.inv(corr_mat)

    EW, EV = np.linalg.eig(corr_mat)
    EW_inv = np.zeros_like(EW)
    for ii in range(len(EW)):
        if(EW[ii] > 10.**(-10)):
            EW_inv[ii] = 1./EW[ii]
        else:
            EW_inv[ii] = 0.

    corr_mat_pseudo_inv = np.dot(EV, np.dot(np.diag(EW_inv), EV.T))
    var_mat_inv = np.linalg.inv(var_mat)

    return np.dot(f_prime, np.dot(corr_mat_pseudo_inv, f_prime)), np.dot(f_prime, np.dot(var_mat_inv, f_prime)) 
