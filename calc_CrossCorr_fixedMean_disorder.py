#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 16:28:34 2022

@author: tobias
"""

#Possibly TODO: Write as class for this

import numpy as np
from solve_sp_equations_space_disorder import fermi, convolve_spatial_kernel, one_point_source, phi_q_self_consistent
from attractor_with_input import rho_T0, J_w, iter_rho, rel_center


def stat_upto_fourth_order(A_one_point, w_one_point, J_list, g, T, xi, 
                           rho, q, lbda, 
                           N_sample, xs):
    """Calculates the first four cumulants and moments of a binary variable
    in an environment with temperature T and subject to a deterministic 
    chemical potential mu_det_list and a stochastic Gaussian input with mean
    0 and variance 2 * g**2  * q. This average is computed drawing
    N_sample Gaussian random variables."""

    N_x = len(xs)
    Delta_x = 1./N_x
    
    mu_det_list = convolve_spatial_kernel(rho, J_list, Delta_x) + one_point_source(xi, xs, A_one_point, w_one_point) + lbda

    
    m = np.zeros((4,len(mu_det_list)), dtype = float)
    c = np.zeros((4,len(mu_det_list)), dtype = float)
    # Note that the i'th moment (cumulant) is stored in the i-1'th element of
    # the list m (c).
    
    Gauss_list_collection = np.random.normal(0., 1., (N_sample, len(mu_det_list)))
            
    for ii in range(0, N_sample):
        # Gauss_list = np.random.normal(0., 1., len(phi)) #Alternative: Stored list of random numbers and reused them here
        
        mu_total = mu_det_list + Gauss_list_collection[ii,:] * g * np.sqrt(2.*q)
        
        m[0] += fermi(mu_total, T)
        m[1] += fermi(mu_total, T)**2
        m[2] += fermi(mu_total, T)**3
        m[3] += fermi(mu_total, T)**4
            
    m = m/N_sample
    c[0] = m[0]
    c[1] = m[0] - m[1]
    c[2] = m[0] - 3. * m[1] + 2. * m[2]
    c[3] = m[0] - 7. * m[1] + 12.* m[2] - 6. * m[3] 
    
    # contr_dGdq2 = -m[1] + 3.*m[2] - 2.*m[3]
    
    # for ii in range(0,N_x):
    #     print('Erste Momente:', rho[ii], ' ', c[0,ii], ' ', m[0,ii])
    
    return c, m#, contr_dGdq2


def projector(A_one_point, w_one_point, J_list, g, T,  xi, 
              rho, q, lbda, K_eff, 
              xs, N_sample, ell):
     # Note that the i'th moment (cumulant) is stored in the i-1'th element of
    # the list moments (cumulants).
    
    N_x = len(rho)
    
    cumulants, moments = stat_upto_fourth_order(A_one_point, w_one_point, J_list, g, T, xi, 
                           rho, q, lbda, 
                           N_sample, xs)
    
    c2 = cumulants[1,:]
    c3 = cumulants[2,:]
    c4 = cumulants[3,:]    

    
    if(g != 0 and ell !=0 ):
        small_block = np.zeros((2,2), dtype=float)
        print('small_block.shape=', small_block.shape)
        small_block[0,0] = N_x * (g**2/T**2) + (g**4/T**4) * np.sum(c4)
        # small_block[0,0] = N_x * (g**2/T**2) + (g**4/T**4) * np.sum(c4 - 2.*c3 + c2) #Irrtum!!!! War doch richtig, wie es vorher war.
        small_block[0,1] = ell * (g**2/T**2) * np.sum(c3)
        # small_block[0,1] = ell * (g**2/T**2) * np.sum(c3 - c2) #Irrtum!!!! War doch richtig, wie es vorher war.
        small_block[1,0] = small_block[0,1]
        small_block[1,1] = ell**2 * np.sum(c2)
    
        dim_change_mat = np.zeros((N_x,2), dtype=float)
        dim_change_mat[:,0] = (g**2/T**2) * c3
        # dim_change_mat[:,0] = (g**2/T**2) * (c3 - c2) #Irrtum!!!! War doch richtig, wie es vorher war.
        dim_change_mat[:,1] = ell * c2
        
        large_block_sandwiched = np.dot(dim_change_mat.transpose(),np.dot(K_eff/T, dim_change_mat))
        
        kernel = np.linalg.inv(small_block + large_block_sandwiched)
        
        return np.dot(dim_change_mat, np.dot(kernel, dim_change_mat.transpose()))
    elif( g == 0):
        kernel = 1./(np.sum(c2) + np.dot(c2, np.dot(K_eff/T, c2)))
        print("kernel=", kernel)
        return np.outer(c2, c2) * kernel
    elif( ell == 0 and g==0):
        return np.zeros((N_x, N_x), dtype=float)        


def pairwise_corr(A_one_point, w_one_point, J_list, g, T, xi, 
                  rho, q, lbda, K_eff,
                  N_sample, xs, ell):
  
    N_x = len(rho)
    
    cumulants, moments = stat_upto_fourth_order(A_one_point, w_one_point, J_list, g, T, xi, 
                           rho, q, lbda, 
                           N_sample, xs)
    
    c2 = cumulants[1,:]
    
    print('c1[0]=', cumulants[0,0], ' c1[1]=', cumulants[0,1])
    print('c1[0] - c1[0]**2=', cumulants[0,0] - cumulants[0,0]**2, ' c1[1] - c1[1]**2 =', cumulants[0,1] - cumulants[0,1]**2)    
    print('c2[0]=', c2[0], ' c2[1]=', c2[1])
    
    print('m2[0]=', moments[1,0], ' m2[1]=', moments[1,1])
    
    print('K_eff[0,0]=', K_eff[0,0], ' K_eff[0,1]=', K_eff[0,1] )

    direct_contr = np.dot(np.diag(c2), np.dot(K_eff/T, np.diag(c2)))
    
    projector_mat = projector(A_one_point, w_one_point, J_list, g, T,  xi, 
                              rho, q, lbda, K_eff, 
                              xs, N_sample, ell)
    matrix_sandwich = np.identity(N_x) + np.dot(np.diag(c2), K_eff/T)
    
    indirect_contr = - np.dot(matrix_sandwich, np.dot(projector_mat, matrix_sandwich.transpose()))
    
    return direct_contr + indirect_contr, direct_contr, indirect_contr

def Hess_G(J_dist, A_one_point, w_one_point, J_list, g, T, xi, rho, q, lbda, N_sample, xs, ell):
    
    N_x = len(J_list)
    
    cumulants, moments = stat_upto_fourth_order(A_one_point, w_one_point, J_list, g, T, xi, 
                                                rho, q, lbda, 
                                                N_sample, xs)
    
    c2 = cumulants[1,:]
    c3 = cumulants[2,:]
    c4 = cumulants[3,:]
    
    # contr_d2Gdq2 = 
    
    J_dist_inv = np.linalg.inv(J_dist)
    
    HessG = np.zeros((N_x+2,N_x+2), dtype=float)
    
    HessG[:N_x,:N_x] = - T * J_dist_inv + np.diag(c2)
    
    HessG[:N_x,N_x] = (g**2/T**2) * c3
    HessG[N_x,:N_x] = (g**2/T**2) * c3
    
    HessG[:N_x,N_x+1] = ell * c2
    HessG[N_x+1,:N_x] = ell * c2
    
    # HessG[N_x, N_x] = N_x * (g**2/T**2) + (g**4/T**4) * np.sum(c4)
    HessG[N_x, N_x] = N_x * (g**2/T**2) + (g**4/T**4) * np.sum(c4)
    HessG[N_x,N_x+1] = ell * (g**2/T**2) * np.sum(c3)
    HessG[N_x+1,N_x] = HessG[N_x,N_x+1]
    HessG[N_x+1,N_x+1] = ell**2 * np.sum(c2)
    
    return HessG

def corr_brute_force(J_dist, A_one_point, w_one_point, J_list, g, T, xi, rho, q, lbda, N_sample, xs, ell):
    
    N_x = len(J_list)
    
    HessG = Hess_G(J_dist, A_one_point, w_one_point, J_list, g, T, xi, rho, q, lbda, N_sample, xs, ell)
    
    cumulants, moments = stat_upto_fourth_order(A_one_point, w_one_point, J_list, g, T, xi, 
                                                rho, q, lbda, 
                                                N_sample, xs)
    
    c2 = cumulants[1,:]
    c3 = cumulants[2,:]
    c4 = cumulants[3,:]
    
    Sandwich = np.zeros((N_x + 2, N_x), dtype=float)
    Sandwich[:N_x, :] = np.diag(c2)
    
    Sandwich[N_x,:] = (g**2/T**2) * c3
    Sandwich[N_x+1,:] = ell * c2
    
    Sandwich_T = Sandwich.transpose()
    Hess_G_inv = np.linalg.inv(HessG)
    
    print('Hess_G_inv.shape = ', Hess_G_inv.shape)
    print('Sandwich.shape =', Sandwich.shape )
    print('Sandwich_T.shape =', Sandwich_T.shape )
        
    return np.dot(Sandwich_T, np.dot(Hess_G_inv, Sandwich))

def calc_rho_cross_corr_theo(xi, mean_act, g, A_one_point, w_one_point, A_two_point, w_two_point, T, 
                             N_x, N_sample_Gauss, space_form, tol_phi, tol_phi_abs, tol_q, tol_q_abs, eps_relax_phi, eps_relax_rho):
    
    Delta_x = 1./N_x
    
    xs = np.linspace(0.,1.,N_x, endpoint = False)
    
    rho_0 = rho_T0(xs, xi, mean_act)
    J_dist_list =  J_w(xs, w_two_point, A_two_point)

    
    if(g == 0):
    
        rho_final, lbda_final = iter_rho(rho_0, J_dist_list, Delta_x, T, mean_act, 
                                          A_one_point, w_one_point, xi, xs, -0.5*w_two_point)  
        
        var_vec = (1 - rho_final) * rho_final
        
    else:
        phi_final, q_final, rho_final, lbda_final = phi_q_self_consistent(rho_0, mean_act**2, g, 
                                                                          J_dist_list, Delta_x, T, mean_act, 
                                                                          A_one_point, w_one_point, xi, xs, 
                                                                          -0.5 * w_two_point, N_sample_Gauss,
                                                                          tol_phi, tol_phi_abs, tol_q, tol_q_abs,
                                                                          eps_relax_phi, eps_relax_rho)
            
    
        cumulants, moments = stat_upto_fourth_order(A_one_point, w_one_point, J_dist_list, g, T, xi, 
                                                    rho_final, q_final, lbda_final, 
                                                    N_sample_Gauss, xs)
        
        var_vec = cumulants[1]
    
    var_mat = np.diag(var_vec)
        
        
    if(space_form == 'ferro'):
        M = np.outer(np.ones(N_x),np.ones(N_x))
        J_dist = (A_two_point/N_x) * M
            
    elif(space_form == 'rectangle'):
        J_dist = np.zeros((N_x, N_x), dtype=float)
            
        J_dist_row = np.zeros(N_x, dtype=float)
        
        for jj in range(0,N_x):
            del_x = rel_center(xs[0], xs[jj])
            J_dist_row[jj] = (J_w(del_x, w_two_point, A_two_point))/N_x
            
            
        for ii in range(0,N_x):
            J_dist[ii,:] = J_dist_row
            J_dist_row = np.roll(J_dist_row, 1)
            
            
    J_dist_eff = T * np.dot(np.linalg.inv(np.identity(N_x) - np.dot(J_dist/T, var_mat)),J_dist/T)
        
    
    if(g == 0):
        corr_theo_direct = np.dot(var_mat, np.dot(J_dist_eff/T, var_mat))
        
        vec_corr_theo_indirect = np.dot((np.identity(N_x) + np.dot(var_mat, J_dist_eff/T)), var_vec)
        corr_theo_indirect = -np.outer(vec_corr_theo_indirect, vec_corr_theo_indirect)/np.sum(vec_corr_theo_indirect)
        
        corr_theo = np.diag(var_vec)  + corr_theo_direct + corr_theo_indirect
        corr_theo_only_cross = corr_theo - np.diag(np.diag(corr_theo))
    
    else:    
        corr_cross_theo, corr_theo_direct, corr_theo_indirect = pairwise_corr(A_one_point, w_one_point, 
                                                                        J_dist_list, g, T, xi, 
                                                                        rho_final, q_final, lbda_final, J_dist_eff,
                                                                        N_sample_Gauss, xs, 1.)
        
        corr_theo = np.diag(var_vec) + corr_cross_theo
        corr_theo_only_cross = corr_theo - np.diag(np.diag(corr_theo))
    
    
    return xs, J_dist_list, rho_final, lbda_final, np.diag(var_vec), corr_theo, corr_theo_direct, corr_theo_indirect, corr_theo_only_cross
    
    
    
     