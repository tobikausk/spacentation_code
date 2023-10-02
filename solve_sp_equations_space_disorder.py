#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 14:32:15 2022

@author: tobias
"""

import numpy as np
import copy
import parameters_awi as params
from attractor_with_input import iter_rho, convolve_spatial_kernel, one_point_source, fermi
# import multiprocessing as mp
from multiprocessing.pool import ThreadPool
import time


# def sign_mod(x):
#     """"Returns 1, if x>=0, -1, otherwise"""
#     return (1. - np.abs(np.sign(x))) + np.abs(np.sign(x)) * np.sign(x)

# def H(x):
#     return 0.5*(np.sign(x) + 1)

# def H_mod(x):
#     return 0.5*(sign_mod(x) + 1)

# def fermi(mu, T):
#     if(T > 0):
#         return 1./(1. + np.exp(-mu/T))
#     else:
#         return H(mu)
    
# def rel_center(x,x0):
#     x_rel = (x - x0)
#     #Implement periodic boundary conditions:
#     return (H_mod(0.5 - np.abs(x_rel)) *  x_rel 
#             + (1 - H_mod(0.5 - x_rel)) * (x_rel - 1.)
#             + (1 - H_mod(0.5 + x_rel)) * (x_rel + 1.))
    
# def one_point_source(xi,x, A, w):
#     # return A * np.exp(- (xi - x) * (xi - x) /(2.*w*w))
#     Delta_xix = rel_center(xi, x)
#     return A * np.exp(- Delta_xix * Delta_xix /(2.*w*w))
    
# def convolve_spatial_kernel(rho_list, J_list, Delta_x):
#     N = len(rho_list)
#     convolution = np.zeros_like(rho_list)
#     for ii in range(0,N):
#         convolution[ii] = Delta_x * np.dot(np.roll(J_list[::-1], ii), rho_list)
#     return convolution


def gen_rnd_Gauss_mult_thr(mu, sigma, N1, N2, n_proc):
    
    def rnd_mp(rnd_seed):
        return rnd_seed.normal(mu, sigma, (int(N1/n_proc), N2))
    
    # rngs = [np.random.default_rng(n) for n in range(n_proc)]
    rngs = [np.random.default_rng() for n in range(n_proc)]
    
    with ThreadPool(processes=n_proc) as p:
        result = p.map_async(rnd_mp, rngs).get()
    
    return np.array(result).reshape((int(N1/n_proc)*n_proc, N2))


def phi_q_self_consistent(rho_init, q_init, g, J_list, Delta_x, T, f, A_single, w_single, xi, xs, lbda_init, N_sample_disorder, 
                          tol_phi, tol_phi_abs, tol_q, tol_q_abs, eps_phi, eps_rho):
    phi = convolve_spatial_kernel(rho_init, J_list, Delta_x)
    phi_old = phi
    lbda = lbda_init
    
    q = q_init
    q_old = q_init
    
    # plt.figure(1)
    
    iter_step = 0
    
    one_point_list = one_point_source(xi, xs, A_single, w_single)
    
    # Gauss_list_collection = np.random.normal(0., 1., (N_sample_disorder, len(phi)))

    # while((np.sum(np.abs(phi_old - phi))/np.sum(np.abs(phi)) > 0.001  and np.sum(np.abs(phi))/(g * np.sqrt(q)) > 0.0001) 
    #        or (np.abs(q - q_old)/np.abs(q) > 0.001 and g * np.sqrt(q)/np.sum(np.abs(phi)) > 0.0001) or iter_step==0):
    while((np.sum(np.abs(phi_old - phi))/np.sum(np.abs(phi)) > tol_phi  and np.sum(np.abs(phi))/(g * np.sqrt(q)) > tol_phi_abs) 
               or (np.abs(q - q_old)/np.abs(q) > tol_q and g * np.sqrt(q)/np.sum(np.abs(phi)) > tol_q_abs) or iter_step==0):
        
        #Keep random numbers only fixed while searching for right lbda
        # Gauss_list_collection = np.random.normal(0., 1., (N_sample_disorder, len(phi)))
        Gauss_list_collection = gen_rnd_Gauss_mult_thr(0., 1., N_sample_disorder, len(phi), 4)
        # print('Gauss_list_collection,shape = ', Gauss_list_collection.shape)
        
        phi_old = copy.deepcopy(phi)
        q_old = q
        
        lbda_plus = lbda + 10. * (np.max(J_list) + g + A_single) #Search for lbda in range comparable to the other energy terms 
        lbda_minus = lbda - 10. * (np.max(J_list) + g + A_single)
        
        # print("lbda_minus=", lbda_minus)
        # print("lbda_plus=", lbda_plus)
        
        # rho_plus = np.zeros_like(phi)
        # rho_minus = np.zeros_like(phi)
        
        
        
        # for ii in range(0,N_sample_disorder):
        #     # Gauss_list = np.random.normal(0., 1., len(phi))
        #     Gauss_list = Gauss_list_collection[ii,:]
            
        #     # mu_plus = convolve_spatial_kernel(rho, J_list, Delta_x) + one_point_source(xi, xs, A_single, w_single) + lbda_plus
        #     # mu_minus = convolve_spatial_kernel(rho, J_list, Delta_x) + one_point_source(xi, xs, A_single, w_single) + lbda_minus
            
        #     # mu_without_lbda = phi + one_point_source(xi, xs, A_single, w_single) + Gauss_list_collection[ii,:] * g * np.sqrt(2.*q_old)
        #     # mu_without_lbda = phi + one_point_source(xi, xs, A_single, w_single) + Gauss_list * g * np.sqrt(2.*q_old)
        #     mu_without_lbda = phi + one_point_list + Gauss_list * g * np.sqrt(2.*q_old)
            
        #     mu_plus =  mu_without_lbda + lbda_plus
        #     mu_minus = mu_without_lbda + lbda_minus
            
        #     rho_plus += fermi(mu_plus, T)
        #     rho_minus += fermi(mu_minus, T)
        
        # print('N_sample_disorder=', N_sample_disorder)
        # print('np.outer(np.ones(N_sample_disorder), phi + one_point_list).shape = ', np.outer(np.ones(N_sample_disorder), phi + one_point_list).shape)
        # print('Gauss_list_collection.shape = ', Gauss_list_collection.shape)
        
        mu_without_lbda_list = (np.outer(np.ones(N_sample_disorder), phi + one_point_list)
                                + Gauss_list_collection * g * np.sqrt(2.*q_old))
        
        mu_plus_list = mu_without_lbda_list + lbda_plus
        mu_minus_list = mu_without_lbda_list + lbda_minus
    
        rho_plus_list = fermi(mu_plus_list, T)
        rho_minus_list = fermi(mu_minus_list, T)
        
        rho_plus = np.mean(rho_plus_list, axis = 0)
        rho_minus = np.mean(rho_minus_list, axis = 0)
    
        # rho_plus = rho_plus/N_sample_disorder
        # rho_minus = rho_minus/N_sample_disorder
        
        # print('np.max(np.abs(rho_plus))=',np.max(np.abs(rho_plus)),
        #       'np.max(np.abs(rho_plus - rho_plus_np))=', np.max(np.abs(rho_plus - rho_plus_np)),
        #       'np.max(np.abs(rho_minus))=', np.max(np.abs(rho_minus)),
        #       ' np.max(np.abs(rho_minus - rho_minus_np))=', np.max(np.abs(rho_minus - rho_minus_np)))
        
        
        total_act_plus = Delta_x * np.sum(rho_plus)
        total_act_minus = Delta_x * np.sum(rho_minus)
        # print('Weitere aeussere Iteration')
        # print('lbda=', lbda)
        # print('total_act_plus=', total_act_plus)
        # print('total_act_minus=', total_act_minus)
        
        while(np.abs(total_act_minus - f) > 0.0001 and np.abs(total_act_plus - f) > 0.0001):
            lbda_mean = 0.5 * (lbda_minus + lbda_plus)
            
            
            # rho_mean = np.zeros_like(phi)
            # rho_squared_mean = np.zeros_like(phi)
            #############################################################
            # print('Begin disorder average with dedicated function.')
            ############################################################
            start_loop = time.time()
            
            # global calc_mu
            
            # def calc_mu(i):
            #     Gauss_list = Gauss_list_collection[i,:]
            #     # Gauss_list = np.random.normal(0., 1., len(phi))
                
            #     mu_without_lbda = phi + one_point_list + Gauss_list * g * np.sqrt(2.*q_old)
            #     mu_mean = mu_without_lbda + lbda_mean
                
            #     return fermi(mu_mean, T), fermi(mu_mean, T)**2 
            
            # rho_disord_list = np.zeros((N_sample_disorder, len(one_point_list)), dtype=float)
            # rho_squared_disord_list = np.zeros_like(rho_disord_list)
            
            # for ii in range(0, N_sample_disorder):
            #     rho_current, rho_squared_current = calc_mu(ii)
            #     rho_disord_list[ii,:] = rho_current
            #     rho_squared_disord_list[ii,:] = rho_squared_current
            
            ##################################
            ### Alternative computation: #####
            ##################################
            
            # Gauss_list_collection = np.random.normal(0., 1., (N_sample_disorder, len(phi)))
            Gauss_list_collection = gen_rnd_Gauss_mult_thr(0., 1., N_sample_disorder, len(phi), 4)
            
            mu_without_lbda_list = (np.outer(np.ones(N_sample_disorder), phi + one_point_list)
                                    + Gauss_list_collection * g * np.sqrt(2.*q_old))
            mu_mean_list = mu_without_lbda_list + lbda_mean
            
            
            rho_list_np = fermi(mu_mean_list, T)
            rho_squared_list_np = fermi(mu_mean_list, T)**2
            
            rho_mean = np.mean(rho_list_np, axis = 0)
            rho_squared_mean = np.mean(rho_squared_list_np, axis = 0)
            
            # rho_mean_altern = np.mean(rho_disord_list, axis=0)
            # rho_squared_mean_altern = np.mean(rho_squared_disord_list, axis=0)
            # rho_mean = np.mean(rho_disord_list, axis=0)
            # rho_squared_mean = np.mean(rho_squared_disord_list, axis=0)
            
            # print('np.max(np.abs(rho_mean))=',np.max(np.abs(rho_mean)),
            #       'np.max(np.abs(rho_mean - rho_mean_np))=', np.max(np.abs(rho_mean - rho_mean_np)),
            #       'np.max(np.abs(rho_squared_mean))=', np.max(np.abs(rho_squared_mean)),
            #       ' np.max(np.abs(rho_squared_mean - rho_squared_mean_np))=', np.max(np.abs(rho_squared_mean - rho_squared_mean_np)))
            
            
            end_loop = time.time()
            #########################################################################################
            # print('End disorder average with dedicated function. Duration:', end_loop - start_loop)
            #########################################################################################
            
            # print('Begin with parallelized disorder average.')
            # start_parallel = time.time()
            
            # pool = mp.Pool(8)
            
            # results_parallel = pool.map(calc_mu, [jj for jj in range(0, N_sample_disorder)])
            # # results_parallel = pool.map_async(calc_mu, [jj for jj in range(0, N_sample_disorder)])
            
            # end_parallel = time.time()
            # print('len(results_parallel) = ', len(results_parallel))
            # print('len(results_parallel[12]) =', len(results_parallel[12]))
            # print('End parallelized disorder average. Duration:', end_parallel - start_parallel)
            
            # print('Begin disorder average without dedicated function.')
            
            # for ii in range(0, N_sample_disorder):
            #     # Gauss_list = np.random.normal(0., 1., len(phi)) #Alternative: Store list of random numbers and reuse them here
            #     Gauss_list = Gauss_list_collection[ii,:]
            #     #TODO: Ersetze
                
            #     # mu_without_lbda = phi + one_point_source(xi, xs, A_single, w_single) + Gauss_list_collection[ii,:] * g * np.sqrt(2.*q_old)
            #     # mu_without_lbda = phi + one_point_source(xi, xs, A_single, w_single) + Gauss_list * g * np.sqrt(2.*q_old)
            #     mu_without_lbda = phi + one_point_list + Gauss_list * g * np.sqrt(2.*q_old)
            #     mu_mean = mu_without_lbda + lbda_mean
                
            #     rho_mean += fermi(mu_mean, T)
            #     rho_squared_mean += fermi(mu_mean, T)**2
 
            # rho_mean = rho_mean/N_sample_disorder #rho_mean is already averaged over t (i.e. disorder)
            # rho_squared_mean = rho_squared_mean/N_sample_disorder
            
            # print('End disorder average without dedicated function.')
            
 
            
            # print('rho_mean_altern=', rho_mean_altern[112:122], ' rho_mean=', rho_mean[112:122])
            # print('(rho_mean_altern - rho_mean)/rho_mean=', (rho_mean_altern[112:122] - rho_mean[112:122])/rho_mean[112:122])
            # print('rho_squared_mean_altern=', rho_squared_mean_altern[112:122], ' rho_squared_mean=', rho_squared_mean[112:122])
            # print('(rho_squared_mean_altern - rho_squared_mean)/rho_squared_mean=', (rho_squared_mean_altern[112:122] - rho_squared_mean[112:122])/rho_squared_mean[112:122])
            
            total_act_mean = Delta_x * np.sum(rho_mean)
            q_new = Delta_x * np.sum(rho_squared_mean)
            
            # print('total_act_mean=', total_act_mean)            
            
            
            if(np.sign(total_act_mean - f) == np.sign(total_act_minus - f)):
                lbda_minus = lbda_mean
                total_act_minus = total_act_mean
                closest = -1
            elif(np.sign(total_act_mean - f) == np.sign(total_act_plus - f)):
                lbda_plus = lbda_mean
                total_act_plus = total_act_mean
                closest = 1
            else:
                break;
        
            
        
        if(closest == -1):
            lbda = lbda_minus
        else:
            lbda = lbda_plus
        
        
        # rho =  rho_mean
        
        # eps_rho = 0.05
        # eps_phi = 0.05
        # eps_rho = 0.3
        # eps_phi = 0.3
        
        # print('eps_relax_rho=', eps_rho)
        # print('eps_relax_phi=', eps_phi)
        
        q = eps_rho * q_new + (1. - eps_rho) * q_old
        phi = eps_phi * convolve_spatial_kernel(rho_mean, J_list, Delta_x) + (1. - eps_phi) * phi_old 
        # mu = phi + one_point_source(xi, xs, A_single, w_single) + lbda
        # rho = 0.5 * fermi(mu, T) + 0.5 * rho
        
        # plt.plot(xs, rho_old)
        # rho_plus
        iter_step += 1
        print('np.sum(np.abs(phi_old - phi))/np.sum(np.abs(phi)) =', np.sum(np.abs(phi_old - phi))/np.sum(np.abs(phi)))
        print('np.sum(np.abs(q_old - q))/np.sum(np.abs(q)) =', np.sum(np.abs(q_old - q))/np.sum(np.abs(q)))
        
    print('Funktion iter_rho fast vorbei.')        
    return phi, q, rho_mean, lbda



def phi_q_self_consistent_single_rnd_number_draw(rho_init, q_init, g, J_list, Delta_x, T, f, A_single, w_single, xi, xs, lbda_init, N_sample_disorder):
    phi = convolve_spatial_kernel(rho_init, J_list, Delta_x)
    phi_old = phi
    lbda = lbda_init
    
    q = q_init
    q_old = q_init
    
    # plt.figure(1)
    
    iter_step = 0
    
    one_point_list = one_point_source(xi, xs, A_single, w_single)
    
    Gauss_list_collection = np.random.normal(0., 1., (N_sample_disorder, len(phi)))

    while((np.sum(np.abs(phi_old - phi))/np.sum(np.abs(phi)) > 0.001  and np.sum(np.abs(phi))/(g * np.sqrt(q)) > 0.0001) 
           or (np.abs(q - q_old)/np.abs(q) > 0.001 and g * np.sqrt(q)/np.sum(np.abs(phi)) > 0.0001) or iter_step==0):
        
        #Keep random numbers only fixed while searching for right lbda
        # Gauss_list_collection = np.random.normal(0., 1., (N_sample_disorder, len(phi)))
        
        phi_old = copy.deepcopy(phi)
        q_old = q
        
        lbda_plus = lbda + 10. * (np.max(J_list) + g + A_single) #Search for lbda in range comparable to the other energy terms 
        lbda_minus = lbda - 10. * (np.max(J_list) + g + A_single)
        
        print("lbda_minus=", lbda_minus)
        print("lbda_plus=", lbda_plus)
        
        # rho_plus = np.zeros_like(phi)
        # rho_minus = np.zeros_like(phi)
        
        
        
        # for ii in range(0,N_sample_disorder):
        #     # Gauss_list = np.random.normal(0., 1., len(phi))
        #     Gauss_list = Gauss_list_collection[ii,:]
            
        #     # mu_plus = convolve_spatial_kernel(rho, J_list, Delta_x) + one_point_source(xi, xs, A_single, w_single) + lbda_plus
        #     # mu_minus = convolve_spatial_kernel(rho, J_list, Delta_x) + one_point_source(xi, xs, A_single, w_single) + lbda_minus
            
        #     # mu_without_lbda = phi + one_point_source(xi, xs, A_single, w_single) + Gauss_list_collection[ii,:] * g * np.sqrt(2.*q_old)
        #     # mu_without_lbda = phi + one_point_source(xi, xs, A_single, w_single) + Gauss_list * g * np.sqrt(2.*q_old)
        #     mu_without_lbda = phi + one_point_list + Gauss_list * g * np.sqrt(2.*q_old)
            
        #     mu_plus =  mu_without_lbda + lbda_plus
        #     mu_minus = mu_without_lbda + lbda_minus
            
        #     rho_plus += fermi(mu_plus, T)
        #     rho_minus += fermi(mu_minus, T)
        
        mu_without_lbda_list = (np.outer(np.ones(N_sample_disorder), phi + one_point_list)
                                + Gauss_list_collection * g * np.sqrt(2.*q_old))
        
        mu_plus_list = mu_without_lbda_list + lbda_plus
        mu_minus_list = mu_without_lbda_list + lbda_minus
    
        rho_plus_list = fermi(mu_plus_list, T)
        rho_minus_list = fermi(mu_minus_list, T)
        
        rho_plus = np.mean(rho_plus_list, axis = 0)
        rho_minus = np.mean(rho_minus_list, axis = 0)
    
        # rho_plus = rho_plus/N_sample_disorder
        # rho_minus = rho_minus/N_sample_disorder
        
        # print('np.max(np.abs(rho_plus))=',np.max(np.abs(rho_plus)),
        #       'np.max(np.abs(rho_plus - rho_plus_np))=', np.max(np.abs(rho_plus - rho_plus_np)),
        #       'np.max(np.abs(rho_minus))=', np.max(np.abs(rho_minus)),
        #       ' np.max(np.abs(rho_minus - rho_minus_np))=', np.max(np.abs(rho_minus - rho_minus_np)))
        
        
        total_act_plus = Delta_x * np.sum(rho_plus)
        total_act_minus = Delta_x * np.sum(rho_minus)
        print('Weitere aeussere Iteration')
        print('lbda=', lbda)
        print('total_act_plus=', total_act_plus)
        print('total_act_minus=', total_act_minus)
        
        while(np.abs(total_act_minus - f) > 0.0001 and np.abs(total_act_plus - f) > 0.0001):
            lbda_mean = 0.5 * (lbda_minus + lbda_plus)
            
            
            # rho_mean = np.zeros_like(phi)
            # rho_squared_mean = np.zeros_like(phi)
            #############################################################
            # print('Begin disorder average with dedicated function.')
            ############################################################
            start_loop = time.time()
            
            # global calc_mu
            
            # def calc_mu(i):
            #     Gauss_list = Gauss_list_collection[i,:]
            #     # Gauss_list = np.random.normal(0., 1., len(phi))
                
            #     mu_without_lbda = phi + one_point_list + Gauss_list * g * np.sqrt(2.*q_old)
            #     mu_mean = mu_without_lbda + lbda_mean
                
            #     return fermi(mu_mean, T), fermi(mu_mean, T)**2 
            
            # rho_disord_list = np.zeros((N_sample_disorder, len(one_point_list)), dtype=float)
            # rho_squared_disord_list = np.zeros_like(rho_disord_list)
            
            # for ii in range(0, N_sample_disorder):
            #     rho_current, rho_squared_current = calc_mu(ii)
            #     rho_disord_list[ii,:] = rho_current
            #     rho_squared_disord_list[ii,:] = rho_squared_current
            
            ##################################
            ### Alternative computation: #####
            ##################################
            
            # Gauss_list_collection = np.random.normal(0., 1., (N_sample_disorder, len(phi)))
            
            mu_without_lbda_list = (np.outer(np.ones(N_sample_disorder), phi + one_point_list)
                                    + Gauss_list_collection * g * np.sqrt(2.*q_old))
            mu_mean_list = mu_without_lbda_list + lbda_mean
            
            
            rho_list_np = fermi(mu_mean_list, T)
            rho_squared_list_np = fermi(mu_mean_list, T)**2
            
            rho_mean = np.mean(rho_list_np, axis = 0)
            rho_squared_mean = np.mean(rho_squared_list_np, axis = 0)
            
            # rho_mean_altern = np.mean(rho_disord_list, axis=0)
            # rho_squared_mean_altern = np.mean(rho_squared_disord_list, axis=0)
            # rho_mean = np.mean(rho_disord_list, axis=0)
            # rho_squared_mean = np.mean(rho_squared_disord_list, axis=0)
            
            # print('np.max(np.abs(rho_mean))=',np.max(np.abs(rho_mean)),
            #       'np.max(np.abs(rho_mean - rho_mean_np))=', np.max(np.abs(rho_mean - rho_mean_np)),
            #       'np.max(np.abs(rho_squared_mean))=', np.max(np.abs(rho_squared_mean)),
            #       ' np.max(np.abs(rho_squared_mean - rho_squared_mean_np))=', np.max(np.abs(rho_squared_mean - rho_squared_mean_np)))
            
            
            end_loop = time.time()
            #########################################################################################
            # print('End disorder average with dedicated function. Duration:', end_loop - start_loop)
            #########################################################################################
            
            # print('Begin with parallelized disorder average.')
            # start_parallel = time.time()
            
            # pool = mp.Pool(8)
            
            # results_parallel = pool.map(calc_mu, [jj for jj in range(0, N_sample_disorder)])
            # # results_parallel = pool.map_async(calc_mu, [jj for jj in range(0, N_sample_disorder)])
            
            # end_parallel = time.time()
            # print('len(results_parallel) = ', len(results_parallel))
            # print('len(results_parallel[12]) =', len(results_parallel[12]))
            # print('End parallelized disorder average. Duration:', end_parallel - start_parallel)
            
            # print('Begin disorder average without dedicated function.')
            
            # for ii in range(0, N_sample_disorder):
            #     # Gauss_list = np.random.normal(0., 1., len(phi)) #Alternative: Store list of random numbers and reuse them here
            #     Gauss_list = Gauss_list_collection[ii,:]
            #     #TODO: Ersetze
                
            #     # mu_without_lbda = phi + one_point_source(xi, xs, A_single, w_single) + Gauss_list_collection[ii,:] * g * np.sqrt(2.*q_old)
            #     # mu_without_lbda = phi + one_point_source(xi, xs, A_single, w_single) + Gauss_list * g * np.sqrt(2.*q_old)
            #     mu_without_lbda = phi + one_point_list + Gauss_list * g * np.sqrt(2.*q_old)
            #     mu_mean = mu_without_lbda + lbda_mean
                
            #     rho_mean += fermi(mu_mean, T)
            #     rho_squared_mean += fermi(mu_mean, T)**2
 
            # rho_mean = rho_mean/N_sample_disorder #rho_mean is already averaged over t (i.e. disorder)
            # rho_squared_mean = rho_squared_mean/N_sample_disorder
            
            # print('End disorder average without dedicated function.')
            
 
            
            # print('rho_mean_altern=', rho_mean_altern[112:122], ' rho_mean=', rho_mean[112:122])
            # print('(rho_mean_altern - rho_mean)/rho_mean=', (rho_mean_altern[112:122] - rho_mean[112:122])/rho_mean[112:122])
            # print('rho_squared_mean_altern=', rho_squared_mean_altern[112:122], ' rho_squared_mean=', rho_squared_mean[112:122])
            # print('(rho_squared_mean_altern - rho_squared_mean)/rho_squared_mean=', (rho_squared_mean_altern[112:122] - rho_squared_mean[112:122])/rho_squared_mean[112:122])
            
            total_act_mean = Delta_x * np.sum(rho_mean)
            q_new = Delta_x * np.sum(rho_squared_mean)
            
            # print('total_act_mean=', total_act_mean)            
            
            
            if(np.sign(total_act_mean - f) == np.sign(total_act_minus - f)):
                lbda_minus = lbda_mean
                total_act_minus = total_act_mean
                closest = -1
            elif(np.sign(total_act_mean - f) == np.sign(total_act_plus - f)):
                lbda_plus = lbda_mean
                total_act_plus = total_act_mean
                closest = 1
            else:
                break;
        
            
        
        if(closest == -1):
            lbda = lbda_minus
        else:
            lbda = lbda_plus
        
        
        # rho =  rho_mean
        
        # eps_rho = 0.1
        # eps_phi = 0.1
        eps_rho = 0.3
        eps_phi = 0.3
        
        q = eps_rho * q_new + (1. - eps_rho) * q_old
        phi = eps_phi * convolve_spatial_kernel(rho_mean, J_list, Delta_x) + (1. - eps_phi) * phi_old 
        # mu = phi + one_point_source(xi, xs, A_single, w_single) + lbda
        # rho = 0.5 * fermi(mu, T) + 0.5 * rho
        
        # plt.plot(xs, rho_old)
        # rho_plus
        iter_step += 1
        print('np.sum(np.abs(phi_old - phi))/np.sum(np.abs(phi)) =', np.sum(np.abs(phi_old - phi))/np.sum(np.abs(phi)))
        print('np.sum(np.abs(q_old - q))/np.sum(np.abs(q)) =', np.sum(np.abs(q_old - q))/np.sum(np.abs(q)))
        
    print('Funktion iter_rho fast vorbei.')        
    return phi, q, rho_mean, lbda



def phi_q_self_consistent_old(rho_init, q_init, g, J_list, Delta_x, T, f, A_single, w_single, xi, xs, lbda_init, N_sample_disorder):
    phi = convolve_spatial_kernel(rho_init, J_list, Delta_x)
    phi_old = phi
    lbda = lbda_init
    
    q = q_init
    q_old = q_init
    
    # plt.figure(1)
    
    iter_step = 0
    
    # Gauss_list_collection = np.random.normal(0., 1., (N_sample_disorder, len(phi)))

    while((np.sum(np.abs(phi_old - phi))/np.sum(np.abs(phi)) > 0.001  and np.sum(np.abs(phi))/(g * np.sqrt(q)) > 0.0001) 
           or (np.abs(q - q_old)/np.abs(q) > 0.001 and g * np.sqrt(q)/np.sum(np.abs(phi)) > 0.0001) or iter_step==0):
        
        #Keep random numbers only fixed while searching for right lbda
        Gauss_list_collection = np.random.normal(0., 1., (N_sample_disorder, len(phi)))
        
        phi_old = copy.deepcopy(phi)
        q_old = q
        
        lbda_plus = lbda + 10. * (np.max(J_list) + g + A_single) #Search for lbda in range comparable to the other energy terms 
        lbda_minus = lbda - 10. * (np.max(J_list) + g + A_single)
        
        print("lbda_minus=", lbda_minus)
        print("lbda_plus=", lbda_plus)
        
        rho_plus = np.zeros_like(phi)
        rho_minus = np.zeros_like(phi)
        
        one_point_list = one_point_source(xi, xs, A_single, w_single)
        
        for ii in range(0,N_sample_disorder):
            Gauss_list = np.random.normal(0., 1., len(phi))
            # Gauss_list = Gauss_list_collection[ii,:]
            
            # mu_plus = convolve_spatial_kernel(rho, J_list, Delta_x) + one_point_source(xi, xs, A_single, w_single) + lbda_plus
            # mu_minus = convolve_spatial_kernel(rho, J_list, Delta_x) + one_point_source(xi, xs, A_single, w_single) + lbda_minus
            
            # mu_without_lbda = phi + one_point_source(xi, xs, A_single, w_single) + Gauss_list_collection[ii,:] * g * np.sqrt(2.*q_old)
            # mu_without_lbda = phi + one_point_source(xi, xs, A_single, w_single) + Gauss_list * g * np.sqrt(2.*q_old)
            mu_without_lbda = phi + one_point_list + Gauss_list * g * np.sqrt(2.*q_old)
            
            mu_plus =  mu_without_lbda + lbda_plus
            mu_minus = mu_without_lbda + lbda_minus
            
            rho_plus += fermi(mu_plus, T)
            rho_minus += fermi(mu_minus, T)
        
        rho_plus = rho_plus/N_sample_disorder
        rho_minus = rho_minus/N_sample_disorder
        
        total_act_plus = Delta_x * np.sum(rho_plus)
        total_act_minus = Delta_x * np.sum(rho_minus)
        print('Weitere aeussere Iteration')
        print('lbda=', lbda)
        print('total_act_plus=', total_act_plus)
        print('total_act_minus=', total_act_minus)
        
        while(np.abs(total_act_minus - f) > 0.0001 and np.abs(total_act_plus - f) > 0.0001):
            lbda_mean = 0.5 * (lbda_minus + lbda_plus)
            
            
            rho_mean = np.zeros_like(phi)
            rho_squared_mean = np.zeros_like(phi)
            #############################################################
            # print('Begin disorder average with dedicated function.')
            ############################################################
            start_loop = time.time()
            
            # global calc_mu
            
            def calc_mu(i):
                # Gauss_list = Gauss_list_collection[ii,:]
                Gauss_list = np.random.normal(0., 1., len(phi))
                
                mu_without_lbda = phi + one_point_list + Gauss_list * g * np.sqrt(2.*q_old)
                mu_mean = mu_without_lbda + lbda_mean
                
                return fermi(mu_mean, T), fermi(mu_mean, T)**2 
            
            rho_disord_list = np.zeros((N_sample_disorder, len(one_point_list)), dtype=float)
            rho_squared_disord_list = np.zeros_like(rho_disord_list)
            
            for ii in range(0, N_sample_disorder):
                rho_current, rho_squared_current = calc_mu(ii)
                rho_disord_list[ii,:] = rho_current
                rho_squared_disord_list[ii,:] = rho_squared_current
                
            # rho_mean_altern = np.mean(rho_disord_list, axis=0)
            # rho_squared_mean_altern = np.mean(rho_squared_disord_list, axis=0)
            rho_mean = np.mean(rho_disord_list, axis=0)
            rho_squared_mean = np.mean(rho_squared_disord_list, axis=0)
            
            end_loop = time.time()
            #########################################################################################
            # print('End disorder average with dedicated function. Duration:', end_loop - start_loop)
            #########################################################################################
            
            # print('Begin with parallelized disorder average.')
            # start_parallel = time.time()
            
            # pool = mp.Pool(8)
            
            # results_parallel = pool.map(calc_mu, [jj for jj in range(0, N_sample_disorder)])
            # # results_parallel = pool.map_async(calc_mu, [jj for jj in range(0, N_sample_disorder)])
            
            # end_parallel = time.time()
            # print('len(results_parallel) = ', len(results_parallel))
            # print('len(results_parallel[12]) =', len(results_parallel[12]))
            # print('End parallelized disorder average. Duration:', end_parallel - start_parallel)
            
            # print('Begin disorder average without dedicated function.')
            
            # for ii in range(0, N_sample_disorder):
            #     # Gauss_list = np.random.normal(0., 1., len(phi)) #Alternative: Store list of random numbers and reuse them here
            #     Gauss_list = Gauss_list_collection[ii,:]
            #     #TODO: Ersetze
                
            #     # mu_without_lbda = phi + one_point_source(xi, xs, A_single, w_single) + Gauss_list_collection[ii,:] * g * np.sqrt(2.*q_old)
            #     # mu_without_lbda = phi + one_point_source(xi, xs, A_single, w_single) + Gauss_list * g * np.sqrt(2.*q_old)
            #     mu_without_lbda = phi + one_point_list + Gauss_list * g * np.sqrt(2.*q_old)
            #     mu_mean = mu_without_lbda + lbda_mean
                
            #     rho_mean += fermi(mu_mean, T)
            #     rho_squared_mean += fermi(mu_mean, T)**2
 
            # rho_mean = rho_mean/N_sample_disorder #rho_mean is already averaged over t (i.e. disorder)
            # rho_squared_mean = rho_squared_mean/N_sample_disorder
            
            # print('End disorder average without dedicated function.')
            
 
            
            # print('rho_mean_altern=', rho_mean_altern[112:122], ' rho_mean=', rho_mean[112:122])
            # print('(rho_mean_altern - rho_mean)/rho_mean=', (rho_mean_altern[112:122] - rho_mean[112:122])/rho_mean[112:122])
            # print('rho_squared_mean_altern=', rho_squared_mean_altern[112:122], ' rho_squared_mean=', rho_squared_mean[112:122])
            # print('(rho_squared_mean_altern - rho_squared_mean)/rho_squared_mean=', (rho_squared_mean_altern[112:122] - rho_squared_mean[112:122])/rho_squared_mean[112:122])
            
            total_act_mean = Delta_x * np.sum(rho_mean)
            q_new = Delta_x * np.sum(rho_squared_mean)
            
            # print('total_act_mean=', total_act_mean)            
            
            
            if(np.sign(total_act_mean - f) == np.sign(total_act_minus - f)):
                lbda_minus = lbda_mean
                total_act_minus = total_act_mean
                closest = -1
            elif(np.sign(total_act_mean - f) == np.sign(total_act_plus - f)):
                lbda_plus = lbda_mean
                total_act_plus = total_act_mean
                closest = 1
            else:
                break;
        
            
        
        if(closest == -1):
            lbda = lbda_minus
        else:
            lbda = lbda_plus
        
        
        # rho =  rho_mean
        
        # eps_rho = 0.1
        # eps_phi = 0.1
        eps_rho = 0.3
        eps_phi = 0.3
        
        q = eps_rho * q_new + (1. - eps_rho) * q_old
        phi = eps_phi * convolve_spatial_kernel(rho_mean, J_list, Delta_x) + (1. - eps_phi) * phi_old 
        # mu = phi + one_point_source(xi, xs, A_single, w_single) + lbda
        # rho = 0.5 * fermi(mu, T) + 0.5 * rho
        
        # plt.plot(xs, rho_old)
        # rho_plus
        iter_step += 1
        print('np.sum(np.abs(phi_old - phi))/np.sum(np.abs(phi)) =', np.sum(np.abs(phi_old - phi))/np.sum(np.abs(phi)))
        print('np.sum(np.abs(q_old - q))/np.sum(np.abs(q)) =', np.sum(np.abs(q_old - q))/np.sum(np.abs(q)))
        
    print('Funktion iter_rho fast vorbei.')        
    return phi, q, rho_mean, lbda

