#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 15:15:55 2022

@author: tobias
"""

import numpy as np
from matplotlib import pyplot as plt
import copy
# from scipy.optimize import bisect
import parameters_awi as params


def sign_mod(x):
    """"Returns 1, if x>=10**-12, -1, otherwise"""
    x += 10**-12
    return (1. - np.abs(np.sign(x))) + np.abs(np.sign(x)) * np.sign(x) 

def H_mod(x):
    return 0.5*(sign_mod(x) + 1)

def H(x):
    return 0.5*(np.sign(x) + 1)


def one_point_source(xi,x, A, w):
    Delta_xix = rel_center(xi, x)
    return A * np.exp(- Delta_xix * Delta_xix /(2.*w*w))

def deriv_one_point_source(xi, x, A, w):
    Delta_xix = rel_center(xi, x)    
    return - (A * Delta_xix/(w*w)) *  np.exp(- Delta_xix * Delta_xix /(2.*w*w))

def rel_center(x,x0):
    x_rel = (x - x0)
    #Implement periodic boundary conditions:
    return (H_mod(0.5 - np.abs(x_rel)) *  x_rel 
            + (1 - H_mod(0.5 - x_rel)) * (x_rel - 1.)
            + (1 - H_mod(0.5 + x_rel)) * (x_rel + 1.))

def mu_T0(x, x0, w, f):
    x_rel_center_abs = np.abs(rel_center(x,x0))
    
    center = 0.5*w * H_mod(0.5*(f-w) - x_rel_center_abs)
    intermediate = (0.5*f - x_rel_center_abs) * (1. - H_mod(0.5*(f-w) - x_rel_center_abs)) * H(0.5*(f+w) - x_rel_center_abs)
    edge = -0.5*w * (1. - H(0.5*(f+w) - x_rel_center_abs ))   
    
    return  center + intermediate + edge

def rho_T0(x, x0, f):
    x_rel_center = rel_center(x,x0)
    print(0.5 * f - np.abs(x_rel_center))
    return 1. - H_mod(-0.5 * f + np.abs(x_rel_center))



def J_w(x, w, J0):
    return J0 * (H_mod(0.5 * w - np.abs(x)) + H_mod(0.5 * w - np.abs(1. - x))) 

def convolve_spatial_kernel(rho_list, J_list, Delta_x):
    N = len(rho_list)
    convolution = np.zeros_like(rho_list)
    for ii in range(0,N):
        #Explanation +1 next line: For ii=0, we want to map the index -1 to the index 1 (and not 0).
        convolution[ii] = Delta_x * np.dot(np.roll(J_list[::-1], ii+1), rho_list)
    return convolution

def mu_from_rho(x0, xi, A_single, w_single, rho, J_list, Delta_x, lbda):
    return convolve_spatial_kernel(rho, J_list, Delta_x) + one_point_source(xi,x0, A_single, w_single) + lbda


def del_rho_del_xi(rho, w_J, x0, xi, A_single, w_single, lbda, J_list, T, Delta_x, Jw_0):
    w_half_index = round(0.5*w_J/Delta_x)
    rho_shift_minus_half_w = np.roll(rho, -w_half_index)
    rho_shift_plus_half_w = np.roll(rho, w_half_index)
    return deriv_fermi(mu_from_rho(x0, xi, A_single, w_single, rho, J_list, Delta_x, lbda), T) * (
           - Jw_0 * (rho_shift_minus_half_w - rho_shift_plus_half_w) 
           + deriv_one_point_source(xi,x0, A_single, w_single))

def del_mu_del_xi(rho, w_J, x0, xi, A_single, w_single, Delta_x, Jw_0):
    w_half_index = round(0.5*w_J/Delta_x)
    rho_shift_minus_half_w = np.roll(rho, -w_half_index)
    rho_shift_plus_half_w = np.roll(rho, w_half_index)
    indirect = - Jw_0 * (rho_shift_minus_half_w - rho_shift_plus_half_w)
    direct = deriv_one_point_source(xi,x0, A_single, w_single)
    return indirect + direct, indirect, direct 


def Fisher_Info(rho, w_J, x0, xi, A_single, w_single, lbda, J_list, T, Delta_x, Jw_0):
    w_half_index = round(0.5*w_J/Delta_x)
    rho_shift_minus_half_w = np.roll(rho, -w_half_index)
    rho_shift_plus_half_w = np.roll(rho, w_half_index)
    
    prefactor = deriv_fermi(mu_from_rho(x0, xi, A_single, w_single, rho, J_list, Delta_x, lbda), T)/T
    indirect_shift = - Jw_0 * (rho_shift_minus_half_w - rho_shift_plus_half_w)
    direct_shift = deriv_one_point_source(xi,x0, A_single, w_single)
    
    return prefactor * (indirect_shift + direct_shift) **2, prefactor, indirect_shift, direct_shift


def Fisher_Info_numerical(rho, rho_half_step, Delta_x):
    # diff_rho = - np.diff(rho_final)/Delta_x
    diff_rho = - (np.roll(rho_half_step,1) - rho_half_step)/Delta_x
    
    # if(np.abs(diff_rho[0] - 0.5 * (diff_rho[1]+diff_rho[-1]))/np.abs(0.5 * (diff_rho[1]+diff_rho[-1])) > 0.5
    #    and np.abs(diff_rho[1] - diff_rho[-1])/np.abs(0.5 * (diff_rho[1]+diff_rho[-1])) < 0.05 ):
    diff_rho[0] = 0.5 * (diff_rho[1]+diff_rho[-1])
    print(0)
    # if(np.abs(diff_rho[-1] - 0.5 * (diff_rho[-2]+diff_rho[0]))/np.abs(0.5 * (diff_rho[-2]+diff_rho[0])) > 0.5 
    #    and np.abs(diff_rho[-2] - diff_rho[0])/np.abs(0.5 * (diff_rho[-2]+diff_rho[0])) < 0.05):
    #     diff_rho[-1] = 0.5 * (diff_rho[-2]+diff_rho[0])
    #     print(-1)
    
    for ii in range(1,len(diff_rho)-1):
        if(np.abs(diff_rho[ii] - 0.5 * (diff_rho[ii+1]+diff_rho[ii-1]))/np.abs(0.5 * (diff_rho[ii+1]+diff_rho[ii-1])) > 0.5
           and np.abs(diff_rho[ii+1] - diff_rho[ii-1])/np.abs(0.5 * (diff_rho[ii+1]+diff_rho[ii-1])) < 0.05 ):
            diff_rho[ii] = 0.5 * (diff_rho[ii+1]+diff_rho[ii-1])
            print(ii)

    return diff_rho**2/(rho * (1. - rho))


def Fisher_Info_integrated(rho, w_J, x0, xi, A_single, w_single, lbda, J_list, T, Delta_x, Jw_0):
    I_list, prefac_list, indirect_list, direct_list = Fisher_Info(rho, w_J, x0, xi, A_single, w_single, lbda, J_list, T, Delta_x, Jw_0)
    
    return Delta_x * np.sum(I_list)

def Fisher_Info_numerical_integrated(rho, rho_half_step, Delta_x):
    I_list = Fisher_Info_numerical(rho, rho_half_step, Delta_x)
    return Delta_x * np.sum(I_list)

def fermi(mu, T):
    if(T > 0):
        return 1./(1. + np.exp(-mu/T))
    else:
        return H(mu)
    
def deriv_fermi(mu, T):
    return 1./(4.*T*np.cosh(mu/(2.*T))**2)
    

def iter_rho(rho_init, J_list, Delta_x, T, f, A_single, w_single, xi, xs, lbda_init):
    rho = rho_init
    rho_old = rho_init
    lbda = lbda_init
    
    iter_step = 0
    while(np.sum((rho_old - rho)**2) > 0.0000001 or iter_step==0):
        rho_old = copy.deepcopy(rho)
        
        lbda_plus = lbda + 10. * (np.max(J_list) + A_single)
        lbda_minus = lbda - 10. * (np.max(J_list) + A_single)
        
        mu_plus = convolve_spatial_kernel(rho, J_list, Delta_x) + one_point_source(xi, xs, A_single, w_single) + lbda_plus
        mu_minus = convolve_spatial_kernel(rho, J_list, Delta_x) + one_point_source(xi, xs, A_single, w_single) + lbda_minus
        
        rho_plus = fermi(mu_plus, T)
        rho_minus = fermi(mu_minus, T)
        
        total_act_plus = Delta_x * np.sum(rho_plus)
        total_act_minus = Delta_x * np.sum(rho_minus)
        print('Weitere aeussere Iteration')
        print('lbda=', lbda)

        
        while(np.abs(total_act_minus - f) > 0.0000001 and np.abs(total_act_plus - f) > 0.0000001):
            lbda_mean = 0.5 * (lbda_minus + lbda_plus)
            
            mu_mean = convolve_spatial_kernel(rho, J_list, Delta_x) + one_point_source(xi, xs, A_single, w_single) + lbda_mean
            
            rho_mean = fermi(mu_mean, T)
            
            total_act_mean = Delta_x * np.sum(rho_mean)
                    
            if(np.sign(total_act_mean - f) == np.sign(total_act_minus - f)):
                lbda_minus = lbda_mean
                total_act_minus = total_act_mean
                closest = -1
            else:
                lbda_plus = lbda_mean
                total_act_plus = total_act_mean
                closest = 1
        
        if(closest == -1):
            lbda = lbda_minus
        else:
            lbda = lbda_plus
        
        mu = convolve_spatial_kernel(rho, J_list, Delta_x) + one_point_source(xi, xs, A_single, w_single) + lbda
        rho = 0.1 * fermi(mu, T) + 0.9 * rho
        
        total_act = Delta_x * np.sum(rho)
        
        print('total_act=', total_act)
        
        # plt.plot(xs, rho_old)
        
        iter_step += 1
        print('np.sum(np.abs(rho_old - rho))/np.sum(np.abs(rho)) =', np.sum(np.abs(rho_old - rho))/np.sum(np.abs(rho)))
    
    total_act = Delta_x * np.sum(rho)
    print('total_act=', total_act)
    print('Funktion iter_rho fast vorbei.')        
    return rho, lbda
