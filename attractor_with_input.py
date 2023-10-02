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
    # return A * np.exp(- (xi - x) * (xi - x) /(2.*w*w))
    Delta_xix = rel_center(xi, x)
    return A * np.exp(- Delta_xix * Delta_xix /(2.*w*w))

def deriv_one_point_source(xi, x, A, w):
    Delta_xix = rel_center(xi, x)    
    # return - (A * (xi - x)/(w*w)) *  np.exp(- (xi - x) * (xi - x) /(2.*w*w))
    return - (A * Delta_xix/(w*w)) *  np.exp(- Delta_xix * Delta_xix /(2.*w*w))

# def phi(r, phi0, DeltaPhi, w):
# # def phi(r):
#     # return params.phi0 + params.DeltaPhi * np.exp(- (r/params.w)**2 )
#     return phi0 + DeltaPhi * np.exp(- (r/w)**2 )

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
    # return H_mod(0.5 * f - np.abs(x_rel_center))
    return 1. - H_mod(-0.5 * f + np.abs(x_rel_center))


# def J_w(x, w, J0):
#     return J0 * (H(0.5 * w - np.abs(x)) + H(0.5 * w - np.abs(1. - x))) 

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
    
    # return (deriv_fermi(mu_from_rho(x0, xi, A_single, w_single, rho, J_list, Delta_x, lbda), T)/T) * (
    #        - Jw_0 * (rho_shift_minus_half_w - rho_shift_plus_half_w) 
    #        + deriv_one_point_source(xi,x0, A_single, w_single))**2
    
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
    
    # plt.figure(4)
    
    iter_step = 0
    # while(np.sum(np.abs(rho_old - rho))/np.sum(np.abs(rho)) > 0.0000001 or iter_step==0):
    while(np.sum((rho_old - rho)**2) > 0.0000001 or iter_step==0):
        rho_old = copy.deepcopy(rho)
        
        lbda_plus = lbda + 10. * (np.max(J_list) + A_single)
        lbda_minus = lbda - 10. * (np.max(J_list) + A_single)
        
        # lbda_plus = lbda + 1.
        # lbda_minus = lbda - 1.
        
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
            
            # print('total_act_mean=', total_act_mean)
            
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
        # rho = 0.5 * fermi(mu, T) + 0.5 * rho
        
        total_act = Delta_x * np.sum(rho)
        
        print('total_act=', total_act)
        
        # plt.plot(xs, rho_old)
        
        iter_step += 1
        print('np.sum(np.abs(rho_old - rho))/np.sum(np.abs(rho)) =', np.sum(np.abs(rho_old - rho))/np.sum(np.abs(rho)))
    
    total_act = Delta_x * np.sum(rho)
    print('total_act=', total_act)
    print('Funktion iter_rho fast vorbei.')        
    return rho, lbda

# def iter_rho(rho_init, J_list, Delta_x, T, f, lbda_init):
    # rho = rho_init
    # rho_old = rho_init
    # lbda = lbda_init
    
    # plt.figure(4)
    
    # iter_step = 0
    # while(np.sum(np.abs(rho_old - rho))/np.sum(np.abs(rho)) > 0.001 or iter_step==0):
    #     rho_old = copy.deepcopy(rho)
        
    #     lbda_plus = lbda + 0.1
    #     lbda_minus = lbda - 0.1
        
    #     mu_plus = convolve_spatial_kernel(rho, J_list, Delta_x) + lbda_plus
    #     mu_minus = convolve_spatial_kernel(rho, J_list, Delta_x) + lbda_minus
        
    #     rho_plus = fermi(mu_plus, T)
    #     rho_minus = fermi(mu_minus, T)
        
    #     total_act_plus = Delta_x * np.sum(rho_plus)
    #     total_act_minus = Delta_x * np.sum(rho_minus)
    #     print('Weitere aeussere Iteration')
    #     print('lbda=', lbda)

        
    #     while(np.abs(total_act_minus - f) > 0.001 and np.abs(total_act_plus - f) > 0.001):
    #         lbda_mean = 0.5 * (lbda_minus + lbda_plus)
            
    #         mu_mean = convolve_spatial_kernel(rho, J_list, Delta_x) + lbda_mean
            
    #         rho_mean = fermi(mu_mean, T)
            
    #         total_act_mean = Delta_x * np.sum(rho_mean)
            
    #         print('Bisection')
            
    #         if(np.sign(total_act_mean - f) == np.sign(total_act_minus - f)):
    #             lbda_minus = lbda_mean
    #             total_act_minus = total_act_mean
    #             closest = -1
    #         else:
    #             lbda_plus = lbda_mean
    #             total_act_plus = total_act_mean
    #             closest = 1
        
    #     if(closest == -1):
    #         lbda = lbda_minus
    #     else:
    #         lbda = lbda_plus
        
    #     mu = convolve_spatial_kernel(rho, J_list, Delta_x) + lbda
    #     rho = 0.5 * fermi(mu, T) + 0.5 * rho
        
    #     plt.plot(xs, rho_old)
        
    #     iter_step += 1
    #     print('np.sum(np.abs(rho_old - rho))/np.sum(np.abs(rho)) =', np.sum(np.abs(rho_old - rho))/np.sum(np.abs(rho)))
        
    # print('Funktion iter_rho fast vorbei.')        
    # return rho, lbda
    
 
       
        
    
    

# # plt.figure(9)
# # plt.plot(connect_widths, Fisher_Info_integrated_list, label='ana')
# # plt.plot(connect_widths, Fisher_Info_numerical_integrated_list, label='num')
# # plt.xlabel('Connection width')
# # plt.ylabel('Fisher information')
# # plt.legend()

# plt.figure(10)
# plt.plot(mean_activities, lbda_list)
# plt.xlabel('Mean activity')
# plt.ylabel('lambda')

# plt.figure(99)
# plt.plot(mean_activities, Fisher_Info_integrated_list, label='ana, A_one_point='+str(params.A_single))
# plt.plot(mean_activities, Fisher_Info_numerical_integrated_list, label='num, A_one_point='+str(params.A_single))
# plt.xlabel('Mean activity')
# plt.ylabel('Fisher information')
# plt.legend()

# plt.figure(999)
# plt.plot(mean_activities, Fisher_Info_integrated_list, label='ana, A_one_point='+str(params.A_single))
# # plt.plot(mean_activities, Fisher_Info_numerical_integrated_list, label='num, A_one_point='+str(params.A_single))
# plt.xlabel('Mean activity')
# plt.ylabel('Fisher information')
# plt.legend()

# # # connect_strengths = [0., 0.7, 1.5]
# # # connect_strengths = [1.]
# # # connect_widths = [0.05, 0.1, 0.2, 0.25, 0.29]
# # # connect_widths = [0.05, 0.2, 0.29]
# # connect_widths = [0.05, 0.1, 0.2, 0.25, 0.3, 0.35]
# # # connect_widths = [0.25, 0.3, 0.35]
# # # connect_strengths = [0., 0.1, 0.3] #, 0.7, 1., 1.2]
# # Fisher_Info_integrated_list = np.zeros_like(connect_widths)
# # Fisher_Info_numerical_integrated_list = np.zeros_like(connect_widths)
# # lbda_list = np.zeros_like(connect_widths)

# # for ii, w_two_point in enumerate(connect_widths):
# #     J_list =  J_w(xs, w_two_point, params.Jspace_0)
    
# #     mu_0 = mu_T0(xs, params.x0, w_two_point, params.f)
# #     rho_0 = rho_T0(xs, params.x0, params.f)
    
# #     mu_0_half_step = mu_T0(xs_half_step, params.x0, w_two_point, params.f)
# #     rho_0_half_step = rho_T0(xs_half_step, params.x0, params.f)
    
# #     # plt.figure(1)
# #     # plt.plot(xs, mu_0)
    
# #     # plt.figure(2)
# #     # plt.plot(xs, rho_0)
    
# #     plt.figure(2)
# #     plt.plot(xs, fermi(mu_0,0.), '*')
# #     plt.plot(xs, rho_0)
    
# #     plt.figure(3)
# #     plt.plot(xs,convolve_spatial_kernel(rho_0, J_list, params.Delta_x)-0.5*w_two_point,'*')
# #     plt.plot(xs, mu_0)
    
# #     rho_final, lbda_final = iter_rho(rho_0, 
# #                                      J_list, params.Delta_x, params.T, params.f, params.A_single, params.w_tc, params.xi, 
# #                                      xs, -0.5*w_two_point)
# #     rho_final_half_step, lbda_final_half_step = iter_rho(rho_0_half_step, 
# #                                                          J_list, params.Delta_x, params.T, params.f, params.A_single, params.w_tc, params.xi, 
# #                                                          xs_half_step, -0.5*w_two_point)
    
    
# #     mu_final_half_step = mu_from_rho(xs, params.xi, params.A_single, params.w_tc, rho_final_half_step, J_list, params.Delta_x, lbda_final_half_step)
# #     mu_final_without_source_half_step = mu_from_rho(xs_half_step, params.xi, 0., params.w_tc, rho_final_half_step, J_list, params.Delta_x, lbda_final_half_step)
# #     mu_final_without_source = mu_from_rho(xs, params.xi, 0., params.w_tc, rho_final, J_list, params.Delta_x, lbda_final_half_step)
    
# #     lbda_list[ii] = lbda_final
    
    
# #     rho_final_check = fermi(mu_from_rho(xs, params.xi, params.A_single, params.w_tc, rho_final, J_list, params.Delta_x, lbda_final), params.T)
    
    
# #     plt.figure(5)
# #     plt.plot(xs, rho_0)
# #     plt.plot(xs, fermi(convolve_spatial_kernel(rho_final, J_list, params.Delta_x)
# #                        + one_point_source(params.xi, xs, params.A_single, params.w_tc)+lbda_final, params.T), '*')
# #     plt.plot(xs, rho_final)
# #     plt.title('figure 5')
    
# #     # deriv_rho_final = np.diff(rho_final_half_step)/Delta_x
# #     # deriv_rho_final = np.hstack((deriv_rho_final, [(rho_final_half_step[0] - rho_final_half_step[-1])/Delta_x ]))
    
    
# #     deriv_rho_final = (np.roll(rho_final_half_step,1) - rho_final_half_step)/params.Delta_x
# #     deriv_mu_final = (np.roll(mu_final_half_step,1) - mu_final_half_step)/params.Delta_x
    
# #     # deriv_mu_indirect_final_num = (np.roll(mu_final_without_source_half_step,1) - mu_final_without_source_half_step)/Delta_x
# #     deriv_mu_indirect_final_num = (np.roll(mu_final_without_source,1) - mu_final_without_source)/params.Delta_x
# #     deriv_mu_direct_final_num = (np.roll(one_point_source(params.xi, xs_half_step, params.A_single, params.w_tc),1) 
# #                                    - one_point_source(params.xi, xs_half_step, params.A_single, params.w_tc))/params.Delta_x
    
# #     diff_rho_ana = del_rho_del_xi(rho_final, w_two_point, xs, params.xi, params.A_single, params.w_tc, 
# #                                   lbda_final, J_list, params.T, params.Delta_x, params.Jspace_0)
# #     diff_mu_ana, diff_mu_indirect_ana, diff_mu_direct_ana = del_mu_del_xi(rho_final, w_two_point, xs, params.xi, 
# #                                                                           params.A_single, params.w_tc, params.Delta_x, params.Jspace_0)

# #     Fisher_Info_ana, prefactor_ana, indirect_shift_ana, direct_shift_ana = Fisher_Info(rho_final, w_two_point, xs, params.xi, params.A_single, params.w_tc,
# #                                                                                         lbda_final, J_list, params.T, params.Delta_x, params.Jspace_0)
# #     # Fisher_Info_num, prefactor_num, indirect_shift_num, direct_shift_num = Fisher_Info_numerical(rho_final, w_two_point, xs, xi, A_single, w_tc, lbda_final, J_list, T, Delta_x, Jspace_0)
# #     Fisher_Info_num = Fisher_Info_numerical(rho_final, rho_final_half_step, params.Delta_x)
    
# #     Fisher_Info_integrated_list[ii] = Fisher_Info_integrated(rho_final, w_two_point, xs, params.xi, params.A_single, params.w_tc, 
# #                                                              lbda_final, J_list, params.T, params.Delta_x, params.Jspace_0)
# #     Fisher_Info_numerical_integrated_list[ii] = Fisher_Info_numerical_integrated(rho_final, rho_final_half_step, params.Delta_x)
    
# #     if(np.abs(w_two_point - 0.05) < 0.0001 or np.abs(w_two_point - 0.2) < 0.0001 or np.abs(w_two_point - 0.35) < 0.0001):
    
# #         plt.figure(6)
# #         plt.plot(xs, deriv_rho_final**2, label='num, w_two_point='+str(w_two_point))
# #         plt.plot(xs, diff_rho_ana**2, '.', label='ana, w_two_point='+str(w_two_point))
# #         plt.title('deriv rho squared')
# #         plt.legend()
        
        
# #         plt.figure(7)
# #         plt.plot(xs, -deriv_rho_final, label='num, w_two_point='+str(w_two_point))
# #         plt.plot(xs, diff_rho_ana, '.', label='ana, w_two_point='+str(w_two_point))
# #         plt.title('deriv rho')
# #         plt.legend()
        
# #         plt.figure(77)
# #         plt.plot(xs, deriv_mu_final, label='num, w_two_point='+str(w_two_point))
# #         plt.plot(xs, diff_mu_ana, '.', label='ana, w_two_point='+str(w_two_point))
# #         plt.title('deriv mu')
# #         plt.xlim(0.2,0.4)
# #         plt.legend()
        
# #         plt.figure(87)
# #         plt.plot(xs, deriv_mu_indirect_final_num, label='num, w_two_point='+str(w_two_point))
# #         plt.plot(xs, diff_mu_indirect_ana, '.', label='ana, w_two_point='+str(w_two_point))
# #         plt.title('indirect part of deriv of mu')
# #         plt.xlim(0.2,0.4)
# #         plt.legend()
        
# #         plt.figure(97)
# #         plt.plot(xs, deriv_mu_direct_final_num, label='num, w_two_point='+str(w_two_point))
# #         plt.plot(xs, diff_mu_direct_ana, '.', label='ana, w_two_point='+str(w_two_point))
# #         plt.title('deriv mu')
# #         plt.xlim(0.2,0.4)
# #         plt.legend()
                
            
            
# #             # plt.figure(8)
# #             # plt.plot(xs, Fisher_Info_ana, label='Fisher Information')
            
# #         # plt.figure(111)
# #         # plt.plot(xs, prefactor_ana, label='prefactor ana w_two_point='+str(w_two_point))
# #         # plt.plot(xs[1:], prefactor_num, '*', label='prefactor num w_two_point='+str(w_two_point))
# #         # plt.legend()
        
# #         # plt.figure(211)
# #         # plt.plot(xs, direct_shift_ana**2, label='direct squared ana w_two_point='+str(w_two_point))
# #         # plt.plot(xs[1:], direct_shift_num**2, '*', label='direct squared num w_two_point='+str(w_two_point))
# #         # plt.legend()
        
# #         # plt.figure(311)
# #         # plt.plot(xs, 2. * direct_shift_ana * indirect_shift_ana, label='2x direct * indirect  ana w_two_point='+str(w_two_point))
# #         # plt.plot(xs[1:], 2. * direct_shift_num * indirect_shift_num, '*', label='2x direct indirect squared num w_two_point='+str(w_two_point))
# #         # plt.legend()
        
# #         plt.figure(11)
# #         plt.plot(xs, Fisher_Info_ana, label='ana w_two_point='+str(w_two_point))
# #         plt.plot(xs, Fisher_Info_num, label='num w_two_point='+str(w_two_point))
# #         plt.xlabel('x')
# #         plt.ylabel('Local Fisher Information')
# #         # plt.xlim(0.2,0.4)
# #         plt.legend()
            
# #         plt.figure(12)
# #         plt.plot(xs, rho_final, label='w_two_point='+str(w_two_point))
# #         plt.plot(xs, rho_final_check, '*', label='check for w_two_point='+str(w_two_point))
# #         plt.xlabel('x')
# #         plt.ylabel(r'$\rho_{0}$')
# #         plt.legend()
        
        
    
    

# # # plt.figure(9)
# # # plt.plot(connect_widths, Fisher_Info_integrated_list, label='ana')
# # # plt.plot(connect_widths, Fisher_Info_numerical_integrated_list, label='num')
# # # plt.xlabel('Connection width')
# # # plt.ylabel('Fisher information')
# # # plt.legend()

# # plt.figure(10)
# # plt.plot(connect_widths, lbda_list)
# # plt.xlabel('Connection width')
# # plt.ylabel('lambda')

# # plt.figure(99)
# # plt.plot(connect_widths, Fisher_Info_integrated_list, label='ana, w_tc='+str(params.w_tc))
# # plt.plot(connect_widths, Fisher_Info_numerical_integrated_list, label='num, w_tc='+str(params.w_tc))
# # plt.xlabel('Connection width')
# # plt.ylabel('Fisher information')
# # plt.legend()



# ##############################################################################
# ### Variation of connection strength #########################################
# ##############################################################################

# # # connect_strengths = [0., 0.7, 1.5]
# # # connect_strengths = [1.]
# # connect_strengths = [0., 0.1, 0.3, 0.7, 1., 1.2, 1.5]
# # # connect_strengths = [0., 0.1, 0.3] #, 0.7, 1., 1.2]
# # Fisher_Info_list = np.zeros_like(connect_strengths)
# # lbda_list = np.zeros_like(connect_strengths)

# # for ii, Jspace_0 in enumerate(connect_strengths):
# #     J_list = J_w(xs, w_J, Jspace_0)
    
# #     mu_0 = mu_T0(xs, x0, w_J, f)
# #     rho_0 = rho_T0(xs, x0, f)
    
# #     plt.figure(1)
# #     plt.plot(xs, mu_0)
    
# #     # plt.figure(2)
# #     plt.plot(xs, rho_0)
    
# #     plt.figure(2)
# #     plt.plot(xs, fermi(mu_0,0.), '*')
# #     plt.plot(xs, rho_0)
    
# #     plt.figure(3)
# #     plt.plot(xs,convolve_spatial_kernel(rho_0, J_list, Delta_x)-0.5*w_J,'*')
# #     plt.plot(xs, mu_0)
    
# #     rho_final, lbda_final = iter_rho(rho_0, J_list, Delta_x, T, f, A_single, w_tc, xi, xs, -0.5*w_J)
    
# #     lbda_list[ii] = lbda_final
    
# #     plt.figure(5)
# #     plt.plot(xs, rho_0)
# #     plt.plot(xs, fermi(convolve_spatial_kernel(rho_final, J_list, Delta_x)+ one_point_source(xi, xs, A_single, w_tc)+lbda_final, T), '*')
# #     plt.plot(xs, rho_final)
# #     plt.title('figure 5')
    
# #     deriv_rho_final = np.diff(rho_final)/Delta_x
    
# #     diff_rho_ana = del_rho_del_xi(rho_final, w_J, xs, xi, A_single, w_tc, lbda_final, J_list, T, Delta_x, Jspace_0)
    
# #     if(np.abs(Jspace_0 - 0.0) < 0.0001 or np.abs(Jspace_0 - 0.7) < 0.0001 or np.abs(Jspace_0 - 1.5) < 0.0001):
    
# #         plt.figure(6)
# #         plt.plot(xs[1:], deriv_rho_final**2, label='num')
# #         plt.plot(xs, diff_rho_ana**2, '*', label='ana')
# #         plt.legend()
        
# #         plt.figure(7)
# #         plt.plot(xs[1:], deriv_rho_final, label='num')
# #         plt.plot(xs, diff_rho_ana, '*', label='ana')
# #         plt.legend()
        
# #         Fisher_Info_ana = Fisher_Info(rho_final, w_J, xs, xi, A_single, w_tc, lbda_final, J_list, T, Delta_x, Jspace_0)
        
# #         plt.figure(8)
# #         plt.plot(xs, Fisher_Info_ana, label='Fisher Information')
    
# #         plt.figure(11)
# #         plt.plot(xs, Fisher_Info_ana, label='J0_space='+str(connect_strengths[ii]))
# #         plt.xlabel('x')
# #         plt.ylabel('Local Fisher Information')
# #         plt.xlim(0.2,0.4)
# #         plt.legend()
        
# #         plt.figure(12)
# #         plt.plot(xs, rho_final, label='J0_space='+str(connect_strengths[ii]))
# #         plt.xlabel('x')
# #         plt.ylabel(r'$\rho_{0}$')
# #         plt.legend()
    
# #     Fisher_Info_list[ii] = Fisher_Info_integrated(rho_final, w_J, xs, xi, A_single, w_tc, lbda_final, J_list, T, Delta_x, Jspace_0)

# # plt.figure(9)
# # plt.plot(connect_strengths, Fisher_Info_list)
# # plt.xlabel('Connection strength')
# # plt.ylabel('Fisher information')

# # plt.figure(10)
# # plt.plot(connect_strengths, lbda_list)
# # plt.xlabel('Connection strength')
# # plt.ylabel('lambda')