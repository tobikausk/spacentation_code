#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 00:55:25 2022

@author: tobias
"""
# %%
# from load_and_plot_activities_deprecated import load_PairwCorr_error, load_MeanAct_error
from load_parameters import load_parameters
from load_and_plot_activities import load_PairwCorr_error, load_MeanAct_error
# import parameters_awi as params
from matplotlib import pyplot as plt

from attractor_with_input import J_w, rel_center, iter_rho, rho_T0, mu_from_rho
from solve_sp_equations_space_disorder import phi_q_self_consistent
from calc_CrossCorr_fixedMean_disorder import pairwise_corr, stat_upto_fourth_order, Hess_G, calc_rho_cross_corr_theo
from calc_FisherInfo import FisherInfo_from_corr
import visualization_prx as vis

import numpy as np

import os
import shutil

# %%
#TODO: Build in test, if parameters belonging to one time_suffix_list really
# differ in only one variable.

#Input-driven:
time_suffix_list = ['2022_07_11_20_31_46', '2022_07_11_20_32_19', '2022_07_11_20_32_50',
                    '2022_07_11_20_33_22', '2022_07_11_20_33_52', '2022_07_11_20_34_22']

#Nearly autonomous:
# time_suffix_list = ['2022_07_21_17_08_49', '2022_07_21_17_09_46', '2022_07_21_17_10_36', 
#                     '2022_07_21_17_11_12', '2022_07_21_17_11_57', '2022_07_21_17_12_34']

tol_phi = 0.001
tol_phi_abs = 0.001
tol_q = 0.001
tol_q_abs = 0.001
eps_relax_phi = 0.1
eps_relax_rho = 0.1

datapath_theo = 'data/Theory/Comparison_to_MC/'

figurepath = 'figures/FisherInfo_vs_g/'

foldername = time_suffix_list[0] + '/'
figurepath_specific = os.path.join(figurepath, foldername)
os.makedirs(figurepath_specific, exist_ok=True)
shutil.copy('data/MC_sim/Fisher_Information/' + time_suffix_list[0] + '/Parameters_attractor_' + time_suffix_list[0] + '.txt', 
                figurepath_specific)

calc_theory_newly = True
# calc_theory_newly = False

def round_to_read(x):
    if(np.round(x) == x):
        return int(x)
    else: 
        return x

N_sample_Gauss = 10000


N_theo_points = 100

FisherInfo_num_mean_list = np.zeros(len(time_suffix_list))
FisherInfo_num_rms_list = np.zeros(len(time_suffix_list))

FisherInfo_theo_fine_list = np.zeros(N_theo_points)

FisherInfo_theo_var_fine_list = np.zeros(N_theo_points)
FisherInfo_theo_crosscorr_direct_fine_list = np.zeros(N_theo_points)
FisherInfo_theo_crosscorr_indirect_fine_list = np.zeros(N_theo_points)

g_b_list = np.zeros(len(time_suffix_list))

for kk, time_suffix in enumerate(time_suffix_list):
    datapath = 'data/MC_sim/Fisher_Information/' + time_suffix + '/'

    (T, A_two_point_b, w_two_point, A_one_point_b, 
    w_one_point, g_b, xi, mean_act, connect_shape, conserved,
    N_pf, N_thermal, N_measure, N_MC_runs, N_J, N_r) = load_parameters(datapath, True)
    
    if(connect_shape == 0):
        space_form = 'Gauss'
    elif(connect_shape == 1):
        space_form = 'rectangle'
    elif(connect_shape == 2):
        space_form = 'ferro'

    xs_num = np.linspace(0.,1., N_pf, endpoint = False)
        
    corr_collection = np.zeros((N_J, N_pf, N_pf), dtype=float)
    corr_only_cross_collection = np.zeros_like(corr_collection)
    mean_act_collection = np.zeros((N_J, N_pf), dtype=float)
    
    for seed_no in range(0, N_J):
    
        corr = load_PairwCorr_error(datapath, space_form, 'MC', round_to_read(T), 
                                                      mean_act, round_to_read(A_one_point_b), w_one_point, 
                                                      round_to_read(A_two_point_b), w_two_point, 
                                                      round_to_read(g_b), 
                                                      xi, N_pf, seed_no, 'connected')
     
        
        mean_activities = load_MeanAct_error(datapath, space_form, 'MC', round_to_read(T), 
                                                      mean_act, round_to_read(A_one_point_b), w_one_point, 
                                                      round_to_read(A_two_point_b), w_two_point, 
                                                      round_to_read(g_b), 
                                                      xi, N_pf, seed_no)
        
        corr_collection[seed_no, :, :] = corr
        mean_act_collection[seed_no, :] = mean_activities
        corr_only_cross_collection[seed_no, :, :] = corr - np.diag(np.diag(corr))
     
    
    corr_mean = np.mean(corr_collection, axis=0)
    corr_rms = np.sqrt(np.var(corr_collection, axis = 0)/N_J) # Estimate for the error on the mean
    corr_mean_only_cross = np.mean(corr_only_cross_collection, axis=0)
    act_mean = np.mean(mean_act_collection, axis=0)
    
    FisherInfo_num = FisherInfo_from_corr(corr_mean, xs_num, xi, A_one_point_b, w_one_point, T)
    FisherInfo_num_rms = FisherInfo_from_corr(corr_rms, xs_num, xi, A_one_point_b, w_one_point, T)
    
    FisherInfo_num_mean_list[kk] = FisherInfo_num/N_pf
    FisherInfo_num_rms_list[kk] = FisherInfo_num_rms/N_pf
    
    g_b_list[kk] = g_b
    
    
if(calc_theory_newly == True):
    g_b_fine_list = np.linspace(0., np.max(g_b_list), N_theo_points, endpoint=False)
        
    for ii, g in enumerate(g_b_fine_list):
            (xs_theo, J_dist_list, rho_final, lbda_final, var_mat, corr_theo, corr_theo_direct, corr_theo_indirect, 
             corr_theo_only_cross) = calc_rho_cross_corr_theo(xi, mean_act, g, A_one_point_b, w_one_point, 
                                        A_two_point_b, w_two_point, T, N_pf, N_sample_Gauss, space_form,
                                        tol_phi, tol_phi_abs, tol_q, tol_q_abs, eps_relax_phi, eps_relax_rho)
                                                              
            FisherInfo_theo = FisherInfo_from_corr(corr_theo, xs_theo, xi, A_one_point_b, w_one_point, T)/N_pf
            
            FisherInfo_var_theo = FisherInfo_from_corr(var_mat, xs_theo, xi, A_one_point_b, w_one_point, T)/N_pf
            FisherInfo_crosscorr_direct_theo = FisherInfo_from_corr(corr_theo_direct, xs_theo, xi, A_one_point_b, w_one_point, T)/N_pf
            FisherInfo_crosscorr_indirect_theo = FisherInfo_from_corr(corr_theo_indirect, xs_theo, xi, A_one_point_b, w_one_point, T)/N_pf
            
            FisherInfo_theo_fine_list[ii] = FisherInfo_theo
            
            FisherInfo_theo_var_fine_list[ii] = FisherInfo_var_theo
            FisherInfo_theo_crosscorr_direct_fine_list[ii] = FisherInfo_crosscorr_direct_theo
            FisherInfo_theo_crosscorr_indirect_fine_list[ii] = FisherInfo_crosscorr_indirect_theo
            
            
    with open( datapath_theo + 'FisherInfo_theo_list_T=' + str(T) + '_f=' + str(mean_act) 
                        + '_A_one_point=' + str(A_one_point_b) +  '_w_one_point=' + str(w_one_point) 
                        + '_A_two_point=' + str(A_two_point_b) + '_w_two_point=' + str(w_two_point) 
                        + '_xi='+str(xi) + '_N_x=' + str(N_pf), 'wb') as f:        
        
        np.save(f, g_b_fine_list)
        np.save(f, FisherInfo_theo_fine_list)
        np.save(f, FisherInfo_theo_var_fine_list)
        np.save(f, FisherInfo_theo_crosscorr_direct_fine_list)
        np.save(f, FisherInfo_theo_crosscorr_indirect_fine_list)

else:
    with open( datapath_theo + 'FisherInfo_theo_list_T=' + str(T) + '_f=' + str(mean_act) 
                        + '_A_one_point=' + str(A_one_point_b) +  '_w_one_point=' + str(w_one_point) 
                        + '_A_two_point=' + str(A_two_point_b) + '_w_two_point=' + str(w_two_point) 
                        + '_xi='+str(xi) + '_N_x=' + str(N_pf), 'rb') as f:

        g_b_fine_list = np.load(f)
        FisherInfo_theo_fine_list = np.load(f)                   
        FisherInfo_theo_var_fine_list = np.load(f)
        FisherInfo_theo_crosscorr_direct_fine_list = np.load(f)
        FisherInfo_theo_crosscorr_indirect_fine_list = np.load(f)

# %%
################
### Plotting ###
################

#create instance of visualization, this sets all the matplotlib rcParams
visualization = vis.visualization()  

# set standard figsize to one column according  to prx guidelines
visualization.set_SCI_1column_fig_style(ratio=vis.panel_wh_ratio)  # here, also another width/height ratio can be entered, eg: 3.

# create fig, here only one panel. Use constrained_layout to get a figure with nice boundaries fitting to the plots
fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

ax.set_title(r'\textbf{b}',loc='left')
# plt.figure(1045)
plt.plot(g_b_fine_list/T, FisherInfo_theo_fine_list, color = 'black', label='theory')
plt.errorbar(g_b_list/T, FisherInfo_num_mean_list, yerr = FisherInfo_num_rms_list, ls='none', color = 'gray', label='Monte Carlo')
plt.xlabel(r'disorder strength $g$')
plt.ylabel(r'$\frac{{\cal I}_{n} (\xi)}{N}$', rotation = 0, labelpad=10)
plt.legend()

plt.savefig(figurepath_specific + 'FisherInfo_vs_g_' + time_suffix_list[0] + '.eps')
plt.savefig(figurepath_specific + 'FisherInfo_vs_g_' + time_suffix_list[0] + '.pdf')
plt.savefig(figurepath_specific + 'FisherInfo_vs_g_' + time_suffix_list[0] + '.jpg')

#create instance of visualization, this sets all the matplotlib rcParams
visualization = vis.visualization()  

# set standard figsize to one column according  to prx guidelines
visualization.set_SCI_1column_fig_style(ratio=vis.panel_wh_ratio)  # here, also another width/height ratio can be entered, eg: 3.

# create fig, here only one panel. Use constrained_layout to get a figure with nice boundaries fitting to the plots
fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

ax.set_title(r'\textbf{b}',loc='left')

# plt.figure(1046)
plt.plot(g_b_fine_list/T, FisherInfo_theo_fine_list, color = 'black', label='total')
plt.plot(g_b_fine_list/T, FisherInfo_theo_var_fine_list, '--', color = 'black', label='variance')
plt.plot(g_b_fine_list/T, FisherInfo_theo_crosscorr_direct_fine_list, '--', color = 'gray', label='cov., direct')
plt.plot(g_b_fine_list/T, FisherInfo_theo_crosscorr_indirect_fine_list, color = 'gray', label='cov., indirect')
plt.xlabel(r'disorder strength $g$')
plt.ylabel(r'$\frac{{\cal I}_{n} (\xi)}{N}$', rotation = 0, labelpad=10)
plt.legend()

plt.savefig(figurepath_specific + 'FisherInfo_vs_g_different_contributions_' + time_suffix_list[0] + '.eps', bbox_inches = 'tight')
plt.savefig(figurepath_specific + 'FisherInfo_vs_g_different_contributions_' + time_suffix_list[0] + '.pdf', bbox_inches = 'tight')
plt.savefig(figurepath_specific + 'FisherInfo_vs_g_different_contributions_' + time_suffix_list[0] + '.jpg', bbox_inches = 'tight')


#create instance of visualization, this sets all the matplotlib rcParams
visualization = vis.visualization()  

# set standard figsize to one column according  to prx guidelines
visualization.set_SCI_1column_fig_style(ratio=vis.panel_wh_ratio)  # here, also another width/height ratio can be entered, eg: 3.

# create fig, here only one panel. Use constrained_layout to get a figure with nice boundaries fitting to the plots
fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

# plt.figure(1046)
plt.plot(g_b_fine_list/T, FisherInfo_theo_fine_list, color = 'black', label='theory')
plt.errorbar(g_b_list/T, FisherInfo_num_mean_list, yerr = FisherInfo_num_rms_list, ls='none', color = 'gray', label='Monte Carlo')
plt.plot(g_b_fine_list/T, FisherInfo_theo_var_fine_list, '--', color = 'blue')
plt.plot(g_b_fine_list/T, FisherInfo_theo_crosscorr_direct_fine_list, '--', color = 'green')
plt.plot(g_b_fine_list/T, FisherInfo_theo_crosscorr_indirect_fine_list, color = 'gray')
plt.xlabel(r'disorder strength $g$')
plt.ylabel(r'Fisher info. per neuron $\frac{{\cal I}_{n} (\xi)}{N}$', rotation = 90, labelpad=10)
plt.legend()

plt.savefig(figurepath_specific + 'FisherInfo_vs_g_for_talk_' + time_suffix_list[0] + '.eps', bbox_inches = 'tight')
plt.savefig(figurepath_specific + 'FisherInfo_vs_g_for_talk_' + time_suffix_list[0] + '.pdf', bbox_inches = 'tight')
plt.savefig(figurepath_specific + 'FisherInfo_vs_g_for_talk_' + time_suffix_list[0] + '.jpg', bbox_inches = 'tight')

# %%
