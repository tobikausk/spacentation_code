#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 19:05:26 2022

@author: tobias
"""
# %%
import numpy as np
from calc_CrossCorr_fixedMean_disorder import calc_rho_cross_corr_theo
from calc_FisherInfo import FisherInfo_from_corr, FisherInfo_from_corr_compare_noise_corr
# import plotfuncs_PRL as pf
import visualization_prx as vis
from matplotlib import pyplot as plt


N_theo_points = 80 #200
xi = 0.2
mean_act = 0.15 # 0.2 
A_one_point = 0.2  
w_one_point = 0.07 
w_two_point = 0.1 
# T = 0.05 # 0.2 
T = 1. 
g_disorder = 0

N_x = 100
N_sample_Gauss = 100000

tol_phi = 0.001
tol_phi_abs = 0.001
tol_q = 0.0001
tol_q_abs = 0.0001
eps_relax_phi = 0.1
eps_relax_rho = 0.1

figurepath = 'figures/Only_local_recurrency/'

space_form = 'rectangle'

A_two_point_max = 40. # Temperature set to 1.0 (instead of 0.05)
r_scale_max = 1.

calc_newly = True
# calc_newly = False

savefig = True


    

if(calc_newly == True):

    r_scale_list = np.linspace(0., r_scale_max, N_theo_points, endpoint=False)
    A_two_point_list = r_scale_list * A_two_point_max
    FisherInfo_theo_fine_list = np.zeros_like(r_scale_list)
    FisherInfo_theo_fine_uncorr_list = np.zeros_like(r_scale_list)

    FisherInfo_crosscov_theo_fine_list = np.zeros_like(r_scale_list)

    rho_collection = np.zeros((N_theo_points, N_x)) 



    for ii in range(0, N_theo_points):    
        (xs_theo, J_dist_list, rho_final, lbda_final, variance, corr_theo, corr_theo_direct, corr_theo_indirect, 
                corr_theo_only_cross) = calc_rho_cross_corr_theo(xi, mean_act, g_disorder, A_one_point, w_one_point, 
                                            A_two_point_list[ii], w_two_point, T, N_x, N_sample_Gauss, space_form,
                                            tol_phi, tol_phi_abs, tol_q, tol_q_abs, eps_relax_phi, eps_relax_rho)
            
                                                            
        FisherInfo_theo = FisherInfo_from_corr(corr_theo, xs_theo, xi, A_one_point, w_one_point, T)/N_x
        FisherInfo_crosscov_theo = FisherInfo_from_corr(corr_theo_direct, xs_theo, xi, A_one_point, w_one_point, T)/N_x

        FisherInfo_theo_alternative_formula, FisherInfo_theo_uncorr = FisherInfo_from_corr_compare_noise_corr(corr_theo, variance,
                                                                                                xs_theo, xi, A_one_point, w_one_point, T)

        FisherInfo_theo_fine_list[ii] = FisherInfo_theo
        FisherInfo_theo_fine_uncorr_list[ii] = FisherInfo_theo_uncorr/N_x

        FisherInfo_crosscov_theo_fine_list[ii] = FisherInfo_crosscov_theo

        rho_collection[ii, :] = rho_final

        with open( figurepath + 'FisherInfo_theo_list_T=' + str(T) + '_f=' + str(mean_act) 
                                + '_A_one_point=' + str(A_one_point) +  '_w_one_point=' + str(w_one_point) 
                                + '_A_two_point_max=' + str(A_two_point_max) + '_w_two_point=' + str(w_two_point)
                                + '_g_disorder='+ str(g_disorder) + '_xi='+str(xi) + '_N_x=' + str(N_x) 
                                + '_N_sample_Gauss=' + str(N_sample_Gauss), 'wb') as f:        
            
            np.save(f, r_scale_list)
            np.save(f, A_two_point_list)
            np.save(f, FisherInfo_theo_fine_list)
            np.save(f, FisherInfo_theo_fine_uncorr_list)
            np.save(f, xs_theo)
            np.save(f, rho_collection)
    
   
        
else:
    with open( figurepath + 'FisherInfo_theo_list_T=' + str(T) + '_f=' + str(mean_act) 
                            + '_A_one_point=' + str(A_one_point) +  '_w_one_point=' + str(w_one_point) 
                            + '_A_two_point_max=' + str(A_two_point_max) + '_w_two_point=' + str(w_two_point)
                            + '_g_disorder='+ str(g_disorder) + '_xi='+str(xi) + '_N_x=' + str(N_x)
                            + '_N_sample_Gauss=' + str(N_sample_Gauss), 'rb') as f:
    
        
        r_scale_list = np.load(f)
        A_two_point_list = np.load(f)
        FisherInfo_theo_fine_list = np.load(f)
        FisherInfo_theo_fine_uncorr_list = np.load(f)
        xs_theo = np.load(f)
        rho_collection = np.load(f)
        
 # %%       

############################
### Figures in APS style ###
############################


#create instance of visualization, this sets all the matplotlib rcParams
visualization = vis.visualization()  

# set standard figsize to one column according  to prx guidelines
visualization.set_SCI_1column_fig_style(ratio=vis.panel_wh_ratio)  # here, also another width/height ratio can be entered, eg: 3.

# create fig, here only one panel. Use constrained_layout to get a figure with nice boundaries fitting to the plots
fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

ax.set_title(r'\textbf{b}',loc='left')

plt.plot(A_two_point_list/T, FisherInfo_theo_fine_list, color = 'black', label = 'with cov')
plt.plot(A_two_point_list/T, FisherInfo_theo_fine_uncorr_list, '--', color = 'black', label = 'without cov')
plt.xlabel(r'strength of recurrency $K_{\mathrm{rec}}$')
plt.ylabel(r'$\frac{{\cal I}_{n}\left(\xi\right)}{N}$', rotation = 0)
plt.legend()

figurename = 'FisherInfo_for_scaling_recurrency_A_two_point_max=' + str(A_two_point_max) + '_A_one_point=' + str(A_one_point) + '_N_sample_Gauss=' + str(N_sample_Gauss)

if(savefig == True):
    plt.savefig(figurepath + figurename + '.pdf', bbox_inches="tight")
    plt.savefig(figurepath + figurename + '.eps', bbox_inches="tight")
    plt.savefig(figurepath + figurename + '.jpg', bbox_inches="tight")

#create instance of visualization, this sets all the matplotlib rcParams
visualization = vis.visualization()  

# set standard figsize to one column according  to prx guidelines
visualization.set_SCI_1column_fig_style(ratio=vis.panel_wh_ratio)  # here, also another width/height ratio can be entered, eg: 3.

# create fig, here only one panel. Use constrained_layout to get a figure with nice boundaries fitting to the plots
fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

plt.plot(A_two_point_list/T, FisherInfo_theo_fine_list, color = 'black', label = 'with cov')
plt.plot(A_two_point_list/T, FisherInfo_theo_fine_uncorr_list, '--', color = 'black', label = 'without cov')
plt.xlabel('strength of recurrency')
plt.ylabel(r'Fisher info per neuron $\frac{{\cal I}_{n}\left(\xi\right)}{N}$', rotation = 90)
plt.legend()

figurename = 'FisherInfo_for_talk_scaling_recurrency_A_two_point_max=' + str(A_two_point_max) + '_A_one_point=' + str(A_one_point) + '_N_sample_Gauss=' + str(N_sample_Gauss)

if(savefig == True):
    plt.savefig(figurepath + figurename + '.pdf', bbox_inches="tight")
    plt.savefig(figurepath + figurename + '.eps', bbox_inches="tight")
    plt.savefig(figurepath + figurename + '.jpg', bbox_inches="tight")


#create instance of visualization, this sets all the matplotlib rcParams
visualization = vis.visualization()  

# set standard figsize to one column according  to prx guidelines
visualization.set_SCI_1column_fig_style(ratio=vis.panel_wh_ratio)  # here, also another width/height ratio can be entered, eg: 3.

# create fig, here only one panel. Use constrained_layout to get a figure with nice boundaries fitting to the plots
fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

ax.set_title(r'\textbf{a}',loc='left')

color_list = ['black', 'blue', 'cyan', 'gray', 'green']

for kk, ii in enumerate([0, int(0.6*N_theo_points),int(0.7*N_theo_points), int(0.8*N_theo_points), int(0.9*N_theo_points)]):
    plt.plot(xs_theo, rho_collection[ii, :], 
             label=r'$K_{\mathrm{rec}}=$' + str(np.round(A_two_point_list[ii]/T )), color = color_list[kk])

plt.xlabel(r'place-field center $x$')
plt.ylabel(r'Activity $\rho(x)$')

plt.legend()


figurename = 'Tuning_curve_var_A_two_point_max=' + str(A_two_point_max) + 'A_one_point=' + str(A_one_point) + '_N_sample_Gauss=' + str(N_sample_Gauss)

if(savefig == True):
    plt.savefig(figurepath + figurename + '.pdf', bbox_inches="tight")
    plt.savefig(figurepath + figurename + '.eps', bbox_inches="tight")
    plt.savefig(figurepath + figurename + '.jpg', bbox_inches="tight")

# %%
