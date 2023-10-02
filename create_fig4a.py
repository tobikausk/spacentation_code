#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 19:05:26 2022

@author: tobias
"""
# %%
import numpy as np
from calc_CrossCorr_fixedMean_disorder import calc_rho_cross_corr_theo
from calc_FisherInfo import FisherInfo_from_corr
# import plotfuncs_PRL as pf
import visualization_prx as vis
from matplotlib import pyplot as plt


N_theo_points = 80 #200
xi = 0.2
mean_act = 0.15 # 0.2 
A_one_point = 0.01 # 0.45 
w_one_point = 0.07 
w_two_point = 0.1 
T = 0.05 # 0.2 

N_x = 100
# N_sample_Gauss = 1000
# N_sample_Gauss = 30000
N_sample_Gauss = 100000

tol_phi = 0.001
tol_phi_abs = 0.001
tol_q = 0.0001
tol_q_abs = 0.0001
# eps_relax_phi = 0.3
# eps_relax_rho = 0.3
eps_relax_phi = 0.1
eps_relax_rho = 0.1

figurepath = 'figures/Scale_synapses/'

space_form = 'rectangle'

g_max = 0.4 #  2.    # Maximal plotted value therefore 0.16 (= 0.4 * 0.4)
A_two_point_max = 1. # Maximal plotted value therefore 0.4 (= 1 * 0.4)
r_scale_max = 0.4

calc_newly = True
# calc_newly = False


    

if(calc_newly == True):

    with open( figurepath + 'Parameters_scaling_synapses', 'wb') as f:        
    
        np.save(f, g_max)
        np.save(f, N_theo_points)
        np.save(f, xi)
        np.save(f, mean_act) 
        np.save(f, A_one_point) 
        np.save(f, w_one_point) 
        np.save(f, A_two_point_max)
        np.save(f, T) 
        np.save(f, N_x)
        np.save(f, N_sample_Gauss)

    r_scale_list = np.linspace(0., r_scale_max, N_theo_points, endpoint=False)
    g_list = r_scale_list * g_max
    A_two_point_list = r_scale_list * A_two_point_max

    FisherInfo_theo_fine_list = np.zeros_like(r_scale_list)


    for ii in range(0, N_theo_points):    
        (xs_theo, J_dist_list, rho_final, lbda_final, variance, corr_theo, corr_theo_direct, corr_theo_indirect, 
            corr_theo_only_cross) = calc_rho_cross_corr_theo(xi, mean_act, g_list[ii], A_one_point, w_one_point, 
                                        A_two_point_list[ii], w_two_point, T, N_x, N_sample_Gauss, space_form,
                                        tol_phi, tol_phi_abs, tol_q, tol_q_abs, eps_relax_phi, eps_relax_rho)
        
                                                        
        FisherInfo_theo = FisherInfo_from_corr(corr_theo, xs_theo, xi, A_one_point, w_one_point, T)/N_x
        FisherInfo_theo_fine_list[ii] = FisherInfo_theo

    with open( figurepath + 'FisherInfo_theo_list_T=' + str(T) + '_f=' + str(mean_act) 
                            + '_A_one_point=' + str(A_one_point) +  '_w_one_point=' + str(w_one_point) 
                            + '_A_two_point_max=' + str(A_two_point_max) + '_w_two_point=' + str(w_two_point)
                            + '_g_max='+ str(g_max) + '_xi='+str(xi) + '_N_x=' + str(N_x) 
                            + '_N_sample_Gauss=' + str(N_sample_Gauss), 'wb') as f:        
        
        np.save(f, r_scale_list)
        np.save(f, g_list)
        np.save(f, A_two_point_list)
        np.save(f, FisherInfo_theo_fine_list)
        
else:
    with open( figurepath + 'FisherInfo_theo_list_T=' + str(T) + '_f=' + str(mean_act) 
                            + '_A_one_point=' + str(A_one_point) +  '_w_one_point=' + str(w_one_point) 
                            + '_A_two_point_max=' + str(A_two_point_max) + '_w_two_point=' + str(w_two_point)
                            + '_g_max='+ str(g_max) + '_xi='+str(xi) + '_N_x=' + str(N_x)
                            + '_N_sample_Gauss=' + str(N_sample_Gauss), 'rb') as f:
    
        
        r_scale_list = np.load(f)
        g_list = np.load(f)
        A_two_point_list = np.load(f)
        FisherInfo_theo_fine_list = np.load(f)
        
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

ax.set_title(r'\textbf{a}',loc='left')

plt.plot(r_scale_list/r_scale_max, FisherInfo_theo_fine_list, color = 'black')
plt.xlabel(r'scaling factor $r$')
plt.ylabel(r'$\frac{{\cal I}_{n}\left(\xi\right)}{N}$', rotation = 0)
plt.yticks([0.068, 0.069, 0.070])

figurename = 'FisherInfo_for_scaled_synapses_A_one_point=' + str(A_one_point) + '_N_sample_Gauss=' + str(N_sample_Gauss)

plt.savefig(figurepath + figurename + '.pdf')
plt.savefig(figurepath + figurename + '.eps')
plt.savefig(figurepath + figurename + '.jpg')

#create instance of visualization, this sets all the matplotlib rcParams
visualization = vis.visualization()  

# set standard figsize to one column according  to prx guidelines
visualization.set_SCI_1column_fig_style(ratio=vis.panel_wh_ratio)  # here, also another width/height ratio can be entered, eg: 3.

# create fig, here only one panel. Use constrained_layout to get a figure with nice boundaries fitting to the plots
fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

plt.plot(r_scale_list/r_scale_max, FisherInfo_theo_fine_list, color = 'black')
plt.xlabel('scaling factor')
plt.ylabel(r'Fisher info per neuron $\frac{{\cal I}_{n}\left(\xi\right)}{N}$', rotation = 90)
plt.yticks([0.068, 0.069, 0.070])

figurename = 'FisherInfo_for_talk_for_scaled_synapses_A_one_point=' + str(A_one_point) + '_N_sample_Gauss=' + str(N_sample_Gauss)

plt.savefig(figurepath + figurename + '.pdf')
plt.savefig(figurepath + figurename + '.eps')
plt.savefig(figurepath + figurename + '.jpg')


    

# %%
