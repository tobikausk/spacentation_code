# %%
import numpy as np
from matplotlib import pyplot as plt
import visualization_prx as vis
#Specific dependencies:
from load_full_sim_data import load_full_sim_data, load_parameters_without_xi, load_xis, generate_file_name
from linear_regression import comp_cov_train_test, lin_reg_with_res, lin_reg_regularized
from import_xis import load_xis
from bias_lin_est import compute_bias, compute_deriv

# %%
datapath = 'data/MC_sim/Linear_readout/'
datapath_theo = 'data/Theory/'
figurepath = 'figures/Linear_readout/'

# time_stamp = '2023_02_23_14_41_05'

# time_stamp_list = ['2023_03_06_23_13_54', '2023_03_06_23_14_19', '2023_03_06_23_15_29', 
#                    '2023_03_07_19_51_13', '2023_03_07_19_59_19', '2023_03_07_16_38_47']

time_stamp = '2023_09_25_22_28_02'

time_stamp_list = ['2023_09_25_22_28_02', '2023_09_29_17_09_59']


time_stamp_list_lumped = '' # Puts all elements of time_stamp_list into one long string
for time_stamp in time_stamp_list:
    time_stamp_list_lumped += time_stamp + '_'

##############################################
### L2 regularization of linear regression ###
##############################################
lbda_list = [0.01] 

xi_min = 0.4
xi_max = 0.6
N_xi = 1000

res_rnd_guess = (xi_max - xi_min)**2 / 12.

load_and_compute_newly = True

g_disord_list = np.zeros(len(time_stamp_list), dtype = float)

res_var_reg_list = np.zeros((len(lbda_list), len(g_disord_list)), dtype = float)
res_var_train_reg_list = np.zeros((len(lbda_list), len(g_disord_list)), dtype = float)

train_portion = 0.8

if(load_and_compute_newly == True):
    for kk, time_stamp in enumerate(time_stamp_list):

        (T, A_two_point, w_two_point, A_one_point,
        w_one_point, g_disord, mean_act, connect_shape,
        N_pf, N_thermal, N_time, N_MC_runs, N_xi, xi_min_load, xi_max_load, N_r) = load_parameters_without_xi(datapath + time_stamp + '/')


        xi_list = load_xis(datapath + time_stamp + '/')

        if(len(xi_list) < N_xi):
            raise Exception('Not enough data!')
        
        if(xi_min != xi_min_load or xi_max != xi_max_load):
            raise Exception('xi_min and xi_max do not match!')

        xi_list = xi_list[:N_xi]

        pop_act_sample = np.zeros((N_pf, N_time, N_xi))
        
        pf_centers = np.linspace(0, 1, int(N_pf), endpoint=False)

        for ii in range(0, N_xi):
            filename = generate_file_name(xi_list[ii], T, A_two_point, w_two_point, A_one_point,
                                            w_one_point, g_disord, mean_act, N_pf)

            pop_act_sample[:,:, ii] = load_full_sim_data(datapath + time_stamp + '/' + filename)

            if(ii < 10):
                pop_act_mean = np.mean(pop_act_sample[:,:, ii], axis = 1)

                plt.figure(1+ii)
                plt.plot(pf_centers, pop_act_mean)
                plt.title('Activity for xi=' + str(xi_list[ii]))


        (cov_act_train, cov_act_xi_train, var_xi_train, cov_act_test, cov_act_xi_test,
         var_xi_test) = comp_cov_train_test(pop_act_sample, xi_list, train_portion)


        for ii, lbda in enumerate(lbda_list):
            J_ridge, resid_train_ridge, resid_test_ridge = lin_reg_regularized(cov_act_train, cov_act_xi_train, var_xi_train,
                                                        cov_act_test,  cov_act_xi_test,  var_xi_test, lbda = lbda)

            res_var_reg_list[ii, kk] = resid_test_ridge
            res_var_train_reg_list[ii, kk] = resid_train_ridge

        g_disord_list[kk] = g_disord
        

    ##################################################
    ### Store results in a file for later analysis ###
    ##################################################

    with open( datapath_theo + 'Linear_readout/Residual_lin_readout_' + time_stamp_list_lumped, 'wb') as f:        
        
        np.save(f, g_disord_list)
        np.save(f, res_var_reg_list)
        np.save(f, res_var_train_reg_list)

else:
    with open( datapath_theo + 'Linear_readout/Residual_lin_readout_' + time_stamp_list_lumped, 'rb') as f:

        g_disord_list = np.load(f)
        res_var_reg_list = np.load(f)
        res_var_train_reg_list = np.load(f)


# %%
####################################################
### Load Fisher information curve for comparison ###
####################################################

xi_Fisher = 0.2

# Load in case it was not done before because the results of the linear regression were loaded from a file
(T, A_two_point, w_two_point, A_one_point,
        w_one_point, g_disord, mean_act, connect_shape,
        N_pf, N_thermal, N_time, N_MC_runs, N_xi, xi_min_load, xi_max_load, N_r) = load_parameters_without_xi(datapath + time_stamp + '/')

with open( datapath_theo + 'Comparison_to_MC/FisherInfo_theo_list_T=' + str(T) + '_f=' + str(mean_act) 
                        + '_A_one_point=' + str(A_one_point) +  '_w_one_point=' + str(w_one_point) 
                        + '_A_two_point=' + str(A_two_point) + '_w_two_point=' + str(w_two_point) 
                        + '_xi=' + str(xi_Fisher) + '_N_x=' + str(N_pf), 'rb') as f:

        g_b_fine_list = np.load(f)
        FisherInfo_theo_per_neur_list = np.load(f)                   
        FisherInfo_theo_var_per_neur_list = np.load(f)
        FisherInfo_theo_crosscorr_direct_per_neur_list = np.load(f)
        FisherInfo_theo_crosscorr_indirect_per_neur_list = np.load(f)


FisherInfo_theo_total_list = FisherInfo_theo_per_neur_list * N_pf 

# %%

#create instance of visualization, this sets all the matplotlib rcParams
visualization = vis.visualization()  

# set standard figsize to one column according  to prx guidelines
visualization.set_SCI_1column_fig_style(ratio=vis.panel_wh_ratio)  # here, also another width/height ratio can be entered, eg: 3.

# create fig, here only one panel. Use constrained_layout to get a figure with nice boundaries fitting to the plots
fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)


g_disord_list_argsort = np.argsort(g_disord_list)
g_disord_sorted_list = g_disord_list[g_disord_list_argsort]

res_var_reg_sorted_list = res_var_reg_list[:,g_disord_list_argsort]

for lbda in lbda_list:
    # plt.plot(g_disord_list, 1./res_var_reg_list[lbda_list.index(lbda),:], '.', color = 'black')
    plt.plot(g_disord_sorted_list/T, 1./res_var_reg_sorted_list[lbda_list.index(lbda),:], '-o', linewidth = .5, color = 'black')
plt.plot(g_b_fine_list/T, np.ones_like(g_b_fine_list)/res_rnd_guess, '--', color = 'black', label = 'random guess')
plt.plot(g_b_fine_list/T, FisherInfo_theo_total_list, color = 'blue', label = r'${\cal I}_{n}\left(\xi\right)$')
plt.xlabel(r'disorder strength $g$')
plt.ylabel(r'$\frac{1}{\mathrm{residual}}$', rotation = 0, labelpad=10)
plt.legend(labelcolor='linecolor')

figurename = ('Inverse_residual_T=' + str(T) + '_A_two_point=' +  str(A_two_point) + '_w_two_point=' + str(w_two_point)
              + '_A_one_point=' + str(A_one_point) + '_w_one_point=' + str(w_one_point) + '_g=' + str(g_disord)
              + '_f=' + str(mean_act) + '_connect_shape=' + str(connect_shape) + '_N_pf=' + str(N_pf) + '_N_thermal=' + str(N_thermal)
                + '_N_time=' + str(N_time) + '_N_MC_runs=' + str(N_MC_runs) + '_N_xi=' + str(N_xi) 
                + '_xi_min=' + str(xi_min_load) + '_xi_max=' + str(xi_max_load)) 

plt.savefig(figurepath + figurename + '.pdf', bbox_inches = 'tight')
plt.savefig(figurepath + figurename + '.eps', bbox_inches = 'tight')
plt.savefig(figurepath + figurename + '.jpg', bbox_inches = 'tight')

# %%
