# spacentation_code
Code accompanying the publication of Kühn and Monasson 2023

Figure 2

Run "create_fig2.py", which uses calc_CrossCorr_fixedMean_disorder.py and calc_FisherInfo.py

#############
#############

Figures 3 and 6

#############################################################################################################################################################################################
### Monte Carlo #############################################################################################################################################################################
#############################################################################################################################################################################################
The program implementing the Monte Carlo (MC) simulation validating our mean-field results for the Fisher information is compiled as follows:

g++ -I ../../../../../eigen/ -L/usr/include -O Fisher_information_attractor.cpp -o Fisher_information_attractor basic_neuronal_functions.cpp population_neuronal_functions_attractor.cpp sim_and_write_av_Fisher_info_attractor.cpp -lm -lgsl -lgslcblas

To run it on linux, use the bash script

run_FisherInfo_numerics_disordered.sh (make it executable with chmod +x run_FisherInfo_numerics_disordered.sh, then type ./run_FisherInfo_numerics_disordered). 
#############################################################################################################################################################################################


When using given data, use the time_suffix_list 
['2022_07_11_20_31_46', '2022_07_11_20_32_19', '2022_07_11_20_32_50', 
'2022_07_11_20_33_22','2022_07_11_20_33_52', '2022_07_11_20_34_22']
to generate the respective panel a
and the time_suffix_list
time_suffix_list = ['2022_07_21_17_08_49', '2022_07_21_17_09_46', '2022_07_21_17_10_36', 
                    '2022_07_21_17_11_12', '2022_07_21_17_11_57', '2022_07_21_17_12_34']
to generate the respective panel b.

Script needs the python files load_parameters.py, load_and_plot_activities.py, attractor_with_input.py, solve_sp_equations_space_disorder.py, calc_CrossCorr_fixedMean_disorder.py and calc_FisherInfo.py

###############
###############

Figure 4a
Run "create_fig4a.py", which uses calc_CrossCorr_fixedMean_disorder.py and calc_FisherInfo.py


Figure 4b

Let Recurrent_feedforward_const_FisherInfo.py run with three different parameters: FisherInfo_wish = 0.01, 0.015 and 0.02. This script uses compute_FisherInfo.py
Then run "create_fig4b"

###############
###############

Figure 5

#############################################################################################################################################################################################
### Monte Carlo #############################################################################################################################################################################
#############################################################################################################################################################################################
In order to generate the MC data to which the linear readout is fit, we use c++ - programs that are compiled by g++ (under linux) as follows:

Compile program that generates time series of full population activity

g++ -I ../../../../../eigen/ -O Generate_raw_states_attractor.cpp -o Generate_raw_states_attractor basic_neuronal_functions.cpp population_neuronal_functions_attractor.cpp sim_and_write_network_states_attractor.cpp -lm -lgsl -lgslcblas

g++  -O generate_xis_and_folder.cpp -o generate_xis_and_folder

To run both programs in linux, use the shell script run_time_dep_pop_act_numerics_disordered.sh (make it executable with chmod +x run_time_dep_pop_act_numerics_disordered.sh, then type ./run_time_dep_pop_act_numerics_disordered). 
##############################################################################################################################################################################################

Then run "create_fig5.py" to generate the figure. Script needs load_full_sim_data.py, linear_regression.py and bias_lin_est.py
