# spacentation_code
Code accompanying the publication of Kühn and Monasson 2023.

The python and c++ scripts provided in this repository enables to generate the data shown in the figures of the present publication. Please use the conda environment provided by spacentation.yml to run the python scripts. For the cpp programs, only standard libraries are used, notably eigen. To run the cpp-files, linux bash-scripts are provided. For different operating systems, this solution would have to be adapted.

The corresponding scripts, ordered by figures are as follows:

**Figure 2**

Run "create_fig2.py", which uses calc_CrossCorr_fixedMean_disorder.py and calc_FisherInfo.py


**Figures 3 and 6**

To exactly reproduce these figures use the data from the Monte-Carlo simulations provided in the respective folder. 

Alternatively regenerate the Monte-Carlo data using the code provided in the cpp folder, as described below. By running the resulting cpp-scripts, there are folders created in the data folder named with the respective time stamps of the simulation, which have to be entered manually into the python script to compare them with the theory and plot the respective figures. Note that because the simulations are stochastic, the results will of course not look exactly the same as in the paper in this case.

*Monte-Carlo simulations* 

The program implementing the Monte Carlo (MC) simulation validating our mean-field results for the Fisher information is compiled as follows:

g++ -I *location_of_eigen-library* -L/usr/include -O Fisher_information_attractor.cpp -o Fisher_information_attractor basic_neuronal_functions.cpp population_neuronal_functions_attractor.cpp sim_and_write_av_Fisher_info_attractor.cpp -lm -lgsl -lgslcblas

To run it on linux, use the bash script

run_FisherInfo_numerics_disordered.sh (make it executable with chmod +x run_FisherInfo_numerics_disordered.sh, then type ./run_FisherInfo_numerics_disordered). 

*Pre-computed data*

When using given data, use the time_suffix_list

['2022_07_11_20_31_46', '2022_07_11_20_32_19', '2022_07_11_20_32_50', 
'2022_07_11_20_33_22','2022_07_11_20_33_52', '2022_07_11_20_34_22']

to generate the respective panel a
and the time_suffix_list

time_suffix_list = ['2022_07_21_17_08_49', '2022_07_21_17_09_46', '2022_07_21_17_10_36', 
                    '2022_07_21_17_11_12', '2022_07_21_17_11_57', '2022_07_21_17_12_34']

to generate the respective panel b.

*Comparison with theory*

Run create_fig_3_and_6.py to eventually create the figures.

Dependencies: load_parameters.py, load_and_plot_activities.py, attractor_with_input.py, solve_sp_equations_space_disorder.py, calc_CrossCorr_fixedMean_disorder.py and calc_FisherInfo.py


**Figure 4a**

Run "create_fig4a.py" (takes several hours on a standard laptop if calc_theory_newly == True). 

Dependencies: calc_CrossCorr_fixedMean_disorder.py and calc_FisherInfo.py


**Figure 4b**

Run Recurrent_feedforward_const_FisherInfo.py with three different parameters: FisherInfo_wish = 0.01, 0.015 and 0.02., adapting A_one_point_list, A_two_point_minus_list and A_two_point_plus_list as indicated by the commented lines (takes several hours each time on a standard laptop). Once the data is generated, this does not have to be done again, of course.

Dependency: compute_FisherInfo.py

Then run create_fig4b.py.


**Figure 5**

Here again, to reproduce the figures from the manuscript exactly, use the data provided in the data folder.

To regenerate Monte-Carlo simulations, use the code in the cpp folder.

*Monte-Carlo simulations*

In order to generate the MC data to which the linear readout is fit, we use c++ - programs that are compiled by g++ (under linux) as follows:

Program that generates time series of full population activity:

g++ -I *location_of_eigen-library* -O Generate_raw_states_attractor.cpp -o Generate_raw_states_attractor basic_neuronal_functions.cpp population_neuronal_functions_attractor.cpp sim_and_write_network_states_attractor.cpp -lm -lgsl -lgslcblas

Program that generates the samples of positions to be inferred later:

g++  -O generate_xis_and_folder.cpp -o generate_xis_and_folder

To run both programs in linux, use the shell script run_time_dep_pop_act_numerics_disordered.sh (make it executable with chmod +x run_time_dep_pop_act_numerics_disordered.sh, then type ./run_time_dep_pop_act_numerics_disordered). 

*Pre-computed data*

Used the time stamps as given in the python script create_fig5.py.

*Compare with theory*

Then run create_fig5.py to generate the figure. 

Dependencies: load_full_sim_data.py, linear_regression.py and bias_lin_est.py
