#ifndef POPULATION_NEURONAL_FUNCTIONS_H
#define POPULATION_NEURONAL_FUNCTIONS_H



#include "PI_Definition.h"
#include "Vector_Definitions.h"
#include "basic_neuronal_functions.h"
#include <iostream>
#include <chrono>
#include <Eigen/Dense>

using std:: cout;
using std:: endl;
using namespace std::chrono;
using namespace Eigen;

void Q_unn_canon_ens(double &Q_final, double &E_final, VectorXd V_list, Vector_double r_list, double xi, Vector_double psi_list, Matrix_double J, double A_one_point, double w_one_point, double T);

void Energy_attr(double &E_final, Vector_int V_list, Vector_double r_list, double xi, Matrix_double J, double A_one_point, double w_one_point);

double Delta_Energy_attr(VectorXd &V_list_eig, Vector_double &psi_list, MatrixXd &J_eig, int i_switch_1, int i_switch_2, int Ising_like);

double Delta_Energy_attr_nonconserved(VectorXd &V_list_eig, Vector_double &psi_list, MatrixXd &J_eig, int i_switch, int Ising_like);

//Functions for Monte-Carlo computations:
void Metro_update_attractor(VectorXd &sites_eig, double &Q_unnorm_sites, double &E_sites, Vector_double &psi_list, MatrixXd &J_eig, double xi, 
double A_one_point, double w_one_point, double T, int Ising_like);

void Metro_update_attractor_nonconserved(VectorXd &sites_eig, double &Q_unnorm_sites, double &E_sites, Vector_double &psi_list, MatrixXd &J_eig, double xi, 
double A_one_point, double w_one_point, double T, int Ising_like);

void compute_time_dep_population_activity_MC(Matrix_double &pop_act, Vector_double r_list, Matrix_double J, 
double xi, double A_one_point, double w_one_point, int total_act, double T, int N_MC_runs, int N_thermal, int N_measure);

void compute_Activity_Fisher_Info_MC(double &S_mean_final, double &S_raw_var_final, Vector_double &mean_act, Matrix_double &corr_raw, Matrix_double &corr, Vector_double r_list, 
Matrix_double J, double xi, double A_one_point, double w_one_point, int total_act, double T, int N_MC_runs, int N_thermal, int N_measure, int Ising_like);

void compute_Activity_Fisher_Info_MC_nonconserved(double &S_mean_final, double &S_raw_var_final, Vector_double &mean_act, Matrix_double &corr_raw, Matrix_double &corr, Vector_double r_list, 
Matrix_double J, double xi, double A_one_point, double w_one_point, int total_act, double T, int N_MC_runs, int N_thermal, int N_measure, int Ising_like);


#endif


