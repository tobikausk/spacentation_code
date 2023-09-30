#ifndef BASIC_NEURONAL_FUNCTIONS_H
#define BASIC_NEURONAL_FUNCTIONS_H



#include "PI_Definition.h"
#include "Vector_Definitions.h"
#include <math.h>
#include <iostream>
#include <random>
#include <Eigen/Dense>

using std:: cout;
using std:: endl;
using std::random_device;
using std::mt19937;
using std::normal_distribution;
using std::uniform_real_distribution;
using std::string;

using namespace Eigen;

double H(double x);

double xrel(double x, double x0);

void connect_dist_Gauss(double &J_dist, double J0, double sigma, double r1, double r2);

void connect_dist_rect(double &J_dist, double J0, double w, double r1, double r2);

void one_point_source(double &psi, double x, double x0, double A, double w);

void one_point_source_prime(double &psi_prime, double x, double x0, double A, double w);

void tc(double &p, double x, double x0, double w, double tc_min, double Delta_tc);
	// Calculates the value of a tuning curve with shape exp^(-(x-x0)^2)
	// The tuning curve has its maximum at x0, its width is specified by w_cos, its minimal (maximal) value is tc_min (tc_min+Delta_tc)
	
void tc_prime(double &pp, double x, double x0, double w, double Delta_tc);

void one_point_source(double &psi, double x, double x0, double A, double w);

void generate_spin_state(Vector_int &spin_state, int k, int Ising_like);

void Susc(double &S, Vector_int V, Vector_double r, double xi, double w, double tc_min, double Delta_tc);

void Susc_attractor(double &S, VectorXd V, Vector_double r, double xi, double A_one_point, double w_one_point, double T);

void generate_J(Matrix_double &J, double g, unsigned int seed);

void connect_dist_Gauss(double &J_dist, double J0, double sigma, double r1, double r2);

void connect_dist_rect(double &J_dist, double J0, double w, double r1, double r2);

void generate_J_dist(Matrix_double &J_dist, Vector_double r_list, double J0, double sigma, int shape);

void generate_pf_centers(Vector_double &r_list);

double rho_sample(double xi);

#endif


