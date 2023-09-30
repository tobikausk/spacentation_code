#include "basic_neuronal_functions.h"


double H(double x){
	if(x>0){
		return 1.;	
	}
	else{
		return 0.;
	}
}

double xrel(double x, double x0){
	//Computes the distance between the points x and x0 considering periodic boundary conditions.
	if(abs(x-x0) <= 0.5){
		return x - x0;
	}else if( x - x0 > 0.5){
		// x -> x-1
		return x- x0 -1.;
	}else if(x - x0 < -0.5){ 
		// x0 -> x0 - 1
		return x - x0 +1;
	}
	return 0;
}

void connect_dist_Gauss(double &J_dist, double J0, double sigma, double r1, double r2){
	double Delta_r = xrel(r1, r2);
	J_dist = J0 * exp(-Delta_r*Delta_r/(sigma * sigma));
}

void connect_dist_rect(double &J_dist, double J0, double w, double r1, double r2){
	double Delta_r = xrel(r1, r2);
	if(abs(Delta_r) <= .5*w){
		J_dist = J0;
	}else{
		J_dist = 0.;
	}
}

void one_point_source(double &psi, double x, double x0, double A, double w){
	double d = xrel(x,x0);
	psi = A * exp(- d *d /(2.*w*w));
}

void one_point_source_prime(double &psi_prime, double x, double x0, double A, double w){
	double d = xrel(x,x0);
	psi_prime = - (A * d/(w*w)) *  exp(- d * d /(2.*w*w));
}


void tc(double &p, double x, double x0, double w, double tc_min, double Delta_tc){
	// Calculates the value of a tuning curve with shape exp^(-(x-x0)^2)
	// The tuning curve has its maximum at x0, its width is specified by w_cos, its minimal (maximal) value is tc_min (tc_min+Delta_tc)
	if(tc_min + Delta_tc > 1.){
		throw "Maximum of tuing curve larger than one.";	
	}
	p = tc_min + Delta_tc * exp(-(x - x0)*(x - x0)/(w*w));
}

void tc_prime(double &pp, double x, double x0, double w, double Delta_tc){
	// Calculates the value of the derivative of a tuning curve with shape exp^(-(x-x0)^2)
	// The tuning curve has its maximum at x0, its width is specified by w_cos, its minimal (maximal) value is tc_min (tc_min+Delta_tc)
	pp = -2.*Delta_tc * (x - x0)/(w*w) * exp(-(x - x0)*(x - x0)/(w*w));
}

void generate_spin_state(Vector_int &spin_state, int k, int Ising_like){
	int N = spin_state.size();        
	for(int i=0; i < N; i++){
		if(Ising_like == 1){
			spin_state[i] = int(k/pow(2,i))%2 * 2 -1;
		}
		else{
			spin_state[i] = int(k/pow(2,i))%2;
		}
	}
}


void Susc(double &S, Vector_int V, Vector_double r, double xi, double w, double tc_min, double Delta_tc){
	int N_pf = r.size();
	
	double p;
	double pp;
	double S_sum = 0;
	for(int i=0; i<N_pf; i++){
		tc(p, xi, r[i], w, tc_min, Delta_tc);
		tc_prime(pp,xi, r[i], w, Delta_tc);
		S_sum += V[i]*pp/((1-V[i])/2. + V[i] * p); 
	}
	S = S_sum;
}


void Susc_attractor(double &S, VectorXd V, Vector_double r, double xi, double A_one_point, double w_one_point, double T){
	int N_pf = r.size();
	
	double psi;
	double psi_prime;
	double S_sum = 0;
	for(int i=0; i<N_pf; i++){
		one_point_source_prime(psi_prime, xi, r[i], A_one_point, w_one_point);

		S_sum += V(i)*psi_prime/T; 
	}
	S = S_sum;
}



void generate_J(Matrix_double &J, double g, unsigned int seed){// Important: seed unsigned int, not int!
	int N = J.size();
	
	cout << "Seed in generate_J:" << seed << endl;
	
	std::mt19937 gen{seed};
 
    	std::normal_distribution<> d{0.,g/sqrt(N)};
 	
    	
	for (int ii=0; ii<N; ii++){
		J[ii].resize(N);
		for(int jj=0; jj<N; jj++){
			J[ii][jj] = d(gen);
		}
	}
	    	
}

void generate_J_dist(Matrix_double &J_dist, Vector_double r_list, double J0, double sigma, int shape){
	int N = J_dist.size();
	
	double current_weight;
	
	
	for (int ii=0; ii<N; ii++){
		J_dist[ii].resize(N);
		for(int jj=0; jj<N; jj++){
			if(shape == 0){
				connect_dist_Gauss(current_weight, J0, sigma, r_list[ii], r_list[jj]);
			}else if(shape == 1){
				connect_dist_rect(current_weight, J0, sigma, r_list[ii], r_list[jj]);
			}
			else if(shape == 2){
				current_weight = -J0;
			}
			J_dist[ii][jj] = current_weight/(2.*N); 
		}
	}
	
	for(int ii=0; ii<N; ii++){
		J_dist[ii][ii] = 0; //Make sure that there no self-connections (autapses)
	}	
}

void generate_pf_centers(Vector_double &r_list){
	random_device rd;
	mt19937 mt(rd());
	uniform_real_distribution<double> dist_uniuni(0, 1.0);
	
	int N = r_list.size();
	  
	for(int i=0; i<N; i++){
		r_list[i] = dist_uniuni(mt);
	}
}

double rho_sample(double xi){
	double rho;
	rho = 1.;
	return rho;
}

