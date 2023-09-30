//#include "population_neuronal_functions.h"
#include "population_neuronal_functions_attractor.h"

void Q_unn_canon_ens(double &Q_final, double &E_final, VectorXd V_list, Vector_double r_list, double xi, Vector_double psi_list, Matrix_double J, double A_one_point, double w_one_point, double T){
	int N_pf = V_list.size();
	double psi;
	double potential;
	
	
	E_final = 0; 
	for(int i=0; i<N_pf; i++){
		one_point_source(psi, xi, r_list[i], A_one_point, w_one_point);
		E_final += - psi * V_list(i);
		for(int j=0; j<i; j++){
			E_final += -J[i][j] * V_list(i) * V_list(j);
		}
		for(int j=i+1; j<N_pf; j++){
			E_final += -J[i][j] * V_list(i) * V_list(j);
		}
	}
	Q_final = exp(-E_final/T);
}

void Energy_attr(double &E_final, Vector_int V_list, Vector_double r_list, double xi, Matrix_double J, double A_one_point, double w_one_point){
	int N_pf = V_list.size();
	double psi;
	double potential;
	E_final = 0; 
	for(int i=0; i<N_pf; i++){
		one_point_source(psi, xi, r_list[i], A_one_point, w_one_point);
		E_final += - psi * V_list[i];
		for(int j=0; j<i; j++){
			E_final += -J[i][j] * V_list[i] * V_list[j];
		}
		for(int j=i+1; j<N_pf; j++){
			E_final += -J[i][j] * V_list[i] * V_list[j];
		}
	}
}


double Delta_Energy_attr(VectorXd &V_list_eig, Vector_double &psi_list, MatrixXd &J_eig, int i_switch_1, int i_switch_2, int Ising_like){
		int N_pf = V_list_eig.size();
		
		double Delta_E;
		if(Ising_like == 1){
			Delta_E = 2.* (psi_list[i_switch_1] *  V_list_eig(i_switch_1) + psi_list[i_switch_2] *  V_list_eig(i_switch_2));
		}
		else if(Ising_like == 0){
			Delta_E = psi_list[i_switch_1] *  (2.* V_list_eig(i_switch_1) - 1.) + psi_list[i_switch_2] *  (2. * V_list_eig(i_switch_2) - 1.);
		}
		
		
		double Delta_E_for = Delta_E;
		for(int i=0; i<N_pf; i++){
			if(i!=i_switch_1 and i != i_switch_2){
				if(Ising_like == 1){
					Delta_E_for += 2. * J_eig(i,i_switch_1) * V_list_eig(i_switch_1) * V_list_eig(i);
					Delta_E_for += 2. * J_eig(i,i_switch_2) * V_list_eig(i_switch_2) * V_list_eig(i); 
				}
				else if(Ising_like == 0){
					Delta_E_for += J_eig(i,i_switch_1) * (2.* V_list_eig(i_switch_1) - 1.) * V_list_eig(i);
					Delta_E_for += J_eig(i,i_switch_2) * (2.* V_list_eig(i_switch_2) - 1.) * V_list_eig(i);
				}
			}
		}
		
		return Delta_E_for;
}

double Delta_Energy_attr_nonconserved(VectorXd &V_list_eig, Vector_double &psi_list, MatrixXd &J_eig, int i_switch, int Ising_like){
		int N_pf = V_list_eig.size();
		
		double Delta_E;
		if(Ising_like == 1){
			Delta_E = 2.* psi_list[i_switch] *  V_list_eig(i_switch);
		}
		else if(Ising_like == 0){
			Delta_E = psi_list[i_switch] *  (2.* V_list_eig(i_switch) - 1.);
		}

		double Delta_E_for = Delta_E;
		for(int i=0; i<N_pf; i++){
			if(i!=i_switch){
				if(Ising_like == 1){
					Delta_E_for += 2. * J_eig(i,i_switch) * V_list_eig(i_switch) * V_list_eig(i); 
				}
				else if(Ising_like == 0){
					Delta_E_for += J_eig(i,i_switch) * (2.* V_list_eig(i_switch) - 1.) * V_list_eig(i);
				}
			}
		}
		
		return Delta_E_for;
}

void Metro_update_attractor(VectorXd &sites_eig, double &Q_unnorm_sites, double &E_sites, Vector_double &psi_list, MatrixXd &J_eig, double xi, 
double A_one_point, double w_one_point, double T, int Ising_like){	
	
	random_device rd;
    	mt19937 mt(rd());
	int N_pf = sites_eig.size();

    	std::uniform_int_distribution<int> dist_uni_N_pf(0,N_pf-1); //Distribution  i n c l u d e s  limits.      
	std::uniform_real_distribution<double> dist_uniuni(0,1);

	double Delta_E_direct;
	
	int i_flip_1 = dist_uni_N_pf(mt); //Randomly choose the first spin to be flipped
	int i_flip_2 = dist_uni_N_pf(mt);
	while((Ising_like == 1 && sites_eig(i_flip_1)*sites_eig(i_flip_2)==1) || (Ising_like == 0 && (2*sites_eig(i_flip_1) -1) * (2*sites_eig(i_flip_2)-1) == 1 )){//Newly draw the second index until one hits a spin with
	//inverse direction in order to maintain total magnetization
		i_flip_2 = dist_uni_N_pf(mt);
	}
	
		
	double Gewichtung_direct;
	
	Delta_E_direct = Delta_Energy_attr(sites_eig, psi_list, J_eig, i_flip_1, i_flip_2, Ising_like);
	Gewichtung_direct =  exp(-Delta_E_direct/T);
	
	double Zufallszahl;
	if(Gewichtung_direct >=1){
		if(Ising_like == 1){
			sites_eig(i_flip_1) = -sites_eig(i_flip_1);
			sites_eig(i_flip_2) = -sites_eig(i_flip_2);
        }
        else if(Ising_like == 0){
			sites_eig(i_flip_1) = 1 - sites_eig(i_flip_1);
			sites_eig(i_flip_2) = 1 - sites_eig(i_flip_2);
		}
        
        
		
	E_sites = E_sites + Delta_E_direct;
	Q_unnorm_sites = Q_unnorm_sites * Gewichtung_direct;		
		
	}else{        
        Zufallszahl = dist_uniuni(mt);
        if(Gewichtung_direct > Zufallszahl){
			if(Ising_like == 1){
				sites_eig(i_flip_1) = -sites_eig(i_flip_1);
				sites_eig(i_flip_2) = -sites_eig(i_flip_2);
			}
			else if(Ising_like == 0){
				sites_eig(i_flip_1) = 1 - sites_eig(i_flip_1);
				sites_eig(i_flip_2) = 1 - sites_eig(i_flip_2);
			}
		    
			E_sites = E_sites + Delta_E_direct;
			Q_unnorm_sites = Q_unnorm_sites * Gewichtung_direct;
		}
	}
}

void Metro_update_attractor_nonconserved(VectorXd &sites_eig, double &Q_unnorm_sites, double &E_sites, Vector_double &psi_list, MatrixXd &J_eig, double xi, 
double A_one_point, double w_one_point, double T, int Ising_like){	
	
	random_device rd;
    	mt19937 mt(rd());
	int N_pf = sites_eig.size();

    	std::uniform_int_distribution<int> dist_uni_N_pf(0,N_pf-1); //Distribution  i n c l u d e s  limits.      
	std::uniform_real_distribution<double> dist_uniuni(0,1);

	double Delta_E_direct;
	
	int i_flip = dist_uni_N_pf(mt); //Bestimme den zu flippenden Spin
	
	double Gewichtung_direct;
	
	Delta_E_direct = Delta_Energy_attr_nonconserved(sites_eig, psi_list, J_eig, i_flip, Ising_like);
	Gewichtung_direct =  exp(-Delta_E_direct/T);
	
	double Zufallszahl;
	if(Gewichtung_direct >=1){
		if(Ising_like == 1){
			sites_eig(i_flip) = -sites_eig(i_flip);
        }
        else if(Ising_like == 0){
			sites_eig(i_flip) = 1 - sites_eig(i_flip);
		}
        
		E_sites = E_sites + Delta_E_direct;
		Q_unnorm_sites = Q_unnorm_sites * Gewichtung_direct;
		
		
	}else{        
        Zufallszahl = dist_uniuni(mt);
        if(Gewichtung_direct > Zufallszahl){
			if(Ising_like == 1){
				sites_eig(i_flip) = -sites_eig(i_flip);
			}
			else if(Ising_like == 0){
				sites_eig(i_flip) = 1 - sites_eig(i_flip);
			}
			
			E_sites = E_sites + Delta_E_direct;
			Q_unnorm_sites = Q_unnorm_sites * Gewichtung_direct;
		}
	}
	
	
}

void compute_time_dep_population_activity_MC(Matrix_double &pop_act, Vector_double r_list, Matrix_double J, 
double xi, double A_one_point, double w_one_point, int total_act, double T, int N_MC_runs, int N_thermal, int N_measure){

	// Draws sample from the probability distribution of the network and returns the activity of the entire population during the measurement time.
	
	int N_pf = r_list.size();
    
	for(int i=0; i < N_pf; i++){
		pop_act[i].resize(N_measure); 
	}
	
	Vector_double psi_list(N_pf);
    double psi_current;
	
	for(int i=0; i < N_pf; i++){
		one_point_source(psi_current, xi, r_list[i], A_one_point, w_one_point);
		psi_list[i] = psi_current;
	}
    
	cout << "All fine before defining J_eig in outer MC function." << endl;
	    
	cout << "N_pf = " << N_pf << endl;
	cout << "J.size()=" << J.size() << endl;
	cout << "J[0].size()=" << J[0].size() << endl;
	    
    
    MatrixXd J_eig(N_pf, N_pf);
    for(int i=0; i<N_pf; i++){
		for (int j = 0; j < N_pf; j++){
			J_eig(i,j) = J[i][j] + J[j][i]; //Note that J_eig is symmetrized - in contrast to J!
		}
	}
	
	cout << "Definition of J_eig successfull" << endl;	
    
    double Q_unnorm_current;
    double E_current;
    
    double total_act_measured;
    
    for(int i_run = 0; i_run < N_MC_runs; i_run++){
		random_device rd;
		mt19937 generate_perm(rd());
	    
	    VectorXd V_current_eig(N_pf);
		
		Vector_int permute_indices(N_pf);
		for (int i = 0; i < N_pf; i++){
			permute_indices[i] = i;
		}
		
		std::shuffle(permute_indices.begin(), permute_indices.end(),generate_perm);
		
		for(int i=0; i < total_act; i++){ //Create first state with average activity 2f-1
			V_current_eig(permute_indices[i]) = 1;
		}
		for(int i=total_act; i < N_pf; i++){ 
			V_current_eig(permute_indices[i]) = 0;
		}
	    	
	    Q_unn_canon_ens(Q_unnorm_current, E_current, V_current_eig, r_list, xi, psi_list, J, A_one_point, w_one_point, T);

	    int Ising_like = 0;
	    for(int i=0; i<N_thermal; i++){
			Metro_update_attractor(V_current_eig, Q_unnorm_current, E_current, psi_list, J_eig, xi, A_one_point, w_one_point, T, Ising_like);	
	    }
	    /////////////////////////////
	    /// Start measurements //////
	    /////////////////////////////
	    
	    for(int i=0; i<N_measure; i++){
			Metro_update_attractor(V_current_eig, Q_unnorm_current, E_current, psi_list, J_eig, xi, A_one_point, w_one_point, T, Ising_like);
	    	for(int j=0; j<N_pf; j++){
				pop_act[j][i] = V_current_eig(j);
			}	
	    }	
	}
}


void compute_Activity_Fisher_Info_MC(double &S_mean_final, double &S_raw_var_final, Vector_double &mean_act, Matrix_double &corr_raw, Matrix_double &corr, 
Vector_double r_list, Matrix_double J, double xi, double A_one_point, double w_one_point, int total_act, double T, int N_MC_runs, int N_thermal, int N_measure, int Ising_like){

	// Computes the Fisher information (more precisely, the first two moments of S) and the mean activity of every neuron according to the tuning curve and the 
	// part of the connectivity depending on space, by Monte-Carlo sampling.
	//TODO first: Adjust list of arguments, adapt to computation of (unnormalized) probability in "Boltzmann style".
	
	int N_pf = r_list.size();
    Vector_double psi_list(N_pf);
    double psi_current;
	for(int i=0; i < N_pf; i++){
		one_point_source(psi_current, xi, r_list[i], A_one_point, w_one_point);
		psi_list[i] = psi_current;
	}
    
    cout << "All fine before defining J_eig in outer MC function." << endl;
    
    cout << "N_pf = " << N_pf << endl;
    cout << "J.size()=" << J.size() << endl;
    cout << "J[0].size()=" << J[0].size() << endl;
    
    MatrixXd J_eig(N_pf, N_pf);
    for(int i=0; i<N_pf; i++){
		for (int j = 0; j < N_pf; j++){
			J_eig(i,j) = J[i][j] + J[j][i]; //Note that J_eig is symmetrized - in contrast to J!
			//cout << J_eig(i,j) << endl;
			//cout << J[i][j] + J[j][i] << endl;
		}
	}
	
	cout << "Definition of J_eig successfull" << endl;	
    
    
    double sum_S = 0.;
    double sum_S_av = 0.;
    double sum_S_squared = 0.;
    double sum_S_squared_av = 0.;
    	
    Vector_double mean_act_trial(N_pf); //Helping variable to record activity for every MC run.
    Matrix_double corr_raw_trial(N_pf);
    for(int jj=0; jj<N_pf; jj++){    //Is determined by adding up the contributions from mean_act_trial (and dividing by the number of runs).
		mean_act[jj] = 0;
		corr_raw[jj].resize(N_pf);
		for(int kk=0; kk<N_pf; kk++){
			corr_raw[jj][kk] = 0;
		}
	}
    	
    double Q_unnorm_current;
    double E_current;
    double Susc_current;
    
    double total_act_measured;
    
    for(int i_run = 0; i_run < N_MC_runs; i_run++){
		random_device rd;
		mt19937 generate_perm(rd());
	    
	    VectorXd V_current_eig(N_pf);
		
		Vector_int permute_indices(N_pf);
		for (int i = 0; i < N_pf; i++){
			permute_indices[i] = i;
		}
		
		std::shuffle(permute_indices.begin(), permute_indices.end(),generate_perm);
		
		for(int i=0; i < total_act; i++){ //Create first state with average activity 2f-1
			V_current_eig(permute_indices[i]) = 1;
		}
		for(int i=total_act; i < N_pf; i++){ 
			if(Ising_like == 1){
				V_current_eig(permute_indices[i]) = -1;
			}
			else if(Ising_like == 0){
				V_current_eig(permute_indices[i]) = 0;
			}
		}
	    	
	    Q_unn_canon_ens(Q_unnorm_current, E_current, V_current_eig, r_list, xi, psi_list, J, A_one_point, w_one_point, T);
	    
	    for(int i=0; i<N_thermal; i++){
			if(i % 1000 == 0){
				total_act_measured = 0;
				for (int jj = 0; jj < N_pf; jj++){
					total_act_measured += V_current_eig(jj);
				}
				total_act_measured = total_act_measured/N_pf;
	    		}
		
	    	Metro_update_attractor(V_current_eig, Q_unnorm_current, E_current, psi_list, J_eig, xi, A_one_point, w_one_point, T, Ising_like);
	    }
	    /////////////////////////////
	    /// Start measurements //////
	    /////////////////////////////
	    sum_S = 0;
	    sum_S_squared = 0;
	    	
	    for(int jj=0; jj<N_pf; jj++){
	    	mean_act_trial[jj] = 0;
	    	corr_raw_trial[jj].resize(N_pf);
			for(int kk=0; kk<N_pf; kk++){
				corr_raw_trial[jj][kk] = 0;
			}	
	    }
	    	
	    for(int i=0; i<N_measure; i++){
			Metro_update_attractor(V_current_eig, Q_unnorm_current, E_current, psi_list, J_eig, xi, A_one_point, w_one_point, T, Ising_like);
	    	Susc_attractor(Susc_current, V_current_eig, r_list, xi, A_one_point, w_one_point, T);
	    		
	    	sum_S += Susc_current;
			sum_S_squared += Susc_current*Susc_current;
			for(int jj=0; jj<N_pf; jj++){
				mean_act_trial[jj] += V_current_eig(jj);
				corr_raw_trial[jj].resize(N_pf);
				for(int kk=0; kk<N_pf; kk++){
					corr_raw_trial[jj][kk] += V_current_eig(jj) * V_current_eig(kk); //Alternatively: (\sum_j mean_act_trial_jj)**2
				}
			}	
	    }
	    	   
	    sum_S = sum_S/N_measure;
	    sum_S_squared = sum_S_squared/N_measure;
	    for(int jj=0; jj<N_pf; jj++){
	    	mean_act_trial[jj] = mean_act_trial[jj]/N_measure;
	    	mean_act[jj] += mean_act_trial[jj];
	    	for(int kk=0; kk<N_pf; kk++){
				corr_raw_trial[jj][kk] = corr_raw_trial[jj][kk]/N_measure;
				corr_raw[jj][kk] += corr_raw_trial[jj][kk];
			}	
	    }
	    sum_S_av += sum_S;
	    sum_S_squared_av += sum_S_squared;
	    	
    }
    S_mean_final = sum_S_av/N_MC_runs;
	S_raw_var_final = sum_S_squared_av/N_MC_runs;
	for(int jj=0; jj < N_pf; jj++){
		mean_act[jj] = mean_act[jj]/N_MC_runs;
		for(int kk=0; kk < N_pf; kk++){
			corr_raw[jj][kk] = corr_raw[jj][kk]/N_MC_runs;
		}	
	}
	for(int jj=0; jj < N_pf; jj++){
		corr[jj].resize(N_pf);
		for(int kk=0; kk < N_pf; kk++){
			corr[jj][kk] = corr_raw[jj][kk] - mean_act[jj] * mean_act[kk];
		}	
	}	
}

void compute_Activity_Fisher_Info_MC_nonconserved(double &S_mean_final, double &S_raw_var_final, Vector_double &mean_act, Matrix_double &corr_raw, Matrix_double &corr, 
Vector_double r_list, Matrix_double J, double xi, double A_one_point, double w_one_point, int total_act, double T, int N_MC_runs, int N_thermal, int N_measure, int Ising_like){

	// Computes the Fisher information (more precisely, the first two moments of S) and the mean activity of every neuron according to the tuning curve and the 
	// part of the connectivity depending on space, by Monte-Carlo sampling.
	
	
	int N_pf = r_list.size();
    Vector_double psi_list(N_pf);
    double psi_current;
	for(int i=0; i < N_pf; i++){
		one_point_source(psi_current, xi, r_list[i], A_one_point, w_one_point);
		psi_list[i] = psi_current;
	}
    
    cout << "All fine before defining J_eig in outer MC function." << endl;
    
    cout << "N_pf = " << N_pf << endl;
    cout << "J.size()=" << J.size() << endl;
    cout << "J[0].size()=" << J[0].size() << endl;
    
    
    MatrixXd J_eig(N_pf, N_pf);
    for(int i=0; i<N_pf; i++){
		for (int j = 0; j < N_pf; j++){
			J_eig(i,j) = J[i][j] + J[j][i]; //Note that J_eig is symmetrized - in contrast to J!
		}
	}
	
	cout << "Definition of J_eig successfull" << endl;	
    
    double sum_S = 0.;
    double sum_S_av = 0.;
    double sum_S_squared = 0.;
    double sum_S_squared_av = 0.;
    	
    Vector_double mean_act_trial(N_pf); //Helping variable to record activity for every MC run.
    Matrix_double corr_raw_trial(N_pf);
    for(int jj=0; jj<N_pf; jj++){    //Is determined by adding up the contributions from mean_act_trial (and dividing by the number of runs).
		mean_act[jj] = 0;
		corr_raw[jj].resize(N_pf);
		for(int kk=0; kk<N_pf; kk++){
			corr_raw[jj][kk] = 0;
		}
	}
    	
    double Q_unnorm_current;
    double E_current;
    double Susc_current;
    
    double total_act_measured;
    
    for(int i_run = 0; i_run < N_MC_runs; i_run++){
		random_device rd;
		mt19937 generate_perm(rd());
			
		VectorXd V_current_eig(N_pf);
			
		Vector_int permute_indices(N_pf);
		for (int i = 0; i < N_pf; i++){
			permute_indices[i] = i;
		}
		
		std::shuffle(permute_indices.begin(), permute_indices.end(),generate_perm);
			
		for(int i=0; i < total_act; i++){ //Create first state with average activity 2f-1
			V_current_eig(permute_indices[i]) = 1;
		}
		for(int i=total_act; i < N_pf; i++){ 
			if(Ising_like == 1){
				V_current_eig(permute_indices[i]) = -1;
			}
			else if(Ising_like == 0){
				V_current_eig(permute_indices[i]) = 0;
			}
		}
		
	    	
		Q_unn_canon_ens(Q_unnorm_current, E_current, V_current_eig, r_list, xi, psi_list, J, A_one_point, w_one_point, T);
			
		for(int i=0; i<N_thermal; i++){
			Metro_update_attractor_nonconserved(V_current_eig, Q_unnorm_current, E_current, psi_list, J_eig, xi, A_one_point, w_one_point, T, Ising_like);	
		}
		/////////////////////////////
		/// Start measurements //////
		/////////////////////////////
		sum_S = 0;
		sum_S_squared = 0;
	    	
		for(int jj=0; jj<N_pf; jj++){
			mean_act_trial[jj] = 0;
			corr_raw_trial[jj].resize(N_pf);
			for(int kk=0; kk<N_pf; kk++){
				corr_raw_trial[jj][kk] = 0;
			}	
		}
	    	
		for(int i=0; i<N_measure; i++){
			Metro_update_attractor_nonconserved(V_current_eig, Q_unnorm_current, E_current, psi_list, J_eig, xi, A_one_point, w_one_point, T, Ising_like);
			Susc_attractor(Susc_current, V_current_eig, r_list, xi, A_one_point, w_one_point, T);
						
			sum_S += Susc_current;
			sum_S_squared += Susc_current*Susc_current;
			for(int jj=0; jj<N_pf; jj++){
				mean_act_trial[jj] += V_current_eig(jj);
				corr_raw_trial[jj].resize(N_pf);
				for(int kk=0; kk<N_pf; kk++){
					corr_raw_trial[jj][kk] += V_current_eig(jj) * V_current_eig(kk); //Alternatively: (\sum_j mean_act_trial_jj)**2
				}
			}	
	    }
	    	   
		sum_S = sum_S/N_measure;
		sum_S_squared = sum_S_squared/N_measure;
		for(int jj=0; jj<N_pf; jj++){
			mean_act_trial[jj] = mean_act_trial[jj]/N_measure;
			mean_act[jj] += mean_act_trial[jj];
			for(int kk=0; kk<N_pf; kk++){
				corr_raw_trial[jj][kk] = corr_raw_trial[jj][kk]/N_measure;
				corr_raw[jj][kk] += corr_raw_trial[jj][kk];
			}	
		}
	    	
		sum_S_av += sum_S;
		sum_S_squared_av += sum_S_squared;
    }
    	
	S_mean_final = sum_S_av/N_MC_runs;
	S_raw_var_final = sum_S_squared_av/N_MC_runs;
	for(int jj=0; jj < N_pf; jj++){
		mean_act[jj] = mean_act[jj]/N_MC_runs;
		for(int kk=0; kk < N_pf; kk++){
			corr_raw[jj][kk] = corr_raw[jj][kk]/N_MC_runs;
		}	
	}
	for(int jj=0; jj < N_pf; jj++){
		corr[jj].resize(N_pf);
		for(int kk=0; kk < N_pf; kk++){
			corr[jj][kk] = corr_raw[jj][kk] - mean_act[jj] * mean_act[kk];
		}	
	}	
}

