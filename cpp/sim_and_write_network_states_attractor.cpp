#include "sim_and_write_network_states_attractor.h"


void sim_and_write_network_states(double A_one_point, double w_one_point, double A_two_point, double w_two_point, double f, double g, double xi, double T, int shape,
int N_pf, int N_thermal, int N_measure, int N_MC_runs, string datapath, string time_suffix, unsigned int seed_disorder){
	//The string datapath is assumed to already include the link to the subfolder named "time_suffix"
  	
  	int Total_act = f * N_pf; //Note: this is in the spirit of the bit-like convention
  	
  	// File vector with centers of place fields:
	Vector_double r_list(N_pf);
	for(int ii=0; ii<N_pf; ii++){
		r_list[ii] = (1.*ii)/N_pf;
	}	
	
	
	Matrix_double J_total(N_pf); //Consist of the sum of the distance-dependent and the disordered part, is therefore newly composed for every disorder realization
	for(int ii=0; ii < N_pf; ii++){
		J_total[ii].resize(N_pf);
	}
	Matrix_double J_dist(N_pf); //Code does not compile if object is only created in case that there is really a spatial connectivity drawn 	

	generate_J_dist(J_dist, r_list, A_two_point, w_two_point, shape);
	
	//Pick the seed_number's random number from the list stored in timesuffix/seeds.txt
	
	//Tensor_double J_disorder_collection(N_J);
	Matrix_double J_disorder(N_pf);
	generate_J(J_disorder, g, seed_disorder);
	
	
	for(int jj=0; jj < N_pf; jj++){
		for(int kk=0; kk < N_pf; kk++){
				//J_total[jj][kk] = J_disorder_collection[ii][jj][kk] + J_dist[jj][kk];
			J_total[jj][kk] = J_disorder[jj][kk] + J_dist[jj][kk];
		}
	}

    Matrix_double pop_act_MC(N_pf);
	
	compute_time_dep_population_activity_MC(pop_act_MC, r_list, J_total, xi, A_one_point, w_one_point, 
									        Total_act, T, N_MC_runs, N_thermal, N_measure);
	
	
	  //////////////////////////////////////////////////
	 /// Define strings for later use in file names ///
  	//////////////////////////////////////////////////

  	ofstream Ausgabe_annealed;
  	string dateiname_string_annealed;
  	const char *dateiname_annealed;
  	
  	string string_N_pf;
  	ostringstream  umwandeln_N_pf;
  	umwandeln_N_pf << N_pf;
  	string_N_pf = umwandeln_N_pf.str();
  	
  	ostringstream umwandeln_T;
  	umwandeln_T << T;
  	string string_T = umwandeln_T.str();
  	
  	ostringstream umwandeln_A_one_point_b;
  	umwandeln_A_one_point_b << A_one_point; //A_one_point is assumed to be already scaled according to the bit-like convention
  	string string_A_one_point_b = umwandeln_A_one_point_b.str();
  	
	ostringstream umwandeln_w_one_point;
  	umwandeln_w_one_point << w_one_point;
  	string string_w_one_point = umwandeln_w_one_point.str();  	

	ostringstream umwandeln_w_two_point;
  	umwandeln_w_two_point << w_two_point;  //For file names resort to bitlike convention
  	string string_w_two_point = umwandeln_w_two_point.str();
  	
  	ostringstream umwandeln_A_two_point_b;
  	umwandeln_A_two_point_b << A_two_point;
  	string string_A_two_point_b = umwandeln_A_two_point_b.str();
  	
  	ostringstream umwandeln_g_b;
	umwandeln_g_b << g;
  	string string_g_b = umwandeln_g_b.str();

	ostringstream umwandeln_xi;
	umwandeln_xi << xi;
	string string_xi = umwandeln_xi.str();
	
	ostringstream umwandeln_f;
	umwandeln_f << f;
	string string_f = umwandeln_f.str();
	
	//The string datapath is assumed to already include the link to the subfolder named "time_suffix"
	
	string file_names_suffix = "";
		  	
	file_names_suffix += "T=" + string_T + "_f=" + string_f + "_A_one_point_b=" + string_A_one_point_b + "_w_one_point=" + string_w_one_point; 
	file_names_suffix += "_A_two_point_b=" + string_A_two_point_b + "_w_two_point=" + string_w_two_point; 
	file_names_suffix += "_g_b=" + string_g_b + "_xi="+ string_xi + "_N_pf=" + string_N_pf + ".txt";
		
  	  /////////////////////////////////////////////////
	 /// Create files for moments spin sum via MC ////
	/////////////////////////////////////////////////
  	
  	/// Write mean neural activity (including error from disorder average)
  	
  	ofstream Ausgabe_activity_MC;
  	string dateiname_string_activity_MC;
  	const char *dateiname_activity_MC;
  	
  	dateiname_string_activity_MC = datapath + "MeanAct_MC_for_rectangle_and_" + file_names_suffix;
  	
  	dateiname_activity_MC = dateiname_string_activity_MC.c_str();
  	Ausgabe_activity_MC.open(dateiname_activity_MC);
  	
  	for(int i=0; i<N_pf-1; i++){
        for(int j=0; j<N_measure; j++){
            Ausgabe_activity_MC << pop_act_MC[i][j] << " ";
        }
        Ausgabe_activity_MC << endl;
  	}

    for(int j=0; j<N_measure; j++){
        Ausgabe_activity_MC << pop_act_MC[N_pf - 1][j] << " ";
    } 
  	
  	
  	Ausgabe_activity_MC.close();
  	cout << "Ausgabe_activity_MC wieder geschlossen." << endl;

}
