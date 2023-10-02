#include "sim_and_write_av_Fisher_info_attractor.h"

unsigned int read_seed_disorder(string datapath, int seed_number){
	
	string file_seedlist = "seeds.txt";
	string dateiname_string_seed_list;
	dateiname_string_seed_list = datapath + file_seedlist;
	
	cout << "Dateiname der seed-Zahlen beim Lesen:" << dateiname_string_seed_list  << endl; 
	
	const char *dateiname_seed_list;
	dateiname_seed_list = dateiname_string_seed_list.c_str();
	
	std::ifstream read_seed_list(dateiname_seed_list);	
	unsigned int rnd_seed; //Important that it is an unsigned int! Most of the random numbers are too large for int.
	
	int count_lines = 0;
	if(read_seed_list.is_open()){
		while(read_seed_list >> rnd_seed && count_lines < seed_number){ // seed_number picks out the seed_number+1'th random number
			
			count_lines++;	
		}
		read_seed_list.close();
		cout << "seed_number=" << seed_number << " should equal "<< count_lines << "'th random number:" << rnd_seed << endl;		
	}
	else{
		cout << "Unable to open file with seed list." << endl;
	}
	
	cout << "In read_seed_disorder, rnd_seed finally is " << rnd_seed << endl;
	return rnd_seed;
}	

void sim_and_write_av_Fisher_information(double A_one_point, double w_one_point, double A_two_point, double w_two_point, double f, double g, double xi, double T, int shape, int N_pf, 
int N_thermal, int N_measure, int N_MC_runs, int Ising_like, int conserved, string datapath, string time_suffix, int seed_number){
	//The string datapath is assumed to already include the link to the subfolder named "time_suffix"
  	
  	
  	int Total_act = f * N_pf; //Note: this is in the spirit of the bit-like convention
  	
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
	unsigned int seed_disorder = read_seed_disorder(datapath, seed_number); //Important that it is an unsigned int! Most of the random numbers are too large for int!!!
	
	//Tensor_double J_disorder_collection(N_J);
	Matrix_double J_disorder(N_pf);
	generate_J(J_disorder, g, seed_disorder);
		
	
	double S_mean_MC;
	double S_raw_var_MC;
	double S_var_MC;
	
	Vector_double pop_act_MC(N_pf);
	
	//Second-order statistics:	
	Matrix_double corr_raw_MC(N_pf);
	Matrix_double corr_MC(N_pf);
	
	for(int jj=0; jj < N_pf; jj++){
		for(int kk=0; kk < N_pf; kk++){
			J_total[jj][kk] = J_disorder[jj][kk] + J_dist[jj][kk];
		}
	}
	
	if(conserved == 1){  
		compute_Activity_Fisher_Info_MC(S_mean_MC, S_raw_var_MC, pop_act_MC, corr_raw_MC, corr_MC, r_list, J_total, xi, A_one_point, w_one_point, 
											Total_act, T, N_MC_runs, N_thermal, N_measure, Ising_like);
	}
	else{
		compute_Activity_Fisher_Info_MC_nonconserved(S_mean_MC, S_raw_var_MC, pop_act_MC, corr_raw_MC, corr_MC, r_list, J_total, xi, A_one_point, w_one_point, 
											Total_act, T, N_MC_runs, N_thermal, N_measure, Ising_like);
	}
	S_var_MC = S_raw_var_MC - S_mean_MC * S_mean_MC;
		

	  //////////////////////////////////////////////////
	 /// Define strings for later use in file names ///
  	//////////////////////////////////////////////////
  	
  	string string_N_pf;
  	ostringstream  umwandeln_N_pf;
  	umwandeln_N_pf << N_pf;
  	string_N_pf = umwandeln_N_pf.str();
  	
  	ostringstream umwandeln_T;
  	umwandeln_T << T;
  	string string_T = umwandeln_T.str();
  	
  	ostringstream umwandeln_A_one_point_b;
  	if(Ising_like == 1){
		umwandeln_A_one_point_b << 2.*A_one_point; //For file names resort to bitlike convention
  	}
  	else if(Ising_like == 0){
		umwandeln_A_one_point_b << A_one_point; //A_one_point is assumed to be already scaled according to the bit-like convention
	}
  	string string_A_one_point_b = umwandeln_A_one_point_b.str();
  	
	ostringstream umwandeln_w_one_point;
  	umwandeln_w_one_point << w_one_point;
  	string string_w_one_point = umwandeln_w_one_point.str();  	

	ostringstream umwandeln_w_two_point;
  	umwandeln_w_two_point << w_two_point;  //For file names resort to bitlike convention
  	string string_w_two_point = umwandeln_w_two_point.str();
  	
  	ostringstream umwandeln_A_two_point_b;
  	if(Ising_like == 1){
		umwandeln_A_two_point_b << 4.*A_two_point;
  	}
  	else if(Ising_like == 0){
		umwandeln_A_two_point_b << A_two_point;
	}
  	string string_A_two_point_b = umwandeln_A_two_point_b.str();
  	
  	ostringstream umwandeln_g_b;
	if(Ising_like == 1){
		umwandeln_g_b << 4.*g; //Account here for the fact that the C++-code is written for Ising variables, whereas in the theory mostly bitlike variables are assumed
  	}
  	else if(Ising_like == 0){
		umwandeln_g_b << g;
	}
  	string string_g_b = umwandeln_g_b.str();

	ostringstream umwandeln_xi;
	umwandeln_xi << xi;
	string string_xi = umwandeln_xi.str();
	
	ostringstream umwandeln_f;
	umwandeln_f << f;
	string string_f = umwandeln_f.str();
	
	ostringstream umwandeln_seed_number;
	umwandeln_seed_number << seed_number;
	string string_seed_number = umwandeln_seed_number.str();
	
	//The string datapath is assumed to already include the link to the subfolder named "time_suffix"
	
	string file_names_suffix = "";
	
	if(conserved == 0){
		file_names_suffix += "nonconserved_meanact_"; 
	}
	  	
	file_names_suffix += "T=" + string_T + "_f=" + string_f + "_A_one_point_b=" + string_A_one_point_b + "_w_one_point=" + string_w_one_point; 
	file_names_suffix += "_A_two_point_b=" + string_A_two_point_b + "_w_two_point=" + string_w_two_point; 
	file_names_suffix += "_g_b=" + string_g_b + "_xi="+ string_xi + "_N_pf=" + string_N_pf; 
	file_names_suffix += "_seed_number=" + string_seed_number + ".txt";
	
	
  	  /////////////////////////////////////////////////
	 /// Create files for moments spin sum via MC ////
	/////////////////////////////////////////////////
  	ofstream Ausgabe_moments_MC;
  	string dateiname_string_moments_MC;
  	const char *dateiname_moments_MC;
  	
  	if(shape == 1){
  		dateiname_string_moments_MC = datapath + "Fisher_info_MC_for_rectangle_and_" + file_names_suffix;
  	}else if(shape == 2){
		dateiname_string_moments_MC = datapath + "Fisher_info_MC_ferro_and_" + file_names_suffix;
	}
  	
  	dateiname_moments_MC = dateiname_string_moments_MC.c_str();
  	Ausgabe_moments_MC.open(dateiname_moments_MC);
  
  	Ausgabe_moments_MC << S_mean_MC/N_pf << " " <<  S_var_MC/N_pf;
  	
  	Ausgabe_moments_MC.close();
  	cout << "Ausgabe_moments_MC closed again." << endl;
  	
  	/// Write mean neural activity (including error from disorder average)
  	
  	ofstream Ausgabe_activity_MC;
  	string dateiname_string_activity_MC;
  	const char *dateiname_activity_MC;
  	
  	if(shape == 1){
  		dateiname_string_activity_MC = datapath + "MeanAct_MC_for_rectangle_and_" + file_names_suffix;
  	}else if(shape == 2){
		dateiname_string_activity_MC = datapath + "MeanAct_MC_ferro_and_" + file_names_suffix;
	}
  	
  	dateiname_activity_MC = dateiname_string_activity_MC.c_str();
  	Ausgabe_activity_MC.open(dateiname_activity_MC);
  	
  	for(int i=0; i<N_pf-1; i++){
  		Ausgabe_activity_MC << pop_act_MC[i] << endl;
  	}
  	Ausgabe_activity_MC << pop_act_MC[N_pf - 1];
  	
  	Ausgabe_activity_MC.close();
  	cout << "Ausgabe_activity_MC closed again." << endl;
  	
  	// Write pairwise correlations (including error from disorder average)
	ofstream Ausgabe_corr_MC;
	string dateiname_string_corr_MC;
	const char *dateiname_corr_MC;
	  	
	if (shape == 1){
		dateiname_string_corr_MC = datapath + "Correlation_connected_MC_for_rectangle_and_" + file_names_suffix;
	}else if(shape == 2){
		dateiname_string_corr_MC = datapath + "Correlation_connected_MC_ferro_and_" + file_names_suffix;
	}
	
	cout << "datapath: " << datapath << endl;
	cout << "dateiname_string_corr_MC: " << dateiname_string_corr_MC << endl;
		
	dateiname_corr_MC = dateiname_string_corr_MC.c_str();
	Ausgabe_corr_MC.open(dateiname_corr_MC);
		
	for (int i = 0; i < N_pf; i++){
		for (int j = 0; j < N_pf; j++){
			Ausgabe_corr_MC << corr_MC[i][j] << " ";
		}
		Ausgabe_corr_MC << endl;
	}
	
	Ausgabe_corr_MC.close();
	cout << "Ausgabe_corr_MC closed again." << endl;
	
	// Write raw pairwise correlations (including error from disorder average)
	ofstream Ausgabe_corr_raw_MC;
	string dateiname_string_corr_raw_MC;
	const char *dateiname_corr_raw_MC;
	  	
	if(shape == 1){
		dateiname_string_corr_MC = datapath + "Correlation_raw_MC_for_rectangle_and_" + file_names_suffix;
	}else if(shape == 2){
		dateiname_string_corr_MC = datapath + "Correlation_raw_MC_ferro_and_" + file_names_suffix;
	}
		
	dateiname_corr_raw_MC = dateiname_string_corr_raw_MC.c_str();
	Ausgabe_corr_raw_MC.open(dateiname_corr_raw_MC);
		
	for (int i = 0; i < N_pf; i++){
		for (int j = 0; j < N_pf; j++){
			Ausgabe_corr_raw_MC << corr_raw_MC[i][j] << " ";
		}
		Ausgabe_corr_raw_MC << endl;
	}
	
	Ausgabe_corr_raw_MC.close();
	cout << "Ausgabe_corr_raw_MC closed again." << endl; 			
}
