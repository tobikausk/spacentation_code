#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <vector>
#include<iostream>
#include<fstream>
#include <sstream>
#include<random>
#include <time.h>

#include "basic_neuronal_functions.h"
#include "population_neuronal_functions_attractor.h"
#include "sim_and_write_av_Fisher_info_attractor.h"


using namespace std;



int main(int argc, char** argv){
  	
	
	cout << "You have entered " << argc
         << " arguments:" << "\n";
  
    for (int i = 0; i < argc; i++)
        cout << argv[i] << "\n";
        
	int seed_number;
	
	seed_number = strtol(argv[1], nullptr, 0);
	
    cout << "The seed_number is " << seed_number << endl;
    
    
    ifstream read_time_stamp("time_stamp.txt");
	
	string time_suffix;
	
	if(read_time_stamp.is_open()){
		
		while(read_time_stamp >> time_suffix){ 	
		}
		read_time_stamp.close();
				
	}
	else{
		cout << "Unable to open file with time stamp." << endl;
	}
	
	cout << "Current time_suffix: " << time_suffix << endl;
	

	double T;
	double A_one_point_I;
	double w_one_point;
	double w_two_point;
	double A_two_point_I;
	double A_one_point_b; //Store the parameters also in the "bit-like" convention.
	double A_two_point_b;
	
	
	Vector_double f_list(0);
	
	double xi;
	double g_I;
	double g_b;
	
	int conserved; //If equal to 1, the total activity is conserved (two spin flips in opposing directions at every MC update), otherwise not (one spin flip at every MC update)
	
	int N_pf;
	int connect_shape; //0 means Gaussian connectivity, 1 means rectangular connectivity, 2 means constant connectivity (ferromagnet)
	
	int N_J;
	int N_r;
	
	//Parameters for Monte-Carlo simulations:
	int N_thermal;
	int N_measure;
	int N_MC_runs;
	
	string datapath;
	
	string line;
	ifstream file("Parameters_attractor.txt", ios_base::in);
	
		
	string parameter_names_list[17] = {"T", "A_two_point_b", "w_two_point", "A_one_point_b", "w_one_point", "g_b", "xi", 
		"f", "connect_shape", "conserved", "N_pf", "N_thermal", "N_measure", "N_MC_runs", "N_J", "N_r", "datapath"};
	
	string temp;
	string parameter_name;
	double parameter_value;
	
	// Variables set at run time:
	int Ising_like;
	Ising_like = 0;
	
	int i=0;
	while(getline(file, line)){
		stringstream sts;
		sts << line;
		
		sts >> parameter_name;
		if(parameter_name != parameter_names_list[i]){
			cerr << "Parameter file is not formatted as expected." << endl;
			cout << "parameter_name: " << parameter_name << " parameter_names_list[i]: " << parameter_names_list[i] << endl;
		}
		//int parameter_count = 0;
		if(parameter_name != "datapath"){	
			while(!sts.eof()){
				sts >> temp;
				if(stringstream(temp) >> parameter_value){
					//cout << "parameter_value=" << parameter_value << endl;
					if(parameter_name == "T"){
						T = parameter_value;
					}
					else if(parameter_name == "A_two_point_b"){
						if(Ising_like == 1){
							A_two_point_I = parameter_value/4.; //Adapt it to the Ising convention
						}
						A_two_point_b = parameter_value;
						
					}
					else if(parameter_name == "w_two_point"){
						w_two_point = parameter_value;
					}
					else if(parameter_name == "A_one_point_b"){
						if(Ising_like == 1){
							A_one_point_I = parameter_value/2.;
						}
						A_one_point_b = parameter_value;
					}
					else if(parameter_name == "w_one_point"){
						w_one_point = parameter_value;
					}
					else if(parameter_name == "g_b"){
						if(Ising_like == 1){
							g_I = parameter_value/4.;
						}
						g_b = parameter_value;
					}
					else if(parameter_name == "xi"){
						xi = parameter_value;
					}
					else if(parameter_name == "f"){
						f_list.push_back(parameter_value);
					}
					else if(parameter_name == "connect_shape"){
						connect_shape = parameter_value;	
					}
					else if(parameter_name == "conserved"){
						conserved = parameter_value;
					}
					else if(parameter_name == "N_pf"){
						N_pf = parameter_value;
					}
					else if(parameter_name == "N_J"){
						N_J = parameter_value;
					}
					else if(parameter_name == "N_r"){
						N_r = parameter_value;
					}
					else if(parameter_name == "N_thermal"){
						N_thermal = parameter_value;
					}
					else if(parameter_name == "N_measure"){
						N_measure = parameter_value;
					}
					else if(parameter_name == "N_MC_runs"){
						N_MC_runs = parameter_value;
					}
				}
			}
		}
		else{
			while(!sts.eof()){
				sts >> datapath;
			}
		}
		i++;		
	}
	cout << "T=" << T << " A_two_point_b=" << A_two_point_b << " w_two_point=" << w_two_point << endl;
	cout << "A_one_point_b=" << A_one_point_b << " w_one_point=" << w_one_point << " g_b=" << g_b << " xi=" << xi << " connect_shape=" << connect_shape << " conserved=" << conserved << endl;
	cout << "f=";
	for(int i=0; i < f_list.size(); i++){
		cout << f_list[i] << ", ";
	}
	cout << endl;
	cout << "N_pf=" << N_pf << " N_thermal=" << N_thermal << " N_measure=" << N_measure << " N_MC_runs=" << N_MC_runs << " N_J=" << N_J << " N_r=" << N_r << " datapath=" << datapath <<endl; 
	
  	
  	//Now really start doing something:
  	for(int ii=0; ii<f_list.size(); ii++){  	
		if(Ising_like == 1){
			sim_and_write_av_Fisher_information(A_one_point_I, w_one_point, A_two_point_I, w_two_point, f_list[ii], g_I, xi, T, 
			connect_shape, N_pf, N_thermal, N_measure, N_MC_runs, Ising_like, conserved, datapath + time_suffix + '/', time_suffix, seed_number);
		}
		else if(Ising_like == 0){
			sim_and_write_av_Fisher_information(A_one_point_b, w_one_point, A_two_point_b, w_two_point, f_list[ii], g_b, xi, T, 
			connect_shape, N_pf, N_thermal, N_measure, N_MC_runs, Ising_like, conserved, datapath + time_suffix + '/', time_suffix, seed_number);
		}
  	}
	return 0;
}
