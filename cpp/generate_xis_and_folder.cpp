#include <fstream>
#include <iostream>
#include <ctime>

#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

int main(){
	// Create suffix containing current time and create folder bearing this as name
	// to store therein a copy of the parameter file, a text file with the seeds to generate the disordered part of the interaction
	// and the generated data
	
	/////////////////////////////////////////////////////////
	// Read parameters to be read from the parameter file ///
	/////////////////////////////////////////////////////////
	int N_xi;
	double xi_min;
	double xi_max;
	string datapath;
	
	string line;
	ifstream file("Parameters_attractor_without_xi.txt", ios_base::in);
			
	string parameter_names_list[15] = {"T", "A_two_point_b", "w_two_point", "A_one_point_b", "w_one_point", "g_b", 
		"f", "connect_shape", "N_pf", "N_thermal", "N_measure", "N_MC_runs", "N_xi", "N_r", "datapath"};
	
	string temp;
	string parameter_name;
	double parameter_value;
		
	int i=0;
	while(getline(file, line)){
		stringstream sts;
		sts << line;
		
		sts >> parameter_name;
		if(parameter_name != parameter_names_list[i]){
			cerr << "Parameter file is not formatted as expected." << endl;
			cout << "parameter_name: " << parameter_name << " parameter_names_list[i]: " << parameter_names_list[i] << endl;
		}
		if(parameter_name != "datapath"){	
			while(!sts.eof()){
				sts >> temp;
				if(stringstream(temp) >> parameter_value){
					if(parameter_name == "N_xi"){
						N_xi = parameter_value;
					}
					else if(parameter_name == "xi_min"){
						xi_min = parameter_value;
					}
					else if(parameter_name == "xi_max"){
						xi_max = parameter_value;
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
	/////////////////////////////////////////
	///// End reading from parameter file ///
	/////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////////
	//// Generate subfolder of data folder named by the current time ///
	////////////////////////////////////////////////////////////////////
	
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];
	char time_suffix[80];

	time (&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(time_suffix,sizeof(time_suffix),"%Y_%m_%d_%H_%M_%S",timeinfo);

	std::cout << time_suffix << std::endl;
	
	string time_suffix_string = "";
	time_suffix_string = string(time_suffix);
	
	cout << "time_suffix_string=" << time_suffix_string << endl;
	
	ofstream write_stime_stamp;
	write_stime_stamp.open("time_stamp.txt");
	write_stime_stamp << time_suffix;
	write_stime_stamp.close();
	
	char folder_loc[100];
	
	
	
	cout << "datapath=" << datapath << endl;
	
	strcpy(folder_loc, datapath.c_str());
	strcat(folder_loc, time_suffix);
	
	cout << "folder_loc = " << folder_loc << endl;
	
	if (mkdir(folder_loc, 0777) == -1){
		cerr << "Error :  " << strerror(errno) << endl;
	}
	else{
		cout << "Directory created" << endl;
	}
	
	////////////////////////////////////////////////////////////////////
	//// Finished creation of subfolder ////////////////////////////////
	////////////////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////////
	//// Copy parameter file into this subfolder, //////////////////////
	//// create list of seeds and store it in a text file int this /////
	//// subfolder as well. ////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////
	
	string params_name_goal = "/Parameters_attractor_without_xi_";
	string txt_ending = ".txt";
	
	//Copy parameter file:
	ifstream source("Parameters_attractor_without_xi.txt", std::ios::binary);
	ofstream goal(datapath + time_suffix + params_name_goal + time_suffix + txt_ending, std::ios::binary);
	
	goal << source.rdbuf();
	
	//Create and write xi-file:
	string folder_loc_string = "";
	folder_loc_string = string(folder_loc);
	
	
	string file_xilist = "/xis.txt";
	string dateiname_string_xi_list;
	dateiname_string_xi_list = folder_loc_string + file_xilist;
	
	cout << dateiname_string_xi_list  << endl; 
		
	random_device rd;
	ofstream write_xi_list;
    	mt19937 mt(rd());
	uniform_real_distribution<double> dist_uniuni(0, 1.0);

    // Generate N_xi uniformly distributed random numbers between 0 and 1
	write_xi_list.open(dateiname_string_xi_list.c_str());
	for(int i=0; i < N_xi; i++){
		write_xi_list << xi_min + (xi_max - xi_min) * dist_uniuni(mt) << endl;  
	}
	write_xi_list.close();
	
	return 0;
}
