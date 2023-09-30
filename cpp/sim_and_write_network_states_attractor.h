#ifndef SIM_AND_WRITE_AV_FISHER_INFO
#define SIM_AND_WRITE_AV_FISHER_INFO



#include "Vector_Definitions.h"
#include "population_neuronal_functions_attractor.h"
#include <fstream>
#include <sstream>

using std::cout;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::uniform_real_distribution;

unsigned int read_seed_disorder(string datapath, int seed_number);

void sim_and_write_network_states(double A_one_point, double w_one_point, double A_two_point, double w_two_point, double f, double g, double xi, double T, int shape,
int N_pf, int N_thermal, int N_measure, int N_MC_runs, string datapath, string time_suffix, unsigned int seed_disorder);

#endif
