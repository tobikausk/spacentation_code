import numpy as np
import os

def load_xis(datapath):
    f = open(datapath + 'xis.txt', "r")
    liste = []

    for line in f:
        liste.append(line.split())
    f.close()

    xis = np.zeros(len(liste))
    N_xi = len(liste)
    for i in range(N_xi):
        xis[i] = liste[i][0]
    
    return xis

def load_parameters_without_xi(dp):
    files = os.listdir(dp)
    
    right_file = [f for f in files if ( f.find('Parameters_attractor_without_xi_')!=-1 )]
    print(right_file[0])
    
    f = open(dp + right_file[0], "r")
    
    liste = []

    for line in f:
        liste.append(line.split())
    
    f.close()
    
    print(liste)
    
    labels_expected = ["T", "A_two_point_b", "w_two_point", "A_one_point_b", 
                           "w_one_point", "g_b", "f", "connect_shape",
                           "N_pf", "N_thermal", "N_measure", "N_MC_runs", "N_xi", 
                            "xi_min",  "xi_max", "N_r",
                           "datapath"]
    parameter_list = np.zeros(len(labels_expected)-1) #Exclude last element because it indicates the datapath
    
    for ii in range(0, len(liste)-1): #Exclude last element because it indicates the datapath
        if(liste[ii][0] != labels_expected[ii]):
            print("Parameters in parameter file not in expected order.")
        else:
            parameter_list[ii] = float(liste[ii][1])
            if(parameter_list[ii] == np.round(parameter_list[ii])):
                parameter_list[ii] = int(parameter_list[ii])  #Needed for correctly addressing the names of the files with
                                                              #data generated numerically.
            
    print("parameter_list=", parameter_list)
    
    [T, A_two_point_b, w_two_point, A_one_point_b, 
         w_one_point, g_b, f, connect_shape,
         N_pf, N_thermal, N_measure, N_MC_runs, N_xi, x_min, x_max, N_r] = parameter_list
    
    
    return (T, A_two_point_b, w_two_point, A_one_point_b, 
         w_one_point, g_b, f, int(connect_shape),
         int(N_pf), int(N_thermal), int(N_measure), int(N_MC_runs), int(N_xi), x_min, x_max, int(N_r))

def generate_file_name(xi, T, A_two_point, w_two_point, A_one_point, 
                        w_one_point, g_disord, f, N_pf):
    
    for param in [xi, T, A_two_point, w_two_point, A_one_point, w_one_point, g_disord, f, N_pf]:
        if(param == np.round(param)):
            param = int(param)
            print('Integer parameter detected: ', param)
    if(T != np.round(T)):
        filename = 'MeanAct_MC_for_rectangle_and_T=' +str(T)
    else:
        filename = 'MeanAct_MC_for_rectangle_and_T=' +str(int(T))
    if(f != np.round(f)):
        filename += '_f=' + str(f)
    else:
        filename += '_f=' + str(int(f))
    if(A_one_point != np.round(A_one_point)):
        filename += '_A_one_point_b=' + str(A_one_point)
    else:
        filename += '_A_one_point_b=' + str(int(A_one_point))
    if(w_one_point != np.round(w_one_point)):
        filename += '_w_one_point=' + str(w_one_point) 
    else:
        filename += '_w_one_point=' + str(int(w_one_point))
    if(A_two_point != np.round(A_two_point)):    
        filename += '_A_two_point_b=' + str(A_two_point)
    else:
        filename += '_A_two_point_b=' + str(int(A_two_point))
    if(w_two_point != np.round(w_two_point)):
        filename += '_w_two_point=' + str(w_two_point)
    else:
        filename += '_w_two_point=' + str(int(w_two_point))
    if(g_disord != np.round(g_disord)):
        filename += '_g_b=' + str(g_disord)
        print('g_b not an integer.')
    else:
        filename += '_g_b=' + str(int(g_disord))
    if(xi != np.round(xi)):
        filename += '_xi=' + str(xi) 
    else:
        filename += '_xi=' + str(int(xi))
    filename += '_N_pf=' + str(N_pf) + '.txt'

    print('In the function: filename=', filename)

    return filename

def load_full_sim_data(datapath):
    """Load the data located in the (text) file datapath.
    It contains the (binary) activity of a network, where in every row
    the time series for a single neuron is stored.
    Returns a matrix of shape (n_neurons, n_time_steps)
    containing the activity of the network."""
    data = np.loadtxt(datapath)
    #Read data line by line and write into an array accordingly
    
    f = open(datapath, "r")
    liste = []

    for line in f:
        liste.append(line.split())
    f.close()
    
    N_neurons = len(liste)
    N_time_steps = len(liste[0])
    print('N_neurons=', N_neurons)
    print('N_time_steps=', N_time_steps)
    activity_list = np.zeros((N_neurons, N_time_steps), dtype=int)

    for i in range(N_neurons):
        for j in range(N_time_steps):
            activity_list[i, j] = liste[i][j]

    return activity_list
