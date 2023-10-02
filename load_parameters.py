#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 18:51:29 2022

@author: tobias
"""

import numpy as np
import os

def load_parameters(dp, determine_if_conserved):
    files = os.listdir(dp)
    
    
    right_file = [f for f in files if ( f.find('Parameters_attractor')!=-1 )]
    
    print(right_file[0])
    
    f = open(dp + right_file[0], "r")
    
    liste = []

    for line in f:
        liste.append(line.split())
    
    f.close()
    
    print(liste)
    
    if(determine_if_conserved == False):
        labels_expected = ["T", "A_two_point_b", "w_two_point", "A_one_point_b", 
                           "w_one_point", "g_b", "xi", "f", "connect_shape", 
                           "N_pf", "N_thermal", "N_measure", "N_MC_runs", "N_J", "N_r",
                           "datapath"]
    else:
        labels_expected = ["T", "A_two_point_b", "w_two_point", "A_one_point_b", 
                           "w_one_point", "g_b", "xi", "f", "connect_shape", "conserved",
                           "N_pf", "N_thermal", "N_measure", "N_MC_runs", "N_J", "N_r",
                           "datapath"]
    parameter_list = np.zeros(len(labels_expected)-1) #Exlude last element because it indicates the datapath
    
    for ii in range(0, len(liste)-1): #Exclude last element because it indicates the datapath
        if(liste[ii][0] != labels_expected[ii]):
            print("Parameters in parameter file not in expected order.")
        else:
            parameter_list[ii] = float(liste[ii][1])
            if(parameter_list[ii] == np.round(parameter_list[ii])):
                parameter_list[ii] = int(parameter_list[ii])  #Needed for correctly addressing the names of the files with
                                                              #data generated numerically.
            
    print("parameter_list=", parameter_list)
    
    if(determine_if_conserved == False):
        [T, A_two_point_b, w_two_point, A_one_point_b, 
         w_one_point, g_b, xi, f, connect_shape, 
         N_pf, N_thermal, N_measure, N_MC_runs, N_J, N_r] = parameter_list
        
        return (T, A_two_point_b, w_two_point, A_one_point_b, 
         w_one_point, g_b, xi, f, int(connect_shape),
         int(N_pf), int(N_thermal), int(N_measure), int(N_MC_runs), int(N_J), int(N_r))
        
    else:
        [T, A_two_point_b, w_two_point, A_one_point_b, 
         w_one_point, g_b, xi, f, connect_shape, conserved,
         N_pf, N_thermal, N_measure, N_MC_runs, N_J, N_r] = parameter_list
    
    
        return (T, A_two_point_b, w_two_point, A_one_point_b, 
         w_one_point, g_b, xi, f, int(connect_shape), int(conserved),
         int(N_pf), int(N_thermal), int(N_measure), int(N_MC_runs), int(N_J), int(N_r))