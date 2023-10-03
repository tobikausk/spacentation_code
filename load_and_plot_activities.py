#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 18:51:29 2022

@author: tobias
"""

import numpy as np
import os
from matplotlib import pyplot as plt

def load_Moments_Susc_error(dp, space_form, method, T, freq, A_one_point_b, w_one_point, A_two_point_b, w_two_point, g_b, xi, N_pf, seed_no):
    files = os.listdir(dp)
    
    
    right_file = [f for f in files if ( #( f.find('T='+str(T)+'_')!=-1 or f.find('T='+str(T)+'.')!=-1)
                                                   f.find('T='+str(T)+'_')!=-1
                                                # and ( f.find('f='+str(freq)+'_')!=-1 or f.find('f='+str(freq)+'.')!=-1)
                                                and f.find('f='+str(freq)+'_')!=-1
                                                # and ( f.find('A_one_point_b='+str(A_one_point_b)+'_')!=-1 or f.find('A_one_point_b='+str(A_one_point_b)+'.')!=-1)
                                                and f.find('A_one_point_b='+str(A_one_point_b)+'_')!=-1
                                                # and ( f.find('w_one_point='+str(w_one_point)+'_')!=-1 or f.find('w_one_point='+str(w_one_point)+'.')!=-1)
                                                and f.find('w_one_point='+str(w_one_point)+'_')!=-1
                                                # and ( f.find('A_two_point_b='+str(A_two_point_b)+'_')!=-1 or f.find('A_two_point_b='+str(A_two_point_b)+'.')!=-1)
                                                and ( f.find('A_two_point_b='+str(A_two_point_b)+'_')!=-1)
                                                # and ( f.find('w_two_point='+str(w_two_point)+'_')!=-1 or f.find('w_two_point='+str(w_two_point)+'.')!=-1)
                                                and ( f.find('w_two_point='+str(w_two_point)+'_')!=-1)
                                                # and ( f.find('g_b='+str(g_b)+'_')!=-1 or f.find('g_b='+str(g_b)+'.')!=-1)  
                                                and f.find('g_b='+str(g_b)+'_')!=-1
                                                # and ( f.find('xi='+str(xi)+'_')!=-1 or f.find('xi='+str(xi)+'.')!=-1)
                                                and f.find('xi='+str(xi)+'_')!=-1
                                                # and  ( f.find('N_pf='+str(N_pf)+'_')!=-1 or f.find('N_pf='+str(N_pf)+'.')!=-1)
                                                and f.find('N_pf='+str(N_pf)+'_')!=-1
                                                and f.find('seed_number='+str(seed_no) + '.') !=-1
                                                and f.find(method)!=-1                                              
                                                and f.find(space_form)!=-1
                                                 
                                                and (f.find('MeanAct') !=-1) and os.stat(dp + f).st_size != 0)]
    
    print(right_file[0])
    
    f = open(dp + right_file[0], "r")
    
    liste = []

    for line in f:
        liste.append(float(line))
    
    f.close()
    
    return np.array(liste)


def load_MeanAct_error(dp, space_form, method, T, freq, A_one_point_b, w_one_point, A_two_point_b, w_two_point, g_b, xi, N_pf, seed_no):
    files = os.listdir(dp)
    
    right_file = [f for f in files if ( #( f.find('T='+str(T)+'_')!=-1 or f.find('T='+str(T)+'.')!=-1)
                                                      f.find('T='+str(T)+'_')!=-1
                                                # and ( f.find('f='+str(freq)+'_')!=-1 or f.find('f='+str(freq)+'.')!=-1)
                                                and  f.find('f='+str(freq)+'_')!=-1
                                                # and ( f.find('A_one_point_b='+str(A_one_point_b)+'_')!=-1 or f.find('A_one_point_b='+str(A_one_point_b)+'.')!=-1)
                                                and  f.find('A_one_point_b='+str(A_one_point_b)+'_')!=-1
                                                # and ( f.find('w_one_point='+str(w_one_point)+'_')!=-1 or f.find('w_one_point='+str(w_one_point)+'.')!=-1)
                                                and f.find('w_one_point='+str(w_one_point)+'_')!=-1
                                                # and ( f.find('A_two_point_b='+str(A_two_point_b)+'_')!=-1 or f.find('A_two_point_b='+str(A_two_point_b)+'.')!=-1)
                                                and f.find('A_two_point_b='+str(A_two_point_b)+'_')!=-1
                                                # and ( f.find('g_b='+str(g_b)+'_')!=-1 or f.find('g_b='+str(g_b)+'.')!=-1)
                                                and f.find('g_b='+str(g_b)+'_')!=-1
                                                # and ( f.find('w_two_point='+str(w_two_point)+'_')!=-1 or f.find('w_two_point='+str(w_two_point)+'.')!=-1)
                                                and f.find('w_two_point='+str(w_two_point)+'_')!=-1
                                                # and ( f.find('xi='+str(xi)+'_')!=-1 or f.find('xi='+str(xi)+'.')!=-1)
                                                and f.find('xi='+str(xi)+'_')!=-1
                                                # and ( f.find('N_pf='+str(N_pf)+'_')!=-1 or f.find('N_pf='+str(N_pf)+'.')!=-1)
                                                and f.find('N_pf='+str(N_pf)+'_')!=-1
                                                and f.find('seed_number='+str(seed_no) + '.') !=-1                                                
                                                and f.find(method)!=-1                                              
                                                and f.find(space_form)!=-1
                                                and (f.find('MeanAct') !=-1) and os.stat(dp + f).st_size != 0)]
    
    
    f = open(dp + right_file[0], "r")
    
    liste = []

    for line in f:
        liste.append(line.split())
    
    f.close()
    
    # print(liste)
    
    activities_mean = np.zeros(len(liste))
    # activities_rms = np.zeros_like(activities_mean)
    
    for ii in range(0, len(liste)):
        for jj in range(0,len(liste[ii])):
            # print(liste[ii][jj])
            activities_mean[ii] = liste[ii][0]
            # activities_rms[ii] = liste[ii][1]
            
    # return activities_mean, activities_rms
    return activities_mean

def load_PairwCorr_error(dp, space_form, method, T, freq, A_one_point_b, w_one_point, A_two_point_b, w_two_point, g_b, xi, N_pf, seed_no, corr_type):
    files = os.listdir(dp)
    
  
    
    if(corr_type == 'connected'):
        right_file = [f for f in files if ( #( f.find('T='+str(T)+'_')!=-1 or f.find('T='+str(T)+'.')!=-1)
                                                    f.find('T='+str(T)+'_')!=-1
                                                    # and ( f.find('f='+str(freq)+'_')!=-1 or f.find('f='+str(freq)+'.')!=-1)
                                                    and  f.find('f='+str(freq)+'_')!=-1
                                                    # and ( f.find('A_one_point_b='+str(A_one_point_b)+'_')!=-1 or f.find('A_one_point_b='+str(A_one_point_b)+'.')!=-1)
                                                    and f.find('A_one_point_b='+str(A_one_point_b)+'_')!=-1
                                                    # and ( f.find('w_one_point='+str(w_one_point)+'_')!=-1 or f.find('w_one_point='+str(w_one_point)+'.')!=-1)
                                                    and  f.find('w_one_point='+str(w_one_point)+'_')!=-1 
                                                    # and ( f.find('A_two_point_b='+str(A_two_point_b)+'_')!=-1 or f.find('A_two_point_b='+str(A_two_point_b)+'.')!=-1)
                                                    and ( f.find('A_two_point_b='+str(A_two_point_b)+'_')!=-1)
                                                    # and ( f.find('g_b='+str(g_b)+'_')!=-1 or f.find('g_b='+str(g_b)+'.')!=-1)
                                                    and f.find('g_b='+str(g_b)+'_')!=-1
                                                    # and ( f.find('w_two_point='+str(w_two_point)+'_')!=-1 or f.find('w_two_point='+str(w_two_point)+'.')!=-1)
                                                    and f.find('w_two_point='+str(w_two_point)+'_')!=-1
                                                    # and ( f.find('xi='+str(xi)+'_')!=-1 or f.find('xi='+str(xi)+'.')!=-1)
                                                    and f.find('xi='+str(xi)+'_')!=-1
                                                    # and ( f.find('N_pf='+str(N_pf)+'_')!=-1 or f.find('N_pf='+str(N_pf)+'.')!=-1)
                                                    and f.find('N_pf='+str(N_pf)+'_')!=-1
                                                    # and ( f.find('N_J='+str(N_J)+'_')!=-1 or f.find('N_J='+str(N_J)+'.')!=-1)
                                                    and f.find('seed_number='+str(seed_no) + '.')  !=-1                                                
                                                    and f.find(method)!=-1                                              
                                                    and f.find(space_form)!=-1
                                                    and f.find('raw') ==-1
                                                    and (f.find('Correlation') !=-1) 
                                                    and os.stat(dp + f).st_size != 0)]
    elif(corr_type == 'raw'):
        right_file = [f for f in files if ( #( f.find('T='+str(T)+'_')!=-1 or f.find('T='+str(T)+'.')!=-1)
                                                f.find('T='+str(T)+'_')!=-1
                                                    # and ( f.find('f='+str(freq)+'_')!=-1 or f.find('f='+str(freq)+'.')!=-1)
                                                    and f.find('f='+str(freq)+'_')!=-1
                                                    # and ( f.find('A_one_point_b='+str(A_one_point_b)+'_')!=-1 or f.find('A_one_point_b='+str(A_one_point_b)+'.')!=-1)
                                                    and f.find('A_one_point_b='+str(A_one_point_b)+'_')!=-1
                                                    # and ( f.find('w_one_point='+str(w_one_point)+'_')!=-1 or f.find('w_one_point='+str(w_one_point)+'.')!=-1)
                                                    and f.find('w_one_point='+str(w_one_point)+'_')!=-1
                                                    # and ( f.find('A_two_point_b='+str(A_two_point_b)+'_')!=-1 or f.find('A_two_point_b='+str(A_two_point_b)+'.')!=-1)
                                                    and f.find('A_two_point_b='+str(A_two_point_b)+'_')!=-1
                                                    # and ( f.find('g_b='+str(g_b)+'_')!=-1 or f.find('g_b='+str(g_b)+'.')!=-1)
                                                    and f.find('g_b='+str(g_b)+'_')!=-1
                                                    # and ( f.find('w_two_point='+str(w_two_point)+'_')!=-1 or f.find('w_two_point='+str(w_two_point)+'.')!=-1)
                                                    and f.find('w_two_point='+str(w_two_point)+'_')!=-1
                                                    # and ( f.find('xi='+str(xi)+'_')!=-1 or f.find('xi='+str(xi)+'.')!=-1)
                                                    and f.find('xi='+str(xi)+'_')!=-1
                                                    # and ( f.find('N_pf='+str(N_pf)+'_')!=-1 or f.find('N_pf='+str(N_pf)+'.')!=-1)
                                                    and f.find('N_pf='+str(N_pf)+'_')!=-1
                                                    # and ( f.find('N_J='+str(N_J)+'_')!=-1 or f.find('N_J='+str(N_J)+'.')!=-1)
                                                    and f.find('seed_number='+str(seed_no) + '.') !=-1                                                   
                                                    and f.find(method)!=-1                                              
                                                    and f.find(space_form)!=-1
                                                    and f.find('raw') !=-1
                                                    and (f.find('Correlation') !=-1) and os.stat(dp + f).st_size != 0)]
    
    print("The right file is", right_file[0], "that was the right file.")
    # print("The right file is", right_file[0], "that was the right file for the second time.")
    
    f = open(dp + right_file[0], "r")
    
    liste = []

    for line in f:
        liste.append(line.split())
    
    f.close()
    
    # print(liste)
    
    # N_pf_vertical = int(round(len(liste)/2. - 1.))
    # N_pf_horizontal = len(liste[0])
    
    # print('N_pf_vertical=', N_pf_vertical, ' N_pf_horizontal=', N_pf_horizontal, ' N_pf=', N_pf)
    
    correlations_mean = np.zeros((N_pf,N_pf), dtype=float)
    # correlations_rms = np.zeros_like(correlations_mean)
    
    for ii in range(0, N_pf):
        for jj in range(0,N_pf):
            # print(liste[ii][jj])
            # print(liste[ii+N_pf+1][jj])
            correlations_mean[ii][jj] = liste[ii][jj]
            # correlations_rms[ii][jj] = liste[ii+N_pf+1][jj]
            
    # return correlations_mean, correlations_rms
    return correlations_mean


