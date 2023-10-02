# %%
import numpy as np
from scipy.optimize import bisect
from compute_FisherInfo import compute_FisherInfo_from_scratch
from matplotlib import pyplot as plt

# Parameters:
T = 0.2
# A_two_point_start = 1.
w_two_point = 0.1
# A_one_point_start = 0.45
w_one_point = 0.07
# g_b = 1.
g = 1.
xi = 0.2
mean_act = 0.15
N_x = 100

N_sample_Gauss = 10000

FisherInfo_wish = 0.015

tol_phi = 0.001
tol_phi_abs = 0.0001
tol_q = 0.001
tol_q_abs = 0.0001
eps_phi = 0.1
eps_rho = 0.1

def find_Atwopoint(A_one_point, FisherInfo_goal, Atwopoint_minus, Atwopoint_plus,
                   w_one_point, w_two_point, xi, mean_act, g, T, N_x, N_sample_Gauss,
                   tol_phi, tol_phi_abs, tol_q, tol_q_abs, eps_phi, eps_rho):
    """
    Find A_two_point that gives FisherInfo_goal
    """
    def match_FisherInfo(Atwopoint):
        print("Atwopoint=", Atwopoint)
        FisherInfo = compute_FisherInfo_from_scratch(A_one_point, w_one_point, Atwopoint, w_two_point, xi, mean_act, 
                                               g, T, N_x, N_sample_Gauss,
                                               tol_phi, tol_phi_abs, tol_q, tol_q_abs, eps_phi, eps_rho)[0]/N_x
        print("Fisher info=", FisherInfo)
        return FisherInfo - FisherInfo_goal
    
    A_two_point = bisect(match_FisherInfo, Atwopoint_minus, Atwopoint_plus)
    return A_two_point

# A_one_point_list = np.array([0.019, 0.02, 0.022, 0.024, 0.026]) # For FisherInfo_wish = 0.02
# A_one_point_list = np.array([0.012, 0.014, 0.016, 0.018]) # For FisherInfo_wish = 0.01
A_one_point_list = np.array([0.016, 0.018, 0.02, 0.022]) # For FisherInfo_wish = 0.015
# A_two_point_minus_list = np.array([8., 5., 5., 1., -1.]) # For FisherInfo_wish = 0.02
# A_two_point_minus_list = np.array([-1., -1., -1., -1.]) # For FisherInfo_wish = 0.01
A_two_point_minus_list = np.array([9., 5., 1., -1.]) # For FisherInfo_wish = 0.015
# A_two_point_plus_list = np.array([16., 14., 12., 10.,8.]) # For FisherInfo_wish = 0.02
# A_two_point_plus_list = np.array([14., 12., 10.,8.]) # For FisherInfo_wish = 0.02
A_two_point_plus_list = np.array([16., 12., 10.,8.]) # For FisherInfo_wish = 0.015

A_two_point_list = np.zeros_like(A_one_point_list)

for ii, A_one_point in enumerate(A_one_point_list):
    A_two_point_list[ii] = find_Atwopoint(A_one_point, FisherInfo_wish, A_two_point_minus_list[ii], A_two_point_plus_list[ii],
                                         w_one_point, w_two_point, xi, mean_act, g, T, N_x, N_sample_Gauss,
                                         tol_phi, tol_phi_abs, tol_q, tol_q_abs, eps_phi, eps_rho)


plt.figure(1)
plt.plot(A_one_point_list, A_two_point_list)
# %%
datapath = 'data/Theory/FisherInfo_vs_Aone_Atwo/'

with open( datapath + 'const_longer_FisherInfowish=' + str(FisherInfo_wish) + '_T=' + str(T) + '_f=' + str(mean_act) 
                        + '_w_one_point=' + str(w_one_point) + '_w_two_point=' + str(w_two_point) 
                        + '_xi='+str(xi) + '_N_x=' + str(N_x), 'wb') as f:        
        
    np.save(f, FisherInfo_wish )
    np.save(f, A_one_point_list)
    np.save(f, A_two_point_list)

with open( datapath + 'const_longer_FisherInfowish=' + str(FisherInfo_wish) + '_T=' + str(T) + '_f=' + str(mean_act) 
                        + '_w_one_point=' + str(w_one_point) + '_w_two_point=' + str(w_two_point) 
                        + '_xi='+str(xi) + '_N_x=' + str(N_x), 'rb') as f:

        FisherInfo_wish_load = np.load(f)
        A_one_point_list_load = np.load(f)
        A_two_point_list_load = np.load(f)


# %%
