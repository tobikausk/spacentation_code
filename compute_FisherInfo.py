import numpy as np
from attractor_with_input import J_w, rho_T0, rel_center
from solve_sp_equations_space_disorder import phi_q_self_consistent
from calc_CrossCorr_fixedMean_disorder import pairwise_corr, stat_upto_fourth_order
from calc_FisherInfo import FisherInfo_from_corr

def compute_FisherInfo_from_scratch(A_one_point, w_one_point, A_two_point, w_two_point, xi, mean_act, g, T, N_x,
                                    N_sample_Gauss, tol_phi, tol_phi_abs, tol_q, tol_q_abs, eps_phi, eps_rho):
    """
    Compute Fisher information given the basic parameters
    """
    ###################################################################
    ### Set derived quantities not requiring elaborate computations ###
    ###################################################################
    Delta_x = 1./N_x
    xs_theo = np.linspace(0.,1.,N_x, endpoint = False)
    rho_0 = rho_T0(xs_theo, xi, mean_act)
    
    J_dist_list =  J_w(xs_theo, w_two_point, A_two_point)
    J_dist = np.zeros((N_x, N_x), dtype=float)
    J_dist_row = np.zeros(N_x, dtype = float)    
        
    for ii in range(0,N_x):
        del_x = rel_center(xs_theo[0], xs_theo[ii])
        J_dist_row[ii] = (J_w(del_x, w_two_point, A_two_point))/N_x
             
    for ii in range(0,N_x):
        J_dist[ii,:] = J_dist_row
        J_dist_row = np.roll(J_dist_row, 1)
    ###################################################################

    (phi_final, q_final, rho_final, 
    lbda_final) = phi_q_self_consistent(rho_0, mean_act**2, g, 
                                        J_dist_list, Delta_x, T, mean_act, 
                                        A_one_point, w_one_point, xi, xs_theo, 
                                        -0.5 * w_two_point, N_sample_Gauss,
                                        tol_phi, tol_phi_abs, tol_q, tol_q_abs, eps_phi, eps_rho)
                                                     
    cumulants, moments = stat_upto_fourth_order(A_one_point, w_one_point, J_dist_list, g, T, xi, 
                                                rho_final, q_final, lbda_final, 
                                                N_sample_Gauss, xs_theo)
                
    var_vec = cumulants[1]
    var_mat = np.diag(var_vec)
                                
    J_dist_eff = T * np.dot(np.linalg.inv(np.identity(N_x) - np.dot(J_dist/T, var_mat)),J_dist/T)

    (cov, cov_direct, cov_indirect) = pairwise_corr(A_one_point, w_one_point, 
                                        J_dist_list, g, T, xi, 
                                        rho_final, q_final, lbda_final, J_dist_eff,
                                        N_sample_Gauss, xs_theo, 1.)
            
    corr_theo = np.diag(var_vec) + cov
                
    FisherInfo = FisherInfo_from_corr(corr_theo, xs_theo, xi, A_one_point, w_one_point, T)
    FisherInfo_from_var = FisherInfo_from_corr(np.diag(var_vec), xs_theo, xi, A_one_point, w_one_point, T)
    FisherInfo_from_cov_direct = FisherInfo_from_corr(cov_direct, xs_theo, xi, A_one_point, w_one_point, T)
    FisherInfo_from_cov_indirect = FisherInfo_from_corr(cov_indirect, xs_theo, xi, A_one_point, w_one_point, T)

    return FisherInfo, FisherInfo_from_var, FisherInfo_from_cov_direct, FisherInfo_from_cov_indirect