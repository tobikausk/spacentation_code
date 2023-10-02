import numpy as np

def compute_bias(J, n, xi_list):
    """ Compute bias for linear readout with the readout weights J, the input n
    and the output xi."""
    N_neurons, N_rep, N_sample = n.shape
    if(N_neurons != len(J)):
        raise ValueError("Number of neurons assumed in J and present in n do not match.")
    
    xi_predict_list = np.mean(np.sum(J[:, None, None] * n, axis=0), axis = 0)
    bias_list = xi_predict_list - xi_list
    return bias_list, xi_predict_list

def compute_deriv(xi, bias):
    """ Compute the derivative of the bias with respect to the sample xi."""

    argsort_xi = np.argsort(xi)
    xi_sorted = xi[argsort_xi]
    bias_sorted = bias[argsort_xi]

    return xi_sorted, bias_sorted, np.diff(bias_sorted)/np.diff(xi_sorted)