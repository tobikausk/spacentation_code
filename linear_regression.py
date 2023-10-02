import numpy as np

def comp_cov_train_test(pop_act_sample, xi_sample, train_ratio):
    N_pf, N_time, N_xi = pop_act_sample.shape
    N_train = int(N_xi * train_ratio)
    N_test = N_xi - N_train

    xi_train = xi_sample[:N_train]
    xi_test = xi_sample[N_train:]
    pop_act_train = pop_act_sample[:,:, :N_train]
    pop_act_test = pop_act_sample[:,:, N_train:]

    xi_repeat_train = np.repeat(xi_train, N_time)
    xi_repeat_test = np.repeat(xi_test, N_time)

    pop_act_flat_train =  np.zeros((N_pf, N_time * N_train))
    pop_act_flat_test =  np.zeros((N_pf, N_time * N_test))

    for jj in range(N_pf):
        pop_act_flat_train[jj,:] = pop_act_train[jj,:,:].T.flatten()
        pop_act_flat_test[jj,:] = pop_act_test[jj,:,:].T.flatten()


    cov_act_train = np.dot(pop_act_flat_train, pop_act_flat_train.T)/(N_time * N_train) #Raw (unconnected) covariance is needed!!!
    cov_act_test = np.dot(pop_act_flat_test, pop_act_flat_test.T)/(N_time * N_test)

    cov_act_xi_train = np.zeros(N_pf)
    cov_act_xi_test = np.zeros(N_pf)
    
    for ii in range(N_pf):
        cov_act_xi_train[ii] = np.mean(pop_act_flat_train[ii,:] * xi_repeat_train)
        cov_act_xi_test[ii] = np.mean(pop_act_flat_test[ii,:] * xi_repeat_test)

    var_xi_train = np.mean(xi_train**2)
    var_xi_test = np.mean(xi_test**2)

    return cov_act_train, cov_act_xi_train, var_xi_train, cov_act_test, cov_act_xi_test, var_xi_test

def lin_reg_with_res(cov_act_train, cov_act_xi_train, var_xi_train, cov_act_test, cov_act_xi_test, var_xi_test):
    N_pf = cov_act_train.shape[0]
    cov_act_train_inv = np.linalg.inv(cov_act_train)
    J = np.dot(cov_act_train_inv, cov_act_xi_train)
    resid_train = np.dot(J, np.dot(cov_act_train,J)) - 2. * np.dot(J, cov_act_xi_train) + var_xi_train
    resid_test = np.dot(J, np.dot(cov_act_test,J)) - 2. * np.dot(J, cov_act_xi_test) + var_xi_test
    return J, resid_train, resid_test

def lin_reg_regularized(cov_act_train, cov_act_xi_train, var_xi_train, cov_act_test, cov_act_xi_test, var_xi_test, lbda):
    N_pf = cov_act_train.shape[0]
    cov_act_train_inv = np.linalg.inv(cov_act_train + lbda * np.eye(N_pf))
    J = np.dot(cov_act_train_inv, cov_act_xi_train)
    resid_train = np.dot(J, np.dot(cov_act_train,J)) - 2. * np.dot(J, cov_act_xi_train) + var_xi_train
    resid_test = np.dot(J, np.dot(cov_act_test,J)) - 2. * np.dot(J, cov_act_xi_test) + var_xi_test
    return J, resid_train, resid_test