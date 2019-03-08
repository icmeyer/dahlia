import numpy as np
from scipy.linalg import expm


def solve_matrix(matrix, conc, flux, steps, time):
    dt = time/steps
    conc_over_time = np.zeros([matrix.shape[0], steps])
    conc_over_time[:, 0] = conc
    for i in range(steps-1):
        conc = matrix_exp_solve(matrix, dt, conc)
        conc_over_time[:, i+1] = conc

    return conc_over_time
        
    
def matrix_exp_solve(matrix, dt, conc):
    # print('exponent evaluation')
    # print(matrix)
    # print(np.exp(-matrix*dt))
    # print(conc)
    return np.matmul(expm(-matrix*dt),conc)


    
    

