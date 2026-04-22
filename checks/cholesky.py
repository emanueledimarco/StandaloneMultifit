import numpy as np
from scipy.linalg import cholesky, solve_triangular


#invcov = np.loadtxt("matrici.csv", delimiter=",")
invcov = np.loadtxt("matrice_pulse_cov_ridotta_16x16.csv", delimiter=",")
#print(np.tril(invcov, k=-1))
invcov += np.tril(invcov, k=-1).T

print(invcov)
L = cholesky(invcov, lower=True)

print(L)
