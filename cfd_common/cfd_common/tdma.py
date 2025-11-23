# cfd_tdma/solver.py
import numpy as np


def TDMA(start, end, AW, AE, AP, SOR):
    n = len(AP)
    X = np.zeros(n)
    c_prime = np.zeros(n)
    d_prime = np.zeros(n)

    c_prime[start] = AE[start] / AP[start]
    d_prime[start] = SOR[start] / AP[start]
    
    for i in range(start + 1, end + 1):
        denom = AP[i] - AW[i] * c_prime[i - 1]
        c_prime[i] = AE[i] / denom
        d_prime[i] = (SOR[i] - AW[i] * d_prime[i - 1]) / denom

    X[end] = d_prime[end]
    for i in range(end - 1, start - 1, -1):
        X[i] = d_prime[i] - c_prime[i] * X[i + 1]

    return X


def LLTDMA(i_start, i_end, j_start, j_end, AW, AE, AS, AN, AP, SOR, n_iter=2):
    IDIM, JDIM = AP.shape
    SOR_x = np.zeros(IDIM)
    SOR_y = np.zeros(JDIM)
    X = np.zeros((IDIM, JDIM))

    for _ in range(n_iter):
        # TDMA: x-direction
        for j in range(j_start, j_end + 1):
            for i in range(i_start, i_end + 1):
                neighbors = AS[i][j] * X[i][j - 1] + AN[i][j] * X[i][j + 1]
                SOR_x[i] = SOR[i][j] + neighbors
            X[:, j] = TDMA(i_start, i_end, -AW[:, j], -AE[:, j], AP[:, j], SOR_x)

        # TDMA: y-direction
        for i in range(i_start, i_end + 1):
            for j in range(j_start, j_end + 1):
                neighbors = AW[i][j] * X[i - 1][j] + AE[i][j] * X[i + 1][j]
                SOR_y[j] = SOR[i][j] + neighbors
            X[i, :] = TDMA(j_start, j_end, -AS[i, :], -AN[i, :], AP[i, :], SOR_y)

    return X
