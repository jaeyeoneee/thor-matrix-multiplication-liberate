import numpy as np

def lower_diagonal_vector(A, k):
    """
    
    A : a matrix
    k : k-th lower diagonal vector

    """

    a, b = A.shape
    max_v = max(a, b)

    Lk = np.zeros(max_v)
    for t in range(max_v):
        i = (t + k) % a
        j = t % b
        Lk[t] = A[i, j]
    
    return Lk

def mask(n, c, ell):
    s = n*c
    mu_l0 = np.zeros(s)
    mu_l1 = np.zeros(s)
    mu_l2 = np.zeros(s)
    mu_l3 = np.zeros(s)

    ell_mod_c = ell % c

    for i in range(s):
        r = i // n
        local_idx = i % n

        is_front = (local_idx < (n - ell))
        
        if (ell_mod_c <= r < c):
            if is_front:
                mu_l1[i] = 1.0
            else:
                mu_l3[i] = 1.0
        else:
            if is_front:
                mu_l0[i] = 1.0
            else:
                mu_l2[i] = 1.0

    return mu_l0, mu_l1, mu_l2, mu_l3