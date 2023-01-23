import numpy as np
import galois
import scipy as sp

GF = galois.GF(2)

def initial_factor(A):
    n, _ = A.shape

    N1 = [
        k for k in range(n)
        if A[k, k] != sum((A[k, j] for j in range(n) if j != k), start=GF(0))
    ]

    N2 = [
        (i, j) for i in range(n) for j in range(n)
        if i < j and A[i, j] != 0
    ]

    E = GF.Zeros((n, len(N1) + len(N2)))

    for k, i in enumerate(N1):
        E[i, k] = 1

    for k, (i, j) in enumerate(N2):
        E[i, len(N1) + k] = 1
        E[j, len(N1) + k] = 1
    
    return E

def improve_factor(B, target):
    n, m = B.shape
    nullsp = B.null_space()

    print(f"{m = }")
    if m == target:
        return B

    # for i in range(nullity):
    #     if (nullsp[i, :] == 1).all():
    #         continue
    #     else:
    #         y = nullsp[i, :]
    #         break
    # else:
    #     raise RuntimeError("couldn't find low-weight vector in the nullspace")
    print(f"{B = }")
    print(f"{target = }")
    y = nullsp[0, :]

    if np.sum(y) == 1:
        y = np.append(y, 1)
        B = np.hstack((B, GF.Zeros((n, 1))))

    for a in range(m):
        for b in range(a):
            if y[a] + y[b] == 1:
                break
        else:
            continue
        break

    z = B[:, a] + B[:, b]
    B = B + np.outer(z, y)
    B = np.delete(B, [a, b], 1)
    return B
    
def find_factor(A):
    n, _ = A.shape
    rank = np.linalg.matrix_rank(A)
    delta = int(all(A[i, i] == 0 for i in range(n)))
    target = rank + delta

    B = initial_factor(A)
    while True:
        Bp = improve_factor(B, target)
        if B.shape == Bp.shape:
            return B
        else:
            B = Bp
