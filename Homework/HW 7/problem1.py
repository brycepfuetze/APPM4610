import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


# my first attempt because my I thought 'maybe I can do this without sparse
# matrices'. I was wrong, this is not sufficient.
def makeSys(a, b, h, alpha, beta, r):
    N = int((b - a) / h)
    A = np.zeros((N + 1, N + 1))
    S = np.zeros(N + 1)
    
    A[0, 0] = 1
    A[N, N] = 1
    S[0] = alpha
    S[N] = beta

    x = np.linspace(a, b, N + 1)

    for i in range(1, N):
        A[i, i - 1] = 1 / h**2
        A[i, i] = -2 / h**2
        A[i, i + 1] = 1 / h**2
        S[i] = r(a + i * h)

    return A, S, x

a = 0
b = 1
alpha = 1
beta = np.exp(2)
r = lambda x: 4 * np.exp(2 * x)
h = np.logspace(-1, -4, 4)

sol = lambda x: np.exp(2 * x)
err = np.zeros(len(h))

for i in range(len(h)):
    A, S, x = makeSys(a, b, h[i], alpha, beta, r)
    y = spsolve(A, S)
    ysol = sol(x)
    print('h = ', h[i])
    err[i] = np.linalg.norm(y - ysol)
    print('Error: ', err[i])
    inv = np.linalg.inv(A)
    K = np.linalg.norm(A) * np.linalg.norm(inv)
    print('Condition number K: ', K)
