import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import plotly.express as px
import plotly.graph_objects as go


# use my first attempt at problem 1...
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
        A[i, i - 1] = (x[i] * np.exp(x[i]) + np.exp(x[i])) / (h**2) + (2*x[i] * np.exp(x[i]) + np.exp(x[i])) / (2*h)
        A[i, i] = -2 * (x[i] * np.exp(x[i]) + np.exp(x[i])) / (h**2)
        A[i, i + 1] = (x[i] * np.exp(x[i]) + np.exp(x[i])) / (h**2) + (2*x[i] * np.exp(x[i]) + np.exp(x[i])) / (2*h)
        S[i] = r(a + i * h)

    return A, S, x

a = 0
b = 1
alpha = 0
beta = 0
r = lambda x: np.exp(x) * (-9*np.pi**2 * (x+1) * np.sin(3* np.pi * x) + 3 * np.pi * (x+1) * np.cos(3*np.pi * x) + 3*np.pi * np.cos(3*np.pi * x))
h = np.logspace(-1, -4, 10)

sol = lambda x: np.sin(3 * np.pi * x)
err = np.zeros(len(h))

for i in range(len(h)):
    A, S, x = makeSys(a, b, h[i], alpha, beta, r)
    y = np.linalg.solve(A, S)
    ysol = sol(x)
    print('h = ', h[i])
    err[i] = np.linalg.norm(y - ysol)
    print('Error: ', err[i])
    inv = np.linalg.inv(A)
    K = np.linalg.norm(A) * np.linalg.norm(inv)
    print('Condition number K: ', K)


fig = go.Figure()
fig.add_trace(go.Scatter(x=h, y=err, mode='lines+markers', name='Forward/Backward Difference'))

# Add titles and labels
fig.update_layout(title='Convergence of Derived Difference Method',
                  xaxis_title='h',
                  yaxis_title='||error||',
                  xaxis_type='log',
                  yaxis_type='log')
fig.show()

fig = go.Figure()
fig.add_trace(go.Scatter(x=x, y=y, mode='lines+markers', name='Approx'))
fig.add_trace(go.Scatter(x=x, y=ysol, mode='lines+markers', name='Exact'))
# Add titles and labels
fig.update_layout(title='Boundary Condition Method Comparison',
                  xaxis_title='x',
                  yaxis_title='y')
fig.show()