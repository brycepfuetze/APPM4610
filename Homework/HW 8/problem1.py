import numpy as np
import plotly as px 
import plotly.graph_objects as go
import scipy.linalg as linalg

def delta(j,k):
    if j == k:
        return 1
    else:
        return 0

def aphi(j,k):
    return np.pi**2 * (j+1) * (k+1) * delta(j,k) + 2 * np.sin((j+1)* np.pi / 2) * np.sin((k+1) * np.pi / 2)

def makeSchrodingerApprox(N):

    A = np.zeros((N, N))
    S = np.zeros(N)
    x = np.linspace(0, 1, N)

    
    for k in range(N):
        S[k] = - np.sqrt(2) / ((k+1) * np.pi) * (np.cos((k+1) * np.pi) - 1)
        for j in range(N):
            A[j,k] = aphi(j,k)
        
    print(A)
    approx = linalg.solve(A, S)
    return approx, x
    

# try with N = 5 and N = 35
# and create our y values using the basis functions sqrt(2) * sin((k+1) * pi * x)
N = 5
approx5, x5 = makeSchrodingerApprox(N)
y5 = np.zeros((N, len(x5)))
yfive = np.zeros(len(x5))
for k in range(N):
    y5[k] = np.sqrt(2) * np.sin((k+1) * np.pi * x5)
    yfive = yfive + approx5[k] * y5[k]



N = 35
approx35, x35 = makeSchrodingerApprox(N)
y35 = np.zeros((N, len(x35)))
y3five = np.zeros(len(x35))
for k in range(N):
    y35[k] = np.sqrt(2) * np.sin((k+1) * np.pi * x35)
    y3five = y3five + approx35[k] * y35[k]




# plot
fig = go.Figure()
fig.add_trace(go.Scatter(x=x5, y=yfive, mode='lines+markers', name='N=5'))
fig.add_trace(go.Scatter(x=x35, y=y3five, mode='lines+markers', name='N=35'))

fig.update_layout(title='Schrodinger Approximation',
                  xaxis_title='x',
                  yaxis_title='Approximation',
                  xaxis_type='linear',
                  yaxis_type='linear')
fig.show()