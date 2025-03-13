import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import scipy.sparse as sp
import plotly.express as px


# attempt number 2 stealing the sparse code from class. I tried it on my
# own, I do not understand enough about how Gillman used sparse matrices
def eval_pqr(x):

    p = 0*x
    q = 0*x
    r = 4* np.exp(2*x)
     
    return(p,q,r)  


def makeSys(a, b, h, alpha, beta, r):
    N = int((b - a) / h)
    x = np.linspace(a, b, N + 1)

# evaluate coefficients of differential equation
    p,q,r = eval_pqr(x)
     
 
# create the finite difference matrix     
    Matypp = 1/h**2*(sp.diags(2*np.ones(N-1)) -sp.diags(np.ones(N-2),-1) - 
           sp.diags(np.ones(N-2),1))
           
    Matyp = 1/(2*h)*(sp.diags(np.ones(N-2),1)-sp.diags(np.ones(N-2),-1))
     
    A = Matypp +sp.diags(p[1:N],0)@Matyp + sp.diags(q[1:N])

    

# create the right hand side rhs: (N-1) in size
    rhs = -r[1:N]
#  update with boundary data   
    rhs[0] = rhs[0] + (1/h**2-1/(2*h)*-p[1])*alpha
    rhs[N-2] = rhs[N-2] + (1/h**2+1/(2*h)*-p[N-1])*beta
     
# solve for the approximate solution

    sol = sp.linalg.spsolve(A,rhs)
     
    yapp = np.zeros(N+1)
    yapp[0] = alpha
    for j in range(1,N):
        yapp[j] = sol[j-1]
    yapp[N] = beta    


    return yapp, x

a = 0
b = 1
alpha = 1
beta = np.exp(2)
r = lambda x: 4 * np.exp(2 * x)
h = np.logspace(-1, -5, 5)

sol = lambda x: np.exp(2 * x)
err = np.zeros(len(h))

for i in range(len(h)):
    yapp, x = makeSys(a, b, h[i], alpha, beta, r)
    ysol = sol(x)
    err[i] = np.linalg.norm(yapp - ysol)
    print(err[i])


fig = px.line(x = h, y = err, log_y = True,log_x=True, title='Error for various h: Finite Difference', labels={'x': 'h', 'y':'||error||'})
fig.update_layout(
    xaxis=dict(exponentformat='none'),
    yaxis=dict(exponentformat='none')
)
fig.show()