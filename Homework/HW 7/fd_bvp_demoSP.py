import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse import csc_matrix

def driver():

# this demo code considers the boundary value problem
# y'' = p(x)y'+q(x)y + r(x)
#y(a) = alpha, y(b) = beta

# boundary data
     a = 1
     b = 2
     alpha = 1
     beta = 2

     # step size
     h = 0.1
     N = int((b-a)/h)
     
     x = np.linspace(a,b,N+1)
     
     yapp = make_FDmatDir(x,h,N,alpha,beta)
          
#  exact solution 
     c2 = 1/70*(8-12*np.sin(np.log(2))-4*np.cos(np.log(2)))
     c1 = 11/10-c2
     y = lambda x: c1*x+c2/(x**2)-3/10*np.sin(np.log(x))-1/10*np.cos(np.log(x))
             
     yex = y(x)
    
     plt.plot(x,yapp,label = 'FD aprox')
     plt.plot(x,yex,label = 'Exact')
     plt.xlabel('x')
     plt.legend(loc = 'upper left')
     plt.show()
     
     err = np.zeros(N+1)
     for j in range(0,N+1):
          err[j] = abs(yapp[j]-yex[j])
          
     print('err = ', err)
          
     plt.plot(x,err,label = 'FD aprox')
     plt.xlabel('x')
     plt.xlabel('absolute error')
     plt.legend(loc = 'upper left')
     plt.show()
          
     return
     
def eval_pqr(x):

     p = -2/x
     q = 2/(x**2)
     r = np.sin(np.log(x))/(x**2)
     
     return(p,q,r)     



def make_FDmatDir(x,h,N,alpha,beta):

# evaluate coefficients of differential equation
     (p,q,r) = eval_pqr(x)
     
 
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

     return yapp
     
driver()     
