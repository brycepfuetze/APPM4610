import numpy as np
import math
import matplotlib.pyplot as plt

def driver():

# this demo code considers the boundary value problem
# y'' = p(x)y'+q(x)y + r(x)
# y'(a) = alpha, y'(b) = beta
# with both first and second order boundary condition approx. 

# boundary data
     a = 0
     b = 1
     alpha = math.pi
     beta = -math.pi

     # step size
     h = 0.1
     N = int((b-a)/h)
     
     x = np.linspace(a,b,N+1)
     
     yapp = make_FDmatNeu1st(x,h,N,alpha,beta)
          
#  exact solution 
     y = lambda x: np.sin(math.pi*x)
             
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

     p = -np.zeros(len(x))
     q = -2*np.ones(len(x))
     r = math.pi**2*np.sin(math.pi*x)+2*np.sin(math.pi*x)
     
     return(p,q,r)     



def make_FDmatNeu1st(x,h,N,alpha,beta):

# evaluate coefficients of differential equation
     (p,q,r) = eval_pqr(x)
 
# create the finite difference matrix     
     Matypp = 1/h**2*(np.diag(2*np.ones(N+1)) -np.diag(np.ones(N),-1) - 
           np.diag(np.ones(N),1))
           
     Matyp = 1/(2*h)*(np.diag(np.ones(N),1)-np.diag(np.ones(N),-1))
     
     B = np.matmul(np.diag(p,0),Matyp)
     
     A = Matypp + B+ np.diag(q)

# correct the entries to enforce the boundary condition
     A[0,0] = 1/h
     A[0,1] = -1/h
     A[N,N-1] = 1/h
     A[N,N] = -1/h
    
     A = -A
     
# create the right hand side rhs: (N) in size
     rhs = r
#  update with boundary data   
     rhs[0] = alpha
     rhs[N] = beta
     
     
# solve for the approximate solution

     Ainv = np.linalg.inv(A)
     yapp = np.matmul(Ainv,rhs)
     
     return yapp
     
def make_FDmatNeu2nd(x,h,N,alpha,beta):

# evaluate coefficients of differential equation
     (p,q,r) = eval_pqr(x)
 
# create the finite difference matrix     
     Matypp = 1/h**2*(np.diag(2*np.ones(N+1)) -np.diag(np.ones(N),-1) - 
           np.diag(np.ones(N),1))
           
     Matyp = 1/(2*h)*(np.diag(np.ones(N),1)-np.diag(np.ones(N),-1))
     
     B = np.matmul(np.diag(p,0),Matyp)
     
     A = Matypp + B+ np.diag(q)

# correct the entries to enforce the boundary condition
     A[0,0] = 1/h
     A[0,1] = -1/h
     A[N,N-1] = 1/h
     A[N,N] = -1/h
    
     A = -A
     
# create the right hand side rhs: (N) in size
     rhs = r
#  update with boundary data   
     rhs[0] = alpha
     rhs[N] = beta
     
     
# solve for the approximate solution

     Ainv = np.linalg.inv(A)
     yapp = np.matmul(Ainv,rhs)
     
     return yapp
     
driver()     
