import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.sparse import csc_matrix
import scipy.sparse as sp
import plotly.express as px
import plotly.graph_objects as go



def driver():

# this code considers the following boundary value
# problem 
# - d/dx( k(x) du/dx) +q(x)u(x) = f(x) for x in (a,b)
# u(a) = u(b) = 0

# exact solution
    u = lambda x: x**2-x

    a = 0
    b = 1 
    
     
# N = number of nodes +1;
    N = 32
# space between nodes
    h = (b-a)/N
     
    xh = np.linspace(a,b,N+1)

    sol = np.zeros((N-1,3))
    for i in range(1,4):
        A = make_Matrix(xh,h,N, i)
        rhs = make_rhs(xh,h,N)
     
# solve for the approximate solution 
# at the interior nodes
        sol[:,i-1] = sp.linalg.spsolve(A,rhs)

# # create the vector with the approximations at the 
# # nodes     
#     uapp = np.zeros(N+1)
     
#     for j in range(1,N):
#          uapp[j] = sol[j-1]
         
#     uex = u(xh)
     
#     abserr = np.linalg.norm(uapp-uex)
#     print(abserr)
     
#     plt.plot(xh,uapp,label = 'FEM aprox')
#     plt.plot(xh,uex,label = 'Exact')
#     plt.xlabel('x')
#     plt.legend(loc = 'upper left')
#     plt.show()

#     err = np.zeros(N+1)
#     for j in range(0,N+1):
#          err[j] = abs(uapp[j]-uex[j])
          
          
#     plt.plot(xh,err,label = 'FEM aprox')
#     plt.xlabel('x')
#     plt.xlabel('absolute error')
#     plt.legend(loc = 'upper left')
#     plt.show()
     # plot each solution
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=xh, y=sol[:,0], mode='lines+markers', name='epsilon = 0.1'))
    fig.add_trace(go.Scatter(x=xh, y=sol[:,1], mode='lines+markers', name='epsilon = 0.25'))
    fig.add_trace(go.Scatter(x=xh, y=sol[:,2], mode='lines+markers', name='epsilon = 1'))
    fig.update_layout(title='FEM Approximation',
                      xaxis_title='x',
                      yaxis_title='Approximation',
                      xaxis_type='linear',
                      yaxis_type='linear')
    fig.show()
    
     
    return

     
def eval_k(i):
    k = 0.1
    return k
      
def eval_q(x):
      q = 1
      return q        

def eval_f(x):
     f = 1
     return f

      
def eval_stiffD(x,xj):
       # evaluates integrand for diagonal of 
       # stiffness matrix
     val = eval_q(x)*(x-xj)**2
     return val        

def eval_stiffO(x,xj,xk):
      # evaluates integrand of off-diag of 
      # stiffness matrix
     val = eval_q(x)*(x-xj)*(xk-x)
     return val        


def eval_rhsInt1(x,xk):
     # evaluate the integrand for the rhs 
     # left side of hat
     val = eval_f(x)*(x-xk)
     return val
     
def eval_rhsInt2(x,xk):
     # evaluate the integrand for the rhs 
     # right side of hat
     val = eval_f(x)*(xk-x)
     return val
          
     
def make_Matrix(x,h,N,i):
     match i:
        case 1:
            k=0.1
        case 2:
            k=0.25
        case 3:
            k = 1
     # make the entries of the mass matrix
     d_main = np.zeros(N-1)
     d_sub = np.zeros(N-2)
     d_sup = np.zeros(N-2)
     
     # make the diagonal entries
     for j in range(0,N-1):
         a = x[j]
         b = x[j+1]
         c = x[j+2]
         tmp1 = quad(eval_k,a,b)
         tmp2 = quad(eval_k,b,c)
         tmp = tmp1[0] + tmp2[0]
         d_main[j] = (1/h**2)*tmp * k
     
     for j in range(0,N-2):
        b = x[j+1]
        c = x[j+2]     
        tmp1 = quad(eval_k,b,c)
        tmp2 = quad(eval_k,b,c) 
        d_sub[j] = -(1/h**2)*tmp1[0] * k
        d_sup[j] = -(1/h**2)*tmp2[0] * k
        
# make entries for the stiffness matrix 
     for j in range(0,N-1):
         a = x[j]
         b = x[j+1]
         c = x[j+2]
         tmp1 = quad(eval_stiffD,a,b,args = (a,))
         tmp2 = quad(eval_stiffD,b,c,args = (c,))
         tmp = tmp1[0] + tmp2[0]
         d_main[j] = d_main[j]+(1/h**2)*tmp

     for j in range(0,N-2):
        a = x[j]
        b = x[j+1]
        c = x[j+2]     
        tmp1 = quad(eval_stiffO,a,b, args = (a,b,))
        tmp2 = quad(eval_stiffO,b,c, args = (b,c,) ) 
        d_sub[j] = d_sub[j]+(1/h**2)*tmp1[0]
        d_sup[j] = d_sup[j]+(1/h**2)*tmp2[0]

# list of all the data
     data = [d_sub, d_main, d_sup]   
 # which diagonal each vector goes into
     diags = [-1,0,1]   
 # create the matrix                
     A = sp.diags(data,diags,format='csc') 
     
     return(A)

def make_rhs(x,h,N):

     rhs = np.zeros(N-1)
     
     for j in range(0,N-1):
        a = x[j]
        b = x[j+1]
        c = x[j+2]     
        tmp1 = quad(eval_rhsInt1,a,b,args = (a,))
        tmp2 = quad(eval_rhsInt2,b,c,args = (c,))
        rhs[j] = (1/h)*tmp1[0] + (1/h)*tmp2[0]
     
     return(rhs)
     
driver()     
