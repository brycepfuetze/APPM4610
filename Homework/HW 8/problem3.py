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
     u = lambda x: np.sin(3 * np.pi * x)
     sol = lambda x: np.sin(3 * np.pi * x)

     a = 0
     b = 1 
     
# N = number of nodes +1;
     N = np.logspace(1,4,100)
     err = np.zeros(len(N))
     h = np.zeros(len(N))
# space between nodes
     for i in range(len(N)):
        h[i] = (b-a)/N[i]
     
        xh = np.linspace(a,b,int(N[i])+1)

        A = make_Matrix(xh,h[i],int(N[i]))

        rhs = make_rhs(xh,h[i],int(N[i]))
        yapprox = sp.linalg.spsolve(A,rhs)
        yapprox = np.pad(yapprox, (1, 1), 'constant', constant_values=(0,0))

        ysol = sol(xh)
        print('h = ', h[i])
        err[i] = np.linalg.norm(yapprox - ysol)
        print('Error: ', err[i])
        # inv = np.linalg.inv(A)
        # K = np.linalg.norm(A) * np.linalg.norm(inv)
        # print('Condition number K: ', K)


     fig = go.Figure()
     fig.add_trace(go.Scatter(x=h, y=err, mode='lines+markers', name='Forward/Backward Difference'))

    # Add titles and labels
     fig.update_layout(title='Convergence of Finite Element Method',
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


     
# solve for the approximate solution 
# at the interior nodes
     sol = sp.linalg.spsolve(A,rhs)

# create the vector with the approximations at the 
# nodes     
     uapp = np.zeros(N+1)
     
     for j in range(1,N):
         uapp[j] = sol[j-1]
         
     uex = u(xh)
     
     abserr = np.linalg.norm(uapp-uex)
     print(abserr)
     
     plt.plot(xh,uapp,label = 'FEM aprox')
     plt.plot(xh,uex,label = 'Exact')
     plt.xlabel('x')
     plt.legend(loc = 'upper left')
     plt.show()

     err = np.zeros(N+1)
     for j in range(0,N+1):
          err[j] = abs(uapp[j]-uex[j])
          
          
     plt.plot(xh,err,label = 'FEM aprox')
     plt.xlabel('x')
     plt.xlabel('absolute error')
     plt.legend(loc = 'upper left')
     plt.show()
     
     
     return

     
def eval_k(x):
     k = (1+x) * np.exp(x)
     return k
      
def eval_q(x):
      q = 0
      return q        

def eval_f(x):
     f = -(-9*np.pi**2 * (x+1) * np.sin(3* np.pi * x) + 3 * np.pi * (x+1) * np.cos(3*np.pi * x) + 3*np.pi * np.cos(3*np.pi * x))*np.exp(x)
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
          
     
def make_Matrix(x,h,N):

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
         d_main[j] = (1/h**2)*tmp
     
     for j in range(0,N-2):
        b = x[j+1]
        c = x[j+2]     
        tmp1 = quad(eval_k,b,c)
        tmp2 = quad(eval_k,b,c) 
        d_sub[j] = -(1/h**2)*tmp1[0]
        d_sup[j] = -(1/h**2)*tmp2[0]
        
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
