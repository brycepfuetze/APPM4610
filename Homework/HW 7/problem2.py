import numpy as np
import math
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

# FROM CODE SUPPLIED ON CANVAS

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
     # remove offset
     yapp = yapp - yapp[0]
          
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

     h = [0.1,0.05, 0.025, 0.0125, 0.00625]
     err1 = np.zeros(len(h))
     err2 = np.zeros(len(h))
     for i in range(len(h)):
        N = int((b-a)/h[i])
        x = np.linspace(a,b,N+1)
        yapp1 = make_FDmatNeu1st(x,h[i],N,alpha,beta)
        yapp2 = make_FDmatNeu2nd(x,h[i],N,alpha,beta)
        yapp1 = yapp1 - yapp1[0]
        yapp2 = yapp2 - yapp2[0]
        yex = y(x)

        err1[i] = np.abs(yex[int(N/2)] - yapp1[int(N/2)])
        err2[i] = np.abs(yex[int(N/2)] - yapp2[int(N/2)])

     # Create the plot
     fig = go.Figure()
     fig.add_trace(go.Scatter(x=h, y=err1, mode='lines+markers', name='Forward/Backward Difference'))
     # Add the second trace
     fig.add_trace(go.Scatter(x=h, y=err2, mode='lines+markers', name='Centered Difference'))
     # Add titles and labels
     fig.update_layout(title='Boundary Condition Method Comparison',
                  xaxis_title='h',
                  yaxis_title='||error||',
                  xaxis_type='log',
                  yaxis_type='log')
     fig.show()


     return
     
def eval_pqr1(x):

     p = -np.zeros(len(x))
     q = -2*np.ones(len(x))
     r = math.pi**2*np.sin(math.pi*x)+2*np.sin(math.pi*x)
     
     return(p,q,r)    

def eval_pqr2(x):

     p = -np.zeros(len(x))
     q = 2*np.ones(len(x))
     r = -math.pi**2*np.sin(math.pi*x)-2*np.sin(math.pi*x)
     
     return(p,q,r)  



def make_FDmatNeu1st(x,h,N,alpha,beta):

# evaluate coefficients of differential equation
     (p,q,r) = eval_pqr1(x)
 
# create the finite difference matrix     
     Matypp = 1/h**2*(np.diag(2*np.ones(N+1)) -np.diag(np.ones(N),-1) - 
           np.diag(np.ones(N),1))
           
     Matyp = 1/(2*h)*(np.diag(np.ones(N),1)-np.diag(np.ones(N),-1))
     
     B = np.matmul(np.diag(p,0),Matyp)
     
     A = Matypp + B+ np.diag(q)

# correct the entries to enforce
# the boundary condition
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
     (p,q,r) = eval_pqr2(x)
 
# create the finite difference matrix     
     Matypp = 1/h**2*(np.diag(2*np.ones(N+1)) -np.diag(np.ones(N),-1) - 
           np.diag(np.ones(N),1))
           
     Matyp = 1/(2*h)*(np.diag(np.ones(N),1)-np.diag(np.ones(N),-1))
     
     B = np.matmul(np.diag(p,0),Matyp)
     
     A = Matypp + B+ np.diag(q)

 # pad the matrix A with zeros on the left, right, top, and bottom
     A = np.pad(A, pad_width=1, mode='constant', constant_values=0)


# correct the entries to 
# enforce the boundary condition
     A[0,0] = -1/(2*h)
     A[0,2] = 1/(2*h)
     A[N+2, N] = 1/(2*h)
     A[N+2, N+2] = -1/(2*h)

    
     A = -A
     
# create the right hand side rhs: (N) in size
     rhs = np.pad(r,pad_width=1, mode='constant', constant_values=0)
#  update with boundary data   
     rhs[0] = alpha
     rhs[N+2] = beta
     
     
# solve for the approximate solution

     Ainv = np.linalg.inv(A)
     yapp = np.matmul(Ainv,rhs)
     # crop approx to remove our invented nodes
     yapp = yapp[1:-1]
     
     return yapp
     
driver()     
