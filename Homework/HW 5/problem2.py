import numpy as np 
import plotly as plt
import plotly.graph_objects as go
from scipy.integrate import solve_ivp

# System equation from Problem 1a
def ueval(t,u1, u2):
     return 3 * u2 - 2*u1 + 6*np.exp(-t)

# Exact solution to system
def eval_exact(t):
    return 2 * np.exp(2*t) - np.exp(t) + np.exp(-t)

# modified from Canvas
def euler_sys(a,b,h,u10,u20):

     N = int((b-a)/h)
     
     u1 = np.zeros(N+1)
     u2 = np.zeros(N+1)
     t = np.zeros(N+1)
     
     u1[0] = u10
     u2[0] = u20

     for jj in range(1, N+1):
        tj = a+(jj-1)*h
        t[jj] = tj+h

        u1[jj] = u1[jj-1] + h * u2[jj-1]
        utmp = ueval(tj, u1[jj-1], u2[jj-1])
        u2[jj] = u2[jj-1] + h * (utmp)

     return (u1,u2,t)

# created with help from Copilot!
def runge_kutta_sys(a, b, h, u10, u20):
    N = int((b - a) / h)
    
    u1 = np.zeros(N + 1)
    u2 = np.zeros(N + 1)
    t = np.zeros(N + 1)
    
    u1[0] = u10
    u2[0] = u20
    t[0] = a
    
    for jj in range(1, N + 1):
        tj = a + (jj - 1) * h
        t[jj] = tj + h
        
        # we're can make this better, but I'll just be doing RK4 twice
        # cause that's easy
        k1_u1 = h * u2[jj - 1]
        k1_u2 = h * ueval(tj, u1[jj - 1], u2[jj - 1])
        
        k2_u1 = h * (u2[jj - 1] + 0.5 * k1_u2)
        k2_u2 = h * ueval(tj + 0.5 * h, u1[jj - 1] + 0.5 * k1_u1, u2[jj - 1] + 0.5 * k1_u2)
        
        k3_u1 = h * (u2[jj - 1] + 0.5 * k2_u2)
        k3_u2 = h * ueval(tj + 0.5 * h, u1[jj - 1] + 0.5 * k2_u1, u2[jj - 1] + 0.5 * k2_u2)
        
        k4_u1 = h * (u2[jj - 1] + k3_u2)
        k4_u2 = h * ueval(tj + h, u1[jj - 1] + k3_u1, u2[jj - 1] + k3_u2)
        
        # use right RK4 weightings for each
        u1[jj] = u1[jj - 1] + (k1_u1 + 2 * k2_u1 + 2 * k3_u1 + k4_u1) / 6
        u2[jj] = u2[jj - 1] + (k1_u2 + 2 * k2_u2 + 2 * k3_u2 + k4_u2) / 6
    
    return u1, u2, t


# set values
h = 0.1 
u10 = 2 # u1 initial
u20 = 2 # u2 initial
a=0
b=1


u1, u2, t = euler_sys(a,b,h,u10, u20)

# print(u1)
# print(u2)

tExact = np.linspace(0,1,11)
yexact = eval_exact(tExact)

# print(yexact)

# Plot!
fig = go.Figure()
fig.add_trace(go.Scatter(x=t, y=u1-yexact, mode='lines+markers', name='Euler Approximation'))
#fig.add_trace(go.Scatter(x=tExact, y=yexact, name='Exact Solution'))

fig.update_layout(title='System Euler Method Error',
                  xaxis_title='t',
                  yaxis_title='y')

#fig.show()

# now use RK4!
rku1, rku2, rkt = runge_kutta_sys(a, b, h, u10, u20)


fig = go.Figure()
fig.add_trace(go.Scatter(x=t, y=u1-yexact, mode='lines+markers', name='Euler Approximation'))

fig.add_trace(go.Scatter(x=rkt, y=rku1-yexact, mode='lines+markers', name='RK Approximation'))
#fig.add_trace(go.Scatter(x=tExact, y=yexact, name='Exact Solution'))


fig.update_layout(title='System Method Comparison',
                  xaxis_title='t',
                  yaxis_title='Error')

#fig.show()

# print(rku1)
# print(rku2)
# print(u1)
# print(u2)
# print(yexact)

# to use the Scipy function we need to redefine our U function to return
# everything in a vector. 

def ueval(t, u):
    u1, u2 = u
    return [u2, 3 * u2 - 2 * u1 + 6 * np.exp(-t)]

# can specify teval to line up with others
#t_eval = np.arange(a, b + h, h)
solution = solve_ivp(ueval, (a, b), [u10, u20]) #, t_eval = t_eval)

print(solution)

print('Final errors: ')
print('Euler: ', np.abs(u1[-1] - eval_exact(b)))
print('RK4: ', np.abs(rku1[-1] - eval_exact(b)))
print('Scipy: ', np.abs(solution.y[0][-1] - eval_exact(b)))