import numpy as np 
import plotly as plt
import plotly.graph_objects as go


def eval_f(t,y):
#     evaluate f(t,y)
     
     f = y**2 / (1+t)
     return f

def eval_exact(t):
#     evaluate f(t,y)
     
     f = -1 / np.log(t+1)
     return f

def RK4(a,b,h,ya):

     N = int((b-a)/h)
     
     yapp = np.zeros(N+1)
     t = np.zeros(N+1)
     
     yapp[0] = ya
     t[0] = a

     for jj in range(1, N+1):
        tj = a+(jj-1)*h
        t[jj] = tj+h
        rk = yapp[jj-1]
        k1 = h*eval_f(tj,rk)
        k2= h*eval_f(tj+h/2,rk+0.5*k1)
        k3 = h*eval_f(tj+h/2,rk+1/2*k2)
        k4 = h*eval_f(tj+h,rk+k3)
        yapp[jj] = rk + 1/6*(k1+2*k2+2*k3+k4)

     return (yapp,t)

def explicit_euler(a,b,h,ya):

     N = int((b-a)/h)
     
     yapp = np.zeros(N+1)
     t = np.zeros(N+1)
     
     yapp[0] = ya
     t[0] = a

     for jj in range(1, N+1):
        tj = a+(jj-1)*h
        t[jj] = tj+h
        ftmp = eval_f(tj,yapp[jj-1])
        yapp[jj] = yapp[jj-1]+h*ftmp

     return (yapp,t)

def modified_euler(a,b,h,ya):

     N = int((b-a)/h)
     
     yapp = np.zeros(N+1)
     t = np.zeros(N+1)
     
     yapp[0] = ya
     t[0] = a

     for jj in range(1, N+1):
        tj = a+(jj-1)*h
        t[jj] = tj+h
        ftmp = eval_f(tj,yapp[jj-1])
        ftmp2 = eval_f(tj + h, yapp[jj-1] + h * ftmp)
        yapp[jj] = yapp[jj-1]+h / 2 * (ftmp + ftmp2)

     return (yapp,t)

y1 = -1 / np.log(2)

hE = 0.025
hME = 0.05
hRK4 = 0.1

t = np.array([1, 1.2, 1.4, 1.6, 1.8, 2])

print(' t   ', 'Explicit Euler', '      Modified Euler', '          RK4')

for b in t:
    euler, tE = explicit_euler(1,b,hE,y1)
    modifiedEuler, tME = modified_euler(1, b, hME, y1)
    RungeKutta, tRK4 = RK4(1, b, hRK4, y1)

    exact = eval_exact(b)

    print(b, ' ', np.abs(euler[-1]-exact),' ', np.abs(modifiedEuler[-1]-exact), ' ', np.abs(RungeKutta[-1]-exact))