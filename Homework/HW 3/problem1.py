import numpy as np
import plotly.express as px
import plotly.graph_objects as go

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


def eval_f(t,y):
    # example 1
    # f = y -t**2+1
    # example 2
    f = (1+t) / (1+y)
    return f

yoft = lambda t: np.sqrt(t**2 + 2*t + 6) - 1
yprime = lambda t, y: (1+t) / (1+y)
ydprime = lambda t: (1+(np.sqrt(t**2 + 2*t + 6) - 1)) - (1+t)*((1+t) / (1+(np.sqrt(t**2 + 2*t + 6) - 1))) / (1+ (np.sqrt(t**2 + 2*t + 6) - 1))**2
ddy = lambda t: -(1+t) / (1+ np.sqrt(t**2 + 2*t + 6) - 1)**2

trace0 = go.Scatter(x=np.linspace(1,2,1000), y=ydprime(np.linspace(1,2,1000)), name='X', yaxis='y1')

data = [trace0]
layout = go.Layout(title='Y`` over [1,2]',
                  height=600,
                  xaxis=dict(
                      title='x'),
                  yaxis=dict(
                          title='y``',
                          anchor='x'),
                  )
fig = go.Figure(data=data, layout=layout)
fig.show()

trace0 = go.Scatter(x=np.linspace(1,2,1000), y=ddy(np.linspace(1,2,1000)), name='X', yaxis='y1')

data = [trace0]
layout = go.Layout(title='d/dy(f(t,y)) over [1,2]',
                  height=600,
                  xaxis=dict(
                      title='x'),
                  yaxis=dict(
                          title='d/dy(f(t,y))',
                          anchor='x'),
                  )
fig = go.Figure(data=data, layout=layout)
fig.show()
L = max(ddy(np.linspace(1,2,1000)))
print('L = ', max(ddy(np.linspace(1,2,1000))))
M = ydprime(2)
print('M = ', ydprime(2))
h = 0.5

errBound1 = h * M / L * (np.exp( L *(1.5-1)) -1)
errBound2 = h * M / L * (np.exp( L *(2-1)) -1)
print('Error bound at t=1.5: ', errBound1)
print('Error bound at t=2: ', errBound2)
print('y(1.5) = ', yoft(1.5))
print('y(2) = ', yoft(2))



a = 1
b = 2
ya = 2
yapp,t = explicit_euler(a,b,h,ya)

print('Explicit Euler approximation at t=1.5: ', yapp[1])
print('Explicit Euler approximation at t=2: ', yapp[2])

print('Error at t=1.5: ', yoft(1.5) - yapp[1])
print('Error at t=1.5: ', yoft(2) - yapp[2])