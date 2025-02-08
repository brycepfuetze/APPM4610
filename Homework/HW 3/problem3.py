import numpy as np
import plotly.express as px
import plotly.graph_objects as go


def Taylor_2nd(a,b,h,ya):
    N = int((b-a)/h)
    yapp = np.zeros(N+1)
    t = np.zeros(N+1)
    yapp[0] = ya
    t[0] = a
    for jj in range(1, N+1):
        tj = a+(jj-1)*h
        t[jj] = tj+h
        ftmp = eval_f(tj,yapp[jj-1])
        ft_tmp = eval_ft(tj,yapp[jj-1])
        yapp[jj] = yapp[jj-1]+h*(ftmp + h/2*ft_tmp)
    return (yapp,t)

def eval_f(t,y):
    # evals the right hand side of the DE
    f = 1/(t**2) - y/t -y**2
    return f

def eval_ft(t,y):
    # evals Df/Dt = total derivative of f wrt t
    # simplify to something stable
    ft = -2*t**-3 + y*t**-2 + (-1/t -2*y)*(t**-2 -y/t -y**2)
    return ft

a = 1
b = 2
h = 0.05
ya = -1

yapp,t = Taylor_2nd(a,b,h,ya)

yactual = lambda t: -1/t

trace0 = go.Scatter(x=t, y=yapp, name='Approximation', yaxis='y1')
trace1 = go.Scatter(x=np.linspace(1,2,1000), y=yactual(np.linspace(1,2,1000)), name='Actual', yaxis='y1')

data = [trace0, trace1]
layout = go.Layout(title='2nd Order Taylor Approximation over [1,2]',
                  height=600,
                  xaxis=dict(
                      title='t'),
                  yaxis=dict(
                          title='y',
                          anchor='x'),
                  )
fig = go.Figure(data=data, layout=layout)
fig.show()

trace0 = go.Scatter(x=t, y=yapp - yactual(t), name='Error', yaxis='y1')

data = [trace0]
layout = go.Layout(title='Error of 2nd Order Taylor Approximation over [1,2]',
                  height=600,
                  xaxis=dict(
                      title='t'),
                  yaxis=dict(
                          title='Error',
                          anchor='x'),
                  )
fig = go.Figure(data=data, layout=layout)
fig.show()

interpolation = np.interp([1.052, 1.555, 1.978], t,yapp)

print('Interpolated y(1.052) = ', interpolation[0])
print('Interpolated y(1.555) = ', interpolation[1])
print('Interpolated y(1.978) = ', interpolation[2])

print('Actual y(1.052) = ', yactual(1.052))
print('Actual y(1.555) = ', yactual(1.555))
print('Actual y(1.978) = ', yactual(1.978))

print('Error at t=1.052: ', yactual(1.052) - interpolation[0])
print('Error at t=1.555: ', yactual(1.555) - interpolation[1])
print('Error at t=1.978: ', yactual(1.978) - interpolation[2])