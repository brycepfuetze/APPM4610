import numpy as np 
import plotly as plt
import plotly.graph_objects as go


def eval_f(t,y):
#     evaluate f(t,y)
     
     f = 1-y
     return f

def diffMethod(a,b,h,ya):
    N = int((b-a)/h)
    yapp = np.zeros(N+1)
    t = np.zeros(N+1)

    yapp[0] = ya[0]
    yapp[1] = ya[1]
    t[0] = a
    t[1] = a+h

    # now time integration with Adams-Bashforth
    for j in range(2,N+1):
        t[j] = a + j*h
        yapp[j] = 4*yapp[j-1] - 3 * yapp[j-2] -2*h* eval_f(t[j], yapp[j-2])

    return yapp, t

# starting y val
y0 = 0
y1 = 1 - np.exp(-0.1)

# set h
h = 0.1

# end points to try
a = 0
b = 1


yapp, t = diffMethod(a,b,h,[y0,y1])

print('t = ', t)
print('yapp = ', yapp)


# Create the plot
fig = go.Figure()
fig.add_trace(go.Scatter(x=t, y=yapp, mode='lines+markers', name='yapp'))

# Add titles and labels
fig.update_layout(title='Approximation with h=0.1',
                  xaxis_title='t',
                  yaxis_title='y')

# Show the plot
fig.show()

# try again with new h and yval
y1 = 1 - np.exp(-0.01)
h = 0.01

yappzero, tzero = diffMethod(a,b,h,[y0,y1])

# Create the plot
fig = go.Figure()
fig.add_trace(go.Scatter(x=t, y=yapp, mode='lines+markers', name='h = 0.1'))
# Add the second trace
fig.add_trace(go.Scatter(x=tzero, y=yappzero, mode='lines+markers', name='h = 0.01'))

# Add titles and labels
fig.update_layout(title='Difference Method Comparison',
                  xaxis_title='t',
                  yaxis_title='y')

# Show the plot
fig.show()