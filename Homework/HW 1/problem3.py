import numpy as np
import plotly.express as px
import plotly.graph_objects as go

f = lambda x: np.exp(x)
C1 = lambda func,x0,h: (func(x0 + h) - 2 * func(x0) + func(x0-h)) / (h**2)
C2 = lambda func,x0,h: (-func(x0 + 2*h) + 16* func(x0+h) - 30*func(x0) + 16*func(x0-h) - func(x0 - 2*h)) / (12*h**2)

h = 0.1*2**(-np.linspace(0,16,17))
x0 = 7*np.pi / 8
actual = f(x0)

center1 = C1(f, x0, h)
center2 = C2(f, x0, h)

print(center1)
print(center2)
print(actual)

err1 = actual - center1
err2 = actual - center2


trace0 = go.Scatter(x=-np.log(h), y=np.log(np.abs(err1)), name='C1', yaxis='y1')
trace1 = go.Scatter(x=-np.log(h), y=np.log(np.abs(err2)), name='C2', yaxis='y1')
#trace2 = go.Scatter(x=np.log(np.arange(0, len(diff1z))), y=np.log(diff1z.flatten()), name='Z', yaxis='y1')

data = [trace0, trace1]#, trace2]
layout = go.Layout(title='Convergence of Centered Difference Equations',
                  height=600,
                  xaxis=dict(
                      title='-Log|h|'),
                  yaxis=dict(
                          title='Log|error|',
                          anchor='x'),
                  )
fig = go.Figure(data=data, layout=layout)
fig.show()

slope = lambda xarr, yarr: (yarr[1] - yarr[0]) / (xarr[1] - xarr[0])
print('Slope of C1 = ', slope(-np.log(h), np.log(np.abs(err1))))
print('Slope of C2 = ', slope(-np.log(h), np.log(np.abs(err2))))