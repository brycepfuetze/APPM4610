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

hVals = 10.0 ** np.arange(-5, -11, -1)

a = 1
b = 1+1e-5
ya = 2

yapps = []
yappsFinal = []
ts = []

for h in hVals:
    yapp,t = explicit_euler(a,b,h,ya)
    yapps.append(yapp)
    ts.append(t)
    yappsFinal.append(yapp[-1])

print (yapps)
print(ts)
print('finals', yappsFinal)

yappsFinalArr = np.array(yappsFinal)

t = 1+1e-5 
yActual = np.sqrt(t**2 + 2*t +6) - 1

error = yappsFinalArr - yActual

print(np.log(hVals))
print(np.log(np.abs(error)))

trace0 = go.Scatter(x=np.log10(hVals), y=np.log10(np.abs(error)), name='Error', yaxis='y1')

data = [trace0]
layout = go.Layout(title='Error of Explicit Euler with Decreasing h',
                  height=600,
                  xaxis=dict(
                      title='Log|h|'),
                  yaxis=dict(
                          title='Log|error|',
                          anchor='x'),
                  )
fig = go.Figure(data=data, layout=layout)
fig.show()
