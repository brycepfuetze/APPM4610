import matplotlib.pyplot as plt
import numpy as np

from matplotlib import cm
from matplotlib.ticker import LinearLocator

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})


# Make data.
kappa = 2.37 # W/ cm K
rho = 2.70 # g/cm^3
c = 0.897 # J / g K
l = 10 # cm
t_end = 20 # s
resolution = 0.01

x = np.arange(0, l, resolution)
t = np.arange(0, t_end, resolution)

x, t = np.meshgrid(x, t)
theta = 1 * np.pi / l

u = np.exp(-kappa * theta **2 *t / (rho * c)) * np.sin(theta * x)

# Plot the surface.
surf = ax.plot_surface(x, t, u, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(0, 1.1)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

ax.set_xlabel('x position [cm]')
ax.set_ylabel('time [s]')
ax.set_zlabel('Temperature [K]')
ax.set_title('Temperature of 10cm bar over 20 seconds')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()