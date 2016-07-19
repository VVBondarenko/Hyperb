#!/usr/bin/python
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')

data = np.genfromtxt('output_LW')
x = data[:,0]
y = data[:,1]
z = data[:,2]

xi = np.linspace(min(x), max(x))
yi = np.linspace(min(y), max(y))

X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z, xi, yi)

surf = ax.plot_surface(X, Y, Z, rstride=5, cstride=5, cmap=cm.jet,
                       linewidth=1, antialiased=True)

ax.set_zlim3d(np.min(Z), np.max(Z))
fig.colorbar(surf)

plt.show()