# Ser paa amplikasjonsfaktoren til theta-regelen for decay ODE

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def A(theta, adt):
	return (1.0 - (1.0-theta)*adt)/(1.0 + theta*adt)

def A_analytic(adt):
	return numpy.exp(-adt)

def theta_analytic(adt):
	return (adt-1+numpy.exp(-adt))/(adt-numpy.exp(-adt)*adt)

N = 100
theta = numpy.linspace(0, 1, N)
adt = numpy.linspace(0, 3, N)
X, Y = numpy.meshgrid(theta, adt)
A_mesh = A(X, Y)
theta_optimal = theta_analytic(adt)

print A_mesh


"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, A_mesh)
plt.hold('on')
ax.plot_surface(X, Y, numpy.zeros((len(theta), len(adt))))
plt.xlabel(r"$\theta$")
plt.show
"""
# or pcolor
fig = plt.figure()
plt.pcolor(X, Y, A_mesh, vmin=-1, vmax=1)

plt.colorbar()
plt.hold('on')
contour = plt.contour(X, Y, A_mesh)
plt.plot(theta_optimal, adt, '--', label='Optimal line')
plt.legend()
plt.xlabel(r'$\theta$')
plt.ylabel(r'$k\Delta t$')
plt.clabel(contour)
plt.show()