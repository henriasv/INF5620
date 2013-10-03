# Simple solver for 2D wave equation

from numpy import asarray, zeros, ones, random, meshgrid, linspace, exp
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D


plt.ion()



def update_ghost_cells(u):
	u[1:-1,0] 		= u[1:-1,2]
	u[1:-1,-1]	= u[1:-1,-3]
	u[0,1:-1]		= u[2,1:-1]
	u[-1,1:-1]	= u[-3,1:-1]

def advance_first_step(u0, u1, q, V, b, Nx, Ny, dx, dt):
	update_ghost_cells(u0)
	u1[1:-1, 1:-1] = \
	dt**2/(4*dx**2)*( \
  	q[2::, 1:-1] * ( u0[2::,1:-1] - u0[1:-1, 1:-1] ) \
	+ q[1:-1, 2::] * ( u0[1:-1,2::] - u0[1:-1,1:-1] ) \
	+ q[1:-1, 0:-2] * ( u0[1:-1,0:-2] - u0[1:-1,1:-1] ) \
	+ q[1:-1, 1:-1] * ( u0[2::,1:-1] + u0[1:-1, 2::] + u0[1:-1,0:-2] -4 * u0[1:-1,1:-1] + u0[0:-2,1:-1]) \
	+ q[0:-2,1:-1] * (-u0[1:-1,1:-1]+u0[0:-2,1:-1])) \
	+ u0[1:-1,1:-1] + V[1:-1, 1:-1] *dt - 0.5*b*V[1:-1,1:-1]*dt**2
	update_ghost_cells(u1)

Nx = 50
Ny = 50

dx = 0.04
dt = 0.02
T = 20

Lx = float((Nx-1)*dx)
Ly = float((Ny-1)*dx)
print "Lx: %.6f, Ly: %.6f" %(Lx, Ly)
q0 = 3 # Mean value of q
b =.2;

u0 = ones((Nx+2, Ny+2))
um1 = ones((Nx+2, Ny+2))
u1 = zeros((Nx+2, Ny+2))

q = q0*random.random((Nx+2, Ny+2))
#q = q0*ones((Nx+2, Ny+2))

# Create figure for plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set initial condition
X, Y = meshgrid(linspace(-dx,Lx+dx, Nx+2), linspace(-dx,Ly+dx, Ny+2))
I = 3*exp(-20*((X-Lx/2)**2+(Y-Ly/2)**2))
I = I-np.mean(I[1:-1,1:-1])
I = ones((Nx+2, Nx+2))*1.01
V = zeros((Nx+2, Nx+2))
#V = 0.3*np.sin(X-Lx/2+Y-Ly/2)
print "Mean(X) = %.3f" % np.mean(X)
print "Sum over initial velocities: %.6f" % np.sum(np.sum(V[1:-1,1:-1])) 
print "Sum over initial u-values: %.6f" % np.sum(np.sum(I[1:-1,1:-1]))
um1[:, :] = I[:,:]


advance_first_step(um1, u0, q, V, b, Nx, Ny, dx, dt)


vmin_ = X.min()
vmax_ = X.max()




myPlot = ax.plot_wireframe(X, Y, u0)
ax.set_zlim3d(vmin_, vmax_)
ax.autoscale_view(tight=None, scalex=True, scaley=True, scalez=False)

for i in range(int(float(T)/dt)):
	u1[1:-1, 1:-1] = \
	dt**2/(dx**2*(dt*b+2))*( \
  	q[2::, 1:-1] * ( u0[2::,1:-1] - u0[1:-1, 1:-1] ) \
	+ q[1:-1, 2::] * ( u0[1:-1,2::] - u0[1:-1,1:-1] ) \
	+ q[1:-1, 0:-2] * ( u0[1:-1,0:-2] - u0[1:-1,1:-1] ) \
	+ q[1:-1, 1:-1] * ( u0[2::,1:-1] + u0[1:-1, 2::] + u0[1:-1,0:-2] - 4 * u0[1:-1,1:-1] + u0[0:-2,1:-1]) \
	+ q[0:-2,1:-1] * (-u0[1:-1,1:-1]+u0[0:-2,1:-1])) \
	+ (1.0/(b*dt +2)) * ((dt*b-2)*um1[1:-1,1:-1] + 4 * u0[1:-1,1:-1])

	um1[:,:]= u0[:,:]
	u0[:,:] = u1[:,:]

	update_ghost_cells(u0);

	plt.draw()
	ax.collections.remove(myPlot)
	myPlot = ax.plot_wireframe(X, Y, u0)
	ax.set_zlim3d(vmin_, vmax_)
	
	print np.sum(np.sum(u0[1:-1, 1:-1]));

raw_input("press enter")

def advance():
	a = 1;


