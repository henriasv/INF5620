# Simple solver for 2D wave equation

from numpy import asarray, zeros, ones, random, meshgrid, linspace, exp
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D


plt.ion()



def update_ghost_cells(u):
	u[:,0] 		= u[:,2]
	u[:,-1]	= u[:,-3]
	u[0,:]		= u[2,:]
	u[-1,:]	= u[-3,:]

def advance_first_step(u0, u1, q, V, b, Nx, Ny, dx, dt):
	u1[1:Nx+1, 1:Ny+1] = \
	dt**2/(4*dx**2)*( \
  	q[2::, 1:Ny+1] * ( u0[2::,1:Ny+1] - u0[1:Nx+1, 1:Ny+1] ) \
	+ q[1:Nx+1, 2::] * ( u0[1:Nx+1,2::] - u0[1:Nx+1,1:Ny+1] ) \
	+ q[1:Nx+1, 0:Nx] * ( u0[1:Nx+1,0:Nx] - u0[1:Nx+1,1:Ny+1] ) \
	+ q[1:Nx+1, 1:Ny+1] * ( u0[2::,1:Ny+1] + u0[1:Nx+1, 2::] + u0[1:Nx+1,0:Ny] -4 * u0[1:Nx+1,1:Ny+1] + u0[0:Nx,1:Ny+1]) \
	+ q[0:Ny,1:Ny+1] * (-u0[1:Nx+1,1:Ny+1]+u0[0:Nx,1:Ny+1])) \
	+ u0[1:Nx+1,1:Ny+1] + V[1:Nx+1, 1:Ny+1] *dt - 2*b*V[1:Nx+1,1:Ny+1]*dt**2

Nx = 100
Ny = 100

dx = 0.1
dt = 0.04
T = 10

Lx = float(Nx*dx)
Ly = float(Ny*dx)

q0 = 1.0 # Mean value of q
b = 0.1;

u0 = ones((Nx+2, Ny+2))
um1 = ones((Nx+2, Ny+2))
u1 = zeros((Nx+2, Ny+2))

#q = q0*random.random((Nx+2, Ny+2))
q = q0*ones((Nx+2, Ny+2))

# Create figure for plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set initial condition
X, Y = meshgrid(linspace(0,Lx, Nx+2), linspace(0,Ly, Ny+2))
I = 0.1*ones((Nx+2, Ny+2))+exp(-(X-Lx/2)**2-(Y-Ly/2)**2)
V = zeros((Nx+2, Ny+2))
um1[:, :] = I[:, :]
advance_first_step(um1, u0, q, V, b, Nx, Ny, dx, dt)
#u0[:, :] = um1[:, :]
vmin_ = um1.min()
vmax_ = um1.max()

update_ghost_cells(u1)
update_ghost_cells(u0)
update_ghost_cells(um1)

myPlot = ax.plot_wireframe(X, Y, u0)
ax.set_zlim3d(vmin_, vmax_)
ax.autoscale_view(tight=None, scalex=True, scaley=True, scalez=False)
#from mayavi import mlab
#s = mlab.surf(u1, warp_scale="auto")
# Timeloop
raw_input("press enter")
for i in range(int(float(T)/dt)):
	u1[1:Nx+1, 1:Ny+1] = \
	dt**2/(dx**2*(dt*b+2))*( \
  	q[2::, 1:Ny+1] * ( u0[2::,1:Ny+1] - u0[1:Nx+1, 1:Ny+1] ) \
	+ q[1:Nx+1, 2::] * ( u0[1:Nx+1,2::] - u0[1:Nx+1,1:Ny+1] ) \
	+ q[1:Nx+1, 0:Nx] * ( u0[1:Nx+1,0:Nx] - u0[1:Nx+1,1:Ny+1] ) \
	+ q[1:Nx+1, 1:Ny+1] * ( u0[2::,1:Ny+1] + u0[1:Nx+1, 2::] + u0[1:Nx+1,0:Ny] -4 * u0[1:Nx+1,1:Ny+1] + u0[0:Nx,1:Ny+1]) \
	+ q[0:Ny,1:Ny+1] * (-u0[1:Nx+1,1:Ny+1]+u0[0:Nx,1:Ny+1]) \
	+ (dx**2)/(dt**2) * ((dt*b-2)*um1[1:Nx+1,1:Ny+1] + 4 * u0[1:Nx+1,1:Ny+1]) \
	)

	#s.mlab_source.scalars = u1[1:Nx+1,1:Ny+1];

	update_ghost_cells(u1);
	update_ghost_cells(u0);

	um1[:,:]= u0[:,:]
	u0[:,:] = u1[:,:]


	plt.draw()
	ax.collections.remove(myPlot)
	myPlot = ax.plot_wireframe(X, Y, u0)
	ax.set_zlim3d(vmin_, vmax_)
	#raw_input("press enter")
	print np.sum(np.sum(u0[1:-2, 1:-2]));



def advance():
	a = 1;
#s = mlab.figure();
#mlab.pipeline.iso_surface(s);

#s.draw()
#mlab.show()
#raw_input("press enter")

# for i in range(Ny+2):
# 	for j in range(Nx+2):
# 		print "%5.2f" % u1[i, j],
# 	print 


