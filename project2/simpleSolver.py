# Simple solver for 2D wave equation

from numpy import asarray, zeros, ones, random, meshgrid, linspace
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
plt.ion()
print "Starting"
Nx = 100
Ny = 100

dx = 0.1
dt = 0.05

Lx = Nx*dx
Ly = Ny*dx

q0 = 0.8 # Mean value of q
b = 0.4

u0 = ones((Nx+2, Ny+2))
um1 = ones((Nx+2, Ny+2))
u1 = zeros((Nx+2, Ny+2))

#q = q0*random.random((Nx+2, Ny+2))
q = q0*ones((Nx+2, Ny+2))


# Create figure for plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));

u0[4::, 3] = 0;


# Advance a timestep
for i in range(10):
	u1[1:Nx+1, 1:Ny+1] = \
	dt**2/(dx**2*(dt*b+2))*( \
  	q[2::, 1:Ny+1] * ( u0[2::,1:Ny+1] - u0[1:Nx+1, 1:Ny+1] ) \
	+ q[1:Nx+1, 2::] * ( u0[1:Nx+1,2::] - u0[1:Nx+1,1:Ny+1] ) \
	+ q[1:Nx+1, 0:Nx] * ( u0[1:Nx+1,0:Nx] - u0[1:Nx+1,1:Ny+1] ) \
	+ q[1:Nx+1, 1:Ny+1] * ( u0[2::,1:Ny+1] + u0[1:Nx+1, 2::] + u0[1:Nx+1,0:Ny] -4 * u0[1:Nx+1,1:Ny+1] + u0[0:Nx,1:Ny+1]) \
	+ q[0:Ny,1:Ny+1] * (-u0[1:Nx+1,1:Ny+1]+u0[0:Nx,1:Ny+1])
	+ (dx**2)/(dt**2)*((dt*b-2)*um1[1:Nx+1,1:Ny+1] + 4 * u0[1:Nx+1,1:Ny+1])
	)
	# Update ghost cells
	u1[:,0] 	= u1[:,2]
	u1[:,Ny+1]	= u1[:,Ny-1]
	u1[0,:]		= u1[2,:]
	u1[Nx+1,:]	= u1[Nx-1,:]

	um1[:,:]= u0[:,:]
	u0[:,:] = u1[:,:]

	#plt.draw()
	#ax.collections.remove(myPlot)
	#myPlot = ax.plot_wireframe(X, Y, u1[1:Nx+1,1:Ny+1])
	

ax.plot_surface(X, Y, u1[1:Nx+1,1:Ny+1])#, rstride=2, cstride=2)
plt.show()
#plt.show()


# for i in range(Ny+2):
# 	for j in range(Nx+2):
# 		print "%5.2f" % u1[i, j],
# 	print 


