# Simple solver for 2D wave equation

from numpy import asarray, zeros, ones, random, meshgrid, linspace, exp
import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import nose.tools as nt


plt.ion()


class user_action:
	def initialize(self, X, Y, initial_data=None, vmin=None, vmax=None):
		pass
	def __call__(self, data2D):
		pass

class wave_plotter(user_action):
	def initialize(self, X, Y, initial_data=None, vmin=None, vmax=None):
		"""
		Create objects to plot on
		"""
		if initial_data == None:
			initial_data = zeros(shape(X));
			if vmin == None:
				self.vmin = -1.0
			if vmax == None:
				self.vmax = 1.0
		else:
			if vmin == None:
				self.vmin = np.min(np.min(initial_data))
			if vmax == None:
				self.vmax = np.max(np.max(initial_data))

		self.X = X; self.Y = Y
		# Create figure for plotting
		self.fig = plt.figure()
		self.ax = self.fig.add_subplot(111, projection='3d')
		self.myPlot = self.ax.plot_wireframe(X, Y, initial_data)
		self.ax.set_zlim3d(self.vmin, self.vmax)
		self.ax.autoscale_view(tight=None, scalex=True, scaley=True, scalez=False)

	def __call__(self, data2D):
		plt.draw()
		self.ax.collections.remove(self.myPlot)
		self.myPlot = self.ax.plot_wireframe(self.X, self.Y, data2D)
		self.ax.set_zlim3d(self.vmin, self.vmax)

class user_action_test_constant(user_action):
	def __init__(self, tolerance=1e-12, constant=0):
		self.tolerance = 1e-12
		self.constant = constant
	def __call__(self, data2D):
		nt.assert_almost_equal(np.max(np.max(data2D)), self.constant, self.tolerance)

def update_ghost_points(u):
	"""
	This function is designed to apply the neumann boundary condition if du/dn = 0 and the spatial integration scheme is centered
	"""
	u[1:-1,0] 		= u[1:-1,2]
	u[1:-1,-1]	= u[1:-1,-3]
	u[0,1:-1]		= u[2,1:-1]
	u[-1,1:-1]	= u[-3,1:-1]

def advance_first_step(u0, u1, q, V, b, dx, dt):
	update_ghost_points(u0)
	u1[1:-1, 1:-1] = \
	dt**2/(4*dx**2)*( \
  	q[2::, 1:-1] * ( u0[2::,1:-1] - u0[1:-1, 1:-1] ) \
	+ q[1:-1, 2::] * ( u0[1:-1,2::] - u0[1:-1,1:-1] ) \
	+ q[1:-1, 0:-2] * ( u0[1:-1,0:-2] - u0[1:-1,1:-1] ) \
	+ q[1:-1, 1:-1] * ( u0[2::,1:-1] + u0[1:-1, 2::] + u0[1:-1,0:-2] -4 * u0[1:-1,1:-1] + u0[0:-2,1:-1]) \
	+ q[0:-2,1:-1] * (-u0[1:-1,1:-1]+u0[0:-2,1:-1])) \
	+ u0[1:-1,1:-1] + V[1:-1, 1:-1] *dt - 0.5*b*V[1:-1,1:-1]*dt**2
	update_ghost_points(u1)

def solver(f_I, f_V, f, f_q, b, L, N, dt, T, user_action=None):
	"""
	Solver will only give quadratic domains. User should in principle not see that the implementation uses ghost cells!
	"""	
	dx = L/(N-1)
	x = linspace(-dx, L+dx, N+2)
	X, Y = meshgrid(x, x)
	um1 = zeros((N+2, N+2))
	u0 = zeros((N+2, N+2))
	u1 = zeros((N+2, N+2))
	um1 = f_I(X, Y, L)
	V = f_V(X, Y, L)
	q = f_q(X, Y, L)

	update_ghost_points(um1)

	advance_first_step(um1, u0, q, V, b, dx, dt)
	user_action.initialize(X, Y, u0)
	for i in range(int(float(T)/dt)):
		advance(um1, u0, u1, q, b, dx, dt)

		if user_action:
			user_action(u0)


def advance(um1, u0, u1, q, b, dx, dt):
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

	update_ghost_points(u0);

if __name__=="__main__":
	N = 100
	L = 1.1
	T = 10
	dx = float(L)/(N-1)
	dt = dx/2.5
	q = ones((N, N))*0.8
	f_I = lambda X, Y, L: 3*exp(-100*((X-L/2)**2+(Y-L/2)**2))
	f_V = lambda X, Y, L: zeros(np.shape(X))
	f = lambda x: zeros(np.shape(X))
	q = lambda X, Y, L: ones(np.shape(X))
	b = 1.2

	myPlotter = wave_plotter()
	solver(f_I, f_V, f, q, b, L, N, dt, T, user_action = myPlotter)


def test_constant_solution():
	c = 1.45
	N = 100
	L = 1.4
	T = 1
	dx = float(L)/(N-1)
	dt = dx/2.01
	q = lambda X, Y, L: zeros(np.shape(X))+0.8
	f_I = lambda X, Y, L: zeros(np.shape(X))+c
	f_V = lambda X, Y, L: zeros(np.shape(X))
	f = lambda x: zeros(np.shape(X))
	b = 0
	tester = user_action_test_constant(constant=c)
	solver(f_I, f_V, f, q, b, L, N, dt, T, user_action = tester)


