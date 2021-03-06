# Simple solver for 2D wave equation

from numpy import asarray, zeros, ones, random, meshgrid, linspace, exp
import numpy as np
import sys
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import nose.tools as nt
import time
from tvtk.api import tvtk
import copy


plt.ion()


class user_action:
	def initialize(self, X, Y, initial_data=None, vmin=None, vmax=None):
		pass
	def __call__(self, data2D, t):
		pass
	def get(self):
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

	def __call__(self, data2D, t):
		plt.draw()
		self.ax.collections.remove(self.myPlot)
		self.myPlot = self.ax.plot_wireframe(self.X, self.Y, data2D)
		self.ax.set_zlim3d(self.vmin, self.vmax)

class wave_plotter_tvtk(user_action):
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
		x_min = np.min(np.min(X))
		y_min = np.min(np.min(Y))
		x_max = np.max(np.max(X))
		y_max = np.max(np.max(Y))
		x_middle = (x_min+x_max)/2
		y_middle = (y_min+y_max)/2
		z_middle = np.mean(np.mean(initial_data))
		L = x_max-x_min;
		diffs = np.shape(X)
		x_diff = diffs[0]
		y_diff = diffs[1]
		z_diff = 1.0#self.vmax-self.vmin

		self.tvtk = tvtk;
		self.sp = tvtk.StructuredPoints(origin = (x_middle, y_middle, z_middle),
										dimensions = (x_diff, y_diff, 1), 
										spacing = (2*L/(x_diff-1), 2*L/(y_diff-1), 100.0))

		self.z = np.transpose(initial_data).flatten()
		
		self.sp.point_data.scalars = self.z
		self.geom_filter = tvtk.ImageDataGeometryFilter(input=self.sp)
		self.warp = tvtk.WarpScalar(input=self.geom_filter.output)
		self.normals = tvtk.PolyDataNormals(input=self.warp.output)

		# The rest of the VTK pipeline.
		self.m = tvtk.PolyDataMapper(input=self.normals.output,
		                        scalar_range=(self.vmin, self.vmax))
		p = tvtk.Property(opacity=0.5, color=(1, 1, 1), representation="s")
		self.a = tvtk.Actor(mapper=self.m, property=p)

		self.ren = tvtk.Renderer(background=(0.0, 0.0, 0.0))
		
		self.ren.add_actor(self.a)

		# Get a nice view.
		self.cam = self.ren.active_camera
		self.cam.azimuth(-50)
		self.cam.roll(90)

		# Create a RenderWindow, add the renderer and set its size.
		self.rw = tvtk.RenderWindow(size=(800, 800))
		self.rw.add_renderer(self.ren)

		# Create the RenderWindowInteractor
		self.rwi = tvtk.RenderWindowInteractor(render_window=self.rw)

		self.rwi.initialize()
		self.ren.reset_camera()
		self.rwi.render()
		#self.rwi.start()
		
	def add_static_background(self, data2D):
		
		self.sp2 = copy.copy(self.sp)
		self.z2 = np.reshape(np.transpose(data2D), (-1,))
		
		self.sp2.point_data.scalars = self.z2
		self.geom_filter2 = tvtk.ImageDataGeometryFilter(input=self.sp2)
		self.warp2 = tvtk.WarpScalar(input=self.geom_filter2.output)
		self.normals2 = tvtk.PolyDataNormals(input=self.warp2.output)

		# The rest of the VTK pipeline.
		self.m2 = tvtk.PolyDataMapper(input=self.normals2.output)#,scalar_range=(self.vmin, self.vmax))
		self.p = tvtk.Property(opacity = 1.0, color=(0.5, 0.5, 0.5), representation="w")
		self.a2 = tvtk.Actor(mapper=self.m2, property=self.p)

		self.ren.add_actor(self.a2)

		self.rwi.initialize()
		self.ren.reset_camera()
		#self.rwi.render()
		


	def __call__(self, data2D, t):
		self.z[:] = np.reshape(np.transpose(data2D), (-1,))
		self.m.scalar_range = min(self.z), max(self.z)
		self.sp.modified()
		self.rwi.render()
		#time.sleep(0.1)

	def get(self):
		self.rwi.start()


class user_action_test_constant(user_action):
	def __init__(self, tolerance=1e-12, constant=0):
		self.tolerance = 1e-12
		self.constant = constant
	def __call__(self, data2D, t):
		nt.assert_almost_equal(np.max(np.max(data2D)), self.constant, self.tolerance)

class user_action_test_symmetric(user_action):
	def __init__(self, tolerance=1e-12, axis='x'):
		self.axis = axis
		self.tolerance = tolerance

	def __call__(self, data2D, t):
		Nx, Ny = np.shape(data2D)
		if self.axis == 'x':
			tmp_left = data2D[0:Nx/2]
			tmp_right = np.fliplr(data2D[:,:])[0:Nx/2]
			tmp_diff = tmp_left-tmp_right;
			nt.assert_almost_equal(np.max(np.max(tmp_diff)), self.tolerance)
		elif self.axis == 'y':
			pass
		else:
			nt.assert_almost_equal(0, 1, 1e-12)

class user_action_convergence_max_error(user_action):
	def __init__(self, tolerance=1e-12):
		self.tolerance = tolerance
		self.exact_function = None
		self.errors = []
	def __call__(self, data2D, t):
		if self.exact_function==None:
			raise Exception('Exact solution function not provided!, Cannot calculate true error')
		diff = data2D-self.exact_function(self.X, self.Y, t)
		self.errors.append(np.max(np.max(np.abs(diff))))

	def initialize(self, X, Y, initial_data=None, vmin=None, vmax=None):
		self.X = X; self.Y = Y

	def set_exact_solution(self, sol):
		self.exact_function = sol

	def get(self):
		return np.asarray(self.errors)

	def reset(self):
		self.errors = []


def update_ghost_points(u):
	"""
	This function is designed to apply the neumann boundary condition if du/dn = 0 and the spatial integration scheme is centered
	"""
	u[1:-1,0] 	= u[1:-1,2]
	u[1:-1,-1]	= u[1:-1,-3]
	u[0,1:-1]	= u[2,1:-1]
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
	t = 0
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
	t = dt
	if not user_action == None:
		user_action.initialize(X[1:-1, 1:-1], Y[1:-1, 1:-1], u0[1:-1, 1:-1])
		try: 
			user_action.add_static_background(-q[1:-1, 1:-1])
		except AttributeError as e:
			print(type(e))
			raise

	for i in range(int(float(T)/dt)):
		advance(um1, u0, u1, q, b, dx, dt)
		t += dt

		if not user_action == None:
			user_action(u0[1:-1, 1:-1], t)


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

def test_constant_solution():
	constant = 1.45
	c = 0.8
	N = 100
	L = 1.4
	T = 1
	dx = float(L)/(N-1)
	dt = dx/2.01
	q = lambda X, Y, L: zeros(np.shape(X))+c**2
	f_I = lambda X, Y, L: zeros(np.shape(X))+constant
	f_V = lambda X, Y, L: zeros(np.shape(X))
	f = lambda x: zeros(np.shape(X))
	b = 0
	tester = user_action_test_constant(constant=constant)
	solver(f_I, f_V, f, q, b, L, N, dt, T, user_action = tester)	

def advance_mistake_1(um1, u0, u1, q, b, dx, dt):
	"""
	(dt*b+2) left out in the first fraction
	Passes  constant test
	"""
	u1[1:-1, 1:-1] = \
	dt**2/(dx**2)*( \
  	q[2::, 1:-1] * ( u0[2::,1:-1] - u0[1:-1, 1:-1] ) \
	+ q[1:-1, 2::] * ( u0[1:-1,2::] - u0[1:-1,1:-1] ) \
	+ q[1:-1, 0:-2] * ( u0[1:-1,0:-2] - u0[1:-1,1:-1] ) \
	+ q[1:-1, 1:-1] * ( u0[2::,1:-1] + u0[1:-1, 2::] + u0[1:-1,0:-2] - 4 * u0[1:-1,1:-1] + u0[0:-2,1:-1]) \
	+ q[0:-2,1:-1] * (-u0[1:-1,1:-1]+u0[0:-2,1:-1])) \
	+ (1.0/(b*dt +2)) * ((dt*b-2)*um1[1:-1,1:-1] + 4 * u0[1:-1,1:-1])

	um1[:,:]= u0[:,:]
	u0[:,:] = u1[:,:]

	update_ghost_points(u0);

def advance_mistake_2(um1, u0, u1, q, b, dx, dt):
	"""
	Wrong indices on u1
	Fails constant test, but passes for c = 0
	"""
	u1[0:-2, 0:-2] = \
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

def advance_mistake_3(um1, u0, u1, q, b, dx, dt):
	pass
def advance_mistake_4(um1, u0, u1, q, b, dx, dt):
	pass
def advance_mistake_5(um1, u0, u1, q, b, dx, dt):
	pass

def plug_wave_solution(direction='x', user_action=wave_plotter()):
	# Courant number
	C = 1.0
	N = 50
	L = 1.01
	T = 10
	c = 1.1
	q = lambda X, Y, L: zeros(np.shape(X))+c**2
	dx = L/(N-1)
	dt = dx*C/c# Formula to get the correct courant number C = c*dt/dx
	
	if (direction == 'x'):
		f_I = lambda X, Y, L: np.logical_and(X>(L/2-2*L/N),X<(L/2+2*L/N))*1.0
	elif (direction == 'y'):
		f_I = lambda X, Y, L: np.logical_and(Y>(L/2-2*L/N),Y<(L/2+2*L/N))*1.0
	else:
		print "Neither x nor y provided as direction for plug wave"

	f_V = lambda X, Y, L: zeros(np.shape(X))
	f = lambda X, Y, L: zeros(np.shape(X))
	b = 0
	solver(f_I, f_V, f, q, b, L, N, dt, T, user_action=user_action)

def test_plug_wave_solution():
	user_action = user_action_test_symmetric(axis='x')
	plug_wave_solution(user_action=user_action)
	user_action = user_action_test_symmetric(axis='y')
	plug_wave_solution(user_action=user_action)


def standing_undamped_waves(user_action = wave_plotter(), h = np.asarray([0.1, 0.05, 0.01, 0.005])):
	C = 0.5
	c = 1.1
	#dx = h
	dt = h*C/c
	L = 1.0
	#Nt = 100
	T = 1.0
	A = 1.1
	m_x = int(2)
	m_y = int(2)
	k_x = m_x*np.pi/L
	k_y = m_y*np.pi/L

	### For testing the solution
	exact_solution = lambda X, Y, t: A*np.cos(k_x*X)*np.cos(k_y*Y)*np.cos(c*np.sqrt(k_x**2+k_y**2)*t)
	try:
		user_action.set_exact_solution(exact_solution)
	except:
		pass

	f_I = lambda X, Y, L: A*np.cos(k_x*X)*np.cos(k_y*Y)
	f_V = lambda X, Y, L: np.zeros(np.shape(X))
	f = lambda X, Y, L: np.zeros(np.shape(X))
	q = lambda X, Y, L: np.zeros(np.shape(X))+c**2
	b = 0.0 # Undamped waves

	if [L/i for i in h] == [int(L/i) for i in h]:
		N = L/h+1
		print N
	else:
		print "one or more h-values are not compatible with the L-value such that L/h is integer"
		sys.exit(1)

	max_errors = []
	for i in range(len(N)):
		solver(f_I, f_V, f, q, b, L, N[i], dt[i], T, user_action)
		try:
			max_errors.append(np.max(user_action.get()))
			user_action.reset()
		except:
			pass
	try:
		plt.plot(np.log10(h), np.log10(max_errors))
		plt.xlabel("$\log_{10} \ h$")
		plt.ylabel("$\log_{10} \ \max(e_{i, j}^n)$")
		p = np.polyfit(np.log10(h), np.log10(max_errors), 1)
		plt.title("Order of true error: " + "%.2f" % p[0])
		plt.show()
		raw_input("press enter")
	except:
		pass

def convergence_standing_undamped_waves():
	user_action = user_action_convergence_max_error()
	standing_undamped_waves(user_action = user_action)

def underwater_hill(B0=1.2, Ba=-0.8, Bmx_scl=0.5, Bmy_scl=0.5, Bs=2.3, b=2.0, shape="gaussian"):
	N = 100
	L = 10.0
	T = 4.0
	g = 9.81
	dx = float(L)/(N-1)
	Bmx = L*Bmx_scl
	Bmy = L*Bmy_scl
	
	c = np.sqrt(g*B0);
	dt = dx/c/2.5
	f_q = lambda X, Y, L: ones(np.shape(X))*c**2
	if (shape == "gaussian"):
		f_q = lambda X, Y, L: (B0+Ba*np.exp(-(((X-Bmx)/Bs)**2+((Y-Bmy)/(b*Bs))**2)))*g

	if (shape == "cosine_hat"):
		f_q = lambda X, Y, L: (B0+(Ba*np.cos(np.pi*(X-Bmx)/2/Bs)*np.cos(np.pi*(Y-Bmy)/2/Bs))*((X-Bmx)**2+(Y-Bmy)**2<=Bs**2))*g

	if (shape == "box"):
		f_q = lambda X, Y, L: (B0+Ba*((Bmx-Bs)<=X)*((Bmx+Bs)>=X)*((Bmy-Bs)<=Y)*((Bmy+Bs)>=Y))*g

	f_I = lambda X, Y, L: 3*exp(-4*(X)**2)
	f_V = lambda X, Y, L: zeros(np.shape(X))
	f = lambda x: zeros(np.shape(X))
	b = 1.1

	myPlotter = wave_plotter_tvtk()
	solver(f_I, f_V, f, f_q, b, L, N, dt, T, user_action = myPlotter)
	myPlotter.get();

def read_command_line():
	if (sys.argv[1] == "plug_wave"):
		try:
			direction = sys.argv[2]
		except:
			direction = 'x'
			print "Direction for plug wave was not provided, set to x"
		plug_wave_solution(direction)

	elif (sys.argv[1] == "standing_undamped_waves"):
		standing_undamped_waves(h = np.asarray([0.02]))

	elif (sys.argv[1] == "convergence_standing_undamped_waves"):
		convergence_standing_undamped_waves()

	elif (sys.argv[1] == "underwater_hill"):
		try:
			shape = sys.argv[2]
			underwater_hill(shape=shape)
		except IndexError as e:
			print type(e)
			underwater_hill(shape = "gaussian")
	else:
		exit(1)

if __name__=="__main__":
	try:
		s = sys.argv[1]
		read_command_line()
	except IndexError as e:
		print type(e)
		#raise
		N = 100
		L = 1.1
		T = 10.0
		dx = float(L)/(N-1)
		
		c = 3.0;
		dt = dx/c/2.5
		f_q = lambda X, Y, L: ones(np.shape(X))*c**2
		f_q = lambda X, Y, L: 0.1*(X/L)**2*c**2+0.1
		f_I = lambda X, Y, L: 3*exp(-200*((X-L/2)**2+(Y-L/2)**2))
		f_V = lambda X, Y, L: zeros(np.shape(X))
		f = lambda x: zeros(np.shape(X))
		b = 2.01

		myPlotter = wave_plotter_tvtk()
		solver(f_I, f_V, f, f_q, b, L, N, dt, T, user_action = myPlotter)
		myPlotter.get();