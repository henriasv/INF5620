import numpy as np
import nose.tools as nt
import matplotlib.pyplot as plt
import sys

# Class for holding and solving a parachute problem
class ParachuteProblem:
	"""
		Class for holding and solving a parachute problem


	"""
	def __init__(	self, 
					m 	= 	85.0, 
					g	=	9.81, 
					C_D = 	1.2,
					A 	= 	0.5,
					V 	= 	0.1, 
					rho =	1.0,
					I 	= 	0,
					Fs	=	lambda t, dt: 0
				):
		self.m 		= 	m
		self.g 		= 	g
		self.C_D 	= 	C_D
		self.A 		= 	A
		self.V 		= 	V
		self.rho 	= 	rho
		self.I 		= 	I
		self.Fs 	= 	Fs

		# Calculate coefficients for the differential equation
		self.rho_b 	= 	m/V
		self.a 		= 	0.5*C_D*rho*A/(self.rho_b*V)
		self.b 		= 	g*(1-rho/self.rho_b)

	def set_initial_condition(self, I):
		self.I = I

	def set_source_function(self, Fs):
		self.Fs = Fs

	def solve(self, T, dt, t0 = None, u0 = None, plot=False):
		dt = float(dt)
		Nt = int(round(T/dt))
		T0 = 0
		T = Nt*dt
		u = np.zeros(Nt+1)

		if not u0==None:
			self.I 	= 	u0[len(u0)-1]
			T0 		= 	t0[len(t0)-1]

		u[0] = self.I;
		t = np.linspace(T0, T+T0, Nt+1)

		a = self.a; b = self.b; Fs = self.Fs;
		for n in xrange(Nt):
			u[n+1] = (u[n] + (b + Fs(t[n], dt))*dt)/(1+a*abs(u[n])*dt)
		
		
		if not u0==None:
			t0 	= 	np.asarray(t0); u0 = np.asarray(u0)
			t 	= 	np.concatenate((t0, t[1::]))
			u 	= 	np.concatenate((u0, u[1::]))


		if plot:
			self.plot(t, u)
		
		return t, u

	def plot(self, t, u):
		import matplotlib.pyplot as plt
		plt.figure()
		plt.plot(t, u)
		plt.xlabel('t')
		plt.ylabel('v(t)')
		plt.title('Velocity of parachuter')
		plt.show()

	def calc_forces(self, t, u, plot = False):
		bouyancy = self.bouyancy(t)
		drag_force = self.drag_force(t, u)
		gravity_force = self.gravity_force(t)
		print bouyancy
		print drag_force
		print gravity_force
		if plot:
			import matplotlib.pyplot as plt
			plt.figure()
			plt.plot(t, drag_force, label="Drag force")
			plt.hold('on')
			plt.plot(t, bouyancy, label='Bouyancy')
			plt.plot(t, gravity_force, label='Gravity')
			plt.plot(t, drag_force+bouyancy+gravity_force, label="Sum")
			plt.title("Forces acting on parachuter")
			plt.legend()
			plt.show()

	def bouyancy(self, t):
		return self.rho*self.g*self.V + 0*t

	def drag_force(self, t, v):
		return -0.5*self.C_D*self.rho*self.A*abs(v)*v

	def gravity_force(self, t):
		return self.m*self.g + 0*t



# Test function for nose testing
def test_ParachuteProblem_solve():
	print "Testing largest deviation from exact solution of discrete equations."
	T	=	10.0
	dt 	= 	0.01
	m 	= 	85.0
	g 	= 	9.81
	C_D =	1.2
	A 	=	0.5 
	V 	=	0.2 
	rho =	1.03 
	I 	=	10
	A_sf= 	2.0 # slope of manufactured solution
	B_sf= 	I
	
	parachuter 		= ParachuteProblem(m, g, C_D, A, V, rho, I)
	source_function = source_function_discrete(parachuter.a, parachuter.b, A_sf, B_sf)
	parachuter.set_source_function(source_function)
	parachuter.set_initial_condition(I)

	t, u	= 	parachuter.solve(T, dt)
	diff 	= 	A_sf*t + B_sf - u
	diff 	= 	max(abs(diff))
	print "Largest error in array %g" % diff
	nt.assert_almost_equal(diff, 0, delta=1e-12)

def test_ParachuteProblem_convergence_rate():
	print "Testing convergence rate for parachute solver"
	T 		= 	1
	dt_s 	= 	[10**i for i in -np.asarray(range(1, 5))]
	m 		= 	85.0
	g 		= 	9.81
	C_D 	=	1.2
	A 		=	0.5 
	V 		=	0.2 
	rho 	=	1.03 
	I 		=	10
	A_sf 	= 	2.0 # slope of manufactured solution
	B_sf 	= 	I
	
	parachuter 		= ParachuteProblem(m, g, C_D, A, V, rho, I)
	source_function = source_function_continous(parachuter.a, parachuter.b, A_sf, B_sf)
	parachuter.set_source_function(source_function)
	parachuter.set_initial_condition(I)

	errors = []
	for dt in dt_s:
		t, u 		= 	parachuter.solve(T, dt)
		e 			= 	A_sf*t + B_sf -u
		l2err 		=	np.sqrt(dt*np.sum(np.power(e,2)))
		errors.append(l2err)

	errors = np.asarray(errors)/T
	plt.plot(np.log10(dt_s), np.log10(errors))
	plt.xlabel('log10(dt)') 
	plt.ylabel('log10(L_2 norm of error)') 
	plt.title('Convergence towards exact solution')

	plt.savefig('report/figures/convergence.pdf')

	p = np.polyfit(np.log10(dt_s), np.log10(errors), 1)
	print "Order of convergence is %g, " % p[0]
	print "based on linear manufactured solution. \n Timesteps tested: ", dt_s
	nt.assert_almost_equal(p[0], 2, delta=0.1)



# Source functions to get a linear solution to the parachute problem
class source_function_discrete:
	def __init__(self, a, b, A, B):
		self.a = a
		self.b = b
		self.A = A
		self.B = B

	def __call__(self, t, dt):
		a = self.a; b = self.b; A = self.A; B = self.B

		return a*abs(A*t+B)*(A*(t+dt)+B) - b + A 

class source_function_continous:
	def __init__(self, a, b, A, B):
		self.a = a; self.b = b; self.A = A; self.B = B

	def __call__(self, t, dt):
		a = self.a; b = self.b; A = self.A; B = self.B

		return a*abs(A*(t+0.5*dt)+B)*(A*(t+0.5*dt)+B) - b + A 


usageString = "Usage: \n Either: >> python skydiving.py <skydiver_free_fall/skydiver_full_dive> \n Or: >> python skydiving.py m g C_D A V rho I T dt"

# Take command line args; Should be called if __name__=="__main__"
def read_command_line():
	"""
		Read parameters from command line
	"""
	if (len(sys.argv) == 10):
		m 	= 	float(sys.argv[1])
		g 	= 	float(sys.argv[2])
		C_D = 	float(sys.argv[3])
		A 	= 	float(sys.argv[4])
		V	= 	float(sys.argv[5])
		rho = 	float(sys.argv[6])
		I 	= 	float(sys.argv[7])
		T 	= 	float(sys.argv[8])
		dt 	= 	float(sys.argv[9])
		return m, g, C_D, A, V, rho, I, T, dt

	elif len(sys.argv) == 2:
		if(sys.argv[1] == "skydiver_free_fall"):
			test = TestCases()
			test.skydiver_free_fall()
		elif(sys.argv[1] == "skydiver_full_dive"):
			test = TestCases()
			test.skydiver_full_dive()
		else:
			print usageString 


class TestCases:
	"""
		Example usages of ParachuteProblem class
	"""
	def __init__(self):
		self.m = 100
		self.g = 9.81
		self.C_D = 1.2
		self.A = 0.5
		self.V = 0.2
		self.rho = 1.03
		self.I = 0

	def skydiver_free_fall(self, t0=None, u0=None, plot=False):
		self.C_D = 1.2
		T = 10
		dt = 0.001

		parachuter = ParachuteProblem(self.m, self.g, self.C_D, self.A, self.V, self.rho, self.I)
		parachuter.set_initial_condition(self.I)
		t, u =  parachuter.solve(T, dt, plot=plot)
		return t, u

	def skydiver_with_parachute(self, t0=None, u0=None, plot=False):
		T = 10
		dt = 0.001
		self.C_D = 1.8
		self.A = 44
		parachuter = ParachuteProblem(self.m, self.g, self.C_D, self.A, self.V, self.rho, self.I)
		t, u = parachuter.solve(T, dt, t0, u0, plot=plot)

		return t, u

	def skydiver_full_dive(self):
		t, u = self.skydiver_free_fall()
		t, u = self.skydiver_with_parachute(t, u, plot=True)


if __name__=="__main__":
	# Do calculations
	if len(sys.argv) == 2:
		read_command_line()
	elif len(sys.argv) == 10:
		m, g, C_D, A, V, rho, I, T, dt = read_command_line()
		parachuter = ParachuteProblem(m, g, C_D, A, V, rho, I)
		t, u = parachuter.solve(T, dt, plot=True)

	else:
		print usageString



