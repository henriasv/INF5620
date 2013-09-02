import numpy as np 

class ParachuteProblem:
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

	def solve(self, T, dt, plot=False):
		dt = float(dt)
		Nt = int(round(T/dt))
		T = Nt*dt
		u = np.zeros(Nt+1)
		u[0] = self.I;
		t = np.linspace(0, T, Nt+1)

		a = self.a; b = self.b; Fs = self.Fs;
		for n in xrange(Nt):
			u[n+1] = (u[n] + (b + Fs(t[n], dt))*dt)/(1+a*abs(u[n])*dt)
		
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