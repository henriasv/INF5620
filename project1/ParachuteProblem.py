import numpy as np 

class ParachuteProblem:
	def __init__(self, m, b, g=9.81, Fs=lambda x, dx: 0):
		self.b = b
		self.g = g
		self.m = m
		self.Fs = Fs
		self.I = 0


	def set_initial_condition(self,I):
		self.I = I

	def solve(self, T, dt):
		dt = float(dt)
		Nt = int(round(T/dt))
		T = Nt*dt
		u = np.zeros(Nt+1)
		u[0] = self.I;
		t = np.linspace(0, T, Nt+1)

		for n in xrange(Nt):
			u[n+1] = (u[n]+(self.g+self.Fs(t[n], dt))*dt)/(1.0+self.b/self.m*abs(u[n])*dt)

		return t, u

	def plot(self, T, dt):
		t, u = self.solve(T, dt)
		import matplotlib.pyplot as plt
		plt.figure()
		plt.plot(t, u)
		plt.xlabel('t')
		plt.ylabel('v(t)')
		plt.title('Velocity of parachuter')
		plt.show()

