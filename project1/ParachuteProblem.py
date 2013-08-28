import numpy as np 

class ParachuteProblem:
	def __init__(self, m, b):
		self.b = b
		self.g = 9.81
		self.m = m


	def set_initial_condition(self,I):
		self.I = I

	def solve(self, T, dt, Fs=0):
		dt = float(dt)
		Nt = int(round(T/dt))
		T = Nt*dt
		u = np.zeros(Nt+1)
		try:
			u[0] = self.I;
		except:
			u[0] = 0;

		t = np.linspace(0, T, Nt+1)

		if (Fs == 0):
			Fs = np.zeros(Nt+1)

		for n in range(0, Nt):
			u[n+1] = (u[n]+(self.g+Fs[n])*dt)/(1.0+self.b/self.m*abs(u[n])*dt)

		return t, u



