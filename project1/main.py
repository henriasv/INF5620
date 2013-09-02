from ParachuteProblem import ParachuteProblem
import matplotlib.pyplot as plt
import sys
# Take command line args
def read_command_line():
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


class TestCases:
	def __init__(self):
		self.m = 85
		self.g = 9.81
		self.C_D = 1.2
		self.A = 0.5
		self.V = 0.2
		self.rho = 1.03
		self.I = 0

	def skydiver_free_fall(self, T0 = 0, I=0):
		self.I = I
		self.m = 100
		T = 10
		dt = 0.001

		parachuter = ParachuteProblem(self.m, self.g, self.C_D, self.A, self.V, self.rho, self.I)
		parachuter.set_initial_condition(self.I)
		t, u =  parachuter.solve(T, dt, plot=True)
		return t + T0, u

	def skydiver_with_parachute(self, u0 = None, t0 = None):
		T = 10
		dt = 0.001

		parachuter = ParachuteProblem(self.m, self.g, self.C_D, self.A, self.V, self.rho, self.I)
		parachuter.set_initial_condition(self.I)
		t, u = parachuter.solve(T, dt, plot=True)
		return t + T0, u

	def skydiver_full_dive(self)
		t, u = self.skydiver_free_fall()
		t, u = self.skydiver_with_parachute()		
if __name__=="__main__":
	# Do calculations
	m, g, C_D, A, V, rho, I, T, dt = read_command_line()
"""
	parachuter = ParachuteProblem(m, g, C_D, A, V, rho)
	parachuter.set_initial_condition(I)
	t, u = parachuter.solve(T, dt, plot=False);
	parachuter.calc_forces(t, u, plot=True)
	"""