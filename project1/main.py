from ParachuteProblem import ParachuteProblem
import matplotlib.pyplot as plt
import numpy as np
import sys
# Take command line args
def read_command_line():
	print sys.argv
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


class TestCases:
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
	else:
		m, g, C_D, A, V, rho, I, T, dt = read_command_line()
"""
	parachuter = ParachuteProblem(m, g, C_D, A, V, rho)
	parachuter.set_initial_condition(I)
	t, u = parachuter.solve(T, dt, plot=False);
	parachuter.calc_forces(t, u, plot=True)
	"""