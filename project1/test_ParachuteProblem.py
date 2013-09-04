# test module for ParachuteProblem
from ParachuteProblem import ParachuteProblem 
import numpy as np
import nose.tools as nt
import matplotlib.pyplot as plt

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