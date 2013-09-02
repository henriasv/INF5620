# test module for ParachuteProblem
from ParachuteProblem import ParachuteProblem 
import numpy as np

def test_ParachuteProblem_solve():
	import nose.tools as nt
	T	=	10.0
	dt 	= 	0.001
	g 	= 	9.81
	b 	= 	1.26
	m 	= 	85.0
	I 	= 	0
	source_function = source_function_def(m, b, g)
	parachuter = ParachuteProblem(m, b, g, source_function)
	parachuter.set_initial_condition(I)
	t, u	= 	parachuter.solve(T, dt)
	diff 	= 	(u-t)
	diff 	= 	max(abs(diff))
	print diff
	nt.assert_almost_equal(diff, 0, delta=1e-13)

def test_ParachuteProblem_convergence_rate():
	T 		= 	1000
	dt_s 	= 	[0.01, 0.1, 1, 10, 100]
	g 		= 	9.81
	b 		= 	10.0
	m 		= 	85.0
	I 		= 	0

	source_function = source_function_def(m, b, g)
	parachuter 		= ParachuteProblem(m, b, g, source_function)
	parachuter.set_initial_condition(I)
	errors 			= np.zeros(len(dt_s)) # Using L2 norm for errors
	
	i = 0
	for dt in dt_s:
		t, u 		= 	parachuter.solve(T, dt)
		e 			= 	t-u
		l2err 		=	np.sqrt(dt*np.sum(e**2))
		errors[i]	= 	l2err
		i += 1
	
	import matplotlib.pyplot as plt
	plt.figure()
	plt.loglog(dt_s, errors)
	plt.show()


class source_function_def:
	def __init__(self, m, b, g):
		self.g = g
		self.b = b
		self.m = m

	def __call__(self, t, dt):
		return self.b/self.m*t*(t+dt) - self.g + 1