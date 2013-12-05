"""
Command line input
argv[1]: degree of elemtents
divisions: x, y, z domain resolution. Domain will be unit anyway
"""
import sys
import numpy
import matplotlib.pyplot as plt

from dolfin import *

def alpha(u):
	return 1 + u*u



def solve_nonlinear_diffusion(degree, divisions, alpha, source_function, exact_function, rho, isPlot=False):
	u_exact = exact_function;
	d = len(divisions)
	domain_type = [UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh]
	mesh = domain_type[d-1](*divisions)
	V = FunctionSpace(mesh, 'Lagrange', degree)
	u0 = Constant(0.0)
	u1 = interpolate(u0, V)

	# misc
	T = 1.0
	t = 0
	h = (1.0/(divisions[0]-1))**2
	print h
	dt = h

	# Variational problem
	u = TrialFunction(V)
	u_k = Function(V)
	v = TestFunction(V)
	a = rho*inner(u, v)*dx + dt*inner(alpha(u_k)*nabla_grad(u), nabla_grad(v))*dx
	L = rho*inner(u1, v)*dx + dt*inner(f, v)*dx

	u = Function(V)
	u = interpolate(u0, V)

	if (isPlot):
		wiz_w = plot(u)
		wiz_w.set_min_max(0, 0.01)

	# Effectively one iteration picard method inside time loop
	while t<T:
		f.t = t;
		for i in range(3):
			solve(a==L, u)
			u_k.assign(u)
		u1.assign(u)
		if (isPlot):
			wiz_w.plot(u)
		t += dt
		if (t > 0.05):		
			u_exact.t = t
			u_e = interpolate(u_exact, V)
			e = u_e.vector().array()-u.vector().array()
			E = numpy.sqrt(numpy.sum(e**2)/u.vector().array().size)
			print t
			print E/h
			break
	return h, E

rho = 1.0
f = Expression("-rho*pow(x[0],3)/3.0 + rho*pow(x[0],2)/2.0 + 8*pow(t,3)*pow(x[0],7)/9.0 - 28*pow(t,3)*pow(x[0],6)/9.0 + 7*pow(t,3)*pow(x[0],5)/2.0 - 5*pow(t,3)*pow(x[0],4)/4.0 + 2*t*x[0] - t", t=0, rho=rho)
u_e = Expression("t*pow(x[0],2)*(0.5-x[0]/3.0)", t=0)
alpha = lambda u: 1+u*u

def read_command_line():
	if len(sys.argv) == 1:
		# run default test
		divisions = [[5*i] for i in range(1, 20, 1)]
		hs = []; errors = []
		for division in divisions:
			h, error = solve_nonlinear_diffusion(1, division, alpha, f,  u_e, rho, isPlot=False)
			hs.append(h); errors.append(error)
		print hs, errors
		hs = numpy.asarray(hs); errors = numpy.asarray(errors)
		plt.figure()
		plt.loglog(hs, errors, '*r');
		plt.xlabel("h")
		plt.ylabel("E")
		plt.show()

if __name__=="__main__":
	read_command_line()