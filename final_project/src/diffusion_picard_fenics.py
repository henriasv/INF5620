"""
Command line input
argv[1]: degree of elemtents
divisions: x, y, z domain resolution. Domain will be unit anyway
"""
import sys
import numpy

from dolfin import *

def alpha(u):
	return 1.0

u_exact = Expression("exp(-pi*pi*t)*cos(pi*x[0])", t=0)

degree = int(sys.argv[1]); print "Degree %d " % degree
divisions = [int(arg) for arg in sys.argv[2:]]
d = len(divisions)
domain_type = [UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh]
mesh = domain_type[d-1](*divisions)
V = FunctionSpace(mesh, 'Lagrange', degree)
u0 = Expression("cos(pi*x[0])")
u1 = interpolate(u0, V)

# misc
T = 1.0
t = 0
h = (1.0/(divisions[0]-1))**2
print h
dt = h
rho = 1.0;

# Variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0.0)
a = rho*inner(u, v)*dx + dt*inner(alpha(u1)*nabla_grad(u), nabla_grad(v))*dx
L = rho*inner(u1, v)*dx + dt*inner(f, v)*dx

u = Function(V)
u = interpolate(u0, V)


wiz_w = plot(u)
wiz_w.set_min_max(-1, 1)

# Effectively one iteration picard method inside time loop
while t<T:
	solve(a==L, u)
	wiz_w.plot(u)
	u1.assign(u)
	t += dt
	if (t > 0.05):		
		u_exact.t = t
		u_e = interpolate(u_exact, V)
		wiz_w.plot(u_e)
		e = u_e.vector().array()-u.vector().array()
		E = numpy.sqrt(numpy.sum(e**2)/u.vector().array().size)
		print t
		print E/h
		break
