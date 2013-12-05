"""
Command line input
argv[1]: degree of elemtents
divisions: x, y, z domain resolution. Domain will be unit anyway
"""
import sys
import numpy
import matplotlib.pyplot as plt

from dolfin import *

class alpha:
	def __init__(self):
		pass
	def __call__(self, u):
		return 1 + u*u

class alpha2(alpha):
	def __init__(self, beta):
		self.beta = beta

	def __call__(self, u):
		return 1 + self.beta*u*u

class User_Action:
	def __init__(self, function=None):
		self.E = None
		self.h = None
		self.t = None

	def act_time_loop(self,t, u,V):
		pass

	def act_after_time_loop(self, t, u, V, h):
		pass

class User_Action_convergence_dt(User_Action):
	def __init__(self, u_e):
		self.u_e = u_e

	def act_after_time_loop(self,t, u, V, h):
		self.h = h
		u_e.t = t
		u_e_in_V = interpolate(u_e, V)
		e = u_e_in_V.vector().array() - u.vector().array()
		self.E = numpy.sqrt(numpy.sum(e**2)/u.vector().array().size)

	def act_time_loop(self, t, u, V):
		return False

class User_Action_different_t(User_Action):
	def __init__(self, u_e):
		self.u_e = u_e

	def act_after_time_loop(self, t, u, V, h):
		self.h = h
		u_e.t = u_e_in_V = interpolate(u_e, V)

class User_Action_plot(User_Action):
	def __init__(self, minmax):
		self.minmax = minmax

	def act_time_loop(self, t, u, V):
		wiz_w = plot(u)
		wiz_w.set_min_max(*self.minmax)
		return False


def solve_nonlinear_diffusion(degree, divisions, alpha, 
								source_function=Constant(0.0), 
								rho=1.0, 
								isPlot=False, 
								T=1.0, 
								user_action=User_Action(), 
								u0 = Constant(0.0)):
	d = len(divisions)
	domain_type = [UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh]
	mesh = domain_type[d-1](*divisions)
	V = FunctionSpace(mesh, 'Lagrange', degree)
	u1 = interpolate(u0, V)
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

	# Effectively N iteration picard method inside time loop
	while t<T:
		f.t = t;
		for i in range(1):
			solve(a==L, u)
			u_k.assign(u)
		u1.assign(u)
		if (isPlot):
			wiz_w.plot(u)
		t += dt
		isBreak = user_action.act_time_loop(t, u, V)
		if isBreak:
			break

	user_action.act_after_time_loop(t, u, V, h)
	return user_action

rho = 1.0
f = Expression("-rho*pow(x[0],3)/3.0 + rho*pow(x[0],2)/2.0 + 8*pow(t,3)*pow(x[0],7)/9.0 - 28*pow(t,3)*pow(x[0],6)/9.0 + 7*pow(t,3)*pow(x[0],5)/2.0 - 5*pow(t,3)*pow(x[0],4)/4.0 + 2*t*x[0] - t", t=0, rho=rho)
u_e = Expression("t*pow(x[0],2)*(0.5-x[0]/3.0)", t=0)
alpha = lambda u: 1+u*u

def read_command_line():
	if len(sys.argv) == 2:
		if sys.argv[1] == "convergence_dt":
		# run default test
			divisions = [[5*i] for i in range(1, 10, 1)]
			hs = []; errors = []
			my_user_action = User_Action_convergence_dt(u_e)

			for division in divisions:
				user_action = solve_nonlinear_diffusion(1, division, alpha, source_function=f, rho=rho, user_action=my_user_action, T = 0.05)
				hs.append(user_action.h); errors.append(user_action.E)
			print hs, errors
			hs = numpy.log10(numpy.asarray(hs)); errors = numpy.log10(numpy.asarray(errors))
			plt.figure()
			p = numpy.polyfit(hs, errors, 1)
			print p
			plt.title(r"Order of error in $\Delta t$: %.2f" % p[0])
			plt.plot(hs, errors, '*r');
			plt.xlabel(r"$\log_{10}(h)$")
			plt.ylabel(r"$\log_{10}(E)$")
			plt.show()

		if (sys.argv[1] == "different_t"):
			divisions = [20, 20]
			t = [0.1, 1, 4, 7]

		if (sys.argv[1] == "gaussian"):
			division = [50, 50]
			u0 = Expression("exp(-0.5/pow(sigma, 2)*(x[0]*x[0]+x[1]*x[1]))", sigma=0.1)
			my_user_action = User_Action_plot((0.0, 0.1))
			user_action = solve_nonlinear_diffusion(1, division, alpha2(100000.0), user_action=my_user_action, T=0.05, u0 = u0)

if __name__=="__main__":
	read_command_line()