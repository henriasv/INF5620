from sympy import *
x, t, rho, dt = symbols('x t rho dt')
def a(u):
	return 1 + u**2

def u_simple(x, t):
	return x**2*(Rational(1,2) - x/3)*t

# Show that u_simple satisfies the BCs
for x_point in 0, 1:
	print 'u_x(%s,t):' % x_point,
	print diff(u_simple(x, t), x).subs(x, x_point).simplify()

u_x(0,t): 0
u_x(1,t): 0
print 'Initial condition:', u_simple(x, 0)
Initial condition: 0
# MMS: full nonlinear problem
u = u_simple(x, t)
f = rho*diff(u, t) - diff(a(u)*diff(u, x), x)
print f.simplify()
-rho*x**3/3 + rho*x**2/2 + 8*t**3*x**7/9 - 28*t**3*x**6/9 +
7*t**3*x**5/2 - 5*t**3*x**4/4 + 2*t*x - t