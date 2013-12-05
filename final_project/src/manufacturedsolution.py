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

print 'Initial condition:', u_simple(x, 0)
# MMS: full nonlinear problem
u = u_simple(x, t)
f = rho*diff(u, t) - diff(a(u)*diff(u, x), x)
print latex(f.simplify())
