from sympy import Symbol, latex

# Definitions for finite difference schemes in R^2\timesR

# u variables are defined as u_ijn, m means minus
u_000 = Symbol('u_{i, j}^{n}')
u_001 = Symbol('u_{i,j}^{n+1}')
u_00m1 = Symbol('u_{i, j}^{n-1}')

u_100 = Symbol('u_{i+1, j}^n')
u_010 = Symbol('u_{i, j+1}^n')
u_m100 = Symbol('u_{i-1, j}^n')
u_0m10 = Symbol('u_{i, j-1}^n')

q_100 = Symbol('q_{i+1, j}^n')
q_010 = Symbol('q_{i, j+1}^n')
q_m100 = Symbol('q_{i-1, j}^n')
q_0m10 = Symbol('q_{i, j-1}^n')

dt = Symbol('\Delta t')
dx = Symbol('\Delta x')

equation = (u_001 -2*u_000 + u_00m1)/(dt**2)
print latex(equation)