
from cs3 import *
from sympy.core.function import AppliedUndef

N = 10

#{"dvdx" = V(i+1) - V(i-1)

# v'' + 2*x*v' - x**2*v = x**2
# f'' + p(x)*f' + q(x)*f = r(x)
# r(i) = x[i]**2
# p(i) = 2*x[i]
# q(i) = -(x[i])**2

# formulation of differential equation using finite difference analysis
# algebra by john geddes
f = s("-r(i)*dx**2 + (1-p(i)*dx/2)*v(i-1) + (-2+q(i)*dx**2)*v(i) + (1+p(i)*dx/2)*v(i+1)")

# build vector of N values on interval (0, 1)
xs = np.linspace(0, 1, N+2)[1:-1]

f = f.subs({s("dx"): np.mean(np.diff(xs))})
print f


funcs = set(map(lambda f: f.func, f.atoms(AppliedUndef)))
locals().update(dict(zip(map(str, funcs), funcs)))

# from wikipedia
# note: function also modifies b[] and d[] params while solving
def TDMASolve(a, b, c, d):
	n = len(d) # n is the numbers of rows, a and c has length n-1
	for i in xrange(n-1):
		d[i+1] -= d[i] * a[i] / b[i]
		b[i+1] -= c[i] * a[i] / b[i]
	for i in reversed(xrange(n-1)):
		d[i] -= d[i+1] * c[i] / b[i+1]
	return [d[i] / b[i] for i in xrange(n)] # return the solution
