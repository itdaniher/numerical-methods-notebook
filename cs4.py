
from cs3 import *
from sympy.core.function import AppliedUndef

# formulation of differential equation using finite difference analysis
f = s("r(i)*dx**2 - (1-p(i)*dx/2)*v(i-1) -(-2+q(i)*dx**2)*v(i) - (1+p(i)*dx/2)*v(i+1)")

print M(f)

# build vector of x values on interval (0, 1)
xs = np.linspace(0,1, 52)[1:-1]

dx = np.mean(np.diff(xs))
print dx

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
