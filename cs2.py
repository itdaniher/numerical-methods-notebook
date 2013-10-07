from casestudyi import *
from matplotlib import pyplot as plt
import numpy.polynomial.polynomial as nppoly
import time

# midpoint method
def midpoint(f, a, b, tol = 0.1, nMax = 100):
	start = time.time()
	n = 0
	c = (a+b)*0.5
	# check provided start points
	if (f(a) > 0 > f(b)) or (f(a) < 0 < f(b)):
		# while f(c) is greater than the tolerance and iteration count is less than norm
		while (abs(f(c)) > tol) & (n < nMax):
			c = (a+b)*0.5
			if sympy.sign(f(c)) == sympy.sign(f(a)):
				a = c
			else:
				b = c
			n += 1
	else:
		c = None
	# zero, iterct, runtime
	return c, n, time.time()-start

# newton's method
def newtons(f, dfdt, xi, tol = 0.1, nMax = 100):
	start = time.time()
	n = 0
	x = xi
	while (abs(f(x)) > tol) & (n < nMax):
		x = x - f(x)/dfdt(x)
		n += 1
	# zero, iterct, runtime
	return x, n, time.time()-start

# linear interpretation
def lerp(f, a, b, tol = 0.1, nMax = 100):
	start = time.time()
	n = 0
	if (f(a) > 0 > f(b)) or (f(a) < 0 < f(b)):
		while (abs(b-a) > tol) & (n < nMax):
			n += 1
			c = b - f(b)*(b-a)/(f(b) - f(a))	
			if sympy.sign(f(c)) == sympy.sign(f(a)):
				a = c
			else:
				b = c
	else:
		c = None
	# zero, iterct, runtime
	return c, n, time.time()-start

# inverse quadratic interpolation
def iqi(f, xs, tol = 0.1, nMax = 100):
	""" function, initial X values, tolerance, max iterations """
	start = time.time()
	n = 0
	while (abs(f(xs[0])) > tol) & (n < nMax):
		p = nppoly.polyfit(xs, map(f, xs), 3)
		guess = map(abs, nppoly.polyroots(p))[1]
		xs.pop()
		xs.insert(0, guess)
		n += 1
	# zero, iterct, runtime
	return xs[1], n, time.time()-start


if __name__ == "__main__":
	# housekeeping, enumeration of constants
	c = { "pduck" : 0.3, # g / cc
		"ph2o": 1.0, # g / cc
		"r": 10, # cm
		"I": 0 # everything we care about is real
		}
	d = s("d")
	vDuckUnderWater = np.pi/3.0*s("(3*r*d**2-d**3)")
	mDuck = s("4/3*pi*r**3*pduck")
	equalibrium = e(vDuckUnderWater*s("ph2o")-mDuck, c)
	# f(p) = 0 / equilibrium evaluated at the root equals zero
	f = sympy.lambdify(d, equalibrium)
	# symbolicly differentiate the "equalibrium" expression wrt d, convert to python function	
	dfdd = sympy.lambdify(d, sympy.diff(equalibrium, d))
	roots = sympy.roots(equalibrium, d, filter="R")
	# find the root that meets physical specifications
	pNumerical = filter(lambda x: 0 < x.subs(c) < c['r'], roots.keys())[0]
	print pNumerical
	tols = np.logspace(-9, -0.1, 200) * pNumerical
	plt.title("execution time vs tolerance")
	plt.xlabel("absolute tolerance")
	plt.ylabel("execution time")
	plt.loglog(tols, [midpoint(f, 0, c["r"], tol)[2] for tol in tols], label="midpoint")
	plt.loglog(tols, [newtons(f, dfdd, c["r"], tol)[2] for tol in tols], label="newtons")
	plt.loglog(tols, [lerp(f, 0, c["r"], tol)[2] for tol in tols], label="lerp")
	plt.loglog(tols, [iqi(f, [c["r"]]*3, tol)[2] for tol in tols], label="iqi")
	plt.legend(loc="best")

	plt.show()
	# evaluate f at points on the interval 0, c['r']
	#plt.plot(np.linspace(0, c['r'], 1000), map(f, np.linspace(0, c['r'], 1000)))
	#plt.show()
