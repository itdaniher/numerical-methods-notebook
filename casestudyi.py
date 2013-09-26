import numpy as np
import sympy
import copy
from matplotlib import pyplot as plt

from sympy.abc import t
from sympy import sympify as s
from sympy import Matrix as M

e = lambda x, c: x.subs(c.items()).evalf()

def feuler(x0, v0, tn, f, c):
	x = x0
	v = v0
	t = 0
	m = s("m")
	dt = s("dt")
	xo = []
	to = []
	while t < tn:
		xo.append(x)
		to.append(t)
		v += e(f(x,v,dt)/m*dt, c)
		x += e(v*dt, c)
		t += e(dt, c)
	return to, xo

def meuler(x0, v0, tn, f, c):
	x = x0
	v = v0
	t = 0
	m = s("m")
	dt = s("dt")
	xo = []
	to = []
	while t < tn:
		xo.append(x)
		to.append(t)
		# get left and right estimates, average
		vi = e(f(x,v,dt)/m*dt, c)
		vf = e(f(x,vi,dt)/m*dt, c)
		v += 0.5*(vi+vf)
		# calculate error constant
		er = e(2.0*(vf-vi)/dt**2, c)
		# evaluate
		x += e(v*dt, c)
		# increment time
		t += e(dt, c)
		#c["dt"] *= e(sympy.sqrt(1.01/(dt*sympy.Abs(vi.norm()-v.norm()))), c)
	return to, xo

if __name__ == "__main__":
	# time without drag
	tf = sympy.solve(1+45*np.sin(35*np.pi*2/180)*t+0.5*(-9.8)*t**2,t)[1]
	# distance without drag
	xf = e(45*np.cos(35*np.pi*2/180)*t, {"t": tf})
	print xf

	# drag equation
	f0 = lambda x, v, dt: M([0, s("m*g")])
	f = lambda x, v, dt: M([0, s("m*g")]) +(s("-1/2*p*Cd*A")*v.norm()**2*v.normalized())
	# initial velocity matrix
	v0 = M([45*np.cos(35*(2*np.pi)/180), # m/s
		45*np.sin(35*(2*np.pi)/180)]) # m/s

	x0 = M([0, 1])

	# constant terms
	c = {
		"Cd": 0.3, # unitless
		"p": 1.3, # Kg/m**3
		"A": 0.0042, # m**2
		"m": 0.145,# kg
		"g": -9.8, # m/s**2
		"dt": .01 # s
	}

	plt.figure()
	ts = np.logspace(-3, 0.1, 8)[::-1]
	ferr = []
	for dt in ts:
		c["dt"] = dt
		err = np.abs(feuler(x0, v0, tf, f, c)[1][-1].norm())
		ferr.append(err)
		#plt.plot(fwd[0], map(lambda x: x[1], fwd[1]), "r-o", label="fwd, "+str(dt)[0:5])
	plt.loglog(ts, ferr, "-o")#, label="fwd, "+str(dt)[0:5])
	merr = []
	for dt in ts:
		c["dt"] = dt
		err = np.abs(meuler(x0, v0, tf, f, c)[1][-1].norm())
		#err = np.abs(meuler(x0, v0, tf, f, c)[1][-1].norm()-xf)/xf
		merr.append(err)
		#plt.plot(mod[0], map(lambda x: x[1], mod[1]), "b-v", label="mod, "+str(dt)[0:5])
	plt.loglog(ts, merr, "-v")#, label="mod, "+str(dt)[0:5])
	plt.xlabel("time (s)")
	plt.ylabel("error")
	plt.title("Error for Heun and Euler integration for varying timesteps - with drag")
	plt.legend(loc="best")
	plt.show()

