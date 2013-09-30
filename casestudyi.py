import numpy as np
import sympy
import rk45

from matplotlib import pyplot as plt

from sympy.abc import t
from sympy import sympify as s
from sympy import Matrix as M

e = lambda x, c: x.subs(c.items()).evalf()

# forward euler numerical integration
def feuler(dfdt, x0, v0, tn, c):
	""" dfdt is a function accepting x, v, and dt, returning a floating point number.
	x0 is the initial position.
	v0 is the initial velocity.
	tn is the final time.
	c is a dictionary of constants."""
	x = x0
	v = v0
	t = 0
	m = c['m']
	dt = c['dt']
	to = []
	xo = []
	vo = []
	while t < tn:
		xo.append(x)
		to.append(t)
		vo.append(v)
		v += dfdt(x,v,dt)*dt
		x += v*dt
		t += dt
	return to, xo, vo

# modified euler numerical integration - trapezoidal approximation
def meuler(dfdt, x0, v0, tn, c):
	""" dfdt is a function accepting x, v, and dt,
		returning either a floating point number or a matrix of floating point numbers.
	x0 is the initial position.
	v0 is the initial velocity.
	tn is the final time.
	c is a dictionary of constants."""
	x = x0
	v = v0
	t = 0
	dt = c['dt']
	xo = []
	to = []
	vo = []
	while t < tn:
		xo.append(x)
		to.append(t)
		vo.append(v)
		# get left and right estimates, average
		vi = dfdt(x,v,dt)*dt
		vf = dfdt(x,vi,dt)*dt
		v += 0.5*(vi+vf)
		# calculate error constant
		er = e(2.0*(vf-vi)/dt**2, c)
		# evaluate
		x += v*dt
		# increment time
		t += dt
		#c["dt"] *= e(sympy.sqrt(1.01/(dt*sympy.Abs(vi.norm()-v.norm()))), c)
	return to, xo, vo

if __name__ == "__main__":
	# constant terms
	c = {
		"Cd": 0.3, # unitless
		"p": 1.3, # Kg/m**3
		"A": 0.0042, # m**2
		"m": 0.145,# kg
		"g": -9.8, # m/s**2
	}

	# time without drag
	tf = sympy.solve(1+45*np.sin(np.deg2rad(35)*2)*t+0.5*(-9.8)*t**2,t)[1]
	# distance without drag
	xf = e(45*np.cos(np.deg2rad(35.)*2)*t, {"t": tf})
	print "distance without drag:",  xf
	print "time without drag:",  tf

	# drag equation
	dfdt = lambda x, v, dt: e(M([0, s("g")]) + (s("-1/2*p*Cd*A")*v.norm()**2*v.normalized()), c)
	v = M(["vx", "vy"])
	print "acceleration:", e(M([0, s("g")]) + (s("-1/2*p*Cd*A")*v.norm()**2*v.normalized()), c)
	# initial velocity matrix
	v0 = M([45.*np.cos(np.deg2rad(35.)*2), 45.*np.sin(np.deg2rad(35.)*2)])
	x0 = M([0, 1])
	c['dt'] = 0.1
	rkres =rk45.rk4(dfdt, x0, v0, tf, c)
	print 'rk45 complete'
	meres = meuler(dfdt, x0, v0, tf, c)
	print 'modeuler complete'
	feres = feuler(dfdt, x0, v0, tf, c)
	print 'fwdeuler complete'
	plt.plot(rkres[0], rkres[2], marker='o', label='rk45')
	plt.plot(meres[0], meres[2], marker='.', label='meuler')
	plt.plot(feres[0], feres[2], marker='v', label='feuler')
	plt.legend(loc='best')
	plt.show()
