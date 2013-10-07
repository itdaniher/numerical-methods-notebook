import numpy as np
import sympy

from sympy.abc import t
from sympy import sympify as s
from sympy import Matrix as M

e = lambda x, c: x.subs(c.items()).evalf()

# forward euler numerical integration
def euler(dfdt, x, v, ec, c):
	""" dfdt is a function accepting x, v, and dt, returning a floating point number.
	x0 is the initial position.
	v0 is the initial velocity.
	tn is the final time.
	c is a dictionary of constants."""
	t = 0
	to = []
	xo = []
	vo = []
	dt = c['dt']
	while ec(t, x, v):
		xo.append(x)
		to.append(t)
		vo.append(v)
		v += dfdt(x,v)*dt
		x += v*dt
		t += dt
	return to, xo, vo

# heun's method numerical integration - trapezoidal approximation
def heun(dfdt, x, v, ec, c):
	""" dfdt is a function accepting x, v, and dt,
		returning either a floating point number or a matrix of floating point numbers.
	x0 is the initial position.
	v0 is the initial velocity.
	tn is the final time.
	c is a dictionary of constants."""
	t = 0
	dt = c['dt']
	tol = c['tol']
	xo = []
	to = []
	vo = []
	while ec(t, x, v):
		xo.append(x)
		to.append(t)
		vo.append(v)
		# get left and right estimates, average
		vi = dfdt(x,v)*dt
		vf = dfdt(x,vi)*dt
		v += 0.5*(vi+vf)
		# evaluate
		x += v*dt
		# increment time
		t += dt
		if tol:
			# modify timestep
			dt *= e(sympy.sqrt(tol/(dt*sympy.Abs(vi.norm()-v.norm()))), c)
	return to, xo, vo


def rk4(dfdt, x, v, ec, c = {'dt': 1./40}):
	dt = c['dt']
	tol = c['tol']
	xo = []
	to = []
	vo = []
	t = 0
	while ec(t, x, v):
		xo.append(x)
		to.append(t)
		vo.append(v)
		x1 = x
		v1 = v
		a1 = dfdt(x1, v1)
	
		x2 = x + 0.5*v1*dt
		v2 = v + 0.5*a1*dt
		a2 = dfdt(x2, v2)

		x3 = x + 0.5*v2*dt
		v3 = v + 0.5*a2*dt
		a3 = dfdt(x3, v3)
		x4 = x + v3*dt
		v4 = v + a3*dt
		a4 = dfdt(x4, v4)
	
		x = x + (dt/6.0)*(v1 + 2*v2 + 2*v3 + v4)
		v = v + (dt/6.0)*(a1 + 2*a2 + 2*a3 + a4)
		t += dt
		if tol:
			dt *= e(sympy.sqrt(tol/(dt*sympy.Abs(v2.norm()-v1.norm()))), c)
	return to, xo, vo

def spring(x, v):
	stiffness = 1
	damping = -0.005
	return -stiffness*x + damping*v
