import sympy

e = lambda x, c: x.subs(c.items()).evalf()

def rk4(a, x, v, ec, c = {'dt': 1./40}):
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
		a1 = a(x1, v1, 0)
	
		x2 = x + 0.5*v1*dt
		v2 = v + 0.5*a1*dt
		a2 = a(x2, v2, dt/2.0)

		x3 = x + 0.5*v2*dt
		v3 = v + 0.5*a2*dt
		a3 = a(x3, v3, dt/2.0)
		x4 = x + v3*dt
		v4 = v + a3*dt
		a4 = a(x4, v4, dt)
	
		x = x + (dt/6.0)*(v1 + 2*v2 + 2*v3 + v4)
		v = v + (dt/6.0)*(a1 + 2*a2 + 2*a3 + a4)
		t += dt
		dt *= e(sympy.sqrt(tol/(dt*sympy.Abs(v2.norm()-v1.norm()))), c)
	return to, xo, vo

def spring(x, v, dt):
	stiffness = 1
	damping = -0.005
	return -stiffness*x + damping*v
