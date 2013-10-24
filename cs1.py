import numpy as np
import sympy
import sys
sys.displayhook = sympy.pprint
import integrators as qd

from matplotlib import pyplot as plt

from sympy.abc import t
from sympy import sympify as s
from sympy import Matrix as M

e = lambda x, c: x.subs(c.items()).evalf()

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
	xf = qd.e(45*np.cos(np.deg2rad(35.)*2)*t, {"t": tf})
	print "distance without drag:",  xf
	print "time without drag:",  tf

	# drag equation
	dfdt = lambda x, v: qd.e(M([0, s("g")]) + (s("-1/2*p*Cd*A")*v.norm()**2*v.normalized()), c)
	v = M(["vx", "vy"])
	print "acceleration:", qd.e(M([0, s("g")]) + (s("-1/2*p*Cd*A")*v.norm()**2*v.normalized()), c)
	# initial velocity matrix
	v0 = M([45.*np.cos(np.deg2rad(35.)*2), 45.*np.sin(np.deg2rad(35.)*2)])
	x0 = M([0, 1])
	c['dt'] = 0.1
	c['tol'] = 0# 0.001*xf # 0.1% of estimated range
	rData = qd.rk4(dfdt, x0, v0, lambda t, x, v: x[1] > 0, c)
	print 'rk45 complete'
	hData = qd.heun(dfdt, x0, v0, lambda t, x, v: x[1] > 0, c)
	print 'modeuler complete'
	eData = qd.euler(dfdt, x0, v0, lambda t, x, v: x[1] > 0, c)
	print 'fwdeuler complete'
	plt.plot(rData[0], rData[1], marker='o', label='RK4')
	plt.plot(hData[0], hData[1], marker='.', label='Heun')
	plt.plot(eData[0], eData[1], marker='v', label='Euler')
	plt.legend(loc='best')
	plt.show()
