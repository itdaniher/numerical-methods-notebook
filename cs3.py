from casestudyi import *

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

	v0 = M([45.*np.cos(np.deg2rad(35.)*2), 45.*np.sin(np.deg2rad(35.)*2)])
	x0 = M([0, 1])
	ps = np.linspace(0, 10, 20)
	c['dt'] = 0.1
	c['tol'] = .5
	ranges = []
	for p in ps:
		c['p'] = p
		# drag equation
		dfdt = lambda x, v, dt: e(M([0, s("g")]) + (s("-1/2*p*Cd*A")*v.norm()**2*v.normalized()), c)
		ranges.append(meuler(dfdt, x0, v0, lambda t, x, v: x[1] > 0, c)[1][-1][0])
	plt.plot(ps, ranges, 'k-o')
	plt.title("Reynolds vs. Range")
	plt.xlabel("Reynolds Number")
	plt.ylabel("Range (m)")
	plt.show()

#v'' + 2xv' − x2**v = x**2 , v(0) = 1, v(1) = 0,

#dvdx = u

#dudx = x**2+x*v-2*x*u


