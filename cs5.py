from cs3 import *

(x, y, z) = s("x, y, z")

f = M([[x**2 + y**2 + z**2 - 9], [x*y*z - 1], [x + y - z**2]])

f += M([[x], [y], [z]])

invdfdt = f.jacobian([x, y, z]).inv()

soln = map(lambda x: map(lambda y: float(y), x), sympy.solve(f, [x,y,z]))

xi = {'x': 0, 'y': 0, 'z': 0}

nSteps = 0

tol = 0.1

while (float(e(f, xi).norm()) > tol) and (nSteps < 100):
	dx = e(invdfdt * f, xi)
	xi['x'] -= dx[0]
	xi['y'] -= dx[1]
	xi['z'] -= dx[2]
	plt.semilogy(nSteps, e(f, xi).norm(), '.')
	nSteps += 1

print e(f, xi), xi, soln

plt.show()
