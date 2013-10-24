from cs3 import *

(x, y, z) = s("x, y, z")

f = M([[x**2 + y**2 + z**2 - 9], [x*y*z - 1], [x + y - z**2]])

f += M([[x], [y], [z]])

dfdt = f.jacobian([x, y, z])
print dfdt
print f

soln = map(lambda x: map(lambda y: float(y), x), sympy.solve(f, [x,y,z]))[1]

xi = {'x': soln[0]*0.99, 'y': soln[1]*0.99, 'z': soln[2]*0.99 }

nSteps = 0

while (float(e(f, xi).norm()) > .1) and (nSteps < 10):
	plt.plot(nSteps, xi['x'], 'o')
	xi['x'] -= e(f.row(0)[0], xi)/e(dfdt.row(0)[0], xi).norm()
	xi['y'] -= e(f.row(0)[1], xi)/e(dfdt.row(0)[1], xi).norm()
	xi['z'] -= e(f.row(0)[2], xi)/e(dfdt.row(0)[2], xi).norm()
	nSteps += 1
	print xi

plt.show()

print e(f, xi), xi
