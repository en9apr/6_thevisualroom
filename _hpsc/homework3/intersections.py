"""
Script to plot g1(x) and g2(x) versus x
and to show where they intersect with black dots

Example use:
	>>> run intersections
"""

import newton
import matplotlib.pyplot as plt
import numpy as np

guesses = np.array([-2.2, -1.6, -0.8, 1.4])
xarray = np.zeros(4)
yarray = np.zeros(4)
i = 0

for xin in guesses:
	(xout, iters) = newton.solve(xin)
	xarray[i] = xout
	yarray[i] = (1-(0.6*xout**2))
	i = i+1

x = np.linspace(-5,5,1000)

g1 = x*np.cos(np.pi*x)

g2 = (1-(0.6*x**2))

plt.figure(1)
plt.clf()
plt.xlim(-5,5)
plt.xlabel('x (-)')
plt.ylabel('y (-)')
plt.title("Two functions and their intersecting points")
#p1 = plt.plot(x, g1, 'b', label=r'$y(x) = xcos(\pix)$')
p1 = plt.plot(x, g1, 'b', label=r'$y(x) = xcos(\pi x)$')
p2 = plt.plot(x, g2, 'r', label=r'$y(x) = 1-(0.6x^2)$')
p3 = plt.plot(xarray, yarray, 'ko', label=r'$intersection$')
plt.legend(loc="best")
#plt.legend([p1, p2, p3],["g1(x) = xcos(pi*x)", "g2(x) = 1-(0.6x**2)", "intersection"])
plt.show() #show the plot on the screen
plt.savefig('intersection.png')   # save figure as .png file
