"""
Script to plot g1(x) and g2(x) versus x
and to show where they intersect with black dots

Example use:
	>>> run intersections
"""

import newton
import matplotlib.pyplot as plt
import numpy as np

guesses = np.array([-4,4])
xarray = np.zeros(2)
yarray = np.zeros(2)
i = 0

for xin in guesses:
	(xout, iters) = newton.solve(xin)
	xarray[i] = xout
	yarray[i] = ((xout-1.0)**4.0)-1.0e-4
	i = i+1

x = np.linspace(-5,5,1000)

g1 = ((x-1.0)**4.0)-1.0e-4
g2 = 0*x

plt.figure(1)
plt.clf()
plt.ylim(-5e-4,5e-4)
plt.xlim(0,2)
plt.xlabel('x (-)')
plt.ylabel('y (-)')
plt.title("Two functions and their intersecting points")
#p1 = plt.plot(x, g1, 'b', label=r'$y(x) = xcos(\pix)$')
p1 = plt.plot(x, g1, 'b', label=r'$y(x) = (x - 1)^4 - 1\times 10^{-4}$')
p2 = plt.plot(x, g2, 'r', label=r'$y(x) = 0$')
p3 = plt.plot(xarray, yarray, 'ko', label=r'$intersection$')
plt.legend(loc="best")
plt.show() #show the plot on the screen
plt.savefig('intersection.png')   # save figure as .png file
