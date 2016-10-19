"""
Module to compute the square root of a number
using Newton's method: f(x) = x**2 - n = 0

Parameters:
	max_iter :- maximum number of iterations for iterative method
	tol :- relative residual tolerance for iterative method
"""
import numpy as np
max_iter = 40
tol = 1.0e-14

def solve(x0, debug=False):
	"""
	Compute the solution to the function using the initial guess & allow debugging

	Inputs:
		fvals_sqrt :- function returning tuple containing values of f(x) and f'(x)
		x0 :- initial guess of the square root
		debug :- prints x0 outside loop and x within loop, default = False
	Outputs:
		x :- computed value of the square root
		iters :- the number of iterations taken 
	Example usage:
		import newton; newton.solve(4., True)
	"""
	iters = 0
	x = x0	
	
	if(debug):
		print "Initial guess: x = %20.15e" %x0

	for i in range(max_iter):

		(fx,dfx) = fvals_sqrt(x)

		if(abs(fx)<tol):			
			break

		x = x - (fx/dfx)
		iters = iters + 1		
		
		if(debug):
			print "At iteration %d x = %20.15e" %(iters,x) 
	return (x, iters)

def fvals_sqrt(x):
	"""
	Return a tuple containing the value of the function & it's derivative at x
	
	Inputs:
		x :- value for which f(x) and f'(x) are to be returned 
	Outputs:
		f :- the value of the function, f(x)
		fp :- the derivative of the function, f'(x)
	Example usage:	
		(f, fp) = fvals_sqrt(10)
	"""
	
	#f = x**2 - 4.0
	#fp = 2.0*x
	f = (x-1.0)**4.0 - 1.0e-4
	fp = 4.0*(x-1.0)**3.0
	return (f, fp)

def test1(debug_solve=False):
	"""
	Test Newton's method using an array of initial gusses
	
	Inputs:
		None
	Outputs:
		Print x, iters, and f(x) to the screen
	Example usage:
		import newton; newton.test1()
	"""
	#for x in [1., 2., 100.]:
	for x in [-2.2, -1.6, -0.8, 1.4]:

		print ""

		(xout, iters) = solve(x, debug=debug_solve)
		print "After %d iterations, the solver returns x = %20.15e" %(iters,xout)
		
		(f, fp) = fvals_sqrt(xout)
		print "The function f(x) = %20.15e" %f

		#assert abs(xout-2.0) < tol, "***Unexpected result: x = %22.15e" %xout

