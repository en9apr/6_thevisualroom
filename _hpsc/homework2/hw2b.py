
"""
Demonstration module for quadratic interpolation.
Update this docstring to describe your code.
Modified by: ** Andrew Roberts **
"""


import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve

def quad_interp(xi,yi):
    """
    Quadratic interpolation.  Compute the coefficients of the polynomial
    interpolating the points (xi[i],yi[i]) for i = 0,1,2.
    Returns c, an array containing the coefficients of
      p(x) = c[0] + c[1]*x + c[2]*x**2.

    """

    # check inputs and print error message if not valid:

    error_message = "xi and yi should have type numpy.ndarray"
    assert (type(xi) is np.ndarray) and (type(yi) is np.ndarray), error_message

    error_message = "xi and yi should have length 3"
    assert len(xi)==3 and len(yi)==3, error_message

    # Set up linear system to interpolate through data points:

    A = np.vstack([np.ones(3), xi, xi**2]).T
    b = yi
    c = solve(A,b)

    return c

def cubic_interp(xi,yi):
    """
    Cubic interpolation.  Compute the coefficients of the polynomial
    interpolating the points (xi[i],yi[i]) for i = 0,1,2,3.
    Returns c, an array containing the coefficients of
      p(x) = c[0] + c[1]*x + c[2]*x**2 +c[3]x**3.

    """

    # check inputs and print error message if not valid:

    error_message = "xi and yi should have type numpy.ndarray"
    assert (type(xi) is np.ndarray) and (type(yi) is np.ndarray), error_message

    error_message = "xi and yi should have length 4"
    assert len(xi)==4 and len(yi)==4, error_message

    # Set up linear system to interpolate through data points:
    
    A = np.vstack([np.ones(4), xi, xi**2, xi**3]).T
    b = yi
    c = solve(A,b)

    return c

def poly_interp(xi,yi):
    
    """Polynomial interpolation.  Compute the coefficients of the polynomial
    interpolating the points (xi[i],yi[i]) for i = 0,1,2... n
    Returns c, an array containing the coefficients of
      p(x) = c[0] + c[1]*x + c[2]*x**2 + c[3]x**3 +...+ c[n]x**n.
    """

    # check inputs and print error message if not valid:

    error_message = "xi and yi should have type numpy.ndarray"
    assert (type(xi) is np.ndarray) and (type(yi) is np.ndarray), error_message

    error_message = "xi and yi should have the same length"
    assert len(xi)== len(yi), error_message

    # Set up linear system to interpolate through data points:

    #stacks them vertically and then takes transpose
    #A = np.empty([len(xi),len(xi)])
    #A[:,0] = 1.
    #for i in range(0,len(xi),1):
    #    for j in range(1,len(xi),1):
    #        A[i,j]=xi[i]**(j)

    # Uses a list comprehension:
    n = len(xi)
    A = np.vstack([xi**j for j in range(n)]).T
    b = yi
    c = solve(A,b)

    return c
    
def plot_quad(xi,yi):
    c = quad_interp(xi,yi)
    x = np.linspace(xi.min()-1, xi.max() + 1, 1000)
    y = c[0]+c[1]*x + c[2]*x**2
    plt.figure(1)
    plt.clf()
    plt.plot(x,y,'b-')
    plt.plot(xi,yi,'ro')
    plt.xlabel('x (-)')
    plt.ylabel('y (-)')
    plt.title("Data points and interpolated polynomial")
    plt.savefig('quadratic.png')

def plot_cubic(xi,yi):
    c = cubic_interp(xi,yi)
    x = np.linspace(xi.min()-1, xi.max() + 1, 1000)
    y = c[0]+c[1]*x + c[2]*x**2 + c[3]*x**3
    plt.figure(1)
    plt.clf()
    plt.plot(x,y,'b-')
    plt.plot(xi,yi,'ro')
    plt.xlabel('x (-)')
    plt.ylabel('y (-)')
    plt.title("Data points and interpolated polynomial")
    plt.savefig('cubic.png')

def plot_poly(xi,yi):
    c = poly_interp(xi,yi)
    x = np.linspace(xi.min()-1, xi.max() + 1, 1000)
    n = len(c)
    y = c[n-1]
    for j in range(n-1, 0, -1):
        y = y*x + c[j-1]

    plt.figure(1)
    plt.clf()
    plt.plot(x,y,'b-')
    plt.plot(xi,yi,'ro')
    plt.xlabel('x (-)')
    plt.ylabel('y (-)')
    plt.title("Data points and interpolated polynomial")
    plt.savefig('poly.png')

def test_quad1():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([-1.,  0.,  2.])
    yi = np.array([ 1., -1.,  7.])
    c = quad_interp(xi,yi)
    c_true = np.array([-1.,  0.,  2.])
    print "c =      ", c
    print "c_true = ", c_true
    # test that all elements have small error:
    assert np.allclose(c, c_true), \
        "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)
        
def test_quad2():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([-10.,  1.,  20.])
    yi = np.array([ 15., -50.,  70.])
    c = quad_interp(xi,yi)
    c_true = np.array([-48.16586922,  -2.24162679,   0.40749601])
    print "c =      ", c
    print "c_true = ", c_true
    # test that all elements have small error:
    assert np.allclose(c, c_true), \
        "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

def test_cubic1():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([-1., 1., 2., 3.])
    yi = np.array([5., 2., 3., 5])
    c = cubic_interp(xi,yi)
    c_true = np.array([ 2.5, -1.41666667, 1., -0.08333333])
    print "c =      ", c
    print "c_true = ", c_true
    # test that all elements have small error:
    assert np.allclose(c, c_true), \
        "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

def test_poly1():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([-1., 1., 2., 3.])
    yi = np.array([5., 2., 3., 5])
    c = cubic_interp(xi,yi)
    c_true = np.array([ 2.5, -1.41666667, 1., -0.08333333])
    print "c =      ", c
    print "c_true = ", c_true
    # test that all elements have small error:
    assert np.allclose(c, c_true), \
        "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

def test_poly2():
    """
    Test code, no return value or exception if test runs properly.
    """
    xi = np.array([-1., 0., 2., 5., 8.])
    yi = np.array([ 1., -1., 7., 15., 25.])
    c = poly_interp(xi,yi)
    c_true = np.array([-1., 1.22777778, 2.51944444, -0.66111111, 0.04722222])
    print "c =      ", c
    print "c_true = ", c_true
    # test that all elements have small error:
    assert np.allclose(c, c_true), \
        "Incorrect result, c = %s, Expected: c = %s" % (c,c_true)

if __name__=="__main__":
    # "main program"
    # the code below is executed only if the module is executed at the command line,
    #    $ python demo2.py
    # or run from within Python, e.g. in IPython with
    #    In[ ]:  run demo2
    # not if the module is imported.
    print "Running test..."
    test_quad1()

