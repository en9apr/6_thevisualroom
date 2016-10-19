"""
Random walk approximate solution to Laplace's equation u_{xx} + u{yy} = 0.

Set demo==True to plot a few random walks interactively.
With demo==False many more walks are used to estimate the solution.

The boundary conditions are set for this test problem by evaluating the
true solution u(x,y) = x^2 - y^2 of Laplace's equation on the
boundary of the domain.

Moreover, the exact solution to the discrete equations
  U_{i-1,j} + U_{i+1,j} + U_{i,j-1} + U_{i,j+1} - 4u_{ij} = 0
with boundary values obtained in this way is easily computed. It is simply
given by evaluating the exact solution at the grid points,
  U_{ij} = x_i^2 - y_j^2
This is because the centered difference approximation to the second
derivative is exact when applied to a quadratic function.

This code implements a random walk on a lattice (rectangular grid) where in
each step the walk goes in one of 4 directions based on the value of a
random number that's uniformly distributed in [0,1].

"""

from pylab import *
from numpy.random import RandomState
import time, sys


# Set plot_walk to:
#   True ==> 10 walks will be done, plotting each for illustration.
#   False ==> Many more walks taken to investigate convergence.

plot_walk = False


# problem description:
ax = 0.
bx = 1.  
ay = 0.4
by = 1.

nx = 19
ny = 11

dx = (bx-ax)/float(nx+1)
dy = (by-ay)/float(ny+1)

debug = False     # to turn on some additional printing



def utrue(x,y):
    """
    Return true solution for a test case where this is known.
    This function should satisfy u_{xx} + u_{yy} = 0.
    """

    utrue = x**2 - y**2
    return utrue


def uboundary(x,y):
    """
    Return u(x,y) if (x,y) is on the boundary.
    """

    if (x-ax)*(x-bx)*(y-ay)*(y-by) != 0:
        print "*** Not on the boundary"
        raise Exception("*** uboundary called at non-boundary point")
        return nan

    # assume we know the true solution and can just evaluate there:
    ub = utrue(x,y)  
    return ub


def plot_initialize(i0,j0):
    """
    Set up plot for demo.
    """

    figure(1,figsize=(7,8))
    clf()
    axes([.1,.3,.8,.6]) # leave room for printing messages
    x = linspace(ax,bx,nx+2)
    y = linspace(ay,by,ny+2)
    X,Y = meshgrid(x,y)  # turn into 2d arrays
    plot(X,Y,'k')
    plot(X.T,Y.T,'k')
    plot([ax+i0*dx],[ay+j0*dy],'ro')
    axis('scaled')
    title("Initial point for random walk")
    draw()
    show()
    time.sleep(1)  # pause to see the initial location
    

def plot_ub(xb,yb,ub):
    """
    Called when we hit the boundary
    """
    plot([xb], [yb], 'ro')
    text(ax, ay-0.2, 'Hit boundary where u = %9.6f' % ub, fontsize=20)
    draw()
    time.sleep(2)


def plot_step(iold, jold, i, j, istep):
    """
    Called each step of the random walk to plot progress
    """

    # plot next segment of walk and new point in red:
    plot([ax+iold*dx, ax+i*dx], [ay+jold*dy, ay+j*dy], 'r',linewidth=2)
    plot([ax+i*dx], [ay+j*dy], 'ro')
    title("After %6i random steps" % (istep+1))
    draw()

    time.sleep(0.05)
    # redraw last segment and point in blue:
    plot([ax+iold*dx, ax+i*dx], [ay+jold*dy, ay+j*dy], 'b',linewidth=2)
    plot([ax+i*dx], [ay+j*dy], 'bo')
    draw()


def random_walk(i0, j0, max_steps):

    """
    Take one random walk starting at (i0,j0) until we reach the boundary or
    exceed max_steps steps.
    Return the value at the boundary point reached, or nan if we failed.
    """

    # starting point:
    i = i0
    j = j0
    if plot_walk:
        plot_initialize(i,j)

    # generate as many random numbers as we could possibly need
    # for this walk, since this is much faster than generating one at a time:
    r = random_generator.uniform(0., 1., size=max_steps)

    if debug:
        print "+++ generated r: ", r
    

    for istep in range(max_steps):
        iold,jold = i,j  # needed for plotting only
    
        # Take the next random step with equal probability in each
        # direction:

        if r[istep] < 0.25:
            i = i-1   # step left
        elif r[istep] < 0.5:
            i = i+1   # step right
        elif r[istep] < 0.75:
            j = j-1   # step down
        else:   
            j = j+1   # step up

        if plot_walk:
            plot_step(iold,jold, i,j, istep)

        # check if we hit the boundary:
        if i*j*(nx+1-i)*(ny+1-j) == 0:
            xb = ax + i*dx
            yb = ay + j*dy
            ub = uboundary(xb, yb)

            if plot_walk:
                plot_ub(xb,yb,ub)
            if debug:
                print "Hit boundary at (%7.3f, %7.3f) after %i7 steps, ub = %7.3f" \
                       % (xb,yb,istep+1,ub)
            break  # end the walk

        if istep==(max_steps-1):
            if debug:
                print "Did not hit boundary after %i steps" % max_steps
            if plot_walk:
                text(0, -0.2, "Did not hit boundary after %i steps" \
                        % max_steps, fontsize=20)
                draw()
                time.sleep(2)
            ub = nan
            
    return ub


def many_walks(i0, j0, max_steps, n_mc):

    ub_sum = 0.   # to accumulate boundary values reached from all walks
    n_success = 0    # to keep track of how many walks reached boundary

    for k in range(n_mc):
        i = i0
        j = j0
        ub = random_walk(i0, j0, max_steps)
        if not isnan(ub):
            # use this result unless walk didn't reach boundary
            ub_sum += ub
            n_success += 1

    u_mc = ub_sum / n_success   # average over successful walks

    return u_mc, n_success


if __name__ == "__main__":
    
    # MAIN PROGRAM
    
    # Try it out from a specific (x0,y0):
    x0 = 0.9
    y0 = 0.6
    
    i0 = round((x0-ax)/dx)
    j0 = round((y0-ay)/dy)

    # shift (x0,y0) to a grid point if it wasn't already:
    x0 = ax + i0*dx
    y0 = ay + j0*dy

    u_true = utrue(x0,y0)

    print "True solution of PDE: u(%7.3f, %7.3f) = %9.5f" % (x0,y0,u_true)
    print "Note: with solution used in demo this is also the solution to the"
    print "      the finite-difference equations on the same grid."

    ans = raw_input("Enter seed for random generator or <return> .... ")
    if ans == "":
        seed = None       # will cause a random seed to be used
    else:
        seed = int(ans)   # convert string returned from raw_input into integer

    random_generator = RandomState(seed)  

    # maximum number of step in each before giving up:
    max_steps = 100*max(nx, ny)

    # initial number of Monte-Carlo walks to take:
    n_mc = 10

    u_mc, n_success = many_walks(i0, j0, max_steps, n_mc)

    error = abs((u_mc - u_true) / u_true)

    print "After %8i random walks, u = %15.9f, rel. error = %15.6e" \
            % (n_success, u_mc, error)
    
    if plot_walk:
        sys.exit()   # quit after a few walks

    # start accumulating totals:
    u_mc_total = u_mc
    n_total = n_success

    for i in range(12):
        u_sum_old = u_mc_total * n_total
        u_mc, n_success = many_walks(i0, j0, max_steps, n_mc)
        u_sum_new = u_mc * n_success
        n_total = n_total + n_success
        u_mc_total = (u_sum_old + u_sum_new) / n_total
        error = abs((u_mc_total - u_true) / u_true)

        print "After %8i random walks, u = %15.9f, rel. error = %15.6e" \
                % (n_total, u_mc_total, error)
        n_mc = 2*n_mc   # double number of trials for next iteration
        
        

