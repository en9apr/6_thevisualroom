def initial_and_boundary_conditions(nx, nt):

    # Initialise data structures
    import numpy as np
    u = np.zeros((nx,nt))

    # Boundary conditions
    u[0,:] = 1.0
    u[nx-1,:] = 0.0

    half = int((nx-1)/2)

    # Initial conditions
    u[0:half,0] = 1.0

    return u

def analytical(sigma, nx, tmax, xmax, c):

    # Increments
    dx = xmax/(nx-1)
    dt =  ((dx*sigma) / c)
    nt = int(sigma*tmax/dt + 1)

    # Initial and Boundary Conditions
    u = initial_and_boundary_conditions(nx, nt)

    # X Loop
    x = np.zeros(nx)
    x = np.linspace(0.0,xmax,nx)

    i = 1
    # Loop - must loop as NumPy cannot be used for dependency reasons
    for n in range(0,nt-1):
        u[i:nx-1, n+1] = u[i-1:nx-2, n]

    return u, x, nt

def lax_friedrichs_convection(sigma, nx, tmax, xmax, c):
    """
    Returns the velocity field and distance for 1D linear convection
    """
    # Increments
    dx = xmax/(nx-1)
    dt =  ((dx * sigma) / c)
    nt = int(tmax/dt + 1)

    # Initial and Boundary Conditions
    u = initial_and_boundary_conditions(nx, nt)

    # X Loop
    x = np.zeros(nx)
    x = np.linspace(0.0,xmax,nx)

    i = 1
    # Loop - must loop as NumPy cannot be used for dependency reasons
    for n in range(0,nt-1):
        u[i:nx-1, n+1] = 0.5*(u[i-1:nx-2, n]+u[i+1:nx, n])-0.5*sigma*(u[i+1:nx, n]-u[i-1:nx-2, n])

    return u, x, nt

def plot(u,x,NT,u_analytical, x_analytical, NT2):
      """
      Plots the 1D velocity field
      """

      import matplotlib.pyplot as plt
      import matplotlib.cm as cm
      fig=plt.figure()
      ax=plt.subplot(111)
      colour=iter(cm.rainbow(np.linspace(0,4,NT)))
      for n in range(0,NT,4):
         c=next(colour)
         ax.plot(x,u[:,n],linestyle='-',c=c,label='n='+str(n))
      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')
      box=ax.get_position()
      fig.set_figwidth(12.0)
      fig.set_figheight(8.0)
      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
      plt.ylim([-0.5,1.5])
      plt.xlim([0.0,5.0])
      plt.xlabel('x (m)')
      plt.ylabel('u (m/s)')
      plt.show()

u00, x00, nt00 = analytical(0.5, 101, 1.0, 5.0, 1.0)

u0, x0, nt0 = lax_friedrichs_convection(0.5, 101, 1.0, 5.0, 1.0)

plot(u0,x0,nt0, u00, x00,nt00)