def initial_and_boundary_conditions_2(xmax, nx, nt):

    # Initialise data structures
    import numpy as np
    u = np.zeros((nx,nt))

    # Boundary conditions
    u[0,:] = 0.0
    u[nx-1,:] = 0.0

    # X Loop
    x = np.zeros(nx)
    x = np.linspace(0.0,xmax,nx)

    quarter = int((nx-1)/4)
    from math import pi as PI
    i=0
    # Initial conditions
    u[i:quarter,0] = 1.0*np.sin(12.0 * PI * (x[i:quarter]))

    return u, x

def analytical_2(sigma, nx, tmax, xmax, c):

    # Increments
    # dt = tmax/(nt-1)
    # dx = (c * dt) / sigma
    dx = xmax/(nx-1)
    dt =  ((dx*sigma) / c)
    nt = int(sigma*tmax/dt + 1)

    # Initial and Boundary Conditions
    u, x = initial_and_boundary_conditions_2(xmax, nx, nt)

    i = 1
    # Loop - must loop as NumPy cannot be used for dependency reasons
    for n in range(0,nt-1):
        u[i:nx-1, n+1] = u[i-1:nx-2, n]

    return u, x, nt

def lax_wendroff_convection_2(sigma, nx, tmax, xmax, c):
    """
    Returns the velocity field and distance for 1D linear convection
    """
    # Increments
    dx = xmax/(nx-1)
    dt =  ((dx * sigma) / c)
    nt = int(tmax/dt + 1)

    # Initial and Boundary Conditions
    u, x = initial_and_boundary_conditions_2(xmax, nx, nt)

    i = 1
    # Loop - must loop as NumPy cannot be used for dependency reasons
    for n in range(0,nt-1):
        u[i:nx-1, n+1] = (u[i:nx-1, n] - 0.5*sigma*(u[i+1:nx, n] - u[i-1:nx-2, n]) +
                         0.5*(sigma**2)*(u[i+1:nx, n] - 2.0*u[i:nx-1, n] + u[i-1:nx-2, n]))

    return u, x, nt

def plot_2(u,x,NT,u_analytical, x_analytical, NT2):
      """
      Plots the 1D velocity field
      """

      import matplotlib.pyplot as plt
      import matplotlib.cm as cm
      fig=plt.figure()
      ax=plt.subplot(111)
      colour=iter(cm.rainbow(np.linspace(0,100,NT)))
      ax.plot(x,u[:,0],linestyle='-',c='b',label='n='+str(0))
      ax.plot(x,u[:,NT-1],linestyle='-',c='r',label='n='+str(NT-1))
      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')
      box=ax.get_position()
      fig.set_figwidth(12.0)
      fig.set_figheight(8.0)
      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
      plt.ylim([-1.5,1.5])
      plt.xlim([0.0,3.0])
      plt.xlabel('x (m)')
      plt.ylabel('u (m/s)')
      plt.show()

u0000, x0000, nt0000 = analytical_2(1.0, 501, 1.0, 5.0, 1.0)
u000, x000, nt000 = lax_wendroff_convection_2(0.5, 501, 1.0, 5.0, 1.0)
plot_2(u000, x000, nt000,u0000, x0000, nt0000)