def navier_stokes_initialisation(niter, r, nx_or_ny, tmax, xmax_or_ymax):
    """
    Returns the velocity field and distance for 2D linear convection
    """
    # Increments:
    nx = ny = nx_or_ny
    xmax = ymax = xmax_or_ymax
    dx = xmax/(nx-1)
    dy = ymax/(ny-1)
    nt = int((tmax / (r*(dx)**2))+1)
    dt = tmax/(nt-1)

    # Initialise data structures:
    import numpy as np
    p = np.zeros((nx,ny))
    u = np.zeros((nx,ny))
    v = np.zeros((nx,ny))

    # linspace is SIMPLER than list comprehensions:
    x = np.linspace(0.0,2.0,nx)
    y = np.linspace(0.0,2.0,ny)

    # Pressure Boundary Conditions:
    p[:, ny-1] = 0.0

    # Velocity Boundary Conditions:
    u[:,ny-1] = 1.0

    return p, x, y, u, v, nx, ny, nt, dx, dy, dt, niter, r

def navier_stokes(rho, nu, niter, r, nx, tmax, xmax):

    (p, x, y, u, v, nx, ny, nt,
    dx, dy, dt, niter, r) = navier_stokes_initialisation(niter, r, nx, tmax, xmax)

    # Increments
    h = dx

    import numpy as np

    # Intermediate copies:
    un = np.zeros((nx, ny))
    vn = np.zeros((nx, ny))
    pm = np.zeros((nx, ny))
    bm = np.zeros((nx, ny)) # bm needs to be exactly zero at the boundaries

    # Loop - use decimal points for all floating point numbers
    for n in range(nt):

        # We know the velocity at i=0, j=0, i=nx-1 and j=ny-1. b is zero at the boundaries.
        bm[1:-1, 1:-1] = ( (rho / (2.0 * h * dt)) * ( u[2:, 1:-1] - u[0:-2, 1:-1]
        + v[1:-1, 2:] - v[1:-1, 0:-2] ) +
        (rho / (4.0*h**2)) * ( (u[2:, 1:-1] - u[0:-2, 1:-1])**2.0 +
        4.0*h*(u[1:-1, 2:] - u[1:-1, 0:-2])*(v[2:, 1:-1] - v[0:-2, 1:-1]) +
        (v[1:-1, 2:] - v[1:-1, 0:-2])**2.0 ) )

        # First points for p. We don't know the pressure at i=0, j=0 and i=nx-1. We DO know the pressure at j=ny-1
        for m in range(niter):

            pm = np.copy(p)
            p[1:-1, 1:-1] = 0.25*( pm[2:, 1:-1] + pm[0:-2, 1:-1] + pm[1:-1, 2:] + pm[1:-1, 0:-2]
            - bm[1:-1, 1:-1]*h**2.0 )

            # Set zero gradient boundary conditions:
            p[0, :] = p[1, :]
            p[:, 0] = p[:, 1]
            p[-1, :] = p[-2, :]

        # First points for u and v. We know the velocity at i=0, j=0, i=nx-1 and j=ny-1.
        # We are simply using the value of pressure here
        un = np.copy(u)
        vn = np.copy(v)

        u[1:-1, 1:-1] = ( un[1:-1, 1:-1] -
        (dt / h) * ( un[1:-1, 1:-1] * ( un[1:-1, 1:-1] - un[0:-2, 1:-1] ) +
        vn[1:-1, 1:-1] * ( un[1:-1, 1:-1] - un[1:-1, 0:-2] ) ) -
        (dt / (2.0 * rho * h)) * ( p[2:, 1:-1] - p[0:-2, 1:-1] ) +
        (dt * nu / h**2.0) * ( un[0:-2, 1:-1] + un[2:, 1:-1] + un[1:-1, 0:-2] + un[1:-1, 2:] -
        4.0 * un[1:-1, 1:-1] ) )

        v[1:-1, 1:-1] = ( vn[1:-1, 1:-1] -
        (dt / h) * ( un[1:-1, 1:-1] * ( vn[1:-1, 1:-1] - vn[0:-2, 1:-1] ) +
        vn[1:-1, 1:-1] * ( vn[1:-1, 1:-1] - vn[1:-1, 0:-2] ) ) -
        (dt / (2.0 * rho * h)) * ( p[1:-1, 2:] - p[1:-1, 0:-2] ) +
        (dt * nu / h**2.0) * ( vn[0:-2, 1:-1] + vn[2:, 1:-1] + vn[1:-1, 0:-2] + vn[1:-1, 2:] -
        4.0 * vn[1:-1, 1:-1] ) )

    return u, v, p, x, y

u33, v33, p33, x33, y33 = navier_stokes(1.0, 0.01, 50, 0.5, 41, 5.0, 2.0)

def vector_contour_2(u, v, p, x, y):
    fig = plt.figure(figsize=(11,9), dpi=100)
    Y,X=np.meshgrid(y,x) #note meshgrid uses y,x not x,y!!!
    plt.pcolor(X,Y,p)
    plt.colorbar()
    #plt.contourf(X,Y,p[:,:],alpha=0.5)    ###plotting the pressure field as a contour
    # plt.colorbar()
    # plt.contour(X,Y,p[:,:])               ###plotting the pressure field outlines
    plt.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2]) ##plotting velocity
    #plt.contour(X,Y,p) ##plotting velocity
    plt.ylim([0,2.0])
    plt.xlim([0,2.0])
    plt.xlabel('X')
    plt.ylabel('Y')

vector_contour_2(u33, v33, p33, x33, y33)