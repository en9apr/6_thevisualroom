def diffusion(nt, nx, tmax, xmax, nu):
   """
   Returns the velocity field and distance for 1D linear convection
   """
   # Increments
   dt = tmax/(nt-1)
   dx = xmax/(nx-1)

   # Initialise data structures
   import numpy as np
   u = np.zeros((nx,nt))
   x = np.zeros(nx)

   # Boundary conditions
   u[0,:] = u[nx-1,:] = 1

   # Initial conditions
   for i in range(1,nx-1):
      if(i > (nx-1)/4 and i < (nx-1)/2):
         u[i,0] = 2
      else:
         u[i,0] = 1

   # Loop
   for n in range(0,nt-1):
      for i in range(0,nx-1):
         u[i,n+1] = u[i,n] + nu*(dt/dx**2.0)*(u[i+1,n]-2.0*u[i,n]+u[i-1,n])

   # X Loop
   for i in range(0,nx):
      x[i] = i*dx

   return u, x

def plot_diffusion(u,x,nt,title):
   """
   Plots the 1D velocity field
   """

   import matplotlib.pyplot as plt
   import matplotlib.cm as cm
   plt.figure()
   colour=iter(cm.rainbow(np.linspace(0,10,nt)))
   for i in range(0,nt,10):
      c=next(colour)
      plt.plot(x,u[:,i],c=c)
   plt.xlabel('x (m)')
   plt.ylabel('u (m/s)')
   plt.ylim([0,3.0])
   plt.title(title)
   plt.show()

u,x = diffusion(151, 51, 0.5, 2.0, 0.1)
plot_diffusion(u,x,151,'Figure 1: nu=0.1, nt=151, nx=51, tmax=0.5s')

u,x = diffusion(151, 51, 0.5, 2.0, 0.242)
plot_diffusion(u,x,151,'Figure 1b: nu=0.242, nt=151, nx=51, tmax=0.5s')

u,x = diffusion(151, 79, 0.5, 2.0, 0.1)
plot_diffusion(u,x,151,'Figure 2: nu=0.1, nt=151, nx=79, tmax=0.5s')

u,x = diffusion(151, 51, 1.217, 2.0, 0.1)
plot_diffusion(u,x,151,'Figure 3: nu=0.1, nt=151, nx=51, tmax=1.217s')