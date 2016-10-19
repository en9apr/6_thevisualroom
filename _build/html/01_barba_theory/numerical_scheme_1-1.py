def diffusion(nt, nx, tmax, xmax, sigma, method):
   """
   Returns the velocity field and distance for 2D diffusion
   """
   # Increments
   dt = tmax/(nt-1)
   dx = xmax/(nx-1)

   # Compute c (given sigma)
   c = sigma * dx / dt

   # Initialise data structures
   import numpy as np
   u = np.zeros((nx,nt))
   x = np.zeros(nx)

   # X Loop
   for i in range(0,nx):
      x[i] = i*dx

   # Boundary conditions
   u[0,:] = u[nx-1,:] = 0

   # Initial conditions
   for i in range(1,nx-1):
      if(0.9<=x[i] and x[i]<=1):
         u[i,0] = 10*(x[i]-0.9)
      elif(1<x[i] and x[i]<=1.1):
         u[i,0] = 10*(1.1-x[i])
      else:
         u[i,0] = 0

   # Loop
   for n in range(0,nt-1):
      for i in range(1,nx-1):
         if(method=='BD'):
            u[i,n+1] = u[i,n]-dt*c*( ( u[i,n]-u[i-1,n] ) /dx )
         elif(method=='CD'):
            u[i,n+1] = u[i,n]-dt*c*( ( u[i+1,n]-u[i-1,n] ) / (2*dx) )

   return u, x

def plot_diffusion(u,x,nt,title):
   """
   Plots the 1D velocity field
   """
   import matplotlib.pyplot as plt
   import matplotlib.cm as cm
   plt.figure()
   colour=iter(cm.rainbow(np.linspace(0,1,nt)))
   for n in range(0,nt,1):
      c=next(colour)
      plt.plot(x,u[:,n],c=c)
      plt.xlabel('x (m)')
      plt.ylabel('u (m/s)')
      plt.title(title)
      plt.show()

u,x = diffusion(6,41, 0.2, 2.0, 0.8, 'CD')
plot_diffusion(u,x,6,'Figure 1: sigma=0.8, Central Differencing in x')

u,x = diffusion(6,41, 0.2, 2.0, 0.8, 'BD')
plot_diffusion(u,x,6,'Figure 1: sigma=0.8, Backward Differencing in x')

u,x = diffusion(6,41, 0.2, 2.0, 1.5, 'BD')
plot_diffusion(u,x,6,'Figure 1: sigma=1.5, Backward Differencing in x')