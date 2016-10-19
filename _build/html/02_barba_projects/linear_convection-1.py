def convection(nt, nx, tmax, xmax, c):
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
      for i in range(1,nx-1):
         u[i,n+1] = u[i,n]-c*(dt/dx)*(u[i,n]-u[i-1,n])

   # X Loop
   for i in range(0,nx):
      x[i] = i*dx

   return u, x

def plot_convection(u,x,nt,title):
   """
   Plots the 1D velocity field
   """

   import matplotlib.pyplot as plt
   plt.figure()
   for i in range(0,nt,10):
      plt.plot(x,u[:,i],'r')
      plt.xlabel('x (m)')
      plt.ylabel('u (m/s)')
      plt.ylim([0,2.2])
      plt.title(title)
      plt.show()

u,x = convection(151, 51, 0.5, 2.0, 0.5)
plot_convection(u,x,151,'Figure 1: c=0.5m/s, nt=151, nx=51, tmax=0.5s')

u,x = convection(151, 1001, 0.5, 2.0, 0.5)
plot_convection(u,x,151,'Figure 2: c=0.5m/s, nt=151, nx=1001, tmax=0.5s')

u,x = convection(151, 51, 2.0, 2.0, 0.5)
plot_convection(u,x,151,'Figure 3: c=0.5m/s, nt=151, nx=51, tmax=2s')