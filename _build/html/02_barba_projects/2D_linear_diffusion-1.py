def diffusion(nt, nx, ny, tmax, xmax, ymax, nu):
   """
   Returns the velocity field and distance for 2D diffusion
   """
   # Increments
   dt = tmax/(nt-1)
   dx = xmax/(nx-1)
   dy = ymax/(ny-1)

   # Initialise data structures
   import numpy as np
   u = np.zeros(((nx,ny,nt)))
   v = np.zeros(((nx,ny,nt)))
   x = np.zeros(nx)
   y = np.zeros(ny)

   # Boundary conditions
   u[0,:,:] = u[nx-1,:,:] = u[:,0,:] = u[:,ny-1,:] = 1
   v[0,:,:] = v[nx-1,:,:] = v[:,0,:] = v[:,ny-1,:] = 1

   # Initial conditions
   u[:,:,:] = v[:,:,:] = 1
   u[(nx-1)/4:(nx-1)/2,(ny-1)/4:(ny-1)/2,0] = 2
   v[(nx-1)/4:(nx-1)/2,(ny-1)/4:(ny-1)/2,0] = 2

   # Loop
   for n in range(0,nt-1):
      for i in range(1,nx-1):
         for j in range(1,ny-1):
            u[i,j,n+1] = ( u[i,j,n]+dt*nu*( ( u[i-1,j,n]-2*u[i,j,n]+u[i+1,j,n] ) /dx**2 +
                                          ( u[i,j-1,n]-2*u[i,j,n]+u[i,j+1,n] ) /dy**2 ) )
            v[i,j,n+1] = ( v[i,j,n]+dt*nu*( ( v[i-1,j,n]-2*v[i,j,n]+v[i+1,j,n] ) /dx**2 +
                                          ( v[i,j-1,n]-2*v[i,j,n]+v[i,j+1,n] ) /dy**2 ) )

   # X Loop
   for i in range(0,nx):
      x[i] = i*dx

   # Y Loop
   for j in range(0,ny):
      y[j] = j*dy

   return u, v, x, y

def plot_3D(u,x,y,time,title,label):
   """
   Plots the 2D velocity field
   """

   import matplotlib.pyplot as plt
   from mpl_toolkits.mplot3d import Axes3D
   fig=plt.figure(figsize=(11,7),dpi=100)
   ax=fig.gca(projection='3d')
   ax.set_xlabel('x (m)')
   ax.set_ylabel('y (m)')
   ax.set_zlabel(label)
   X,Y=np.meshgrid(x,y)
   surf=ax.plot_surface(X,Y,u[:,:,time],rstride=2,cstride=2)
   plt.title(title)
   plt.show()

u,v,x,y = diffusion(151, 51, 51, 0.5, 2.0, 2.0, 0.1)

plot_3D(u,x,y,0,'Figure 1: nu=0.1, nt=151, nx=51, ny=51, t=0sec','u (m/s)')
plot_3D(u,x,y,100,'Figure 2: nu=0.1, nt=151, nx=51, ny=51, t=0.5sec','u (m/s)')
plot_3D(v,x,y,0,'Figure 3: nu=0.1, nt=151, nx=51, ny=51, t=0sec','v (m/s)')
plot_3D(v,x,y,100,'Figure 4: nu=0.1, nt=151, nx=51, ny=51, t=0.5sec','v (m/s)')