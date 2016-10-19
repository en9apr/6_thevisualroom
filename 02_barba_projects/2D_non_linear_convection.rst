====================================
2D First-order Non-Linear Convection
====================================

.. contents::
   :local:

Understand the Problem
======================

* What is the final velocity profile for 2D non-linear convection when the initial conditions are a square wave and the boundary conditions are unity?

* 2D non-linear convection is described as follows:

.. math:: {\partial u \over \partial t} + u {\partial u \over \partial x} + v {\partial u \over \partial y} = 0

.. math:: {\partial v \over \partial t} + u {\partial v \over \partial x} + v {\partial v \over \partial y} = 0

Formulate the Problem
=====================

Input Data
~~~~~~~~~~

* `nt` = 51 (number of temporal points)
* `tmax` = 0.5
* `nx` = 21 (number of x spatial points)
* `nj` = 21 (number of y spatial points)
* `xmax` = 2
* `ymax` = 2

**Initial Conditions:** :math:`t=0`

====================== ========================== ========================= ======================== ===================
x                      i                           y                        j                        u(x,y,t), v(x,y,t)
====================== ========================== ========================= ======================== ===================
:math:`0`              :math:`0`                  :math:`0`                 :math:`0`                :math:`1`
:math:`0 < x \le 0.5`  :math:`0 < i \le 5`        :math:`0 < y \le 0.5`     :math:`0 < j \le 5`      :math:`1`
:math:`0.5 < x \le 1`  :math:`5 < i \le 10`       :math:`0.5 < y \le 1`     :math:`5 < j \le 10`     :math:`2`
:math:`1 < x < 2`      :math:`10 < i < 20`        :math:`1 < y < 2`         :math:`10 < j < 20`      :math:`1`
:math:`2`              :math:`20`                 :math:`2`                 :math:`20`               :math:`1`
====================== ========================== ========================= ======================== ===================

**Boundary Conditions:** :math:`x=0` and :math:`x=2`, :math:`y=0` and :math:`y=2`

========================= ======================== ===================
t                         n                        u(x,y,t), v(x,y,t)
========================= ======================== ===================
:math:`0 \le t \le 0.5`   :math:`0 \le n \le 50`   :math:`1`
========================= ======================== ===================

Output Data
~~~~~~~~~~~

========================= =========================== ========================= =========================
x                         y                           t                         u(x,y,t), v(x,y,t)
========================= =========================== ========================= =========================
:math:`0 \le x \le 2`     :math:`0 \le y \le 2`       :math:`0 \le t \le 0.5`   :math:`?`
========================= =========================== ========================= =========================


Design Algorithm to Solve Problem
=================================

Space-time discretisation
~~~~~~~~~~~~~~~~~~~~~~~~~

* i :math:`\rightarrow` index of grid in x
* j :math:`\rightarrow` index of grid in y
* n :math:`\rightarrow` index of grid in t

Numerical scheme
~~~~~~~~~~~~~~~~

* FD in time
* BD in space

Discrete equation
~~~~~~~~~~~~~~~~~

.. math::

   {{u_{i,j}^{n+1} - u_{i,j}^n} \over {\Delta t}} + 
   u_{i,j}^n {{u_{i,j}^n - u_{i-1,j}^n} \over \Delta x} + 
   v_{i,j}^n {{u_{i,j}^n - u_{i,j-1}^n} \over \Delta y} = 0 

.. math::

   {{v_{i,j}^{n+1} - v_{i,j}^n} \over {\Delta t}} + 
   u_{i,j}^n {{v_{i,j}^n - v_{i-1,j}^n} \over \Delta x} + 
   v_{i,j}^n {{v_{i,j}^n - v_{i,j-1}^n} \over \Delta y} = 0 

Transpose
~~~~~~~~~

.. math::

   u_{i,j}^{n+1} = u_{i,j}^n - \Delta t \left( {u_{i,j}^n \over \Delta x} (u_{i,j}^n - u_{i-1,j}^n)+ 
                   {v_{i,j}^n \over \Delta y} (u_{i,j}^n - u_{i,j-1}^n) \right)

.. math::

   v_{i,j}^{n+1} = v_{i,j}^n - \Delta t \left( {u_{i,j}^n \over \Delta x} (v_{i,j}^n - v_{i-1,j}^n)+ 
                   {v_{i,j}^n \over \Delta y} (v_{i,j}^n - v_{i,j-1}^n) \right)
   
Pseudo-code
~~~~~~~~~~~

::

   #Constants
   nt = 51
   tmax = 0.5
   dt =  tmax/(nt-1) 
   nx =  21
   xmax = 2
   ny = 21
   ymax = 2
   dx = xmax/(nx-1)
   dy = ymax/(ny-1)

   #Boundary Conditions
   #u boundary, bottom and top
   u(0:20,0,0:50)=u(0:20,20,0:50)=1
   #v boundary, bottom and top
   v(0:20,0,0:50)=v(0:20,20,0:50)=1 
   #u boundary left and right
   u(0,0:20,0:50)=u(20,0:20,0:50)=1
   #v boundary left and right
   v(0,0:20,0:50)=v(20,0:20,0:50)=1
   
   #Initial Conditions
   u(:,:,0) = 1
   v(:,:,0) = 1
   u(5:10,5:10,0) = 2
   v(5:10,5:10,0) = 2
   
   #Iteration
   for n between 0 and 49
      for i between 1 and 19
         for j between 1 and 19
             u(i,j,n+1) = u(i,j,n)-dt*[(u(i,j,n)/dx)*(u(i,j,n)-u(i-1,j,n))+
                                       (v(i,j,n)/dy)*(u(i,j,n)-u(i,j-1,n))]
             v(i,j,n+1) = v(i,j,n)-dt*[(u(i,j,n)/dx)*(v(i,j,n)-v(i-1,j,n))+
                                       (v(i,j,n)/dy)*(v(i,j,n)-v(i,j-1,n))]

Implement Algorithm in Python
=============================

.. plot::
   :include-source:

   def non_linear_convection(nt, nx, ny, tmax, xmax, ymax):
      """
      Returns the velocity field and distance for 2D linear convection
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
               u[i,j,n+1] = (u[i,j,n]-dt*((u[i,j,n]/dx)*(u[i,j,n]-u[i-1,j,n])+
                                          (v[i,j,n]/dy)*(u[i,j,n]-u[i,j-1,n])))
               v[i,j,n+1] = (v[i,j,n]-dt*((u[i,j,n]/dx)*(v[i,j,n]-v[i-1,j,n])+
                                          (v[i,j,n]/dy)*(v[i,j,n]-v[i,j-1,n])))

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

   u,v,x,y = non_linear_convection(101, 81, 81, 0.5, 2.0, 2.0)

   plot_3D(u,x,y,0,'Figure 1: c=0.5m/s, nt=101, nx=81, ny=81, t=0sec','u (m/s)')
   plot_3D(u,x,y,100,'Figure 2: c=0.5m/s, nt=101, nx=81, ny=81, t=0.5sec','u (m/s)')
   plot_3D(v,x,y,0,'Figure 3: c=0.5m/s, nt=101, nx=81, ny=81, t=0sec','v (m/s)')
   plot_3D(v,x,y,100,'Figure 4: c=0.5m/s, nt=101, nx=81, ny=81, t=0.5sec','v (m/s)')


Conclusions
===========

Why isn't the square wave maintained?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* As with 1D, the first order backward differencing scheme in space creates false diffusion.

Why does the wave shift?
~~~~~~~~~~~~~~~~~~~~~~~~

* The square wave is being convected by the  wave speed `u` and `v` in 2 dimensions
* Profiles shift by :math:`u \Delta t` and :math:`v \Delta t` - compare Figure 1, 2, 3 and 4 
