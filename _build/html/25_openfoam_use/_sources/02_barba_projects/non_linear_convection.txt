=====================================================================
1D First-order Non-Linear Convection - The Inviscid Burgers' Equation
=====================================================================

.. contents::
   :local:

Understand the Problem
======================

* What is the final velocity profile for 1D non-linear convection when the initial conditions are a square wave and the boundary conditions are constant?

* 1D non-linear convection is described as follows:

.. math:: {\partial u \over \partial t} + u {\partial u \over \partial x} = 0

* This equation is capable of generating discontinuities (shocks)

Formulate the Problem
=====================

* Same as Linear Convection


Design Algorithm to Solve Problem
=================================

Space-time discretisation
~~~~~~~~~~~~~~~~~~~~~~~~~

* Same as Linear Convection

Numerical scheme
~~~~~~~~~~~~~~~~

* Same as Linear Convection

Discrete equation
~~~~~~~~~~~~~~~~~

.. math::

   {{u_i^{n+1} - u_i^n} \over {\Delta t}} + c {{u_i^n - u_{i-1}^n} \over \Delta x}=0 

Transpose
~~~~~~~~~

.. math::

   u_i^{n+1} = u_i^n - c{\Delta t \over \Delta x}(u_i^n - u_{i-1}^n)
   
Pseudo-code
~~~~~~~~~~~

* Very similar to Linear Convection

Implement Algorithm in Python
=============================

* Very similar to Linear Convection

.. plot::

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
            u[i,n+1] = u[i,n]-u[i,n]*(dt/dx)*(u[i,n]-u[i-1,n])

      # X Loop
      for i in range(0,nx):
         x[i] = i*dx

      return u, x

   def plot_convection(u,x,nt,title):
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
         plt.ylim([0,2.2])
         plt.title(title)
         plt.show()

   u,x = convection(151, 51, 0.5, 2.0, 0.5)
   plot_convection(u,x,151,'Figure 1: c=0.5m/s, nt=151, nx=51, tmax=0.5s')

   u,x = convection(151, 302, 0.5, 2.0, 0.5)
   plot_convection(u,x,151,'Figure 2: c=0.5m/s, nt=151, nx=302, tmax=0.5s')

   u,x = convection(151, 51, 2.0, 2.0, 0.5)
   plot_convection(u,x,151,'Figure 3: c=0.5m/s, nt=151, nx=51, tmax=2s')

Conclusions
===========

Why isn't the square wave maintained?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* The first order backward differencing scheme in space still creates false diffusion as before.
* However, due to the non-linearity in the governing equation, if the spatial step is reduced, the solution can develop shocks, see Figure 2.
* Clearly a square wave is not best represented with the inviscid Burgers Equation.

Why does the wave shift to the right?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* The square wave is being convected by the velocity, `u` which is not constant.
* The greatest shift is where the velocity is greatest, see Figure 1

What happens at the wall?
~~~~~~~~~~~~~~~~~~~~~~~~~

* As there is no viscosity, there is a non-physical change in the profile near the wall, see Figure 3.
* Comparing this with the linear example, there is clearly much more numerical diffusion in the non-linear example, perhaps due to the convective term being larger, causing a greater magnitude in numerical diffusion.

