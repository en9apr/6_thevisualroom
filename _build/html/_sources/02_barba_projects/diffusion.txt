====================================================
1D Second-order Linear Diffusion - The Heat Equation
====================================================

.. contents::
   :local:

Understand the Problem
======================

* What is the final temperature profile for 1D diffusion when the initial conditions are a square wave and the boundary conditions are constant?

* 1D diffusion is described as follows:

.. math:: {\partial u \over \partial t} = \nu {\partial^2 u \over \partial x^2}

* Consider looking for a solution of type:

.. math:: u = \hat u e^{i(kx-\omega t)}

* Represents a wave of amplitude :math:`\hat u`, :math:`\omega = 2 \pi f`

* Introducing into PDE, we obtain:

.. math:: i \omega = \nu k^2

* Leading to a solution:

.. math::  u = \hat u e^{ikx}e^{- \nu k^2 t}

where :math:`e^{\nu k^2 t}` is the exponential damping term. So diffusion is an exponentially damped wave.

Note: :math:`\nu > 0` for physical diffusion (if :math:`\nu < 0` would represent an exponentially growing phenomenon, e.g. an explosion or 'the rich get richer' model)

The physics of diffusion are:

* An expotentially damped wave in time
* Isotropic in space - the same in all spatial directions - it does not distinguish between upstream and downstream

The phenomenon of diffusion is isotropic - so the finite difference formula that represents that physics is **central differencing** CD, because CD takes values from upstream and downstream equally.
 
Formulate the Problem
=====================

* Same as Linear Convection


Design Algorithm to Solve Problem
=================================

Space-time discretisation
~~~~~~~~~~~~~~~~~~~~~~~~~

* FD in time
* **CD in space**

Numerical scheme
~~~~~~~~~~~~~~~~

* i :math:`\rightarrow` index of grid in x
* n :math:`\rightarrow` index of grid in t

Discrete equation
~~~~~~~~~~~~~~~~~

.. math::

   {{u_i^{n+1} - u_i^n} \over {\Delta t}} = \nu {{u_{i+1}^n -2u_i^n+ u_{i-1}^n} \over \Delta x^2}

Transpose
~~~~~~~~~

.. math::

   u_i^{n+1} = u_i^n + \nu {\Delta t \over \Delta x^2}(u_{i+1}^n -2u_i^n+ u_{i-1}^n)

   
Pseudo-code
~~~~~~~~~~~

* Very similar to Linear Convection

Implement Algorithm in Python
=============================

* Very similar to Linear Convection

.. plot::

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

Conclusions
===========

Why isn't the square wave maintained?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* The square wave isn't maintained because the system is attempting to reach equilibrium - the rate of change of velocity being equal to the shear force per unit mass. There are no external forces and no convective acceleration terms.

Why does increasing the viscosity, spatial points and time period cause instability?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the viscosity is too large, or if the number of spatial points is too large or if the timestep is too large, then the central differencing method becomes unstable. This is due to the ratio, r. If r is too large, the method becomes unstable:

.. math::

   r = \nu {\Delta t \over (\Delta x)^2} 

