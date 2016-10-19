====================================================
1D First-order Linear Convection - The Wave Equation
====================================================

.. contents::
   :local:

Understand the Problem
======================

* What is the final velocity profile for 1D linear convection when the initial conditions are a square wave and the boundary conditions are constant?

* 1D linear convection is described as follows:

.. math:: {\partial u \over \partial t} + c {\partial u \over \partial x} = 0

Formulate the Problem
=====================

Input Data
~~~~~~~~~~

* `nt` = 51 (number of temporal points)
* `nx` = 21 (number of spatial points)
* `tmax` = 0.5
* `xmax` = 2
* `c` = 1

====================== ========================== ========================= ======================== ===========
x                      i                           t                        n                        u(x,t)
====================== ========================== ========================= ======================== ===========
:math:`0`              :math:`0`                  :math:`0 \le t \le 0.5`   :math:`0 \le i \le 20`   :math:`1`
:math:`0 < x \le 0.5`  :math:`0 < n \le 12.5`     :math:`0`                 :math:`0`                :math:`1`
:math:`0.5 < x \le 1`  :math:`12.5 < n \le 25`    :math:`0`                 :math:`0`                :math:`2`
:math:`1 < x < 2`      :math:`25 < n < 50`        :math:`0`                 :math:`0`                :math:`1`
:math:`2`              :math:`50`                 :math:`0 \le t \le 0.5`   :math:`0 \le i \le 20`   :math:`1`
====================== ========================== ========================= ======================== ===========


Output Data
~~~~~~~~~~~

====================== ========================= =========================
x                      t                         u(x,t)
====================== ========================= =========================
:math:`0 < x < 2`      :math:`0.5`               :math:`?`
====================== ========================= =========================


Design Algorithm to Solve Problem
=================================

Space-time discretisation
~~~~~~~~~~~~~~~~~~~~~~~~~

* i :math:`\rightarrow` index of grid in x
* n :math:`\rightarrow` index of grid in t

Numerical scheme
~~~~~~~~~~~~~~~~

* FD in time
* BD in space

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

::

   #Constants
   nt = 51
   tmax = 0.5
   dt =  tmax/(nt-1) 
   nx =  21
   xmax = 2
   dx = xmax/(nx-1)

   #Boundary Conditions
   for n between 0 and 20
      u(0,n)=1
      u(50,n)=1 
   
   #Initial Conditions
   for i between 1 and 49
      if(12.5 < i < 25)
          u(i,0) = 2
      else
          u(i,0) = 1
   
   #Iteration
   for n between 1 and 20
      u(i,n+1) = u(i,n)-c*(dt/dx)*(u(i,n)-u(i-1,n)
   

Implement Algorithm in Python
=============================

.. plot::
   :include-source:

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

Conclusions
===========

Why isn't the square wave maintained?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* The first order backward differencing scheme in space creates false diffusion.
* If the spatial step is reduced, the error reduces - compare Figure 1 and Figure 2

Why does the wave shift to the right?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* The square wave is being convected by the constant linear wave speed `c`
* For :math:`c > 0` profiles shift to the right by :math:`c \Delta t` - see Figure 2 

What happens at the wall?
~~~~~~~~~~~~~~~~~~~~~~~~~

* As there is no viscosity, there is a non-physical change the the profile near the wall.

