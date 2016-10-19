===================================================================
1D Second-order Non-linear Convection-Diffusion - Burgers' Equation
===================================================================

.. contents::
   :local:

Understand the Problem
======================

* What is the profile for 1D convection-diffusion when the initial conditions are a saw tooth wave and the boundary conditions are periodic?

* How does this compare with the analytical solution?

* 1D convection-diffusion is described as follows:

.. math:: {\partial u \over \partial t} + u {\partial u \over \partial x} = \nu {\partial^2 u \over \partial x^2}

Formulate the Problem
=====================

Input Data
~~~~~~~~~~

**Constants**

* :math:`nt` = 51 (number of temporal points)
* :math:`nx` = 21 (number of spatial points)
* :math:`tmax` = 0.5
* :math:`xmax = 2 \pi`

Different initial and boundary conditions to linear convection:

**Initial Conditions**

.. math:: u_i^{n=0} = -2 \nu {{\partial \phi / \partial x} \over \phi} + 4

where:

.. math:: \phi = exp \left ({{-x^2} \over {4 \nu}} \right ) + exp \left [ -(x-2 \pi)^2 \over {4 \nu} \right ]

.. math:: {{\partial \phi} \over {\partial x}} =
 -{{2x} \over {4 \nu}} exp \left ( {{-x^2} \over {4 \nu}} \right )
 -{2(x-2 \pi) \over {4 \nu}} exp \left [ -(x-2 \pi)^2 \over {4 \nu} \right ]

.. math:: =
 -{{0.5x} \over {\nu}} exp \left ( {{-x^2} \over {4 \nu}} \right )
 -{0.5(x-2 \pi) \over {\nu}} exp \left [ -(x-2 \pi)^2 \over {4 \nu} \right ]

**Boundary Conditions**

Periodic:

.. math:: u_{i=0}^n = u_{i=imax}^n

Output Data
~~~~~~~~~~~

For all :math:`x` and :math:`t`:

.. math:: u(x,t)

Verification
~~~~~~~~~~~~

**Verify the output with the analytical solution**

.. math:: u_i^n = -2 \nu {{\partial \phi / \partial x} \over \phi} + 4

.. math:: \phi = exp \left ({{-(x-4t)^2} \over {4 \nu(t+1)}} \right ) + exp \left [ -(x-4t-2 \pi)^2 \over {4 \nu(t+1)} \right ]

.. math:: {{\partial \phi} \over {\partial x}} =
 -{{2(x-4t)} \over {4 \nu(t+1)}} exp \left ( {{-(x-4t)^2} \over {4 \nu(t+1)}} \right )
 -{{2(x-4t-2 \pi)} \over {4 \nu(t+1)}} exp \left [{-(x-4t-2 \pi)^2} \over {4 \nu(t+1)} \right ]

.. math:: =
 -{{0.5(x-4t)} \over {\nu(t+1)}} exp \left ( {{-(x-4t)^2} \over {4 \nu(t+1)}} \right )
 -{{0.5(x-4t-2 \pi)} \over {\nu(t+1)}} exp \left [{-(x-4t-2 \pi)^2} \over {4 \nu(t+1)} \right ]


Design Algorithm to Solve Problem
=================================

Space-time discretisation
~~~~~~~~~~~~~~~~~~~~~~~~~

* **FD for transient term**
* **BD for convection term**
* **CD for diffusion term**

Numerical scheme
~~~~~~~~~~~~~~~~

* i :math:`\rightarrow` index of grid in x
* n :math:`\rightarrow` index of grid in t

Discrete equation
~~~~~~~~~~~~~~~~~

.. math::

   {{u_i^{n+1} - u_i^n} \over {\Delta t}} + u_i^n {{u_i^n - u_{i-1}^n} \over {\Delta x}} = \nu {{u_{i+1}^n -2u_i^n+ u_{i-1}^n} \over \Delta x^2}

Transpose
~~~~~~~~~

.. math::

   u_i^{n+1} = u_i^n -  u_i^n {\Delta t \over \Delta x} {{(u_i^n - u_{i-1}^n)}} + \nu {\Delta t \over \Delta x^2}(u_{i+1}^n -2u_i^n+ u_{i-1}^n)

   
Pseudo-code
~~~~~~~~~~~

::

   # Constants

     nt = 51
     tmax = 0.5
     dt =  tmax/(nt-1) 
     nx =  21
     xmax = 2
     dx = xmax/(nx-1)
     viscosity = 0.1

   # Range of i is between 0 and nx-1
   # Range of n is between 0 and nt-1

   # This allows the number of points to be nx and nt

   # Periodic Boundary Conditions
   # Create points outside computational domain and set them to their equivalent within the computational domain

     for i between 0 and nx-1
        x(i) = i*dx
        ipos(i) = i+1
        ineg(i) = i-1
   
   # Set Periodicity
   # i:     -1    0,  1,..  nx-2, nx-1, nx 
   # ipos:            start =>    =>    end 
   # ineg:  start =>  =>    end      

     ipos(nx-1) = 0    i.e. nx = 0
     ineg(0) = nx-1    i.e. -1 = nx-1
   
   # Initial Conditions
     for i between 0 and nx-1
        phi = exp( -x(i)^2/(4*vis) ) + exp( -(x(i)-2*pi)^2 / (4*vis) )

        dphi = -(0.5/vis)*exp( -(x^2) / (4*vis) ) -
               (0.5*(x-2*pi) / vis )*exp(-(x-2*pi)^2 / (4*vis) )

        u(i,0) = -2*vis(dphi/phi) + 4
   
   # Analytical Solution (this loop is not time marching, 
     and has no initial conditions, so runs the full range of n)
    
     for n between 0 and nt-1

        t = n*dt
 
        for i between 0 and nx-1
           phi = exp( -(x(i)-4*t)^2/(4*vis*(t+1)) ) + exp( -(x(i)-4*t-2*pi)^2/(4*vis*(t+1)) )

           dphi = -0.5*(x(i)-4*t)^2/(vis*(t+1))*exp( -(x(i)-4*t)^2/(4*vis*(t+1)) )
                  -0.5*(x(i)-4*t-2*pi)^2/(vis*(t+1))*exp( -(x(i)-4*t-2*pi)^2/(4*vis*(t+1)) )

           u_analytical(i,n) = -2*vis(dphi/phi) + 4        

   # Numerical Computation (this loop is time marching, so stops one before the end)

     for n between 1 and nt-2 
         for i between 0 and nx-1
             u(i,n+1) = u(i,n) - u(i,n)*(dt/dx)*(u(i,n)-u(ineg(i),n))+
                        viscosity*(dt/dx^2)*(u(ipos(i),n)-2*u(i,n)+u(ineg(i),n))
           

Implement Algorithm in Python
=============================

* Constants are shown with CAPITAL letters
* Aliases are used for imported functions and constants to shorten formulae

.. plot::
   :include-source:

   from math import pi as PI
   from math import exp as exp



   def analytical_solution(NT, NX, TMAX, XMAX, NU):
      """
      Returns the velocity field and distance for the analytical solution
      """

      # Increments
      DT = TMAX/(NT-1)
      DX = XMAX/(NX-1)

      # Initialise data structures
      import numpy as np
      u_analytical = np.zeros((NX,NT))
      x = np.zeros(NX)
      t = np.zeros(NT)
    
      # Distance
      for i in range(0,NX):
          x[i] = i*DX

      # Analytical Solution 
      for n in range(0,NT):
          t = n*DT
 
          for i in range(0,NX):
              phi = exp( -(x[i]-4*t)**2/(4*NU*(t+1)) ) + exp( -(x[i]-4*t-2*PI)**2/(4*NU*(t+1)) )

              dphi = ( -0.5*(x[i]-4*t)/(NU*(t+1))*exp( -(x[i]-4*t)**2/(4*NU*(t+1)) )
                  -0.5*(x[i]-4*t-2*PI)/(NU*(t+1))*exp( -(x[i]-4*t-2*PI)**2/(4*NU*(t+1)) ) )

              u_analytical[i,n] = -2*NU*(dphi/phi) + 4       

      return u_analytical, x
 
   def convection_diffusion(NT, NX, TMAX, XMAX, NU):
      """
      Returns the velocity field and distance for 1D non-linear convection-diffusion
      """

      # Increments
      DT = TMAX/(NT-1)
      DX = XMAX/(NX-1)

      # Initialise data structures
      import numpy as np
      u = np.zeros((NX,NT))
      u_analytical = np.zeros((NX,NT))
      x = np.zeros(NX)
      t = np.zeros(NT)
      ipos = np.zeros(NX)
      ineg = np.zeros(NX)

      # Periodic boundary conditions
      for i in range(0,NX):
          x[i] = i*DX
          ipos[i] = i+1
          ineg[i] = i-1
  
      ipos[NX-1] = 0
      ineg[0] = NX-1

      # Initial conditions  
      for i in range(0,NX):
          phi = exp( -(x[i]**2)/(4*NU) ) + exp( -(x[i]-2*PI)**2 / (4*NU) )
          dphi = -(0.5*x[i]/NU)*exp( -(x[i]**2) / (4*NU) ) - (0.5*(x[i]-2*PI) / NU )*exp(-(x[i]-2*PI)**2 / (4*NU) )
          u[i,0] = -2*NU*(dphi/phi) + 4      
       
      # Numerical solution
      for n in range(0,NT-1): 
          for i in range(0,NX):
              u[i,n+1] = (u[i,n]-u[i,n]*(DT/DX)*(u[i,n]-u[ineg[i],n])+
                         NU*(DT/DX**2)*(u[ipos[i],n]-2*u[i,n]+u[ineg[i],n]))

      return u, x

   def plot_diffusion(u_analytical,u,x,NT,TITLE):
      """
      Plots the 1D velocity field
      """

      import matplotlib.pyplot as plt
      import matplotlib.cm as cm
      plt.figure()
      ax=plt.subplot(111)
      colour=iter(cm.rainbow(np.linspace(0,20,NT)))   
      for n in range(0,NT,20):
         c=next(colour)
         ax.plot(x,u[:,n],'ko', markerfacecolor='none', alpha=0.5, label='i='+str(n)+' numerical')
         ax.plot(x,u_analytical[:,n],linestyle='-',c=c,label='i='+str(n)+' analytical')
      box=ax.get_position()
      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
      plt.xlabel('x (radians)')
      plt.ylabel('u (m/s)')
      plt.ylim([0,8.0])
      plt.xlim([0,2.0*PI])
      plt.title(TITLE)
      plt.show()

   u,x = convection_diffusion(151, 151, 0.5, 2.0*PI, 0.1)
   u_analytical,x = analytical_solution(151, 151, 0.5, 2.0*PI, 0.1)
   plot_diffusion(u_analytical,u,x,151,'Figure 1: nu=0.1, nt=151, nx=151, tmax=0.5s')

   u,x = convection_diffusion(151, 151, 0.5, 2.0*PI, 0.01)
   u_analytical,x = analytical_solution(151, 151, 0.5, 2.0*PI, 0.01)
   plot_diffusion(u_analytical,u,x,151,'Figure 2: nu=0.01, nt=151, nx=151, tmax=0.5s')

Conclusions
===========

Why doesn't the numerical simulation agree with the analytical solution exactly?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* The numerical solution shows more dissipation through time and space than the analytical solution, despite the fact that the viscosity is the same in both cases (a lot in time, perhaps less in space)
* It is likely that **numerical dissipation** is the cause of the difference between the analytic and numerical solutions
* When physical viscosity is reduced, so reducing physical dissipation, the effect of numerical dissipation is seen more clearly (compare Figure 1 and Figure 2)
* This is caused by the truncation error in the discretisation of the governing equations
* The first order approximations for the transient and convection terms contain numerical diffusion in time and space respectively - probably need to use a higher than first order method for the transient term at least.




