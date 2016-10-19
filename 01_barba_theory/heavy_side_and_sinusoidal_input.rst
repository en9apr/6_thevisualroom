===============================================
The Diffusion of Heaviside and Sinusoidal Input
===============================================

.. contents::
   :local:

New Schemes for Linear Convection
=================================

* Lax-Friedrichs :math:`\Rightarrow` stabilizes CD, 1st order
* Lax-Wendroff :math:`\Rightarrow` stable CD with an additional dissipative term, 2nd order
* Leapfrod :math:`\Rightarrow` 2nd order, 3 level scheme (it is not self-starting)

Heaviside Function Input
========================

The Heaviside Function represents a **shock wave** in velocity and pressure for compressible flow. One of the great challenges in numerical methods is representing sharp gradients.

* The Wave Equation for these tests is described as follows:

.. math:: {\partial u \over \partial t} + c {\partial u \over \partial x} = 0

Upwind Differencing with :math:`\sigma = 1`, :math:`\sigma = 0.5` and :math:`\sigma = 0.25`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`\sigma = 1`

* If the value of :math:`\sigma = 1` then the new velocity simply equals the old velocity. 
* The phase angle is zero, so there is no diffusion or dispersion error

:math:`\sigma = 0.5`

* The jump is diffused by the numerical diffusion arising from the first order truncation error
* The amount of diffusion is increasing as we move through time (as n increases)
* The numerical speed with which the solution moves through the domain is also reduced

:math:`\sigma = 0.25`

* Clearly, the smaller Courant Number reduces the speed at which the solution travels through the domain
* At these very low frequencies, there is not much more added diffusion between :math:`\sigma = 0.50` and :math:`\sigma = 0.25`

.. plot::

	def initial_and_boundary_conditions(nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 1.0
	    u[nx-1,:] = 0.0

	    half = int((nx-1)/2)

	    # Initial conditions      
	    u[0:half,0] = 1.0
	    
	    return u

	def analytical(sigma, nx, tmax, xmax, c):
	    
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u = initial_and_boundary_conditions(nx, nt)

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def upwind_convection(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u = initial_and_boundary_conditions(nx, nt)

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i:nx-1, n]-sigma*(u[i:nx-1, n]-u[i-1:nx-2, n])

	    return u, x, nt

	def plot(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,4,NT)))   
	      for n in range(0,NT,4):
	         c=next(colour)
	         ax.plot(x,u[:,n],linestyle='-',c=c,label='n='+str(n))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-0.5,1.5])
	      plt.xlim([0.0,5.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show()

	u00, x00, nt00 = analytical(1.0, 101, 1.0, 5.0, 1.0)

	u0, x0, nt0 = upwind_convection(1.0, 101, 1.0, 5.0, 1.0)   

	plot(u0,x0,nt0, u00, x00,nt00)   


.. plot::

	def initial_and_boundary_conditions(nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 1.0
	    u[nx-1,:] = 0.0

	    half = int((nx-1)/2)

	    # Initial conditions      
	    u[0:half,0] = 1.0
	    
	    return u

	def analytical(sigma, nx, tmax, xmax, c):
	    
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u = initial_and_boundary_conditions(nx, nt)

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def upwind_convection(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u = initial_and_boundary_conditions(nx, nt)

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i:nx-1, n]-sigma*(u[i:nx-1, n]-u[i-1:nx-2, n])

	    return u, x, nt

	def plot(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,4,NT)))   
	      for n in range(0,NT,4):
	         c=next(colour)
	         ax.plot(x,u[:,n],linestyle='-',c=c,label='n='+str(n))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-0.5,1.5])
	      plt.xlim([0.0,5.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show()

	u00, x00, nt00 = analytical(0.5, 101, 1.0, 5.0, 1.0)

	u0, x0, nt0 = upwind_convection(0.5, 101, 1.0, 5.0, 1.0)   

	plot(u0,x0,nt0, u00, x00,nt00)   


.. plot::

	def initial_and_boundary_conditions(nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 1.0
	    u[nx-1,:] = 0.0

	    half = int((nx-1)/2)

	    # Initial conditions      
	    u[0:half,0] = 1.0
	    
	    return u

	def analytical(sigma, nx, tmax, xmax, c):
	    
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u = initial_and_boundary_conditions(nx, nt)

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def upwind_convection(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u = initial_and_boundary_conditions(nx, nt)

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i:nx-1, n]-sigma*(u[i:nx-1, n]-u[i-1:nx-2, n])

	    return u, x, nt

	def plot(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,4,NT)))   
	      for n in range(0,NT,4):
	         c=next(colour)
	         ax.plot(x,u[:,n],linestyle='-',c=c,label='n='+str(n))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-0.5,1.5])
	      plt.xlim([0.0,5.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show()

	u00, x00, nt00 = analytical(0.25, 101, 1.0, 5.0, 1.0)

	u0, x0, nt0 = upwind_convection(0.25, 101, 1.0, 5.0, 1.0)   

	plot(u0,x0,nt0, u00, x00,nt00)   


Lax-Friedrichs with :math:`\sigma = 0.5`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Numerical dissipation (more than upwind scheme) and odd-even decoupling
* Amount of diffusion is still increasing with increasing n


.. plot::

	def initial_and_boundary_conditions(nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 1.0
	    u[nx-1,:] = 0.0

	    half = int((nx-1)/2)

	    # Initial conditions      
	    u[0:half,0] = 1.0
	    
	    return u

	def analytical(sigma, nx, tmax, xmax, c):
	    
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u = initial_and_boundary_conditions(nx, nt)

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def lax_friedrichs_convection(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u = initial_and_boundary_conditions(nx, nt)

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)

	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = 0.5*(u[i-1:nx-2, n]+u[i+1:nx, n])-0.5*sigma*(u[i+1:nx, n]-u[i-1:nx-2, n])

	    return u, x, nt

	def plot(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,4,NT)))   
	      for n in range(0,NT,4):
	         c=next(colour)
	         ax.plot(x,u[:,n],linestyle='-',c=c,label='n='+str(n))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-0.5,1.5])
	      plt.xlim([0.0,5.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show()

	u00, x00, nt00 = analytical(0.5, 101, 1.0, 5.0, 1.0)

	u0, x0, nt0 = lax_friedrichs_convection(0.5, 101, 1.0, 5.0, 1.0) 

	plot(u0,x0,nt0, u00, x00,nt00)

Lax-Wendroff with :math:`\sigma = 0.5`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* This more accurately represents the step change
* However, there is an oscillatory response

.. plot::

	def initial_and_boundary_conditions(nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 1.0
	    u[nx-1,:] = 0.0

	    half = int((nx-1)/2)

	    # Initial conditions      
	    u[0:half,0] = 1.0
	    
	    return u

	def analytical(sigma, nx, tmax, xmax, c):
	    
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u = initial_and_boundary_conditions(nx, nt)

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def lax_wendroff_convection(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u = initial_and_boundary_conditions(nx, nt)

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)

	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = (u[i:nx-1, n] - 0.5*sigma*(u[i+1:nx, n] - u[i-1:nx-2, n]) +
	                         0.5*(sigma**2)*(u[i+1:nx, n] - 2.0*u[i:nx-1, n] + u[i-1:nx-2, n]))

	    return u, x, nt

	def plot(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,4,NT)))   
	      for n in range(0,NT,4):
	         c=next(colour)
	         ax.plot(x,u[:,n],linestyle='-',c=c,label='n='+str(n))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-0.5,1.5])
	      plt.xlim([0.0,5.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show()

	u00, x00, nt00 = analytical(0.5, 101, 1.0, 5.0, 1.0)

	u0, x0, nt0 = lax_wendroff_convection(0.5, 101, 1.0, 5.0, 1.0)   

	plot(u0,x0,nt0, u00, x00,nt00) 

Leapfrog with :math:`\sigma = 0.5`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* More oscillatory than Lax-Wendroff, more accurate than Lax-Friedrichs


.. plot::

	def initial_and_boundary_conditions(nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 1.0
	    u[nx-1,:] = 0.0

	    half = int((nx-1)/2)

	    # Initial conditions      
	    u[0:half,0] = 1.0
	    
	    return u

	def analytical(sigma, nx, tmax, xmax, c):
	    
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u = initial_and_boundary_conditions(nx, nt)

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def leapfrog_convection(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u = initial_and_boundary_conditions(nx, nt)
	    
	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    # Initialise using Upwind
	    n = 0
	    i = 1
	    u[i:nx-1, n+1] = u[i:nx-1, n]-sigma*(u[i:nx-1, n]-u[i-1:nx-2, n])
	    
	    # Proceed using Leapfrog
	    for n in range(1, nt-1):
	        u[i:nx-1, n+1] = u[i:nx-1, n-1]-sigma*(u[i+1:nx, n]-u[i-1:nx-2, n])
	    
	    return u, x, nt

	def plot(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,4,NT)))   
	      for n in range(0,NT,4):
	         c=next(colour)
	         ax.plot(x,u[:,n],linestyle='-',c=c,label='n='+str(n))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-0.5,1.5])
	      plt.xlim([0.0,5.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show()

	u00, x00, nt00 = analytical(0.5, 101, 1.0, 5.0, 1.0)

	u0, x0, nt0 = leapfrog_convection(0.5, 101, 1.0, 5.0, 1.0)   

	plot(u0,x0,nt0, u00, x00,nt00) 

Results of Test 1: Heaviside Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math:: \sigma = 0.5

.. math:: \Delta x = 0.05

.. math:: n_{steps} = 40

Why does Lax-Friedrichs show step changes in the output?
--------------------------------------------------------

This is a double solution effect:

* :math:`u_i^{n+1}` does  not depend on :math:`u_i^n`
* Shifting the stencil by :math:`i` shows that :math:`u_i^{n+1}` and :math:`u_{i+1}^{n+1}` do not share a single mesh point of their stencils
* This is called "odd-even decoupling" 

Odd-Even Decoupling:
* Solutions on the odd points and even points have different error levels and can't communicate information
* One solution sightly ahead/slightly behind

.. figure:: ../_images/lax_friedrichs_2.png
   :align: center
   :scale: 100%

Why does Lax-Wendroff show good comparison with the analytical solution?
------------------------------------------------------------------------

* Lax-Wendroff is second order, so has reduced numerical diffusion
* However, numerical oscillations occur. More oscillations occur with Leapfrog than Lax-Wendroff


Sinusoidal Input, :math:`k = 4 \pi`
===================================

Travelling sinusoidal wave, 2 periods in a distance of 1m. Corresponding wave number:

.. math:: k = {{2 \pi} \over \lambda}  = 4 \pi

This is the initial condition.


Upwind Differencing with :math:`\sigma = 0.5`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Amplitude is being diffused, effective numerical diffusion after a number of timesteps
* Damped by backward difference method

.. plot::

	def initial_and_boundary_conditions(xmax, nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 0.0
	    u[nx-1,:] = 0.0

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    quarter = int((nx-1)/4)

	    from math import pi as PI
	    i=0
	    # Initial conditions 
	    u[i:quarter,0] = 1.0*np.sin(4 * PI * (x[i:quarter]))
	    
	    return u, x

	def analytical(sigma, nx, tmax, xmax, c):
	    # Increments
	    # dt = tmax/(nt-1)
	    # dx = (c * dt) / sigma
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions(xmax, nx, nt)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def upwind_convection(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions(xmax, nx, nt)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i:nx-1, n]-sigma*(u[i:nx-1, n]-u[i-1:nx-2, n])
	    return u, x, nt

	def plot(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,100,NT)))
	      ax.plot(x,u[:,0],linestyle='-',c='b',label='n='+str(0))
	      ax.plot(x,u[:,NT-1],linestyle='-',c='r',label='n='+str(NT-1))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-1.5,1.5])
	      plt.xlim([0.0,3.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show() 

	u00, x00, nt00 = analytical(0.5, 501, 1.0, 5.0, 1.0)

	u0, x0, nt0 = upwind_convection(0.5, 501, 1.0, 5.0, 1.0)

	plot(u0,x0,nt0, u00, x00,nt00)

Lax-Friedrichs with :math:`\sigma = 0.5`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Substantial numerical diffusion over time

.. plot::

	def initial_and_boundary_conditions(xmax, nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 0.0
	    u[nx-1,:] = 0.0

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    quarter = int((nx-1)/4)

	    from math import pi as PI
	    i=0
	    # Initial conditions 
	    u[i:quarter,0] = 1.0*np.sin(4 * PI * (x[i:quarter]))
	    
	    return u, x

	def analytical(sigma, nx, tmax, xmax, c):
	    # Increments
	    # dt = tmax/(nt-1)
	    # dx = (c * dt) / sigma
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions(xmax, nx, nt)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def lax_friedrichs_convection(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions(xmax, nx, nt)

	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = 0.5*(u[i-1:nx-2, n]+u[i+1:nx, n])-0.5*sigma*(u[i+1:nx, n]-u[i-1:nx-2, n])

	    return u, x, nt

	def plot(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,100,NT)))
	      ax.plot(x,u[:,0],linestyle='-',c='b',label='n='+str(0))
	      ax.plot(x,u[:,NT-1],linestyle='-',c='r',label='n='+str(NT-1))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-1.5,1.5])
	      plt.xlim([0.0,3.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show() 

	u00, x00, nt00 = analytical(0.5, 501, 1.0, 5.0, 1.0)

	u0, x0, nt0 = lax_friedrichs_convection(0.5, 501, 1.0, 5.0, 1.0)

	plot(u0,x0,nt0, u00, x00,nt00)

Lax-Wendroff with :math:`\sigma = 0.5`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Better representation of the wave
* Wiggles at the back of the wave, where there is a non-smooth slope

.. plot::

	def initial_and_boundary_conditions(xmax, nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 0.0
	    u[nx-1,:] = 0.0

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    quarter = int((nx-1)/4)

	    from math import pi as PI
	    i=0
	    # Initial conditions 
	    u[i:quarter,0] = 1.0*np.sin(4 * PI * (x[i:quarter]))
	    
	    return u, x

	def analytical(sigma, nx, tmax, xmax, c):
	    # Increments
	    # dt = tmax/(nt-1)
	    # dx = (c * dt) / sigma
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions(xmax, nx, nt)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def lax_wendroff_convection(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions(xmax, nx, nt)

	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = (u[i:nx-1, n] - 0.5*sigma*(u[i+1:nx, n] - u[i-1:nx-2, n]) +
	                         0.5*(sigma**2)*(u[i+1:nx, n] - 2.0*u[i:nx-1, n] + u[i-1:nx-2, n]))

	    return u, x, nt

	def plot(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,100,NT)))
	      ax.plot(x,u[:,0],linestyle='-',c='b',label='n='+str(0))
	      ax.plot(x,u[:,NT-1],linestyle='-',c='r',label='n='+str(NT-1))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-1.5,1.5])
	      plt.xlim([0.0,3.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show() 

	u00, x00, nt00 = analytical(0.5, 501, 1.0, 5.0, 1.0)

	u0, x0, nt0 = lax_wendroff_convection(0.5, 501, 1.0, 5.0, 1.0)

	plot(u0,x0,nt0, u00, x00,nt00)

Leapfrog with :math:`\sigma = 0.5`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Train of oscillations at the back

.. plot::

	def initial_and_boundary_conditions(xmax, nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 0.0
	    u[nx-1,:] = 0.0

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    quarter = int((nx-1)/4)

	    from math import pi as PI
	    i=0
	    # Initial conditions 
	    u[i:quarter,0] = 1.0*np.sin(4 * PI * (x[i:quarter]))
	    
	    return u, x

	def analytical(sigma, nx, tmax, xmax, c):
	    # Increments
	    # dt = tmax/(nt-1)
	    # dx = (c * dt) / sigma
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions(xmax, nx, nt)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def leapfrog_convection(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions(xmax, nx, nt)
	    
	    # Initialise using Upwind
	    n = 0
	    i = 1
	    u[i:nx-1, n+1] = u[i:nx-1, n]-sigma*(u[i:nx-1, n]-u[i-1:nx-2, n])
	    
	    # Proceed using Leapfrog
	    for n in range(1, nt-1):
	        u[i:nx-1, n+1] = u[i:nx-1, n-1]-sigma*(u[i+1:nx, n]-u[i-1:nx-2, n])
	    
	    return u, x, nt

	def plot(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,100,NT)))
	      ax.plot(x,u[:,0],linestyle='-',c='b',label='n='+str(0))
	      ax.plot(x,u[:,NT-1],linestyle='-',c='r',label='n='+str(NT-1))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-1.5,1.5])
	      plt.xlim([0.0,3.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show() 

	u00, x00, nt00 = analytical(0.5, 501, 1.0, 5.0, 1.0)

	u0, x0, nt0 = leapfrog_convection(0.5, 501, 1.0, 5.0, 1.0)

	plot(u0,x0,nt0, u00, x00,nt00)

Results of Test 2: Sinusoidal Input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math:: \sigma = 0.5

.. math:: \Delta x = 0.01

.. math:: n_{steps} = 200

Implications
------------

* Must avoid 1st order methods for the time propagation of a wave, i.e. Upwind and Lax-Friedrichs

* Lax-Wendroff and Leapfrog much better, discontinuity causes oscillations in wave


Sinusoidal Input, :math:`k = 10 \pi`
====================================

Travelling sinusoidal wave. Corresponding wave number:

.. math:: k =  10 \pi

This is the initial condition.


Upwind Differencing with :math:`\sigma = 0.5`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::

	def initial_and_boundary_conditions_2(xmax, nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 0.0
	    u[nx-1,:] = 0.0

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    quarter = int((nx-1)/4)
	    from math import pi as PI
	    i=0
	    # Initial conditions 
	    u[i:quarter,0] = 1.0*np.sin(12.0 * PI * (x[i:quarter]))
	    
	    return u, x

	def analytical_2(sigma, nx, tmax, xmax, c):
	    
	    # Increments
	    # dt = tmax/(nt-1)
	    # dx = (c * dt) / sigma
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions_2(xmax, nx, nt)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def upwind_convection_2(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions_2(xmax, nx, nt)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i:nx-1, n]-sigma*(u[i:nx-1, n]-u[i-1:nx-2, n])
	    return u, x, nt

	def plot_2(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,100,NT)))
	      ax.plot(x,u[:,0],linestyle='-',c='b',label='n='+str(0))
	      ax.plot(x,u[:,NT-1],linestyle='-',c='r',label='n='+str(NT-1))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-1.5,1.5])
	      plt.xlim([0.0,3.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show() 

	u0000, x0000, nt0000 = analytical_2(1.0, 501, 1.0, 5.0, 1.0)
	u000, x000, nt000 = upwind_convection_2(0.5, 501, 1.0, 5.0, 1.0)
	plot_2(u000, x000, nt000,u0000, x0000, nt0000)


Lax-Friedrichs with :math:`\sigma = 0.5`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::

	def initial_and_boundary_conditions_2(xmax, nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 0.0
	    u[nx-1,:] = 0.0

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    quarter = int((nx-1)/4)
	    from math import pi as PI
	    i=0
	    # Initial conditions 
	    u[i:quarter,0] = 1.0*np.sin(12.0 * PI * (x[i:quarter]))
	    
	    return u, x

	def analytical_2(sigma, nx, tmax, xmax, c):
	    
	    # Increments
	    # dt = tmax/(nt-1)
	    # dx = (c * dt) / sigma
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions_2(xmax, nx, nt)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def lax_friedrichs_convection_2(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions_2(xmax, nx, nt)

	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = 0.5*(u[i-1:nx-2, n]+u[i+1:nx, n])-0.5*sigma*(u[i+1:nx, n]-u[i-1:nx-2, n])

	    return u, x, nt

	def plot_2(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,100,NT)))
	      ax.plot(x,u[:,0],linestyle='-',c='b',label='n='+str(0))
	      ax.plot(x,u[:,NT-1],linestyle='-',c='r',label='n='+str(NT-1))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-1.5,1.5])
	      plt.xlim([0.0,3.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show() 

	u0000, x0000, nt0000 = analytical_2(1.0, 501, 1.0, 5.0, 1.0)
	u000, x000, nt000 = lax_friedrichs_convection_2(0.5, 501, 1.0, 5.0, 1.0)
	plot_2(u000, x000, nt000,u0000, x0000, nt0000)

Lax-Wendroff with :math:`\sigma = 0.5`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::

	def initial_and_boundary_conditions_2(xmax, nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 0.0
	    u[nx-1,:] = 0.0

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    quarter = int((nx-1)/4)
	    from math import pi as PI
	    i=0
	    # Initial conditions 
	    u[i:quarter,0] = 1.0*np.sin(12.0 * PI * (x[i:quarter]))
	    
	    return u, x

	def analytical_2(sigma, nx, tmax, xmax, c):
	    
	    # Increments
	    # dt = tmax/(nt-1)
	    # dx = (c * dt) / sigma
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions_2(xmax, nx, nt)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def lax_wendroff_convection_2(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions_2(xmax, nx, nt)

	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = (u[i:nx-1, n] - 0.5*sigma*(u[i+1:nx, n] - u[i-1:nx-2, n]) +
	                         0.5*(sigma**2)*(u[i+1:nx, n] - 2.0*u[i:nx-1, n] + u[i-1:nx-2, n]))

	    return u, x, nt

	def plot_2(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,100,NT)))
	      ax.plot(x,u[:,0],linestyle='-',c='b',label='n='+str(0))
	      ax.plot(x,u[:,NT-1],linestyle='-',c='r',label='n='+str(NT-1))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-1.5,1.5])
	      plt.xlim([0.0,3.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show() 

	u0000, x0000, nt0000 = analytical_2(1.0, 501, 1.0, 5.0, 1.0)
	u000, x000, nt000 = lax_wendroff_convection_2(0.5, 501, 1.0, 5.0, 1.0)
	plot_2(u000, x000, nt000,u0000, x0000, nt0000)

Leapfrog with :math:`\sigma = 0.5`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::

	def initial_and_boundary_conditions_2(xmax, nx, nt):
	    
	    # Initialise data structures
	    import numpy as np
	    u = np.zeros((nx,nt))
	    
	    # Boundary conditions
	    u[0,:] = 0.0
	    u[nx-1,:] = 0.0

	    # X Loop
	    x = np.zeros(nx)
	    x = np.linspace(0.0,xmax,nx)
	    
	    quarter = int((nx-1)/4)
	    from math import pi as PI
	    i=0
	    # Initial conditions 
	    u[i:quarter,0] = 1.0*np.sin(12.0 * PI * (x[i:quarter]))
	    
	    return u, x

	def analytical_2(sigma, nx, tmax, xmax, c):
	    
	    # Increments
	    # dt = tmax/(nt-1)
	    # dx = (c * dt) / sigma
	    dx = xmax/(nx-1)
	    dt =  ((dx*sigma) / c)
	    nt = int(sigma*tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions_2(xmax, nx, nt)
	    
	    i = 1
	    # Loop - must loop as NumPy cannot be used for dependency reasons
	    for n in range(0,nt-1):
	        u[i:nx-1, n+1] = u[i-1:nx-2, n]
	        
	    return u, x, nt

	def leapfrog_convection_2(sigma, nx, tmax, xmax, c):
	    """
	    Returns the velocity field and distance for 1D linear convection
	    """
	    # Increments
	    dx = xmax/(nx-1)
	    dt =  ((dx * sigma) / c)
	    nt = int(tmax/dt + 1)

	    # Initial and Boundary Conditions
	    u, x = initial_and_boundary_conditions_2(xmax, nx, nt)
	    
	    # Initialise using Upwind
	    n = 0
	    i = 1
	    u[i:nx-1, n+1] = u[i:nx-1, n]-sigma*(u[i:nx-1, n]-u[i-1:nx-2, n])
	    
	    # Proceed using Leapfrog
	    for n in range(1, nt-1):
	        u[i:nx-1, n+1] = u[i:nx-1, n-1]-sigma*(u[i+1:nx, n]-u[i-1:nx-2, n])
	    
	    return u, x, nt

	def plot_2(u,x,NT,u_analytical, x_analytical, NT2):
	      """
	      Plots the 1D velocity field
	      """

	      import matplotlib.pyplot as plt
	      import matplotlib.cm as cm
	      fig=plt.figure()
	      ax=plt.subplot(111)
	      colour=iter(cm.rainbow(np.linspace(0,100,NT)))
	      ax.plot(x,u[:,0],linestyle='-',c='b',label='n='+str(0))
	      ax.plot(x,u[:,NT-1],linestyle='-',c='r',label='n='+str(NT-1))
	      ax.plot(x_analytical,u_analytical[:,NT2-1],linestyle='--',c='k',label='n='+str(NT2-1)+' analytical')   
	      box=ax.get_position()
	      fig.set_figwidth(12.0)
	      fig.set_figheight(8.0)
	      ax.set_position([box.x0, box.y0, box.width*0.7,box.height])
	      ax.legend( bbox_to_anchor=(1.02,1), loc=2)
	      plt.ylim([-1.5,1.5])
	      plt.xlim([0.0,3.0])
	      plt.xlabel('x (m)')
	      plt.ylabel('u (m/s)')
	      plt.show() 

	u0000, x0000, nt0000 = analytical_2(1.0, 501, 1.0, 5.0, 1.0)
	u000, x000, nt000 = leapfrog_convection_2(0.5, 501, 1.0, 5.0, 1.0)
	plot_2(u000, x000, nt000,u0000, x0000, nt0000)



Results of Test 3: High Frequency Sinusoidal Input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math:: \sigma = 0.5

.. math:: \Delta x = 0.01

.. math:: n_{steps} = 200

Implications
------------

* Upwind and Lax Friedrichs are catastrophically dissipative
* Lax Wendroff shows some dissipation and a lag
* Leapfrog shows less dissipation, but more oscillations and still a lag

Summary
=======

1st order schemes:

* Have poor accuracy - gets even worse for solutions with higher frequency, damping is catastrophic

2nd order schemes:

* Provide better accuracy
* Generate numerical oscillations - associated with locations where solution is not smooth. Oscillations are stronger with Leapfrog scheme
* Numerical errors are sensitive to the frequency content of the solution, i.e frequency content of initial condition

We have obtained results from simple 1D models, however they are representative of real flow situations (2D, 3D etc)

