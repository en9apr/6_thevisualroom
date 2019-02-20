==========================
The Marker and Cell Method
==========================

.. contents::
   :local:

Collocated vs. Staggered Grids
------------------------------

* Illustration: consider a 1D Navier-Stokes System

.. math:: {\partial u \over \partial x} = 0
   :label: 1

.. math:: {\partial u \over \partial t} + {\partial u^2 \over \partial x} = 
          - {1 \over \rho} {\partial p \over \partial x} + \nu {\partial^2 u \over \partial x^2}
   :label: 2

* On a standard **collocated grid** we have the following with **FTCS**:

.. math:: {{u_{i+1}^{n+1} - u_{i-1}^{n+1}} \over {2 \Delta x}} = 0
   :label: 3

.. math:: {{u_{i}^{n+1} - u_{i}^{n}} \over {\Delta t}} +
          {({u_{i+1}^{n})^2 - (u_{i-1}^{n})^2} \over {2 \Delta x}} = -
          {1 \over \rho}{ {p_{i+1}^n - p_{i-1}^n} \over {2 \Delta x}} +
          \nu{{u_{i+1}^{n} -2u_i^n + u_{i-1}^{n}} \over {\Delta x^2}}
   :label: 4

.. figure:: ../_images/collocated.png
   :scale: 50%
   :align: center

Apply Pressure Correction Approach
----------------------------------

* Split the N-S Equations into two parts
* **Step One:** Ignore pressure gradient in momentum equation for intermediate velocity (:math:`u^*` is not divergence free). From Equation :eq:`4`:

.. math:: {{u_{i}^{*} - u_{i}^{n}} \over {\Delta t}} +
          {({u_{i+1}^{n})^2 - (u_{i-1}^{n})^2} \over {2 \Delta x}} = -
          \nu{{u_{i+1}^{n} -2u_i^n + u_{i-1}^{n}} \over {\Delta x^2}}
   :label: 5

* **Step Two:** Use intermediate velocity and correct it with the pressure (pressure enforces continuity) to obtain :math:`u^{n+1}`. From Equation :eq:`4`:

.. math:: {{u_{i}^{n+1} - u_{i}^{*}} \over {\Delta t}}  = -
          {1 \over \rho}{ {p_{i+1}^n - p_{i-1}^n} \over {2 \Delta x}}
 
      
.. math:: u_{i}^{n+1} =  u_{i}^{*} - {\Delta t \over \rho}{ {p_{i+1}^n - p_{i-1}^n} \over {2 \Delta x}}
   :label: 6

* **Step Three:** Use pressure correction Equation :eq:`6` to satsfy continuity, i.e. Equation :eq:`3` (just replace i in LHS of :eq:`6` with i+1 or i-1 and correct RHS)

.. math:: {1 \over {2 \Delta x}} \left[ \left({u_{i+1}^{*} - {\Delta t \over \rho}{ {p_{i+2}^n - p_{i}^n} \over {2 \Delta x}}} \right) -
          \left( {u_{i-1}^{*} - {\Delta t \over \rho}{ {p_{i}^n - p_{i-2}^n} \over {2 \Delta x}}} \right) \right] = 0

Re-arranging this gives a form of the Poisson Equation, that ensures continuity:

.. math:: {{p_{i+2}^n - 2p_{i}^n + p_{i-2}^n} \over {4 \Delta x^2}} = 
          {{u_{i+1}^{*} - u_{i-1}^{*}} \over {2 \Delta x}} {\rho \over {\Delta t} }
   :label: 7

This is similar to using FTCS on the divergence of the momentum equation and setting the divergence of velocity to zero, as shown in `the Poisson Equation for pressure <http://www.thevisualroom.com/poisson_for_pressure.html>`_. The form here includes only the first term on the RHS of that complex expression.

Odd-Even Decoupling
-------------------

Note:

* The discretisation of pressure is on a :math:`2 \Delta x` grid.
* The discretisation of velocity is on a :math:`\Delta x` grid.

In the final equation for :math:`p_i^n` if the pressure used is at an **even** point on the grid, the velocity is at an **odd** point, this is called **odd-even decoupling**

.. math:: p_{i}^n = {{p_{i+2}^n + p_{i-2}^n} \over {2}} - {{\Delta x \rho} \over {\Delta t}} ({u_{i+1}^{*} - u_{i-1}^{*}})

**Possible drawback of odd-even decoupling:**

The pressure at point :math:`i` is not influence by the velocity component :math:`u_i^n` and viceversa :math:`\Rightarrow` can result in **high-frequency oscillations**.

* Stencil for :math:`p`: :math:`\quad i-2, \qquad i, \qquad i+2`
* Stencil for :math:`u^*`: :math:`\ \ \qquad i-1, \quad i+1`

Demonstration of Odd-Even Decoupling
------------------------------------

For the `cavity flow <http://nbviewer.ipython.org/github/en9apr/sphinx/blob/master/Navier_Stokes_Cavity_Slices.ipynb>`_, using:

* 41 x 41 mesh
* dx = dy = 0.05
* :math:`\nu` = 0.01
* Re = 200

**Pressure is not monotonic due to odd-even decoupling**

The stencil with the previous pressure correction method was:

* Stencil for :math:`p`: :math:`\quad i-1, \qquad i, \  \qquad i+1`
* Stencil for :math:`u`: :math:`\quad i-1, \qquad \qquad \quad i+1`

So the pressure at :math:`p_i` is **not influenced** by :math:`u_i`

.. plot::

    def navier_stokes_initialisation(niter, r, nx_or_ny, tmax, xmax_or_ymax):
        """
        Returns the velocity field and distance for 2D linear convection
        """
        # Increments:
        nx = ny = nx_or_ny
        xmax = ymax = xmax_or_ymax
        dx = xmax/(nx-1)
        dy = ymax/(ny-1)
        nt = int((tmax / (r*(dx)**2))+1)
        dt = tmax/(nt-1)
        
        # Initialise data structures:
        import numpy as np
        p = np.zeros((nx,ny))
        u = np.zeros((nx,ny))
        v = np.zeros((nx,ny))
        
        # linspace is SIMPLER than list comprehensions:
        x = np.linspace(0.0,2.0,nx)
        y = np.linspace(0.0,2.0,ny)
        
        # Pressure Boundary Conditions:
        p[:, ny-1] = 0.0
        
        # Velocity Boundary Conditions:
        u[:,ny-1] = 1.0
            
        return p, x, y, u, v, nx, ny, nt, dx, dy, dt, niter, r

    def navier_stokes(rho, nu, niter, r, nx, tmax, xmax):
                      
        (p, x, y, u, v, nx, ny, nt, 
        dx, dy, dt, niter, r) = navier_stokes_initialisation(niter, r, nx, tmax, xmax)
        
        # Increments
        h = dx
        
        import numpy as np
        
        # Intermediate copies:
        un = np.zeros((nx, ny))
        vn = np.zeros((nx, ny))
        pm = np.zeros((nx, ny))
        bm = np.zeros((nx, ny)) # bm needs to be exactly zero at the boundaries
            
        # Loop - use decimal points for all floating point numbers
        for n in range(nt):    
            
            # We know the velocity at i=0, j=0, i=nx-1 and j=ny-1. b is zero at the boundaries.  
            bm[1:-1, 1:-1] = ( (rho / (2.0 * h * dt)) * ( u[2:, 1:-1] - u[0:-2, 1:-1] 
            + v[1:-1, 2:] - v[1:-1, 0:-2] ) +
            (rho / (4.0*h**2)) * ( (u[2:, 1:-1] - u[0:-2, 1:-1])**2.0 + 
            4.0*h*(u[1:-1, 2:] - u[1:-1, 0:-2])*(v[2:, 1:-1] - v[0:-2, 1:-1]) + 
            (v[1:-1, 2:] - v[1:-1, 0:-2])**2.0 ) )
            
            # First points for p. We don't know the pressure at i=0, j=0 and i=nx-1. We DO know the pressure at j=ny-1
            for m in range(niter):

                pm = np.copy(p)
                p[1:-1, 1:-1] = 0.25*( pm[2:, 1:-1] + pm[0:-2, 1:-1] + pm[1:-1, 2:] + pm[1:-1, 0:-2]
                - bm[1:-1, 1:-1]*h**2.0 )
                
                # Set zero gradient boundary conditions:
                p[0, :] = p[1, :]
                p[:, 0] = p[:, 1]
                p[-1, :] = p[-2, :]
            
            # First points for u and v. We know the velocity at i=0, j=0, i=nx-1 and j=ny-1.
            # We are simply using the value of pressure here
            un = np.copy(u)
            vn = np.copy(v)
                  
            u[1:-1, 1:-1] = ( un[1:-1, 1:-1] - 
            (dt / h) * ( un[1:-1, 1:-1] * ( un[1:-1, 1:-1] - un[0:-2, 1:-1] ) + 
            vn[1:-1, 1:-1] * ( un[1:-1, 1:-1] - un[1:-1, 0:-2] ) ) - 
            (dt / (2.0 * rho * h)) * ( p[2:, 1:-1] - p[0:-2, 1:-1] ) + 
            (dt * nu / h**2.0) * ( un[0:-2, 1:-1] + un[2:, 1:-1] + un[1:-1, 0:-2] + un[1:-1, 2:] - 
            4.0 * un[1:-1, 1:-1] ) ) 

            v[1:-1, 1:-1] = ( vn[1:-1, 1:-1] - 
            (dt / h) * ( un[1:-1, 1:-1] * ( vn[1:-1, 1:-1] - vn[0:-2, 1:-1] ) + 
            vn[1:-1, 1:-1] * ( vn[1:-1, 1:-1] - vn[1:-1, 0:-2] ) ) - 
            (dt / (2.0 * rho * h)) * ( p[1:-1, 2:] - p[1:-1, 0:-2] ) + 
            (dt * nu / h**2.0) * ( vn[0:-2, 1:-1] + vn[2:, 1:-1] + vn[1:-1, 0:-2] + vn[1:-1, 2:] -
            4.0 * vn[1:-1, 1:-1] ) )          
        
        return u, v, p, x, y

    u33, v33, p33, x33, y33 = navier_stokes(1.0, 0.01, 50, 0.5, 41, 5.0, 2.0)

    def vector_contour_2(u, v, p, x, y):
        fig = plt.figure(figsize=(11,9), dpi=100)
        Y,X=np.meshgrid(y,x) #note meshgrid uses y,x not x,y!!!
        plt.pcolor(X,Y,p)
        plt.colorbar()
        #plt.contourf(X,Y,p[:,:],alpha=0.5)    ###plotting the pressure field as a contour
        # plt.colorbar()
        # plt.contour(X,Y,p[:,:])               ###plotting the pressure field outlines
        plt.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2]) ##plotting velocity
        #plt.contour(X,Y,p) ##plotting velocity
        plt.ylim([0,2.0])
        plt.xlim([0,2.0])
        plt.xlabel('X')
        plt.ylabel('Y')
        
    vector_contour_2(u33, v33, p33, x33, y33)


Solution to Odd-Even Decoupling
-------------------------------

Define velocity and pressure **on separate grids**

* Collocated grid (where odd-even decoupling took place):

.. figure:: ../_images/collocated_2.png
   :scale: 65%
   :align: center

* Staggered grid (to prevent odd-even decoupling):

.. figure:: ../_images/staggered.png
   :scale: 65%
   :align: center

This solution was due to Harlow and Welch (1965)

What effect does this have on the discretised equations?

* Continuity Equation written as:

.. math:: {{u_{i+{1 \over 2}}^{n+1} - u_{i-{1 \over 2}}^{n+1}} \over {\Delta x}} = 0
   :label: 8

* "Fractional step" or "pressure correction"

.. math:: {{u_{i+{1 \over 2}}^{n+1} - u_{i+{1 \over 2}}^{*}} \over {\Delta t}} = -
          {1 \over \rho}{ {p_{i+1}^n - p_{i}^n} \over {\Delta x}} +
   :label: 9

As before, substitute :eq:`9` into :eq:`8`

.. math:: {{p_{i+1}^n - 2p_{i}^n + p_{i-1}^n} \over {\Delta x^2}} = 
          {{u_{i+{1 \over 2}}^{*} - u_{i-{1 \over 2}}^{*}} \over {\Delta x}} {\rho \over {\Delta t} }
   :label: 10

**How do we obtain** :math:`{u_{i+{1 \over 2}}^{*}}` **etc ?**

Answer: **Averaging** We get values at :math:`i+{1 \over 2}` and :math:`i-{1 \over 2}` by averaging adjacent values

This makes the equations **fully coupled** and **eliminates odd-even decoupling**

Extension of Staggered Grid to 2D
---------------------------------

.. figure:: ../_images/2D_staggered.png
   :scale: 75%
   :align: center

This uses three independent grids:

* One for :math:`p_{i,j}` etc
* One for :math:`u_{{i \pm {1 \over 2}}, {j}}` etc
* One for :math:`v_{i,{j \pm {1 \over 2}}}` etc

Equations to be solved:

.. math:: {\partial u \over \partial x}+ {\partial v \over \partial y} = 0
   :label: 11

.. math:: {\partial u \over \partial t} + {\partial u^2 \over \partial x}+ {\partial {uv} \over \partial y} = 
          - {1 \over \rho} {\partial p \over \partial x} + \nu \left( {{\partial^2 u \over \partial x^2} +{ \partial^2 u \over \partial y^2} } \right) = - {\partial \psi \over \partial x} + \nu \left( {{\partial^2 u \over \partial x^2} +{ \partial^2 u \over \partial y^2} } \right)
   :label: 12

.. math:: {\partial v \over \partial t} + {\partial vu \over \partial x}+ {\partial {v^2} \over \partial y} = 
          - {\partial \psi \over \partial y} + \nu \left( {{\partial^2 v \over \partial x^2} +{ \partial^2 v \over \partial y^2} } \right)
   :label: 14

:math:`\psi \Rightarrow` pressure divided by density

Conversion from Conservative Form to Non-Conservative Form (to check equivalence) via product rule:

.. math:: {\partial {uu} \over \partial x} + {\partial {uv} \over \partial y} = 
          u {\partial u \over \partial x} + {\partial u \over \partial x}u +
          {\partial u \over \partial y}v + u{\partial {v} \over \partial y} =
          u \left( {\partial u \over \partial x} + {\partial v \over \partial y} \right)+
          u {\partial u \over \partial x} +  v {\partial u \over \partial y} =
          u {\partial u \over \partial x} +  v {\partial u \over \partial y}

Application of Finite Difference Method to Staggered Grid
---------------------------------------------------------

.. math:: \left . {\partial u \over \partial t} \right |_{i+{1 \over 2},j}^n =
                   {{u_{i+{1 \over 2},j}^{n+1} - u_{i+{1 \over 2},j}^{n}} \over {\Delta t}}


.. math:: \left . {\partial u^2 \over \partial x} \right |_{i+{1 \over 2},j}^n =
                  {{(u_{i,j}^n)^2 - (u_{i+1,j}^{n})^2} \over {\Delta x}}

.. math:: \left . {\partial {uv} \over \partial y} \right |_{i+{1 \over 2},j}^n =
                  { { (u_{i+{1 \over 2},j+{1 \over 2}}^n)(v_{i+{1 \over 2},j+{1 \over 2}}^n)-
                   (u_{i+{1 \over 2},j-{1 \over 2}}^n)(v_{i+{1 \over 2},j-{1 \over 2}}^n)}
                   \over {\Delta y} }

.. math:: \left . {\partial \psi \over \partial x} \right |_{i+{1 \over 2},j}^n =
                  {{(\psi_{i,j}^n) - (\psi_{i+1,j}^{n})} \over {\Delta x}}

.. math:: \left . {\partial^2 u \over \partial x^2} \right |_{i+{1 \over 2},j}^n =
                  {{u_{i+{3 \over 2},j}^n - 2u_{i+{1 \over 2},j}^n + u_{i-{1 \over 2},j}^{n}} \over {\Delta x^2}}

Where is the Velocity Known?
----------------------------

Velocities :math:`u,v` are known only at half mesh points - i.e. at **mid points of vertical and horizontal cell edges**.

* Known at :math:`u_{i \pm {1 \over 2},j}` and :math:`v_{i \pm {1 \over 2},j}` etc

* Unknown at :math:`u_{i,j}` or :math:`u_{i \pm {1 \over 2},j \pm {1 \over 2}}` etc 

You need to average from the half mesh points

How do we get the Velocity at the Corners and Centre of the Cells?
------------------------------------------------------------------

* The velocity is known at the midpoint of the cell edges
* The velocity is unknown at the corners and centres of cells (e.g. points a, b, c and d)

.. figure:: ../_images/2D_staggered_averages.png
   :scale: 65%
   :align: center

So for example, if we wanted the following derivative:

.. math:: \left . {\partial {uv} \over \partial y} \right |_{i+{1 \over 2},j}^n =
                  { { (u_{i+{1 \over 2},j+{1 \over 2}}^n)(v_{i+{1 \over 2},j+{1 \over 2}}^n)-
                   (u_{i+{1 \over 2},j-{1 \over 2}}^n)(v_{i+{1 \over 2},j-{1 \over 2}}^n)}
                   \over {\Delta y} }

At Point a:

.. math:: u_{i+{1 \over 2},j+{1 \over 2}}^n = 
          {1 \over 2}({ u_{i+{1 \over 2},j+{1}}^n + u_{i+{1 \over 2},j}^n })

.. math:: v_{i+{1 \over 2},j+{1 \over 2}}^n = 
          {1 \over 2}({ v_{i+{1},j+{1 \over 2}}^n + v_{i,j+{1 \over 2}}^n })

Similarly for Point b

If we wanted to know this derivative:

.. math:: \left . {\partial u^2 \over \partial x} \right |_{i+{1 \over 2},j}^n =
                  {{(u_{i,j}^n)^2 - (u_{i+1,j}^{n})^2} \over {\Delta x}}

At Point c:

.. math:: u_{i+1,j}^{n} = 
          {1 \over 2}({ u_{i+{3 \over 2},j}^n + u_{i+{1 \over 2},j}^n })

At Point d:

.. math:: u_{i,j}^n =
          {1 \over 2}({ u_{i+{1 \over 2},j}^n + u_{i-{1 \over 2},j}^n })


Historical Notes
----------------

Harlow and Welch (1965) introduced the "Marker and Cell" method:

* It included marker particles to follow free surfaces (now obsolete)
* Current usage of **"Marker and Cell"** method (MAC method) means - **"Projection Method using a Staggered Grid"**
