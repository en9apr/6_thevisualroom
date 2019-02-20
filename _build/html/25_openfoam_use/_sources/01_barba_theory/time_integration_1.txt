===================
The Method of Lines
===================

.. contents::
   :local:

e.g. Convection :math:`\Rightarrow` combinations of space and time schemes show very different behaviour

Using CD in space for the convection term :math:`c \partial u / \partial x`:

 * FD in time :math:`\Rightarrow` unstable
 * BD in time :math:`\Rightarrow` unconditionally stable
 * CD in time :math:`\Rightarrow` conditionally stable

**Question: What criteria should we impose on a time-discretisation method?** Applied to given space-discretisation, so that we obtain a stable and accurate scheme.

General Method
==============

1) Perform a space discretisation :math:`\Rightarrow` system of ODEs
   And express space discretisation in **matrix** form (including BCs)
2) Perform a **spectral** analysis of the matrix :math:`\Rightarrow` eigenvalues
3) Stability conditions on the space discretisation
   
   * The exact time integrated solution of the system of ODEs should not grow unbounded in time
   * The eigenvalues :math:`\Omega` of the space discretisation matrix must have a **non-positive real part**

4) Select a time integration method for the **semi-discretised** system (discretised in space not time)

   :math:`\Rightarrow` stability analysis as a function of eigenvalues :math:`\Omega`

5) **Compatibility** :math:`\Rightarrow` the **stability** region of the time discretisation must include the whole spectrum :math:`\Omega`

The Method of Lines: Analysis of the Space-Discretised Systems
==============================================================

Consider a general conservation law as I.BVP over domain :math:`\Omega` with boundary :math:`\Gamma`:

.. math:: {\partial u \over \partial t} + \nabla \cdot \mathbf{F} = 0

* I.C. :math:`u(\mathbf{x},0) = u^0(\mathbf{x})`, :math:`\mathbf{x} \in \Omega`
* B.C. :math:`u(\mathbf{x},t) = g(\mathbf{x},t)`, :math:`\mathbf{x} \in \Gamma`

In 1D: discretise flux differential operator:

System of ODES:

.. math:: {d u_i \over dt} = -{1 \over \Delta x_i}(f_{i+1/2}^* - f_{i-1/2}^*)

Group all mesh points to form a column vector :math:`\mathbf{u}` and write a system of ODEs

.. math:: {d \mathbf{u} \over dt} = \mathbf{S} \cdot \mathbf{u} + \mathbf{Q}

Where :math:`\mathbf{S}` is a matrix representing space discritisation

And :math:`\mathbf{Q}` is the contribution from BCs

In the FVM formulation, the RHS is the "residual"

This is the "Method of Lines"

Example 1 - Diffusion
=====================

.. math:: u_t = \alpha u_{xx}

CD in space:

.. math:: {d u_i \over dt} = {\alpha \over \Delta x_i}(u_{i+1} - 2u_{i} + u_{i-1})

To get complete matrix :math:`S`, we need BCs

Dirichlet BCs
-------------

Domain :math:`(0, L)` with :math:`N` mesh points (0 and N-1 are the BCs)

.. math:: x_i = i \Delta x

.. math:: \Delta x = {L \over {N-1}}

B.C. :math:`u(0,t) = a` and :math:`u(L,t) = b`

i = 1 (a occurs at i=0)

.. math:: {d u_1 \over dt} = {\alpha \over \Delta x_i}(u_{2} - 2u_{1} + a) 

i = N-2 (b occurs at i=N-1)

.. math:: {d u_{N-2} \over dt} = {\alpha \over \Delta x_i}(b - 2u_{N-2} + u_{N-3})

Finally:

.. math::
   {du \over dt} =
   {\alpha \over {\Delta x^2}}
   \begin{bmatrix}
    -2 & 1 & & &  0 \\
    1 & -2 & 1 & & \\
    \vdots   & \ddots & \ddots & \ddots & \vdots \\
    & &  1 & -2 &  1 \\
    0 & & & 1 & -2
   \end{bmatrix}
   \begin{bmatrix}
   u_1 \\
   u_2 \\
   \vdots \\
   u_{N-3} \\
   u_{N-2}
   \end{bmatrix} +
   \begin{bmatrix}
   a \alpha / \Delta x^2 \\
   0 \\
   \vdots \\
   0 \\
   b \alpha / \Delta x^2 \\
   \end{bmatrix} =
   \mathbf{S} \cdot \mathbf{u} + \mathbf{Q}
 
Neumann BCs
-----------

:math:`{\partial u \over \partial x} = a` at :math:`x = 0` and :math:`u_N = u(L,t) = b`

The Neumann BC needs to be discretised using a one-sided difference

i = 0:

.. math:: {{u_1 - u_0} \over {\Delta x}} = a \Rightarrow u_0 = u_1 - a \Delta x 

i = 1:

.. math:: {d u_1 \over dt} = {\alpha \over \Delta x_i}(u_2 - 2u_1 + u_0) = 
                             {\alpha \over \Delta x_i}(u_2 - 2u_1 + u_1 - a \Delta x) =
                             {\alpha \over \Delta x_i}(u_2 - u_1) - {{a \alpha} \over {\Delta x}}

Finally (we don't know :math:`u_0` anymore, but we have substituted for it at i=1)

.. math::
   {du \over dt} =
   {\alpha \over {\Delta x^2}}
   \begin{bmatrix}
    -1 & 1 & & &  0 \\
    1 & -2 & 1 & & \\
    \vdots   & \ddots & \ddots & \ddots & \vdots \\
    & &  1 & -2 &  1 \\
    0 & & & 1 & -2
   \end{bmatrix}
   \begin{bmatrix}
   u_1 \\
   u_2 \\
   \vdots \\
   u_{N-3} \\
   u_{N-2}
   \end{bmatrix} +
   \begin{bmatrix}
   -a \alpha / \Delta x \\
   0 \\
   \vdots \\
   0 \\
   b \alpha / \Delta x^2 \\
   \end{bmatrix} =
   \mathbf{S} \cdot \mathbf{u} + \mathbf{Q}                          

Importance of BCs for stability
-------------------------------

* The form of :math:`\mathbf{S}` depends on BCs
* The eigenvalues depend on :math:`\mathbf{S}`
* The stability depends on the eigenvalues

**Hence, the stability depends on the BCs**

Periodic BCs
------------

**We must solve all the way to the boundaries, and set any values outside the domain to be equal to the values inside the domain at the opposite end of the domain**

:math:`u_N = u_0` and :math:`u_{-1} = u_{N-1}`

In this case, write an equation for i = 0:

.. math:: {d u_0 \over dt} = {\alpha \over \Delta x^2}(u_1 - 2u_0 + u_{N-1})

For i = N-1:

.. math:: {d u_{N-1} \over dt} = {\alpha \over \Delta x^2}(u_0 - 2u_{N-1} + u_{N-2})
                
.. math::
   {du \over dt} =
   {\alpha \over {\Delta x^2}}
   \begin{bmatrix}
    -2 & 1 & & &  1 \\
    1 & -2 & 1 & & \\
    \vdots   & \ddots & \ddots & \ddots & \vdots \\
    & &  1 & -2 &  1 \\
    1 & & & 1 & -2
   \end{bmatrix}
   \begin{bmatrix}
   u_0 \\
   u_1 \\
   \vdots \\
   u_{N-2} \\
   u_{N-1}
   \end{bmatrix}  =
   \mathbf{S} \cdot \mathbf{u} + \mathbf{Q}     

:math:`\mathbf{S}` similar to the case with Dirichlet BCs, but with 1 in upper right and lower left corners.

Typical for periodic matrices

Example 2 - Convection
======================

.. math:: u_t + a u_x = 0

with :math:`a \gt 0`

Upwind
------

Upwind for interior points:

.. math:: {{d u_i} \over dt} = -{a \over {\Delta x}}(u_i - u_{i-1})

B.C. :math:`u(0,t) = g(t)` at the left boundary

i = 1 (first internal point)

.. math:: {d u_1 \over dt} = -{a \over \Delta x}(u_1 - u_{0})
                           = -{a \over \Delta x}u_1 + {a \over \Delta x} g(t)

i = N-1 (at the right boundary)

.. math:: {d u_{N-1} \over dt} = -{a \over \Delta x}(u_{N-1} - u_{N-2})

Note: no additional boundary condition is needed at the right boundary (because we are using backward differencing, or upwind)

Finally:

.. math::
   {{d \mathbf{u}} \over dt} =
   {a \over {\Delta x}}
   \begin{bmatrix}
    1 & 0 & & &  0 \\
    -1 & 1 & 0 & & \\
    \vdots   & \ddots & \ddots & \ddots & \vdots \\
    & &  -1 & 1 &  0 \\
    0 & & & -1 & 1
   \end{bmatrix}
   \begin{bmatrix}
   u_1 \\
   u_2 \\
   \vdots \\
   u_{N-2} \\
   u_{N-1}
   \end{bmatrix} +
   \begin{bmatrix}
   a g(t) / \Delta x \\
   0 \\
   \vdots \\
   0 \\
   0 \\
   \end{bmatrix} =
   \mathbf{S} \cdot \mathbf{u} + \mathbf{Q}  

**Question: What would happen if we imposed a boundary condition at:** :math:`u_{N-1} = g(t)`

The last equation would be:

.. math:: {d u_{N-2} \over dt} = -{a \over \Delta x}(u_{N-2} - u_{N-3})

This is decoupled from the data at :math:`i = N-1`

No way to ensure that a numerical condition at :math:`u_0` will lead to a downstream value that satisfies :math:`u_{N-1} = g(t)`. Problem is not well-posed

.. figure:: ../_images/downstream.png
   :align: center
   :scale: 70%

**For a > 0**

The value at A, :math:`u_0` must be specified at i=0. The corresponding value at B will be part of the solution, following the characteristic :math:`dx/dt = a`

Any other BC at x=L would be incompatible with :math:`u_0` at A.

Central space discretisation with periodic BCs
----------------------------------------------

Assume:

:math:`u_N = u_0` and :math:`u_{-1} = u_{N-1}`

i = 0:

.. math:: {d u_0 \over dt} = -{a \over {2 \Delta x}}(u_1 - u_{-1}) = -{a \over {2 \Delta x}}(u_1 - u_{N-1})

i = N-1:

.. math:: {d u_{N-1} \over dt} = -{a \over {2 \Delta x}}(u_N - u_{N-2}) = -{a \over {2 \Delta x}}(u_0 - u_{N-2})

Finally:

.. math::
   {{d \mathbf{u}} \over dt} =
   -{a \over {2 \Delta x}}
   \begin{bmatrix}
    0 & 1 & & &  -1 \\
    -1 & 0 & 1 & & \\
    \vdots   & \ddots & \ddots & \ddots & \vdots \\
    & &  -1 & 0 &  1 \\
    1 & & & -1 & 0
   \end{bmatrix}
   \begin{bmatrix}
   u_0 \\
   u_1 \\
   \vdots \\
   u_{N-2} \\
   u_{N-1}
   \end{bmatrix}
