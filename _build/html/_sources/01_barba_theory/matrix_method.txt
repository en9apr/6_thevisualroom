=========================================================
The Matrix Method for Stability of Spatial Discretisation
=========================================================

.. contents::
   :local:

The Eigenvalue Spectrum of :math:`\mathbf{S}`
=============================================

von Neumann stability analysis assumes a periodic solution and a uniform mesh

We need to establish a more general method that includes the influence of the BCs 

Let :math:`\Omega_j` where :math:`j = 1 .... N` be eigenvalues of :math:`\mathbf{S}`:

.. math:: det \left| \mathbf{S} - \Omega I \right| = 0

Associated eigenvectors :math:`V^j`  :math:`\Rightarrow`  

.. math:: \mathbf{S} \cdot V^j(\mathbf{x}) = \Omega_j V^j

:math:`T`: the matrix formed by :math:`V^j` as columns:

.. math:: \mathbf{S} \cdot \mathbf{T} = \mathbf{T} \Omega

where:

.. math::
   \Omega = 
   \begin{bmatrix}
    \Omega_1 &  & & &   \\
    & \Omega_2 & & & \\
    & & \ddots & & \\
    & & & & \Omega_N
   \end{bmatrix}

.. math:: \Omega = \mathbf{T}^{-1} \mathbf{S} \mathbf{T}

Eigenvectors are a complete basis set 

Recall
======

.. math:: {d \mathbf{u} \over dt} = \mathbf{S} \cdot \mathbf{u} + \mathbf{Q}

Model Decomposition
===================

Exact solution :math:`\overline{\mathbf{u}}(t,\mathbf{x})`

.. math:: \overline{\mathbf{u}}(t, \mathbf{x}) = \sum_{j=1}^N \overline{\mathbf{u}_j}(t) \cdot \mathbf{V}^j(\mathbf{x})

First term depends only on time, the second only on space.

.. math:: Q = \sum_{j=1}^N Q_{j} V^j

Model Equation
==============

Leads to **"Modal Equation"** for the time dependent coefficient

.. math:: {\overline{\mathbf{u}_j} \over t} = \Omega_j \overline{\mathbf{u}_j} + Q_j

The exact solution is the sum of the homogeneous solution:

.. math:: {\overline{\mathbf{u}_{jT}} \over t} = \Omega_j \overline{\mathbf{u}_{jT}} \Rightarrow
          \overline{\mathbf{u}_{jT}}(t) = c_{0j} e^{\Omega_j t}

And a particular solution e.g. a solution of the steady state equation:

.. math:: \Omega_j \overline{\mathbf{u}_{jS}} + Q_j = 0 \Rightarrow
          \overline{\mathbf{u}_{jS}} = - {Q_j \over {\Omega_j}}

Modal Solution
==============

.. math:: \overline{\mathbf{u}_j}(t) = c_{0j} e^{\Omega_j t} - {Q_j \over {\Omega_j}}

Where :math:`c_{0j}` are from the I.C.

IC: 

.. math:: \mathbf{u}^0(\mathbf{x}) = u(t=0, \mathbf{x})

In matrix representation:

.. math:: \mathbf{u}^0 =  
   \begin{bmatrix}
    u_0^0 \\
    \vdots \\
    u_i^0 \\
    \vdots \\
    u_{N-1}
   \end{bmatrix} 

And:

.. math:: \mathbf{u}^0 = \sum_{j=1}^N \mathbf{u}_j^0 \mathbf{V}^j

At t=0 we get the coefficients

.. math:: c_{0j} = \mathbf{u}_j^0 + {Q_j \over \Omega_j}

Hence:

.. math:: \overline{\mathbf{u}_j}(t) = \mathbf{u}_j^0 e^{\Omega_j t} - 
                                       {Q_j \over {\Omega_j}} \left( e^{\Omega_j t} - 1  \right)    

**The eigenvalues of the space-discretisation matrix completely determine the stability of the solution**

**S completely determines the behaviour of the solution!**

Stability Condition
===================

For the ODE system:

.. math:: {d \mathbf{u} \over dt} = \mathbf{S} \cdot \mathbf{u} + Q

The exact solution :math:`\overline{\mathbf{u}}(t,\mathbf{x})` must remain bounded

All the modal components must be bounded

Require that exponential terms **do not grow** 

Hence **the real part of the eigenvalues must be negative or zero**

.. math:: Re(\Omega_j) \le 0 \quad \forall j
   :label: 1

:math:`\Omega`- plane: The eigenvalue spectrum has to be restricted to the left half plane

Note that if :eq:`1` is satisfied:

.. math:: \lim_{t \Rightarrow \infty} \mathbf{u}(t) = - \sum_{j=1}^N {Q_j \over \Omega_j} V^j
   :label: 2

this is a solution to the stationary problem:

.. math:: \mathbf{S} \cdot \overline{\mathbf{u}}_S = 0

Want :math:`\Omega` to have large negative real parts, to converge quickly to a steady state solution using an iterative method

.. math:: exp(-Re (\Omega_j)) \Rightarrow 0 

Matrix Method for Stability Analysis
====================================

From the previous examples: the structure of :math:`\mathbf{S}` depends on the BCs and how they are implemented (e.g. one sided difference etc)

We now have a method of investigating the influence of BCs on stability

We can also investigate effects of non-uniform meshes

**More general method than von Neumann analysis**

But, difficult for general boundary conditions to get analytical expressions for the eigenspectrum of :math:`\mathbf{S}`

We use periodic BCs :math:`\Rightarrow V^j = e^{Ik_jx}` (Fourier modes in 1D)

.. math:: x_i = v_i^j = e^{Ik_ji \Delta x} = e^{Ii \phi_j}

where :math:`\phi_j = k_j \Delta x`

Applied to the :math:`\mathbf{S}`:

.. math:: \mathbf{S} \cdot e^{Ik_j i \Delta x} = \Omega (\phi_j)e^{Ik_j i \Delta x}

Example 1 - 1D Diffusion
------------------------

CD scheme - leaving time discretisation undefined, so that we write ODEs:

.. math:: {d u_i \over dt} = {\alpha \over \Delta x^2} (u_{i+1} - 2 u_i + u_{i-1}) = 
                             \mathbf{S} \cdot \mathbf{u}_i

Using eigenfunctions: :math:`V^j (\mathbf{x}) = e^{Ik_j x}`

.. math:: \mathbf{S} \cdot e^{Ik_j i \Delta x} = 
          {\alpha \over \Delta x^2} \left( e^{Ik_j (i+1) \Delta x} - 2 e^{Ik_j i \Delta x} + e^{Ik_j (i-1) \Delta x} \right) =  
          {\alpha \over \Delta x^2}  \left( e^{Ik_j \Delta x} - 2 + e^{-Ik_j \Delta x} \right) e^{Ik_j i \Delta x} =
          {2 \alpha \over \Delta x^2} (cos \phi_j - 1) e^{Ik_j i \Delta x} 

Eigenvalues are real and negative between :math:`{-4 \alpha} / {\Delta x^2}` and :math:`0`

Example 2 - 1D Convection
-------------------------

With 1st order upwind

.. math:: \left( a {du \over dx} \right) = -{a \over {\Delta x}} (u_i - u_{i-1}) = \mathbf{S} \cdot \mathbf{u}_i

.. math:: \mathbf{S} \cdot e^{Ik_j i \Delta x} = 
          -{a \over {\Delta x}}  \left( e^{Ik_j i \Delta x} - e^{Ik_j (i-1) \Delta x} \right) = 
          -{a \over {\Delta x}}  \left( 1 - e^{Ik_j \Delta x} \right) e^{Ik_j i \Delta x} = 
          -{a \over {\Delta x}}  \left( 1 - cos \phi_j + I sin \phi_j \right) e^{Ik_j i \Delta x}

The eigenvalues are complex with negative real part, so stable

.. figure:: ../_images/convection_spectrum.png
   :align: center
   :scale: 70%

Example 3 - 1D Convection
-------------------------

With CD scheme


.. math:: \left( a {du \over dx} \right) = -{a \over {\Delta x}} (u_{i+1} - u_{i-1}) = \mathbf{S} \cdot \mathbf{u}_i

.. math:: \mathbf{S} \cdot e^{Ik_j i \Delta x} = 
          -{a \over {\Delta x}}  \left( e^{Ik_j (i+1) \Delta x} - e^{Ik_j (i-1) \Delta x} \right) = 
          -I {a \over {\Delta x}} sin \phi_j e^{Ik_j i \Delta x}

The eigenvalues are imaginary in range :math:`{{-Ia} / \Delta x}`, :math:`{{Ia} / \Delta x}`


.. figure:: ../_images/convection_spectrum_2.png
   :align: center
   :scale: 70%

Conclusion
==========

* **All 3 examples satisfy the stability condition**

* Also, negative real part contributes :math:`e^{-\left| I Re \Omega_j \right| t }`, this creates damping

* Numerical diffusion
