=============
ODEs and PDEs
=============

.. contents::
   :local:

.. highlight:: latex

How to decide if a PDE is elliptic, hyperbolic or parabolic?
============================================================

.. math::

    a {{\partial^2 u} \over {\partial x^2}} + 
    b {{\partial^2 u} \over {\partial x \partial y}} +
    c {{\partial^2 u} \over {\partial y^2}} = 
    g(x, y)
    
or for a :math:`n \times n` system:

.. math::

    {\partial \vec{U} \over \partial t} + \vec{A} {{\partial \vec{U}} \over {\partial x}}=0
The type depends on the determinant:

* :math:`b^2-4ac>0 \longrightarrow` all n eigenvalues are **real** :math:`\longrightarrow` **Hyperbolic**
* :math:`b^2-4ac=0 \longrightarrow` all n eigenvalues are **real and identical** :math:`\longrightarrow` **Parabolic**
* :math:`b^2-4ac<0 \longrightarrow` all n eigenvalues are **complex** :math:`\longrightarrow` **Elliptic**

What common flow conditions are hyperbolic, elliptic and parabolic?
===================================================================

* **Hyperbolic** :math:`\longrightarrow`  (Steady :math:`\lor` Unsteady) :math:`\land` **Compressible** :math:`\land` **Inviscid** :math:`\land`  **Supersonic** 

* **Parabolic** :math:`\longrightarrow` (Compressible :math:`\lor` Incompressible) :math:`\land` **Viscous Boundary Layer** 

* **Elliptic** :math:`\longrightarrow` (Compressible :math:`\lor` Incompressible) :math:`\land` **Inviscid** :math:`\land` **Steady** :math:`\land` **Subsonic**

* **Mixed** (most cases are mixed) :math:`\longrightarrow` There is a parabolic boundary layer, but parts of the domain away from the boundary are hyperbolic or elliptic

What are some examples of hyperbolic, elliptic and parabolic PDEs?
==================================================================

* Wave Equation :math:`\longrightarrow` **Hyperbolic** :math:`\longrightarrow` :math:`{\partial u \over \partial t} + a {{\partial u} \over {\partial x}}=0`
* Laplace Equation :math:`\longrightarrow` **Elliptic** :math:`\longrightarrow`  :math:`{\partial^2 u \over \partial x^2} + {{\partial^2 u} \over {\partial y^2}}=0`
* Poisson Equation :math:`\longrightarrow` **Elliptic** :math:`\longrightarrow`  :math:`{\partial^2 p \over \partial x^2} + {{\partial^2 p} \over {\partial y^2}}=f(u)`
* Stokes Flow :math:`\longrightarrow` **Parabolic** :math:`\longrightarrow` :math:`{\partial p \over \partial x} = \mu {\partial^2 u \over \partial x^2}`
* Heat Conduction :math:`\longrightarrow` **Parabolic** :math:`\longrightarrow` :math:`{\partial T \over \partial t} = \alpha {\partial^2 T \over \partial x^2}`
* Boundary Layer :math:`\longrightarrow` **Parabolic** :math:`\longrightarrow` :math:`\rho \left( {u {\partial u \over \partial x} + v {\partial u \over \partial y}} \right) = \mu {\partial^2 u \over \partial y^2}`
