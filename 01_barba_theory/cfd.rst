============================
Computational Fluid Dynamics
============================

Meaning of the Navier-Stokes Equations
======================================

.. math:: {\partial \vec V \over \partial t} + \vec V (\nabla \cdot \vec V) = -{\nabla p \over \rho} + \nu \nabla^2 \vec V

	  \text{Unsteady term + Convective term = Pressure gradient term + Viscous term}

* The momentum equation in the Navier Stokes Equations is a 2nd order, non-linear partial differential equation:

  - :math:`\vec V (\nabla \cdot \vec V)` the convective term is the non-linear term
  - :math:`\nu \nabla^2 \vec V` the viscous term is also the 2nd order term

Assumptions:
````````````
* A Newtonian Fluid - shear stress proportional to strain
* Incompressible flow - density constant
* Isothermal flow - temperature constant
* If the fluid is inviscid :math:`\nu = 0`. Then we have the Euler Equations. But it's still non-linear:

.. math:: {\partial \vec V \over \partial t} + \vec V (\nabla \cdot \vec V) = -{\nabla p \over \rho}

What is a solution?
===================

* The velocity field (vector)
* The associated pressure field (scalar)

Why do we Need CFD?
===================

* There are very few known analytical solutions to the Navier Stokes Equations, e.g. when the convective term goes to zero, when there is flow in a pipe, flow between parallel plates and flow between concentric circles.

* Systems may be difficult to test through experimentation - e.g. experiment doesn't allow us to see inside or the instrumentation is limited, or it's dangerous to do the experiment.

* Faster and easier than experiment - CFD allows us to ask "what if" questions about a situation. Don't have to build prototype.

* CFD can be used in animation for films

Components of a CFD Model
=========================

Mathematical Model
``````````````````

* A set of partial differential equations or integral-differential equations

* Associated boundary conditions

* The non linearity is a source of turbulence, vorticity, shock waves, combustion, multi-phase, bubble dynamics, evaporation, condensation. Therefore the model is associated with a target application - e.g. incompressible, inviscid, turbulent, 2D or 3D model

Discretization Method - the major part of CFD
`````````````````````````````````````````````

* Method for approximating the partial differential equations or the integral-differential equations by a system of algebraic equations, i.e. we have a PDE :math:`\mathcal{L} [u(\underline{x})] = f(\underline{x})` and we need to convert it to arithmetic :math:`A\underline{x} = \underline{b}`

* The most important methods for doing this are:

  - Finite Difference (1950s)
  - Finite Element (1960s)
  - Spectral Methods (1970s)
  - Finite Volume (1980s)
  - Boundary Element Method
  - Particle Methods

* Discretization has two aspects:

  - The Geometry (grid or mesh or particles) - gives us a vessel for the solution
  - The Model - all mathematical operators converted into arithmetic operations on grid

Analyse Numerical Scheme
````````````````````````
* All numerical schemes must satisfy certain conditions to be accepted:

  - Consistency
  - Stability
  - Convergence
  - Accuracy

Solve
`````

* Obtain grid/point values of all flow variables

* Two situations:

  - Time dependent :math:`\Rightarrow` ODEs
  - Steady :math:`\Rightarrow` algebraic system of equations

* To solve these equations we require:

  - Time integrators
  - Linear solvers

Post Processing
```````````````

* Visualization
