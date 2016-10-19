=====================================
The Courant-Friedrichs-Lewy Condition
=====================================

.. contents::
   :local:

Comments on the 12 Steps to Navier Stokes
=========================================

* We have really used the simplest possible methods for the 12 Steps to Navier-Stokes

* These methods are unsuitable for practical purposes, because of diffusion e.g. the use of **1st order upwinding (BD in space)** for the convection term in the momentum equation for example.

* We must now develop new methods that are higher order methods.

Stability and Convergence, Modified Differential Equation and Truncation Error
==============================================================================

* For 1D convection, CD in space with FD in time produced a **Modified Differential Equation**, which was not a **convection equation**, but a **convection-diffusion equation** with a negative diffusion coefficient (helping to explain why the scheme was unstable). 

.. math:: \bar{u}_t + c \bar{u}_x = -{{\Delta t} \over 2} c^2 \bar{u}_{xx} + O({\Delta t^2}, \Delta x^2)

Negative diffusion coefficient represents an explosion (unstable)

* In contrast, 1st order upwind for 1D convection:

.. math:: \bar{u}_t + c \bar{u}_x = {{c \Delta x} \over 2} (1+{{c \Delta t} \over {\Delta x}}) \bar{u}_{xx}

Contained numerical diffusion proportional to :math:`\Delta x` (first order) which is why the diffusion is excessive

Stability condition :math:`c \gt 0` and :math:`\sigma = {{c \Delta t} \over \Delta x}` (CFL Condition)

**Then we performed von Neumann analysis for stability**

Summary of Stability
====================

Schemes can have:

* Conditional Stability
* Unconditional Stability
* Unconditional Instability

Implicit and Explicit Schemes
-----------------------------

Implicit schemes are generally unconditionally stable. Implicit methods have no such restriction on time step, but require more computation per timestep (you must solve a linear system of equations)

Explicit schemes lead at best to conditional stability. Conditional stability puts a limit on the timestep - we cannot progress too rapidly in time with the numerical scheme - as it may nnot be able to "trasmit" the information in the solution

**Convection:** The limit on the time step can be a severe requirement, especially for convection dominated equations (c is large compared to other effects). For most explicit schemes:

.. math:: {{c \Delta t} \over {\Delta x}} \le 1
 
**Diffusion:** For diffusion dominated equations, the restriction is less severe, i.e.

.. math:: {{\nu \Delta t} / {\Delta x^2}} \le {1 \over 2}

Physical Interpretation of the CFL Condition
============================================

1st order upwind CFL condition (for most explicit schemes):

.. math:: {{c \Delta t} \over {\Delta x}} \le 1

**Physical Interpretation:** The distance travelled by the solution in one timestep :math:`c \Delta t` must be less than the distance between two mesh points :math:`\Delta x`

.. figure:: ../_images/domain_of_dependence.png
   :align: center
   :scale: 70%

* Hyperbolic PDEs (such as the one for linear convection) have characteristic lines along which the solution travels
* For upwinding (backward differencing in space) the domain of dependence at `n` is from `i-1` to `i`. This is suitable for a positive wave speed `c`
* For forward differencing the domain of dependence at `n` is from `i` to `i+1`. This is suitable for a negative wave speed `c`
* **The solution at the next timestep must be able to include all the physical information that influences the solution from the previous timestep**
* The CFL condition :math:`\sigma \lt 1` ensures that the domain of dependence of the governing equation is **entirely** contained in the domain of dependence of the numerical scheme
* Can extend this to more complex cases where deriving the stability condition is more difficult for more complex numerical schemes.
* Also demonstrates why backward differencing is unstable for a negative wave speed, i.e. if the wave move from right to left, the solution should depend on points `i` and `i+1` not `i` and `i-1`
