===================
Numerical Schemes 2
===================

Definition of Consistency
=========================

* Consistency is a condition on the **numerical scheme**, i.e.

* The scheme must tend to the differential equation when the step in time and space tend to zero.

Definition of Stability
=======================

* Stability is a condition on the **numerical solution**, i.e.

* All errors must **remain bounded** when the iteration process progresses

* For finite values of :math:`\Delta t` and :math:`\Delta x`, the error has to remain bounded when the number of time steps :math:`n` tends to infinity.

* If the error:

.. math:: 
   \bar{\epsilon}_i^n = u_i^n-\bar{u}_i^n

where:

:math:`\qquad u_i^n` = computed solution

:math:`\qquad \bar{u}_i^n` = discretised solution

* Then the stability condition is:

.. math::

  \lim_{n \to \infty} \left \vert \bar{\epsilon}_i^n \right \vert \le \kappa \quad \text{  at fixed  } \Delta t

* **NOTE:**

1) The stability condition is a requirement on the **numerical scheme** only (does not require any condition on the differential equation)
2) Stability does not ensure that the **error** will not become unacceptably large at intermediate timesteps :math:`t^n = n \Delta t` (while remaining bounded)
3) A more general definition of stability later!

Definition of Convergence
=========================

* Convergence is a condition on the **numerical solution**
* The numerical solution must **tend** to the exact solution of the mathematical model, when steps in :math:`t` and :math:`x` tend to zero (i.e. when the mesh is refined).
 
* If the error:

.. math:: 
   \tilde{\epsilon}_i^n = u_i^n-\tilde{u}_i^n

where:

:math:`\qquad u_i^n` = computed solution

:math:`\qquad \tilde{u}_i^n = \tilde{u}_i^n(i\Delta x, n \Delta t)` = solution of differential equation

* **Convergence condition:**

.. math::
   \lim_{\Delta x \to 0, \Delta t \to 0} \left \vert \tilde{\epsilon}_i^n \right \vert = 0

Equivalence Theorem of Lax
==========================

* For a well-posed IVP and a **consistent** discretisation scheme, **stability** is the necessary and sufficient condition for **convergence**

.. figure:: ../_images/consistency_stability_convergence.png
   :scale: 100%
   :align: center



Consistency and the modified differential equation
==================================================

* A consistent scheme is one in which the **truncation error** tends to zero for :math:`\Delta t`, :math:`\Delta x \rightarrow 0`

* e.g. CD in space and FD in time for linear convection:

.. math:: {{u_i^{n+1} - u_i^n} \over {\Delta t}} + c {{u_{i+1}^n - u_{i-1}^n} \over {2 \Delta x}}=0 
   :label: one

* Taylor expansions:

.. math:: u_i^{n+1} = u_i^n + \Delta t \left . {\partial u \over \partial t} \right \vert_i^n + {\Delta t^2 \over 2} \left . {\partial^2 u \over \partial t^2} \right \vert_i^n + {\Delta t^3 \over 6} \left . {\partial^3 u \over \partial t^3} \right \vert_i^n + {\Delta t^4 \over 24} \left . {\partial^4 u \over \partial t^4} \right \vert_i^n + {\Delta t^5 \over 120} \left . {\partial^5 u \over \partial t^5} \right \vert_i^n
   :label: two

.. math:: u_{i+1}^n = u_i^n + \Delta x \left . {\partial u \over \partial x} \right \vert_i^n + {\Delta x^2 \over 2} \left . {\partial^2 u \over \partial x^2} \right \vert_i^n + {\Delta x^3 \over 6} \left . {\partial^3 u \over \partial x^3} \right \vert_i^n + {\Delta x^4 \over 24} \left . {\partial^4 u \over \partial x^4} \right \vert_i^n + {\Delta x^5 \over 120} \left . {\partial^5 u \over \partial x^5} \right \vert_i^n
   :label: three

.. math:: u_{i-1}^n = u_i^n - \Delta x \left . {\partial u \over \partial x} \right \vert_i^n + {\Delta x^2 \over 2} \left . {\partial^2 u \over \partial x^2} \right \vert_i^n - {\Delta x^3 \over 6} \left . {\partial^3 u \over \partial x^3} \right \vert_i^n + {\Delta x^4 \over 24} \left . {\partial^4 u \over \partial x^4} \right \vert_i^n - {\Delta x^5 \over 120} \left . {\partial^5 u \over \partial x^5} \right \vert_i^n
   :label: four

* Re-arranging Equation :eq:`two`:

.. math:: {{u_i^{n+1} - u_i^n} \over {\Delta t}} = \left . {\partial u \over \partial t} \right \vert_i^n + {\Delta t \over 2} \left . {\partial^2 u \over \partial t^2} \right \vert_i^n + O(\Delta t^2)
   :label: five

* Equation :eq:`three` minus Equation :eq:`four` and rearranging:

.. math:: {c \over {2 \Delta x}} {(u_{i+1}^n -  u_{i-1}^n)}  = c \left . {\partial u \over \partial x} \right \vert_i^n + c {\Delta x^2 \over 6} \left . {\partial^3 u \over \partial x^3} \right \vert_i^n + O(\Delta x^4)
   :label: six

* Including the truncation error in :eq:`one`, defined as the difference between the numerical approximation and the differential equation:

.. math:: {{u_i^{n+1} - u_i^n} \over {\Delta t}} + c {{u_{i+1}^n - u_{i-1}^n} \over {2 \Delta x}} - 
          \left ( \left . {\partial u \over \partial t}  \right \vert_i^n  +
          \left . c {\partial u \over \partial x} \right \vert_i^n \right )  =
          {\Delta t \over 2} \left . {\partial^2 u \over \partial t^2} \right \vert_i^n +
          c {\Delta x^2 \over 6} \left . {\partial^3 u \over \partial x^3} \right \vert_i^n +
          O(\Delta t^2) + O(\Delta x^4)

* The RHS of the above equation is the **truncation error** 

.. math::
   
   \epsilon_T =  {\Delta t \over 2} \left . {\partial^2 u \over \partial t^2} \right \vert_i^n +
          c {\Delta x^2 \over 6} \left . {\partial^3 u \over \partial x^3} \right \vert_i^n +
          O(\Delta t^2) + O(\Delta x^4)

* The scheme is consistent because when :math:`\Delta t`, :math:`\Delta x \rightarrow 0`, :math:`\epsilon_T \rightarrow 0`
