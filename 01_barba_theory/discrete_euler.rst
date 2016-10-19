==================================
 Discretising the Euler Equations
==================================

.. contents::
   :local:

Using a Central Difference formula for a **vector** function:

.. math:: {{\partial \mathbf{F}(x_i)} \over \partial x}  =
          {{\mathbf{F}(x_{i+1}) - \mathbf{F}(x_{i-1})} \over 2 \Delta x} =
          {{\mathbf{F}_{i+1} - \mathbf{F}_{i-1}} \over 2 \Delta x}

Simply replaced scalar u with vector :math:`\mathbf{F}`

In all our schemes, similarly replace the scalar derivative:

.. math:: A = {{d F} \over {d u}} 

By the Jacobian matrix:

.. math:: \mathbf{A} = {{d \mathbf{F}} \over {d \mathbf{U}}}

Lax-Friedrichs
==============

FTCS n to n+1, with spatial average for u at n

.. math:: \mathbf{U}_i^{n+1} = {1 \over 2}(\mathbf{U}_{i+1}^n + \mathbf{U}_{i-1}^n) +
                               {\Delta t \over {2 \Delta x}} ({\mathbf{F}_{i+1}^n - \mathbf{F}_{i-1}^n})


Lax-Wendroff
============

* Taylor at n+1, 
* First derivative - replace du/dt with -dF/dx and then CS at i
* Second derivative - replace d/dt (du/dt) with d/dt(-dF/dx),
                  - then replace d/dx (-dF/dt) with d/dx(A dF/dx)
                  - then outer derivative is CS at i+0.5 and CS at i-0.5
                  - inner derivative is CS at i+0.5 and i-0.5
                  - Jacobian is average at i+0.5 and i-0.5


.. math:: \mathbf{U}_i^{n+1} = \mathbf{U}_i^n - {\Delta t \over {2 \Delta x}}({\mathbf{F}_{i+1}^n - \mathbf{F}_{i-1}^n}) +
          {{\Delta t^2} \over {2 \Delta x^2}} \left[ {\mathbf{A}_{i+{1/2}}^n(\mathbf{F}_{i+1}^n - \mathbf{F}_i^n) -
                                                  \mathbf{A}_{i-{1/2}}^n(\mathbf{F}_i^n - \mathbf{F}_{i-1}^n)} \right]

e.g.

.. math:: \mathbf{A}_{i+{1/2}}^n = {1 \over 2}( \mathbf{A}_{i+1}^n + \mathbf{A}_i^n )

.. math:: \mathbf{A}_{i-{1/2}}^n = {1 \over 2}( \mathbf{A}_i^n + \mathbf{A}_{i-1}^n )

Note that :math:`\mathbf{A}` has 9 entries and so evaluating it is expensive 

Plus we have Matrix-Vector multiplication - thus this method is expensive.


Richtmyer
=========

Step 1: Predictor - Lax-Friedrichs: FTCS n to n+1/2, at i=i+1/2, with spatial average for u at n):

.. math:: \mathbf{U}_{i+1/2}^{n+1/2} = {1 \over 2}(\mathbf{U}_{i+1}^n + \mathbf{U}_{i}^n) -
                               {\Delta t \over {2 \Delta x}} ({\mathbf{F}_{i+1}^n - \mathbf{F}_{i}^n})

Step 2 (Corrector - Leapfrog: FTCS n to n+1)

.. math:: \mathbf{U}_{i}^{n+1} = \mathbf{U}_i^n - {\Delta t \over {\Delta x}}({\mathbf{F}_{i+1/2}^{n+1/2} - \mathbf{F}_{i-1/2}^{n+1/2}})

MacCormack
==========

Step 1 (Predictor - FTFS, n to n+1, with :math:`\Delta t`):

.. math:: \mathbf{\tilde{U}}_i^{n+1} = \mathbf{U}_i^n - {{\Delta t} \over {\Delta x}}({\mathbf{F}_{i+1}^n - \mathbf{F}_i^n})

Step 2 (Corrector - FTBS, n+1/2 to n+1, with :math:`\Delta t/2`):

.. math:: \mathbf{U}_i^{n+1} = {1 \over 2} (\mathbf{\tilde{U}}_i^{n+1} + \mathbf{U}_i^n) -
                               {{\Delta t} \over {2 \Delta x}}({\mathbf{\tilde{F}}_i^{n+1} - \mathbf{\tilde{F}}_{i-1}^{n+1}})
