=============================
Tri-Diagonal Matrix Algorithm
=============================

.. contents::
   :local:

Use of the Tri-Diagonal Matrix Algorithm
========================================

The **Tri-Diagonal Matrix Algorithm (TDMA)** or **Thomas Algorithm** is a simplified form of Gaussian elimination that can be used to solve tri-diagonal systems of equations.  

Advantages of the TDMA:

* Less calculations and less storage than Gaussian Elimination
* Cost per unknown is independent of the number of unknowns (good scaling w.r.t. iterative methods)

Disadvantages of the TDMA:

* Round off error is still significant
* May not be suitable for non-linear problems **unless equations can be linearised using the Jacobian**
* Not stable in general, but it is stable if the matrix is **diagonally dominant** or **symmetric positive definite**

Usually direct calculation using TDMA is not used, but instead iterative solvers are used such as:

* Jacobi
* Gauss-Seidel
* SOR
* Conjugate Gradient
* Multi-grid
* Parallel Multi-grid

Matrix Configuration
====================

Tri-diagonal systems for :math:`nx` unknowns may be written as:

.. math:: a_i u_{i-1} + b_i u_i + c_i u_{i+1} = d_i

We know the values at the boundaries (:math:`B`):

.. math:: u_0 = B_0

.. math:: u_{nx-1} = B_{nx-1}

So the matrix looks like this, with known coefficients :math:`a, b, c, d`. The vector :math:`u` is unknown.

.. math::

   \begin{bmatrix}
   1 & 0 & 0 & \cdots & &0 \\
   a_1 & b_1 & c_1 & & & & \\
   0 & a_2 & b_2 & c_2 & & \\
   \vdots &  & \ddots & \ddots & \ddots & \vdots \\
   & &  a_{nx-3} & b_{nx-3} &  c_{nx-3} & 0\\
   & & & a_{nx-2} & b_{nx-2} &  c_{nx-2}\\
   0 & &\dots &0 &0 & 1
   \end{bmatrix}
   \begin{bmatrix}
   u_0 \\
   u_1 \\
   u_2 \\
   \vdots \\
   u_{nx-3} \\
   u_{nx-2} \\
   u_{nx-1}
   \end{bmatrix}
   =
   \begin{bmatrix}
   B_0 \\
   d_1 \\
   d_2 \\
   \vdots \\
   d_{nx-3} \\
   d_{nx-2} \\
   B_{nx-1}
   \end{bmatrix}

Influence of the Boundary Values
================================

Notice the equation on Line 1 can be re-arranged, since we don't need to solve for :math:`u_0`, since it equals :math:`B_0`:

.. math:: a_1 u_{0} + b_1 u_1 + c_1 u_{2} = d_1

.. math:: b_1 u_1 + c_1 u_{2} = d_1 - a_1 u_{0}

Hence:

.. math:: b_1 u_1 + c_1 u_{2} = d_1 - a_1 B_{0}

Similarly with line :math:`nx-2`, we don't need to solve for :math:`u_{nx-1}` since it equals :math:`B_{nx-1}`:

.. math:: a_{nx-2} u_{nx-3} + b_{nx-2} u_{nx-2} + c_{nx-2} u_{nx-1} = d_{nx-2}

.. math:: a_{nx-2} u_{nx-3} + b_{nx-2} u_{nx-2}  = d_{nx-2} - c_{nx-2} u_{nx-1}

Hence:

.. math:: a_{nx-2} u_{nx-3} + b_{nx-2} u_{nx-2}  = d_{nx-2} - c_{nx-2} B_{nx-1}

Re-arrangement of the Matrix
============================

.. math::

   \begin{bmatrix}
    b_1 & c_1 & & &  0 \\
    a_2 & b_2 & c_2 & & \\
    \vdots   & \ddots & \ddots & \ddots & \vdots \\
    & &  a_{nx-3} & b_{nx-3} &  c_{nx-3} \\
    0 & & & a_{nx-2} & b_{nx-2} 
   \end{bmatrix}
   \begin{bmatrix}
   u_1 \\
   u_2 \\
   \vdots \\
   u_{nx-3} \\
   u_{nx-2}
   \end{bmatrix}
   =
   \begin{bmatrix}
   d_1-a_1 B_0 \\
   d_2 \\
   \vdots \\
   d_{nx-3} \\
   d_{nx-2}-c_{nx-2} B_{nx-1} \\
   \end{bmatrix}

Step 1
------

Forward elimination adjusts the upper diagonal & RHS and eliminates lower diagonal:

.. math:: i = 1 \qquad \Rightarrow \qquad c_1^* = {c_1 \over b_1}

.. math:: i = 2,3,... nx-3 \quad \Rightarrow \quad c_i^* = {c_i \over {b_i - c_{i-1}^* a_i}}

.. math:: i = 1 \quad \Rightarrow \quad d_1^* = {d_1 \over b_1}

.. math:: i = 2,3,... nx-2 \quad \Rightarrow \quad d_i^* = {{d_i-d_{i-1}^* a_i} \over {b_i - c_{i-1}^* a_i}}

Step 2
------

The solution is then obtained by back substitution:

.. math:: i = nx-2 \qquad \Rightarrow \qquad u_{nx-2} = d_{nx-2}^*

.. math:: i = nx-3, nx-4,... 1 \quad \Rightarrow \quad u_i = d_i^* - c_i^* u_{i+1}

Alternative Step 1
------------------

Forward elimination adjusts the central diagonal & RHS and elminates the lower diagonal (From Ferziger and Peric)

The line at i=1 is unchanged

.. math:: i = 2,3,... nx-2 \quad \Rightarrow \quad b_i^* = b_i - { {a_i c_{i-1}} \over {b_{i-1}}}

.. math:: i = 2,3,... nx-2 \quad \Rightarrow \quad d_i^* = d_i - {{a_{i} d_{i-1}^* } \over {d_{i-1}}}

Alternative Step 2
------------------

The solution is obtained by back substitution:

Ferziger and Peric omitted this, but it's true:

.. math:: i = nx-2 \qquad \Rightarrow \qquad u_{nx-2} = {d_{nx-2}^* \over b_{nx-2}^*}

Then

.. math:: i = nx-3, nx-4,... 1 \quad \Rightarrow \quad u_i = {{d_i^* - {c_{i}^*}u_{i+1}} \over b_{i}^*}


How does this translate into Python Code?
=========================================

For the Implicit Beam-Warming Method:

.. math:: - {\Delta t \over {4 \Delta x}} \left( A_{i-1}^n u_{i-1}^{n+1} \right) + 
          u_i^{n+1} + {\Delta t \over {4 \Delta x}} \left( A_{i+1}^n u_{i+1}^{n+1} \right) = 
          u_i^n - {1 \over 2} {\Delta t \over \Delta x} (F_{i+1}^n - F_{i-1}^n) + 
          {\Delta t \over 4 \Delta x}(A_{i+1}^n u_{i+1}^n - A_{i-1}^n u_{i-1}^n)

For the inviscid Burgers' Equation:

.. math:: A_i^n = u_i^n

And:

.. math:: F_i^n = {{(u_i^n)}^2 \over 2}

Identifying :math:`a, b, c` and :math:`d`:

.. math:: a_i = - {\Delta t \over {4 \Delta x}} A_{i-1}^n

.. math:: b_i = 1

.. math:: c_i = {\Delta t \over {4 \Delta x}} A_{i+1}^n

.. math:: d_i =  u_i^n - {1 \over 2} {\Delta t \over \Delta x} (F_{i+1}^n - F_{i-1}^n) + 
          {\Delta t \over 4 \Delta x}(A_{i+1}^n u_{i+1}^n - A_{i-1}^n u_{i-1}^n)

Modification of the first and last values of the :math:`d` vector:

.. math:: d_1 =  d_1 - (- {\Delta t \over {4 \Delta x}} A_0^n)B_0

.. math:: d_{nx-2} =  d_{nx-2} - ({\Delta t \over {4 \Delta x}} A_{nx-1}^n)B_{nx-1}

Use a TDMA Solver
=================

This is from Ofan's Blog:

.. code-block:: python
   
    try:
        import numpypy as np # for compatibility with numpy in pypy
    except:
        import numpy as np # if using numpy in cpython
     
    ## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
    def TDMAsolver(a, b, c, d):
        '''
        TDMA solver, a b c d can be NumPy array type or Python list type.
        refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
        '''
        nf = len(a) # number of equations
        ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy the array
        for it in xrange(1, nf):
            mc = ac[it]/bc[it-1]
            bc[it] = bc[it] - mc*cc[it-1]
            dc[it] = dc[it] - mc*dc[it-1]

        xc = ac
        xc[-1] = dc[-1]/bc[-1]

        for il in xrange(nf-2, -1, -1):
            xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

        del bc, cc, dc # delete variables from memory

        return xc 

When implementing this, change the index of the coefficients
============================================================

Change index of coefficients, so that they run from :math:`0` to :math:`nx-3` (i.e. we need two less coefficients than there are in the :math:`u` vector).

.. math::

   \begin{bmatrix}
    b_0 & c_0 & & &  0 \\
    a_1 & b_1 & c_1 & & \\
    \vdots   & \ddots & \ddots & \ddots & \vdots \\
    & &  a_{nx-4} & b_{nx-4} &  c_{nx-4} \\
    0 & & & a_{nx-3} & b_{nx-3} 
   \end{bmatrix}
   \begin{bmatrix}
   u_1 \\
   u_2 \\
   \vdots \\
   u_{nx-3} \\
   u_{nx-2}
   \end{bmatrix}
   =
   \begin{bmatrix}
   d_0-a_0 B_0 \\
   d_1 \\
   \vdots \\
   d_{nx-4} \\
   d_{nx-3}-c_{nx-3} B_{nx-1} \\
   \end{bmatrix}

Introduction of Damping
=======================

The Beam-Warming method may become unstable. Therefore introduce fourth order damping term onto the RHS

.. math:: D_e = -\epsilon_e(\Delta x)^4 {\partial^4 u \over \partial x^4}

i.e. with fourth order CD:

.. math:: D_e = - \epsilon_e(u_{i+2}^n - 4u_{i+1}^n + 6u_i^n - 4u_{i-1}^n + u_{i-2}^n)

For the points :math:`i=0` and :math:`i=1`, the values of :math:`u_{-1}` and :math:`u_{-2}` are unknown, so set them to 1 (as no periodic boundaries are used).

For the points :math:`i=nx-1` and :math:`i=nx-2` the values of :math:`u_{nx}` and :math:`u_{nx+1}` are unknown, so set them equal to 0 (for the same reason)
