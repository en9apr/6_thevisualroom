======================
Python Quick Reference
======================

Python is a high-level open-source language. This page is written with the **use** of Python in mind, as opposed to computer science definitions. iPython is used here because:

* print statement can be avoided (you can just type the variable and get output)
* the output is automatically computed by the iPython interpreter (no need for testoutput blocks)
* it allows reference to variables beyond the code block

.. contents::
   :local:

Modules
=======

Definition
----------

* Modules provide useful functions, such as array operations, plotting and much more.
* We can import the module to allow access to the functions 

.. ipython::
   
   In [1]: # comments in python are denoted by the hash tag

   In [2]: import numpy as np # for matrix operations

   In [3]: import matplotlib.pyplot as plt # for 2D plotting

* We have defined the following aliases:

  =================== ======= =====================
  Module              Alias   Purpose
  =================== ======= =====================
  numpy               np      Matrix Operations
  matplotlib.pyplot   plt     2D Plotting
  =================== ======= =====================

Use
---

* This allows reference to the modules via the alias, e.g.

.. ipython::
   
   In [1]: myarray = np.linspace(0,5,10)

   In [2]: myarray

* If you don't preface the linspace function with `np` python will throw an error:
* To learn the new functions available to you, visit: `Numpy Reference <http://docs.scipy.org/doc/numpy/reference/>`_
* Or for Matlab users: `Numpy for Matlab Users <http://wiki.scipy.org/NumPy_for_Matlab_Users/>`_

Variables
=========

Definition
----------

Python doesn't require explicitly declared variable types like C and Fortran.

.. ipython::

   In [1]: a = 5 # a is an integer 5

   In [1]: b = 'five' # b is a string of the word 'five'

   In [1]: c = 5.0 # c is a floating point 5

Use
---

Use `type` to determine the type of variable:

.. ipython::

   In [1]: type(a)

   In [2]: type(b)

   In [3]: type(c)

Division using integers is in two ways:

* Integer / Integer = Integer

.. ipython::
   
   In [1]: 14 / 3 

* Integer / Float = Float

.. ipython::

   In [1]: 14 / 3.0 

Whitespace in Python
====================

Python uses indents and whitespaces to group statements together. To write a loop in C you might use:

.. code-block:: c

   for (i = 0, i < 5, i++){
       printf("Hi! \n");
   }

Or in Fortran:

.. code-block:: fortran
   
   do i = 1, 4
       print *, "Hi!" 
       print *, " "
   enddo

Python doesn't use curly braces like C or the enddo like Fortran. The Python equivalent is:

.. code-block:: python

   for i in range(5):
       print "Hi! \n"


If you have nested for-loops, there is a further indent for the inner loop:

.. ipython::

   In [1]: for i in range(3):
      ...:     for j in range(3):
      ...:        print i, j
      ...:     print "This statement is within the i-loop but not the j-loop"


Arrays
======

Creation of Arrays
------------------

* NumPy arrays are created from lists

.. ipython::

   In [1]: myvals = np.array([1,2,3,4,5])
   
   In [2]: myvals

   In [3]: type(myvals)

* Or using linspace:

.. ipython::

   In [1]: myvals2 = np.linspace(1,5,5)

   In [2]: myvals2

   In [3]: type(myvals2)

Index
-----

* Python uses a **zero-based index**:
 
  - The first element is `0`
  - The last element is `n-1` where `n` is the number of values in the array

.. ipython::

   In [1]: myvals[0], myvals[4] 

Slicing Arrays
--------------

* The slice is inclusive on the front end and exclusive on the back, so the following command gives us the values of `myvals[0]`, `myvals[1]` and `myvals[2]`
* `myvals[3]` is excluded

.. ipython::

   In [1]: myvals[0:3]

* For a[start:end] THE INDEX DENOTES THE POSITION OF THE SLICE - NOT THE ELEMENT:

.. code-block:: python

    #                  +------+------+---   ---+------+------+
    #                  |  11  |  22  |   33    |  44  |  55  |
    #                  +------+------+---   ---+------+------+
    # Long version:    0      1      2      nx-2   nx-1     nx
    # Short version:                          -2     -1      

.. ipython::

   In [1]: hello = np.array([11, 22, 33, 44, 55])

* So now test it:
 
.. ipython::

   In [1]: hello[0:-1]

   In [2]: hello[1:-2]

   In [3]: hello[:-1]


Slicing Arrays (Short Notation)
-------------------------------

.. ipython::

   In [1]: alpha = ['a', 'b', 'c', 'd', 'e', 'f']

Indexing
~~~~~~~~

Index notation occurs when only one character appears in the square brackets

* The first value:

.. ipython::

   In [1]: alpha[0]

* The last value:

.. ipython::

   In [1]: alpha[-1]


* All the values

.. ipython::

   In [1]: alpha[:]


Slicing
~~~~~~~

The slice notation occurs with a number included to the left or right of the colon

* From the second value onwards:

.. ipython::

   In [1]: alpha[1:]

* From the third value onwards:

.. ipython::

   In [1]: alpha[2:]


* All the values except the last value:

.. ipython::

   In [1]: alpha[0:-1]

   In [2]: alpha[:-1]

* All the values except the last two values:

.. ipython::

   In [1]: alpha[0:-2]

   In [2]: alpha[:-2]

* From the second value to the second to last value:

.. ipython::

   In [1]: alpha[1:-1]

Assigning Array Variables
-------------------------

One of the strange little quirks/features in Python that often confuses people comes up when assigning and comparing arrays of values.

* Create 1D array called `a`:

.. ipython::

   In [1]: a = np.linspace(1,5,5)
   
   In [2]: a

* Make a copy of `a` and call it `b` (this is actually assignment by reference)

.. ipython::

   In [1]: b = a

   In [2]: b

* Now try changing the values in `a`:

.. ipython::

   In [1]: a[2] = 17
 
   In [2]: a

* But this also changed `b`!

.. ipython::

   In [1]: b

**Explanation:** Python created a pointer called `b` that tells us to route it to `a`. This is called **assignment by reference**.

Copying Arrays
--------------

If you want to make a true copy of the array you have to tell Python to copy every element of `a` into a new array: 

* Create an empty array `c` the same length as `a`:

.. ipython::

   In [1]: c = np.empty_like(a)
 
   In [2]: len(c) # tells us how long c is

  
* Copy the values from `a` to `c`:

.. ipython::

   In [1]: c[:] = a[:]

   In [2]: c

* Now change a value in `a`, which doesn't change `c`:

.. ipython::

   In [1]: a[0] = 200

   In [2]: a

   In [3]: c

Array Operations
----------------

Operators
~~~~~~~~~

Addition on a list is concatenation:

.. ipython::

   In [1]: a = [1,2,3]

   In [2]: a + a

Arithmetic operations on an array are element-wise:

.. ipython::

   In [1]: b = np.array([1,2,3])

   In [2]: b + b # element-wise addition

   In [3]: b - b # element-wise subtract

   In [4]: b / b # element-wise divide

   In [5]: b * b # element-wise muliply

   In [6]: b ** 2 # element-wise power

Functions
~~~~~~~~~

Functions on arrays are element-wise:

.. ipython::

   In [2]: np.sin(b) # element-wise sin - use np.sin for this operation not math.sin

Dot Product
~~~~~~~~~~~

.. math::

   \begin{bmatrix}
   1 & 2 & 3 \\
   4 & 5 & 6
   \end{bmatrix} \cdot
   \begin{bmatrix}
   1 & 2 \\
   3 & 4 \\
   5 & 6
   \end{bmatrix} = 
   \begin{bmatrix}
   22 & 28 \\
   49 & 64
   \end{bmatrix}

Dot product on an array:

.. ipython::

   In [1]: horz = np.array([[1,2,3],[4,5,6]])

   In [2]: vert = np.array([[1,2],[3,4],[5,6]])

   In [3]: np.dot(horz,vert)

Array Creation
~~~~~~~~~~~~~~

Lists can be created using `range`:

.. ipython::

   In [1]: range(5)

Arrays can be created using `arange` or `linspace`:

.. ipython::

   In [2]: np.arange(5)

   In [3]: np.linspace(0,4,5)

List Comprehension
~~~~~~~~~~~~~~~~~~

* Imagine a list of three numbers:

.. ipython::

   In [1]: c = [1,2,3]

* The list comprehension:

.. ipython::

   In [2]: cc = [x+y for x,y in zip(c,c)]

   In [3]: cc

* Where `zip` returns a list of tuples, i.e.

.. ipython::

   In [1]: zip(c,c)

* The list comprehension is equivalent to:

.. ipython::

   In [1]: dd = []

   In [2]: for x,y in zip(c,c):
      ...:    dd.append(x+y)

   In [3]: dd

* Or:

.. ipython::

   In [1]: e = np.array(c)

   In [2]: ee = e + e

   In [3]: ee

Which is Faster: Lists or Arrays?
---------------------------------

Create a list and time a list comprehension:

.. ipython::
  
   In [1]: f = range(10000)

   In [2]: timeit ff = [x + y for x,y in zip(f,f)]

Create a NumPy array and time the addition:

.. ipython::

   In [1]: g = np.array(f)

   In [2]: timeit gg = g + g

   In [3]: timeit hh = np.add(g,g)

* NumPy is over 100 times faster than Lists. 
* Not much between `np.add` and `+`.
*  Readability of `+` probably outweighs slight speed penalty.
