=========
 Modules
=========

.. contents::
   :local:

Semantics of Modules
~~~~~~~~~~~~~~~~~~~~

-  Modules are **Abstract Classes**, meaning they contain methods and
   constants, but cannot be instantiated.
-  Modules are libraries that implement useful operations not included
   in the basic code

Syntax of Modules
~~~~~~~~~~~~~~~~~

-  Modules can be accessed via the ``import`` statement
-  Examples of modules are ``math`` and ``random``
-  Don't use modules names as variable names

Math Module
~~~~~~~~~~~

.. testcode::

    import math
    
    print math.ceil(5.6) #Round up
    print math.floor(5.6) #Round down
    
    print math.pow(3,4) #3 to the power of 4
    
    print math.fabs(-5) #absolute value
    
    print math.sqrt(4) #square root
    
    print math.radians(180) #convert to radians
    print math.degrees(3.1415926) # convert to degrees
    
    print math.pi #pi
    print math.e # e^1
    print math.exp(1) #e^1 

    # import only two using capitals for constants

    from math import pi as PI
    from math import exp as exp

    print PI
    print exp(1)

.. testoutput::

    6.0
    5.0
    81.0
    5.0
    2.0
    3.14159265359
    179.99999693
    3.14159265359
    2.71828182846
    2.71828182846
    3.14159265359
    2.71828182846

Numpy Module
~~~~~~~~~~~~

.. testcode::

    import numpy as np
    
    ones=np.ones((5,5))
    print ones
    zeros=np.zeros((5,5))
    print zeros
    print np.sin(ones) 
    print np.linspace(0.0, 2.0, 5)

.. testoutput::
 
    [[ 1.  1.  1.  1.  1.]
    [ 1.  1.  1.  1.  1.]
    [ 1.  1.  1.  1.  1.]
    [ 1.  1.  1.  1.  1.]
    [ 1.  1.  1.  1.  1.]]
    [[ 0.  0.  0.  0.  0.]
    [ 0.  0.  0.  0.  0.]
    [ 0.  0.  0.  0.  0.]
    [ 0.  0.  0.  0.  0.]
    [ 0.  0.  0.  0.  0.]]
    [[ 0.84147098  0.84147098  0.84147098  0.84147098  0.84147098]
    [ 0.84147098  0.84147098  0.84147098  0.84147098  0.84147098]
    [ 0.84147098  0.84147098  0.84147098  0.84147098  0.84147098]
    [ 0.84147098  0.84147098  0.84147098  0.84147098  0.84147098]
    [ 0.84147098  0.84147098  0.84147098  0.84147098  0.84147098]]
    [ 0.   0.5  1.   1.5  2. ]

Random Module
~~~~~~~~~~~~~

.. testcode::

    import random
    
    print random.choice([1,2,3,4,5]) #Random choice from a list
    print random.randint(0,10) #random integer from 0 to 10 (including 10)
    print random.randrange(0,10) #random integer from 0 to 10 (not including 10)
    print random.random() #random float from 0.0 to 1.0 (not including 1.0)

.. testoutput::

    1
    3
    2
    0.622180301476



