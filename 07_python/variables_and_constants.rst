=========================
 Variables and Constants
=========================

.. contents::
   :local:


Variables
---------

-  Placeholders for important values
-  Make program more efficient
-  Make program more understandable

Variable Names
~~~~~~~~~~~~~~

-  Letters, numbers, underscore
-  Starts with letter or underscore
-  Case sensitive

Good Names
^^^^^^^^^^

-  Variable names should be descriptive
-  Python convention is to join multiple words with underscore in
   variable names, e.g. ``favorite_number = 18``

Bad Names
^^^^^^^^^

-  Variables cannot start with a number, e.g. ``1337ninja``
-  Don't use keywords as variables (names that change colour), e.g.
   ``min``, ``max``, ``int``, ``print``, ``import``
-  Avoid too longer a name

More examples of variables
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. testcode::

    my_name = "joe warren"
    print my_name

.. testoutput::

    joe warren


.. testcode::

    my_age = 51
    print my_age

.. testoutput::

    51


.. testcode::

    #the story of the magic pill
    magic_pill = 30
    print my_age - magic_pill
    
    my_granddad = 74
    
    print my_granddad - 2 * magic_pill

.. testoutput::

    21
    14


.. testcode::

    #Temperature example
    
    temp_fahrenheit = 212
    temp_celcius = 5.0 / 9.0 * (temp_fahrenheit - 32)
    print temp_celcius # 0degC = 32degF and 100degC = 212degF

.. testoutput::

    100.0


.. testcode::

    temp_celcius = 100
    temp_fahrenheit = 9.0 / 5.0 * temp_celcius + 32
    print temp_fahrenheit # 0degC = 32degF and 100degC = 212degF

.. testoutput::

    212.0


Constants
---------

Python has no constants, just use capital letters, e.g.

.. testcode::

    ACCELERATION_G = 9.81
    TRUCK_MASS = 10000
    WATER_VISCOSITY = 10**-3
