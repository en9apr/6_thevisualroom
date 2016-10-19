===========
 Functions
===========

.. contents::
   :local:

Purpose of Functions
~~~~~~~~~~~~~~~~~~~~

-  Functions are reusable pieces of programs that take an input and
   produce an output.
-  Functions are **Abstract Methods** in other words they are not associated with an object and cannot be instantiated

Separation of Interface and Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Everything within the box is the implementation
-  Everything outside the box is the interface

.. image:: ../_images/functions_classless_methods_7_0.png


Comments
~~~~~~~~

-  Decide what the function must do and type this as a docstring

.. testcode::
   
   def function_name():
       """ Purpose and Interface of a function """
       return 0

-  Any comments about implementation w.r.t. mathematics can be included
   in the body as an inline comment:

.. testcode::

   a = simple_formula # Implementation comment

Syntax of Functions
~~~~~~~~~~~~~~~~~~~

-  A function definition consists of a header and a body
-  The header includes:

   -  the keyword ``def`` (short for define)
   -  the function name
   -  0 or more parameters enclosed by parentheses ``()``
   -  followed by a colon ``:``

-  The body includes:

   -  a sequence of sentences, all indented by 4 spaces
   -  can ``return`` or ``print`` a value
   -  can include the word ``pass`` if you know you need a function, but
      are not sure about implementation

Deciding Names
~~~~~~~~~~~~~~

-  All valid variable names are valid function names - but a function
   should not have the same name as a variable name - and should not be
   a built-in function
-  Function names must be descriptive, e.g. triangle_area
-  Parameter names must be descriptive, e.g. base and height

Unit Test Functions
~~~~~~~~~~~~~~~~~~~

-  To evaluate a function call, replace the function's parameters in the
   body of the function by their associated values in the call and
   execute the body of the function
-  Do this over a range of values
-  If there is no returned value, ``None`` is returned
-  If you forget the brackets in the function call, it won't throw an
   error, but assume that it's a variable

Global and Local Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Local Variables are only available inside a function
-  Global Variables are available inside and outside a function
-  If you want to modify Global Variables from within functions, they
   must be declared at the start using: global variable\_name


Built-in functions
~~~~~~~~~~~~~~~~~~

.. testcode::

    print int(3.0) #convert to int
    print float("3") #convert to float
    print str(34) #convert to string
    print abs(-5) #absolute value
    print max(3,5,6,7) #maximum of a set of numbers
    print min(-50, -40, -20) #minimum of a set of numbers
    print round(0.123456, 4) #round to 4 d.p.

.. testoutput::

    3
    3.0
    34
    5
    7
    -50
    0.1235

Examples of Functions
~~~~~~~~~~~~~~~~~~~~~

.. testcode::

    # Computes the area of a triangle
    def triangle_area(base, height):    # header
        area = (1.0 / 2) * base * height             # body of function
        return area                 # output         # body of function

.. testcode::

    a1 = triangle_area(3, 8)
    print a1   # you must always know what to expect from a function answer = 12

.. testoutput::

    12.0


.. testcode::

    a1 = triangle_area(4, 7)
    print a1  # output is 14

.. testoutput::

    14.0


.. testcode::

    # Converts farenheit to celcius
    def fahrenheit2celcius(fahrenheit):
        celcius = (5.0 / 9) * (fahrenheit - 32)
        return celcius
    
    # test
    c1 = fahrenheit2celcius(32)
    c2 = fahrenheit2celcius(212)
    print c1,c2

.. testoutput::

    0.0 100.0


.. testcode::

    # Converts fahrenheit to kelvin
    def fahrenheit2kelvin(fahrenheit):
        c = fahrenheit2celcius(fahrenheit) # code is reused - not re-written
        kelvin = c + 273.15
        return kelvin
    
    #test
    k1 = fahrenheit2kelvin(32)
    k2 = fahrenheit2kelvin(212)
    print k1, k2

.. testoutput::

    273.15 373.15


.. testcode::

    # prints hello world
    def hello():     #no inputs 
        print "Hello, world"
        #no outputs
    #test
    hello()
    h = hello()
    print h    # None shows because there is no returned value

.. testoutput::

    Hello, world
    Hello, world
    None

Lambda Functions
~~~~~~~~~~~~~~~~

Python supports the creation of anonymous functions (i.e. functions that are not bound to a name) at runtime, using a construct called "lambda". This is not exactly the same as lambda in functional programming languages, but it is a very powerful concept that's well integrated into Python and is often used in conjunction with typical functional concepts like filter(), map() and reduce().

This piece of code shows the difference between a normal function definition ("f") and a lambda function ("g"):

.. ipython::
   
   In [1]: def f (x): return x**2

   In [2]: print f(8)

   In [3]: g = lambda x: x**2

   In [4]: print g(8)


As you can see, f() and g() do exactly the same and can be used in the same ways. Note that the lambda definition does not include a "return" statement -- it always contains an expression which is returned. Also note that you can put a lambda definition anywhere a function is expected, and you don't have to assign it to a variable at all.


