========================
 Arithmetic Expressions
========================

.. contents::
   :local:

Arithmetic Expressions
----------------------

Arithmetic Expressions are:

-  A number
-  Binary operator applied to two expressions ``** * / + -``
-  Unary operator applied to one expression ``+ -`` (but ``+`` is
   useless as a unary operator)

Expressions are formed by grouping operands and operators via
precedence. A simplified order is given below for common operators:

1. Please (parentheses) ``()``

2. Excuse (exponentiation) ``**`` (this includes roots e.g. :math:`\sqrt 2  = 2^{-0.5}`
  
3. My (multiplication) ``*`` Dear (division) ``/``

4. Aunt (addition) ``+`` Sally (subtraction) ``-``

Multiplication and division have equal precedence, so they are
done in the order they appear. The same applies to addition and
subtraction. You can remember this because Multiplication is the opposite
of Division and Addition is the opposite of Subtraction.

P E MD AS


Concatenate Expressions on Multiple Lines
-----------------------------------------

-  Use brackets ``()`` to implicitly concatenate over multiple lines -
   works for all operators

.. testcode::

   def this_function():
        if(1 == 1 or 
           3 != 43):
            return (0 + 0 - 0 * 0 - 0 / 1 - 0 ** 1 +
                1 + 1 + 1 + 1)

   print this_function()

.. testoutput::

   4

-  Or use backslash ``\``

.. testcode::

   def this_function_2():
        if(1 == 1 or 
           3 != 43):
            logical_statement = 1 == 2 or \
                                2 == 2
            return 0 + 0 - 0 * 0 - 0 / 1 - 0 ** 1 + \
                1 + 1 + 1 + 1

   print this_function_2()

.. testoutput::

    4

Number
------

Types of Number
~~~~~~~~~~~~~~~

-  An integer number - ``int()``
-  A decimal number - ``float()``
-  Use ``type()`` to determine the type of number

.. testcode::

   print type(3), type(3.142)

.. testoutput::
   
   <type 'int'> <type 'float'>

Conversion of Numbers
~~~~~~~~~~~~~~~~~~~~~

-  Convert from ``float`` to ``int`` using ``int()`` - **takes whole
   part of decimal number**

.. testcode::

   print int(3.142), int(-2.8)

.. testoutput::

   3 -2

-  Convert from ``int`` to ``float`` using ``float()`` - **adds decimal
   point**

.. testcode::
  
   print float(3), float(-1)

.. testoutput::

   3.0 -1.0

Floating Point Error
~~~~~~~~~~~~~~~~~~~~

- Floating point numbers are an approximation to decimal numbers. 
- Python floating point numbers can only display 12 decimal digits.
- This is called the "floating point error"

.. testcode::

   print 3.1415926535897932384626433832795028841971

.. testoutput::

   3.14159265359

Arithmetic Operators
--------------------

Types of Arithmetic Operator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  addition ``+``
-  subtraction ``-``
-  multiplication ``*``
-  division ``/``
-  exponentiation ``**``

.. testcode::

   print 1 + 2, 3 - 4, 5 * 6, 2 ** 5

.. testoutput::

   3 -1 30 32

Division in Python 2
~~~~~~~~~~~~~~~~~~~~

If one operand is a decimal (float), the answer is decimal

.. testcode::

   print 1.0 / 3, 5.0 / 2.0, -7 / 3.0

.. testoutput::

   0.333333333333 2.5 -2.33333333333

If both operands are ints, the answer is an int (next lowest integer
after division). Postive integers round "down" negative integers round "up"

.. testcode::

   print 1 / 3, 5 / 2, -7 / 3

.. testoutput::

   0 2 -3

Integer Division
~~~~~~~~~~~~~~~~

The integer division operator ``//`` returns the quotient of two numbers to get whole numbers from floats:

.. testcode::

   print 27.0//6.0

.. testoutput::

   4.0

``//`` is always the same as ``int()`` for +ve but not for -ve numbers:

- Integer division rounds toward -ve infinity for -ve numbers
- Converting float to integer rounds towards zero for -ve numbers

.. testcode::

   print "Same:", 6.0//5.0, int(6.0/5.0)
   print "Different", -6.0//5.0, int(-6.0/5.0)

.. testoutput::

   Same: 1.0 1
   Different -2.0 -1

Division by Zero
~~~~~~~~~~~~~~~~

Also, you should always check for division by zero, e.g.

.. testcode::
   
   numerator = 4
   denominator = 3 * 4 - 12

   if (denominator == 0):
       print "Error: Divide by zero"
   else:
       print numerator / denominator

.. testoutput::

   Error: Divide by zero

Standard Long Division
~~~~~~~~~~~~~~~~~~~~~~

Standard long division yields a quotient and a remainder.

Integer Division
~~~~~~~~~~~~~~~~

``//`` Yields the quotient

Modulus
~~~~~~~

``%`` Yields the remainder

By Definition
~~~~~~~~~~~~~

For any integers ``a`` and ``b``

``a == b * (a // b) + (a % b)``

Use of Modular Arithmetic
~~~~~~~~~~~~~~~~~~~~~~~~~

In Python ``a % b`` always returns an answer between 0 and b (even if a
and or b is negative)

Remainders and modular arithmetic are very useful in games for the
purpose of "wrapping" the canvas, i.e. causing objects that pass of one
side of the canvas to reappear on the opposite side of the canvas.

.. testcode::

    # problem - get the ones digit of a number
    
    num = 49
    tens = num // 10
    ones = num % 10
    print tens, ones
    print 10 * tens + ones

.. testoutput::

    4 9
    49


.. testcode::

    # application - 24 hour clock: what time is 8hrs after 20:00?
    
    hour = 20
    shift = 8
    print (hour + shift) % 24 # modulus has higher precedence than addition

.. testoutput::

    4

.. testcode::

    # application - screen wraparound: what is the position after it has moved through a width?
    
    width = 800
    
    position = 797
    move = 5
    
    position = (position + move) % width # what is left after we divide be screen  width
    print position
    
    # what happens if we moved backwards (move is negative)?
    
    position = 2
    move = -5
    
    position = (position + move) % width # we get back a number in the range 0 to width
    print position

.. testoutput::

    2
    797


.. testcode::

    # how do we convert from int to hours?
    
    hour = 3
    ones = hour % 10
    tens = hour // 10
    print tens, ones, ":00"
    print str(tens), str(ones), ":00" #srt converts int into a string
    print str(tens)+str(ones)+":00" #string formatting

.. testoutput::

    0 3 :00
    0 3 :00
    03:00


.. testcode::

    # division in python is not exact:
    
    print 0.9 % 0.3

.. testoutput::

    5.55111512313e-17


.. testcode::

    #modulus
    
    print 4 % 7 #Returns quotient 4 (as 4/7 = 0 r4)
    print 7 % 7 #Returns remainder 0 (as 7/7 = 1) 
    print 9 % 7 #Returns 2 (as 9/7 = 1 r2)
    
    # So the modulus will return 0 to 7 (not including 7) for any +ve number

.. testoutput::

    4
    0
    2


.. testcode::

    # Modulus with negative numbers:
    # First number - do we move in the +ve or -ve direction?
    # Second number - are we on the +ve or -ve number line?
    
    print 5 % 3   # a % b       returns values from 0 to 2
    print -5 % 3  # b - (a % b) returns values from 2 to 0
    print -5 % -3 # -(a % b)
    print 5 % -3  # -(b - (a % b))

.. testoutput::

    2
    1
    -2
    -1
