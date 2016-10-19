======================
Relational Expressions
======================

.. contents::
   :local:

Relational Operators
~~~~~~~~~~~~~~~~~~~~

The values of two arithmetic expressions can be compared using the
operators:

 - `==`
 - `!=`
 - `<`
 - `>`
 - `<=`
 - `>=`

These return either:

 - `True`
 - `False`

Order of Precedence
~~~~~~~~~~~~~~~~~~~

Now we have looked at all the operators, we can give the order of precedence:

======= ========================= =========================
Order   Name                      Symbol 
======= ========================= =========================
1       Brackets                  ``()``
2       Unary                     ``not - +``
3       Binary Arithmetic 1       ``**``
4       Binary Arithmetic 2       ``* / %``
5       Binary Arithmetic 3       ``+ -``
6       Relational 1              ``< <= => >``
7       Relational 2              ``== !=``
8       Logical And               ``and``
9       Logical Or                ``or`` 
10      Assignment                ``= += -= *= /= %= **=``
======= ========================= =========================

Basically:

- Brackets
- Unary
- Arithmetic 1,2,3
- Relational 1,2
- Logical 1,2
- Assignment 


Examples of Relational Expressions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. testcode::

    # True and False not directly set, but we have comparison operators
    
    # Relational operators
    # >
    # <
    # >=
    # <=
    # ==
    # !=
    
    a = 7 > 3
    print a
    
    x = 5
    y = 5
    b = x > y
    print b
    
    c = "Hello" == 'Hello'  #only the text is compared not the quotes
    print c
    
    d = 20.6 <= 18.3
    print d

.. testoutput::

    True
    False
    True
    False


