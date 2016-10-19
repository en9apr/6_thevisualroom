===================
Logical Expressions
===================

Syntax of Logical Expressions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Can defined constants equal to ``True`` and ``False`` of the type
   ``bool``
-  These constants can be combined to form a Logical Expression via the
   Logical Operators ``and``, ``or`` and ``not``

Semantics of Logical Expressions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``a and b`` is True only if both expressions are True
-  ``a or b`` is True if at least one expression is True
-  ``not a`` is True if the single expression is False

Looked at differently:

.. testcode::

    True and True is True
    True and False is False
    False and True is False
    False and False is False

    True or True is True
    True or False is True
    False or True is True
    False or False is False

    not True is False
    not False is True

Truth Table for NOT
~~~~~~~~~~~~~~~~~~~
==== =====
a    NOT a
==== =====  
T    F
F    T
==== =====

Truth Table for AND
~~~~~~~~~~~~~~~~~~~

==== ==== =======
a    b    a AND b
==== ==== =======
F    F    F      
F    T    F      
T    F    F     
T    T    T      
==== ==== =======

Truth Table for OR
~~~~~~~~~~~~~~~~~~

==== ==== =======
a    b    a OR b  
==== ==== ======= 
F    F    F     
F    T    T
T    F    T
T    T    T
==== ==== =======


Examples of Logical Expressions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. testcode::

    a = True
    b = False
    c = True
    d = False
    
    print a
    print b
    
    print "==="
    
    print not a
    print a and b
    print a or b
    
    print (a and b) or (c and (not d))

.. testoutput::

    True
    False
    ===
    False
    False
    True
    True


.. testcode::

    type(a)

.. testoutput::

    bool



