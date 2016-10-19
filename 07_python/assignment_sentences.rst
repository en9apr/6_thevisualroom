======================
 Assignment Sentences
======================

.. contents::
   :local:

Syntax of Assignment Sentences
==============================

-  Single equals ``=`` is used for assignment to variables
-  Double equals ``==`` is used for testing equality
-  Plus equals ``+=`` is used to add to an existing variable (can also
   be used with all other operators)

Semantics of Assignment Sentences
=================================

Setting two variables equal to each other copies the value stored in one
and assigns it to the other. However, if one of the variables is changed
later, the other one will not change e.g.

.. testcode::

    # Assignment
    # =
    # +=
    # -=
    # /=
    # *=
    # **=
    
    a=2
    b=3
    print "example 1:",a,b
    a=b
    print "example 2:", a,b
    b=4
    print "example 3:",a,b

.. testoutput::

    example 1: 2 3
    example 2: 3 3
    example 3: 3 4


.. testcode::

    my_age = 51
    print my_age

.. testoutput::

    51


.. testcode::

    #birthday  - add one
    my_age = 51 + 1
    print my_age
    
    #better
    my_age = my_age + 1
    print my_age
    
    #better still
    my_age += 1
    print my_age

.. testoutput::

    52
    53
    54


