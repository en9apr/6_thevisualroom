=======================
 Conditional Sentences
=======================

.. contents::
   :local:

Semantics of Conditional Sentences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Conditional Sentences, like functions affect the flow of control
-  Conditional Sentences affect the flow of control by values

Syntax of Conditional Sentences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Conditional Sentences are compound sentences consisting of one or
   more clauses headed by the keywords:

   -  ``if``
   -  ``elif``
   -  ``else``

-  Each ``if`` or ``elif`` clause is followed by a Logical Expression
   and a colon ``:``
-  If the Logical Expression for a clause is ``True`` the body of the
   clause is executed
-  The body must be indented by 4 spaces (or a tab)

Examples of Conditional Sentences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. testcode::

    def greet(friend, money):        # Function has two predicates
        if friend and (money > 20):
            print "Hi!"
            money -= 20
        elif friend:
            print "Hello!"
        else: 
            print "Ha ha!"
            money += 10
        return money
    
    money = 15
    
    money = greet(True, money)
    print "Money:", money
    print ""
    
    money = greet(False, money)
    print "Money:", money
    print ""
    
    money = greet(True, money)
    print "Money:", money
    print ""

.. testoutput::

    Hello!
    Money: 15
    
    Ha ha!
    Money: 25
    
    Hi!
    Money: 5
    


.. testcode::

    # Conditionals Examples
    
    # Return True if year is a leap year, False otherwise
    
    def is_leap_year(year):          # Testing a predicate - convention to use is_function_name
        if(year % 400) == 0:
            return True
        elif (year % 100) == 0:
            return False
        elif (year % 4) == 0:
            return True
        else:
            return False
        
    year = 2010
    leap_year = is_leap_year(year)
    
    if leap_year:
        print year, "is a leap year"
    else:
        print year, "is not a leap year"
        

.. testoutput::

    2010 is not a leap year


