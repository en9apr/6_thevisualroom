====================
 Programming Tips 1
====================

.. contents::
   :local:


Docstrings
~~~~~~~~~~

-  Docstring tells you what the program does (interface) not how it does it.
-  Docstring is surrounded by three quotes, before and after

Comments
~~~~~~~~

-  Comments tells you how program does it (implementation)
-  Comments are preceded by a hash tag

Example of Docstring and Comment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. testcode::

    import math
    
    def area_triangle_sss(side1, side2, side3):
        """
        Returns area of a triangle, given the lengths of its three sides
        """
        
        # Heron's formula
        semi_perim = (side1 + side2 + side3) / 2.0
        return math.sqrt(semi_perim *
                         (semi_perim - side1) *
                         (semi_perim - side2) *
                         (semi_perim - side3))

Example of Always Giving Conditional Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Conditional Output should always be given if no criteria are met, i.e. always use else

.. testcode::

    def favourites(instructor):
        """
        Returns the favourite thing of a given instructor
        """
        
        if instructor == "Joe":
            return "games"
        elif instructor == "Scott":
            return "ties"
        elif instructor == "John":
            return "outdoors"
        else:
            print "Program saw invalid input:", instructor
        
    print favourites("John")
    print favourites("Jeannie")  #Returns None is there is no value
        

.. testoutput::

    outdoors
    Program saw invalid input: Jeannie
    None


