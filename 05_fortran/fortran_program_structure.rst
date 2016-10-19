=================
Program Structure
=================

.. contents::
   :local:

.. highlight:: fortran

Modular Programming
===================

Modules
-------

* Breaks up a task into **small steps**
* Prevents code having to be **re-written**, if we are just using the same again in another program
* Modules usually contain variables, parameters, functions and subroutines that are for the same purpose, e.g. a module for integration might contain the midpoint rule, the trapezoid rule and simpson's rule. This helps to add **meaning** to what is being programmed.

Programs
--------

* Usually used for **testing** the functions and subroutines of a module 
* Programs are often used for **input-output** operations, e.g. read/write to/from files or screen
* The **executable** file (.exe) is usually given the same name as the program

Functions
---------

* Functions return only one value
* Functions can be passed as values to other functions or subroutines

Subroutines
-----------

* Allow multiple returned values by modifying input arguments
* Can't be passed as arguments to other functions or subroutines (although I haven't tried this - I'm not sure what it would mean)

Steps
-----

There might be several steps as the program goes from simple to complex:

1) For **internal** subroutines/functions (for simple programs): 

 * Write a subroutine/function under **contains** as part of the main program. 
 * Use the internal subroutine/function e.g.
   
   - To use a function: ``y = fun_quadratic(a, b)`` 
   - To use a subroutine:  ``call sub_euler(x, y, z)`` 

2) For **external** subroutines/functions (for simple programs): 

 * Write an external subroutine/function under **end program** in the same file or in a different file to the main program (perhaps as the first step before writing a module)
 * Declare a subroutine using the **use** statement within the main program, e.g. ``use sub_euler`` 
 * Declare a function using the **external** attribute in the declaration in the main program, e.g. ``real(kind=8), external :: fun_quadratic``
 * Use the subroutine/function by name as if it were internal

3) To **modularise** a subroutine/function (for complex programs - e.g. to collect similar parts): 

 * Write a module which contains the subroutine/function
 * Declare the subroutine/function's **use** in the main program, e.g. ``use sub_euler, only: first_order``
 * Use the subroutine/function by name as if it were internal
   
4) **Pass a function as an argument** to a subroutine - this may allow the subroutine to be independent of the function passed to it, e.g. integration subroutines apply to all kinds of functions, e.g.

 * Use ``solve(fun_specific, x, y)`` in the main program (i.e. send a specific function such as a quadratic) 
 * The subroutine the function is passed could declare a more general argument list, e.g.

  - ``subroutine solve(fun_general, x, y)`` (i.e. it can recieve any kind of function)
  - ``real(kind=8), external :: fun_general`` (to declare that the function is external)

Probably can't pass subroutines as arguments, as they usually aren't single values.

Program
=======

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           program pro_integrate

               use mod_euler, lname => mod_1
               use mod_runga_kutta, only: fun_quadratic, fun_quartic
               
               implicit none

               interface
                  real function fun_cubic(a, b, c)
                      real, intent(in) :: a, b, c
                  end function fun_cubic
               end interface
                
               ! declaration and initialisation
               ! program body statements

               stop 'message'

           contains
               ! internal subroutines or functions

           end program pro_integrate              
 
     - **Program Structure**

       * ``program pro_integrate`` :math:`\Rightarrow` Start of program called ``pro_integrate``
       * ``use mod_euler, lname => mod_1`` :math:`\Rightarrow` Use external module with rename
       * ``use mod_runga_kutta, only: fun_quadratic, fun_quartic`` :math:`\Rightarrow` Selective use of **external** functions
       * ``implicit none`` :math:`\Rightarrow` Require variable declaration - use AFTER ``use`` statement OR after ``program``, ``module``, ``subroutine``, ``function`` or ``interface`` declaration
       * ``interface`` :math:`\Rightarrow` This is similar to a ``use`` statement for functions contained in external modules (but also requiring variable types and intent)
       * ``stop "message"`` :math:`\Rightarrow` Stops the program cleanly with a message
       * ``contains`` :math:`\Rightarrow` **internal** subroutines or functions (i.e. used and contained within the same file) - no external use is expected. 
       * ``end program pro_integrate`` :math:`\Rightarrow` End of program  

Use Explicit Imports, Not Implicit Imports
------------------------------------------

* Avoid:

::

   use mod_integrator

* Instead say what you are using:

::

   use mod_integrator, only: midpoint


Module
======

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           module mod_runge_kutta

               use mod_quadrature 
  
               implicit none

               private
               public :: fun_quadratic, fun_quartic
               
               interface
                   real function fun_cubic(a, b, c)
                       real, intent(in) :: a, b, c
                   end function fun_cubic
               end interface

               ! declaration and initialisation
               real(kind=8) :: pi
               contains

               ! internal subroutines or functions

            end module mod_runge_kutta
       
     - **Module Structure** - Put **external** subroutines and functions (i.e. external to program) within modules. Then call them via ``use`` statements within program. Modules can also use other modules. 

       * ``module mod_runge_kutta`` :math:`\Rightarrow` Start of module called ``mod_runge_kutta``
       * ``use mod_quadrature`` :math:`\Rightarrow` Use external module
       * ``implicit none`` :math:`\Rightarrow` Require variable declaration
       * ``public`` :math:`\Rightarrow` Visible outside module
       * ``private`` :math:`\Rightarrow` Visible inside module - keeping this empty makes everything is private by default, you then only make public what you want to make public
       * ``interface`` :math:`\Rightarrow` This is similar to ``use`` 
       * ``contains`` :math:`\Rightarrow` **Internal** subroutines or functions - could be public or private
       * ``end module runge_kutta`` :math:`\Rightarrow` End of module  
   * - * ::

           program pro_pi

               use mod_declare, only: pi
               implicit none

               call sub_initialise()            
               print *, "pi = ", pi

           end program pro_pi

       * ::

           subroutine sub_initialise
               
               use mod_declare, only: pi
               implicit none

               pi = acos(-1.d0)

           end subroutine sub_initialise

       * ::

           module mod_declare
                         
               implicit none
               real(kind=8) :: pi
               save

               ! use of pi in a formula

           end module mod_declare

     - **Module Variable** - Declare a **module variable** in a module and initialise it in a subroutine. This allows the program to initialise the module variable only once. 

       Modules that use a module variable must save the value from the initialisation (or it may revert to the default value, whatever that is).

       * ``save`` :math:`\Rightarrow` The declared variable retains it's value from the external initialisation


Subroutine
==========

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           subroutine sub_momentum_equation(nx, x, u, error, nu)
  
               implicit none
            
               integer, intent(in) :: nx
               real(kind=8), allocatable, dimension(:) :: x
               real(kind=8), intent(out) :: u
               integer, intent(inout) :: error
               real(kind=8), optional :: nu
           
               if(present(nu))
                   ! block
                   return
               else
                   ! block
                   return
               end if

           end subroutine sub_momentum_equation
       
     - **Subroutine Structure** - Subroutines differ from functions in that they are allowed to have more than one output. It does this by modifying it's input values (it doesn't actually return a new value). Subroutines may also have no output.

       * ``if(present(nu))`` :math:`\Rightarrow` Detects whether optional variable :math:`\nu` is present
       * ``return`` :math:`\Rightarrow` Forced exit of subroutine (before end statement) - this is useful if no more of the code need


Function
========


.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           real(kind=8) function fun_quadratic(a, b, c, x)
  
               implicit none
               real(kind=8), intent(in) :: a, b, c, x
           
               if(c .ne. 0) then
                   fun_quadratic = a*x**2 + b*x + c
                   return
               else 
                   fun_quadratic = a*x**2 + b*x
                   return
               end if

           end function fun_quadratic
       
     - **Function Structure 1** - Functions differ from subroutines, in that they can only have one output. Functions return a new value, not among the input arguments.

       * ``fun_quadratic`` :math:`\Rightarrow` New variable to be returned by function
       * ``return`` :math:`\Rightarrow` Forced exit of function (before end statement)
   * - ::

           function fun_quadratic(a, b, c, x)
  
               implicit none
               real(kind=8), intent(in) :: a, b, c, x
               real(kind=8) :: fun_quadratic

               ! block

           end function fun_quadratic
       
     - **Function Structure 2:** Notice that the output variable is declared with other variables in this form. Intent is implicit.

       * ``fun_quadratic`` :math:`\Rightarrow` New variable to be returned by function
   * - ::

           recursive function fun_recursion(a, b, c, x) result(y)
  
               implicit none
               real(kind=8), intent(in) :: a, b, c, x
               real(kind=8) :: y

               ! block

           end function fun_quadratic
       
     - **Function Structure 3:** "Recursive Function"

       * ``recursive`` :math:`\Rightarrow` Recursion allowed
       * ``result(y)`` :math:`\Rightarrow` Renames the expected output from ``fun_recursion`` to ``y``
   * - ::

           integer :: fun_increment
           integer :: x

           fun_increment(x) = x + 1
       
     - **Function Structure 4:** "Statement Function" 

       * ``integer :: fun_increment`` :math:`\Rightarrow` Declaration of "statement function" type
       * ``fun_increment(x)`` :math:`\Rightarrow` Takes input of x and adds one to it, assigning the result to ``fun_increment``

Interface
=========

I haven't used interfaces, but here are some definitions.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Code
     - Meaning
   * - ::

           interface

               ! external subroutine or functions 

           end interface

     - **Interface Structure 1** - These are found where a program or a module must use external subroutines or functions (similar to ``use`` statement)
   * - ::

           interface

               ! external subroutines or functions 

               module procedure list

           end interface

     - **Interface Structure 2** 

       * ``! external subroutines or functions`` :math:`\Rightarrow` External subroutines or functions
       * ``module procedure list`` :math:`\Rightarrow` Internal subroutines or functions
   * - ::

           interface operator op

               ! external subroutines or functions 

               module procedure list

           end interface

     - **Interface Structure 3** 

       * ``operator op`` :math:`\Rightarrow` Operator interface
   * - ::

           interface assignment (=)

               ! external subroutines or functions 

               module procedure list

           end interface

     - **Interface Structure 4** 

       * ``assignment (=)`` :math:`\Rightarrow` Conversion interface
