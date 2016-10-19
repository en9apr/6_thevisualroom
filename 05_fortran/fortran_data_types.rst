==========
Data Types
==========

.. contents::
   :local:

.. highlight:: fortran

Primitive Data Types
====================

Declaration and Initialisation
------------------------------

* Fortran does not initialise values, these will be "undefined" unless you give an initial value
* It's best to **separate declaration from initialisation**, because intialisation and declaration together creates a static variable (which may not be what you intended)

.. list-table::
   :header-rows: 1
   :widths: 40 50

   * - Data Type
     - Meaning
   * - ::

             integer :: i, j
             i = 1
             j = 0
     - **Integer** declaration with default 4 bytes and initialisation
   * - ::

             real :: r
             r = 1.0e-7
     - **Real** declaration with default 4 bytes and initialisation (e is for single precision exponent)
   * - ::

             double precision :: d
             d = 2.0d-12
     - **Double precision** declaration with default 8 bytes, equivalent to real(kind=8) and initialisation (d means double precision exponent)
   * - ::

             complex :: z
             z = cmplx(1.0, 1.0)
     - **Complex** declaration with default 8 bytes and initialisation
   * - ::
 
             logical :: b, c
             b = .true.
             c = .false.
     - **Boolean** declaration with default 4 bytes and initialisation
   * - ::

             character :: s, t
             s = "string_1"
             t = 'string_2'
     - **String** declaration with default 1 byte and initialisation (can use single or double quote)

Specifications
~~~~~~~~~~~~~~

These are included in brackets after the data type, e.g.

::

    ! real(specification) :: name_1, name_2

.. list-table::
   :header-rows: 1
   :widths: 40 50

   * - Specification
     - Meaning
   * - ::

           real(kind = 8) :: value
           value = 3.142d20
     - Real declaration with **8 bytes** (double precision)
   * - ::

           character(len = 10) :: s
           s = "string_1"
     - String declaration with length **10 bytes** (10 characters)
   * - ::

           character(len = *) :: s
     - String declaration with length **declared elsewhere**

Specifications are based on bytes and a byte is:

* **One** character
* An integer between **-128 and 127**
* The logical values ``.true.`` and ``.false.``

The meaning of the specifications is in the possible range of values:

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Exponent
     - Range
   * - ::

           integer :: short_integer
     - -2147483648 to 2147483647
   * - ::

           integer(kind=8) :: long_integer
     - -9223372036854775808 to 9223372036854775807
   * - ::

           real :: small_number
     - 1.175E-38 to 3.402E+38
   * - ::

           real(kind=8) :: large_number
     - 2.225D-308 to 1.798D+308

Attributes
~~~~~~~~~~

These are included  after the data type and after the specification, e.g.

::

    ! real(specification), attributes :: name_1, name_2

Quite a lot of the time we use double precision, so I used that is this example

.. list-table::
   :header-rows: 1
   :widths: 40 50

   * - Specification
     - Meaning
   * - ::
 
           real(kind=8), parameter  :: pi = acos(-1.d0)
     - Real **constant** with 4 bytes - it's ok to declare and initialise on the same line here, as it's a constant
   * - ::

           real(kind=8), intent(in) :: x
     - Real variable is expected as an **input**
   * - ::
 
           real(kind=8), intent(out) :: y
     - Real variable is expected as an **output**
   * - ::

           real(kind=8), external :: f
     - Result of function f is real and is evaluated **externally** this can allow function names to be used as arguments to other functions/subroutines, or as the equivalent to **call** for subroutines that are in the same file, but external to a program
   * - ::
 
           real(kind=8), dimension(n) :: array_one
     - Array of **rank** 1 (1D), with **shape** nx1 containing real values and **extent** of possible indices are 1 (default)
   * - ::

           real(kind=8), dimension(10) :: array_two
     - Array of **rank** 1 (1D), with **shape** 10x1 containing real values and **extent** of possible indices are 1 to 10
   * - ::

           real(kind=8), dimension(0:10) :: array_three
     - Array of **rank** 1 (1D), with **shape** 11x1 containing real values and **extent** of possible indices 0 to 10
   * - ::

           real(kind=8), allocatable, dimension(:,:) :: array_three
     - Array of **rank** 2 (2D), with **shape determined elsewhere** containing real values and **extent** of possible indices is 1 (default)
   * - ::

           real(kind=8), dimension(-1:1,3) :: array_four
     - Array of **rank** 2 (2D), with **shape** 3x2 containing real values and **extent** of possible indices is -1 to 1 and 1 to 3

Other attributes I haven't used yet include:

* ``pointer`` - behaves like ``allocatable``, except it uses existing memory not new memory
* ``target`` - the variables to which ``pointer`` is associated
* ``public`` - variables are used externally
* ``private`` - variables are used internally
* ``optional`` - used if the variable is optional on the argument list
* ``save`` - to preserve values between successive calls
* ``intrinsic`` - not sure what this means - **how could you declare an intrinsic variable?**

Substrings
----------

Let's say we had the string ``string = "12345"``

To extract elements, use:

::

    ! string(start : end)

where start and end are the indices of the start and end values.

.. list-table::
   :header-rows: 1
   :widths: 20 20

   * - Code
     - Meaning
   * - ::

           s1 = string(2:5)
     - 2345
   * - ::

           s2 = string(:3)
     - 123
   * - ::

           s3 = string(3:)
     - 345

Initialise a Range of Values to the Same Value
----------------------------------------------

Fortran allows us to **initialise a range of values** to the same value using the ``data`` statement:

::

    ! data nlist /clist/

* ``nlist`` is a list of names of variables, arrays, elements of arrays, substrings, implied do lists
* ``clist`` is either ``c`` or ``r*c`` where ``r`` is the number of values to initialise and ``c`` is the value to give 

.. list-table::
   :header-rows: 1
   :widths: 20 20

   * - Code
     - Meaning
   * - ::
           
           data i, j, k /3*0/
     - Initialise all the values to zero



Array Manipulation
==================

Array Initialisation and Allocation
-----------------------------------

In the Primative Types section, we declared arrays and sometimes gave them shape. We also need to be able to give arrays initial values and if the shape is deferred, specify the shape using ``allocate``.

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

           real(kind=8), dimension(1:10) :: a 
           a = 0.d0
     - **Real array** declaration and intialisation of an array to zero
   * - ::

           real(kind=8):: a(5,5), b(-4,16)
           a = 1.d0
           b = 0.d0
     - **Real array** declaration and intialisation of different shaped arrays in the same line
   * - ::

           a = real(10,15); data a/150*0.0/
     - **Real array** declaration and intialisation of 150 values to zero using Fortran 77 data statement (note semi-colon separates commands). This is a bit dodgy because the type of real values is not specified.
   * - ::

           a = (/100.0d0, 200.0d0, 300.0d0/)
     - **Array constructor**
   * - ::

           a = (/100.0d0, A(1:5, :), 300.0d0/)
     - **Array constructor** with 5 rows of A included in the order (in the same row as 100 and 300)
   * - ::

           a = (/ (j**3), j = 1, m /)
     - **Implied do array constructor** useful if the array contains a formula
   * - ::

           allocate(b(5,5), c(2,1), stat = allocate_status)
           if(allocate_status /= 0) STOP "***Not enough memory ***"
     - **Array allocation** where stat can be used to see if there is enough memory
   * - ::

           deallocate(a)
     - **Array deallocation** to free up memory if the array isn't needed
   * - ::

           x = 1/y + c(2:6,10)
     - **Array expression** (haven't used this one)   
   * - ::

           x = y + z
     - **Array addition** element-by-element addition (also applies to subtraction, multiplication, division)  

Array Extraction
----------------

Let's say we had the array

::

    real(kind=8), dimension(m, n) :: a

To extract elements, use 


::

    ! a(start : end [:stride])

where start and end are the indices of the start and end values and stride is optional, default is stride = 1.

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

           a(:, 2)
     - Read as "All rows, (in the) second column"
   * - ::

           a(m, :)
     - Read as "All columns, (in the) last row"
   * - ::

           a(:10, :10)
     - Leading 10 by 10 submatrix

Reference Section
=================

This section contains methods I haven't used, but they are useful in order to interpret the code of others.

Derived Data Types
------------------

I haven't used these, but, this is how a Dervied Data Type is used:

Definition
~~~~~~~~~~

::

   type person
       character(len=10) :: name
       integer :: age
   end type person

Instantiatation
~~~~~~~~~~~~~~~

::

    type(person) = person_one

Constructor
~~~~~~~~~~~

::
  
    person_one = person("andrew", 24)

How to Access Values
~~~~~~~~~~~~~~~~~~~~

::
 
    name = person_one%name  

Pointers
--------

I haven't used Pointers in Fortran, but this is a list of commands:

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Code
     - Meaning
   * - ::

           real, pointer :: p
     - Pointer declaration
   * - ::
 
           real, pointer :: a(:)
     - Array declaration (with deferred shape)
   * - ::

           real, target :: t
     - Defines the target
   * - ::

           p => t
     - Set pointer p to t
   * - ::

           associated(p, [target])
     - Pointer associated with target
   * - ::

           nullify(p)
     - Associate the pointer with NULL


Fortran 77
----------

To be able to read old Fortran, we might need these definitions. I am filling this table out to decode Ferzinger and Peric's notation (which uses old formatting)

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Data Type
     - Meaning
   * - ::

           integer integer*2 integer*4 integer*8
     - Integer with 2, 4 and 8 bytes (4 is default)
   * - ::
 
           real real*4 real*8 real*16
     - Real with 4, 8 and 16 bytes (4 is default)
   * - ::

           double precision
     - Double precision 8 bytes
   * - ::

           complex complex*8 complex*16 complex*32
     - Complex with 8, 16 and 32 bytes (8 is default)
   * - ::

           logical logical*1 logical*2 logical*4 logical*8
     - Boolean with 1, 2, 4 and 8 bytes (4 is default)
   * - ::

           character character*n
     - String with n bytes (1 is default)
   * - ::
 
           parameter(a, b)
     - List of typeless constants (bad for giving type)
   * - ::

           common a, b, c
     - List of global variables (bad for encapsulation)
   * - ::

           dimension a(n)
     - Array declaration
