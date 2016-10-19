============
Input/Output
============

.. contents::
   :local:

.. highlight:: fortran

Print Statements
================

The print statement is used to send output to monitor - can use single or double quotes.

* The format statement of ``es22.15`` is typical for double precision, which is usually to 15 digits of precision. We have 22 character positions (including dot, exponent and 2 signs), so the largest number that could be printed is (which may depend on what we expect) - but this is nothing to do with machine precision:

::

     most positive = +9.999999999999999E+99
     smallest = +9.999999999999999E-99
     most negative = -9.999999999999999E+99

* Choosing the format of integers may depend on what we expect - so if the maximum number of iterations is 50, we might choose ``i2`` which would cover this range - again this is nothing to do with machine precision of integers:

::

     largest = 99
     smallest = 0

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning
   * - ::

           print '(i10, es22.15)', 2, 1.0d10
           print "(i10, es22.15)", 2, 1.0d10

     - **Print to screen with format**: integer with 10 digits, scientific number with 22 digits and 15 decimal places. Syntax is ``print format, list`` 
   * - ::

           print 10, x, iterations
           10 format("solver returns x=", es22.15," after", i3, " iterations")

           print 11, y, iterations
           11 format('solver returns y=', es22.15,' after', i3, ' iterations')

     - **Print to screen with format** (for extra long lines)
 
   * - ::

           print *, "Velocity v=", v, "m/s"
           print *, 'Velocity u=', v, 'm/s'

     - **List-directed I/O** The asterisk (*) means print in a format compatible with the items in the list. Either double or single quotes may be used. Using undefined format will leave a blank space at the start (because in old Fortran this was column control). It's probably always best to give a format anyway. The default print width for integers is 12 (with i) and it's 25 for doubles (with es).
   * - ::

           print *, ""
           
     - **Print a Blank Line**
   * - ::

           print "(a)", "x = "
           do i = 1, nx, 1
	       print "(11g12.4)", x(i, :)
           end do

     - **Print an array** - 11 columns wide, nx rows long

Opening, Closing, Reading from and Writing to Files
===================================================

* In Fortran each file is associated with a **unit number**, an integer between 1 and 99.
* Some unit numbers are reserved (I usually start from 10): 

   - 5 is keyboard input
   - 6 is keyboard output
 
* ``iostat`` will be non-zero if there is an error - can be used to detect the type of error
* ``err`` then specifies the line number with which to proceed should an error be encountered, perhaps to print to the screen. Error lables can start at 1.

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning
   * - ::

          write(10, '(i10, 2e22.15)', iostat=write_error, err=10) n, velocity, pressure
          
     - **Write list to unit** ``write(unit, format, specifiers) list``
   * - ::

          read(10, '(i10, 2e22.15)', iostat=read_error, err=20) n, velocity, pressure

     - **Read list from unit** ``read(unit, format, specifiers) list``
   * - ::

          open(10, file='data.txt', iostat=open_error, status='unknown', err=30)

     - **Open file** ``open(unit, specifiers)``
   * - ::

          close(10, iostat=close_error, err=40)

     - **Close file** ``close(unit, specifiers)``


Format of the File
==================

Take the following read statement:

::

   read(10, '(i10, 2e22.15)', iostat=read_error, err=20) n, velocity, pressure

The format of the corresponding file must match with comma separated values, i.e.:

::

   1, 200.0, 350.0


Format Statements (print, read, write)
======================================

The key to format statements:

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning

   * - ::

           w

     - Full length
   * - ::

           m

     - Minimum digits
   * - ::

           d

     - Decimal places
   * - ::

           e

     - Exponent length
   * - ::

           n

     - Positions to skip
   * - ::

           c

     - Positions to move
   * - ::

           r

     - Repetitions


.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning

   * - ::

           format = "(F10.3, A, ES14.7)"

     - Format string
   * - ::

           Iw Iw.m

     - Integer form
   * - ::

           Bw.m Ow.m Zw.m

     - Binary, octal, hex integer form
   * - ::

           Fw.d

     - Decimal form real format
   * - ::

           Ew.d 

     - Exponential form (0.12..E-11)
   * - ::

           Ew.dEe
     - Specified exponent length
   * - ::

           ESw.d ESw.dEe

     - Scientific form (1.2...E-10)
   * - ::

           ENw.d ENw.dEe

     - Engineering form (123.4...E-12)
   * - ::

           Gw.d

     - Generalized form
   * - ::

           Gw.dEe

     - Generalized exponent form
   * - ::

           Lw

     - Logical format (T, F)
   * - ::

           A Aw

     - Characters format
   * - ::

           nX
     - Horizontal positioning (skip)
   * - ::

           Tc TLc TRc
     - Move (absolute, left, right)
   * - ::

           r/
     - Vert. positioning (skip lines)
   * - ::

           r(...)

     - Grouping / repetition
   * - ::

           :

     - Format scanning control
   * - ::

           S SP SS

     - Sign control
   * - ::

           BN BZ

     - Blank control (blanks as zeros)

Data-Transfer Specifiers (Read, Write)
======================================

Important specifiers are in **bold**, typical example:

::

   write(10, '(i10, e22.15, e15.6)', iostat=write_error, err=10) n, velocity, pressure

   read(10, '(i10, e22.15, e15.6)', iostat=read_error, err=20) n, velocity, pressure

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning
   * - ::

           iostat=integer_variable

     - **Save iocode (error) to variable**
   * - ::

           err=label

     - **Label to jump to on error**
   * - ::

           advance='yes' 
           advance='no'

     - (Non-)advancing data transfer
   * - ::

           end=label

     - Label to jump to on end of file
   * - ::

           eor=label

     - Label for end of record
   * - ::

           rec=integer_variable


     - Record number to read or write
   * - ::

           size=integer_variable

     - Number of characters read


I/O Specifiers
==============

Open Specifiers
---------------

Important specifiers are in **bold**, typical example:

::

   open(10, file='data.txt', iostat=open_error, status='unknown', err=30)


.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning
   * - ::

           file='filename'

     - **Name of file to open**, e.g. ``file = 'data.txt'`` (if it's in the same directory as the complied code) 
   * - ::

           iostat=integer_variable

     - **Save iocode (error) to variable**, e.g. ``iostat = open_error``
   * - ::

           status='old' 
           status='new' 
           status='replace'
           status='scratch' 
           status='unknown'

     - **Status of input file**, e.g. ``status='unknown'``
   * - ::

           err=label

     - **Label to jump to on error**
   * - ::

           access='sequential' 
           access='direct'

     - Access method
   * - ::

           form='formatted' 
           form='unformatted'

     - Formatted/unformatted I/O
   * - ::

           recl=integer_variable

     - Length of record
   * - ::

           blank='null' 
           blank='zero'

     - Ignore blanks/treat them as 0 
   * - ::

           access='sequential' 
           access='direct' 

     - Access method
   * - ::

           position='asis'
           position='rewind'
           position='append'
           

     - Position, if sequential I/O
   * - ::

           action='read' 
           action='write'
           action='readwrite'

     - Read/write mode
   * - ::

           
           delim='quote' 
           delim='apostrophe'    
           delim='none'

     - Delimiter for char constants
   * - ::

           pad='yes' 
           pad='no'
           
     - Pad with blanks


Close Specifiers
----------------

Important specifiers are in **bold**, typical example (to close the file we opened, with unit ``10``):

::

   close(10, iostat=close_error, err=50)


.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning
   * - ::

           iostat

     - **Save iocode (error) to variable**
   * - ::

           err

     - **Label to jump to on error**
   * - ::

           status='keep' 
           status='delete'

     - Status of closed file

Reference List
==============

This contains commands I haven't used, but maybe useful for translation/understanding.

Getarg, Inquire, Backspace, Endfile, Rewind Functions
-----------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning
   * - ::

           call getarg(2, var)

     - Put 2nd CLI-argument in var
   * - ::

          inquire(unit, spec)

     - Inquiry by unit
   * - ::

          inquire(file=filename, spec)

     - Inquiry by filename
   * - ::

          inquire(iolength=iol) outlist

     - Inquiry by output item list
   * - ::

          backspace(unit, spec)

     - Go back one record
   * - ::

          endfile(unit, spec)

     - Write eof record
   * - ::

          rewind(unit, spec)

     - Jump to beginning of file

Inquire Specifiers
------------------

I haven't used these

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning
   * - ::

           access

     - 
   * - ::

           action

     - 
   * - ::

           blank

     - 
   * - ::

           delim

     - 
   * - ::

           direct

     - 
   * - ::

           exist

     - 
   * - ::

           form

     - 
   * - ::

           formatted

     - 
   * - ::

           iostat

     - 
   * - ::

           name

     - 
   * - ::

           named

     - 
   * - ::

           nextrec

     - 
   * - ::

           number

     - 
   * - ::

           opened

     - 
   * - ::

           pad

     - 
   * - ::

           position

     - 
   * - ::

           read

     - 
   * - ::

           readwrite

     - 
   * - ::

           recl

     - 
   * - ::

           sequential

     - 
   * - ::

           unformatted

     - 
   * - ::

           write

     - 
   * - ::

           iolength

     - 

Backspace-, endfile-, rewind-specifiers
---------------------------------------

I haven't used these

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Code
     - Meaning
   * - ::

           iostat

     - Save iocode (error) to variable
   * - ::

           err

     - Label to jump to on error





























