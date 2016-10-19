! $MYHPSC/homework4/test1.f90
! The purpose is to print a table for the effect of the number
! of integration points on the accuracy of the integration

! Example use: 
! $ gfortran quadrature.f90 test1.f90
! $ ./a.out

program test1

	! To provide input values to print error_table to the screen

    use quadrature, only: trapezoid, error_table

    implicit none
    real(kind=8) :: a,b,int_true
    integer :: nvals(7), i, n

    a = 0.d0
    b = 2.d0
    int_true = (b-a) + (b**4 - a**4) / 4.d0

!    print 10, int_true
! 10 format("true integral: ", es22.14)
!    print *, " "  ! blank line

!	n=5

!	print 11, trapezoid(f,a,b,n)
!11  format("calculated integral: ", es22.14)
!   values of n to test:
    do i=1,7
        nvals(i) = 5 * 2**(i-1)
    enddo

    call error_table(f, a, b, nvals, int_true)

contains

	! To return the the value of the function to integrate, f

    real(kind=8) function f(x)
        implicit none
        real(kind=8), intent(in) :: x 
        
        f = 1.d0 + x**3
    end function f

end program test1
