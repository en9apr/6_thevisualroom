! $MYHPSC/homework4/test2.f90
! The purpose is to print a table for the effect of the number
! of integration points on the accuracy of the integration

! Example use: 
! $ gfortran quadrature.f90 test2.f90
! $ ./a.out

program test2

	! To provide input values to print error_table to the screen

    use quadrature, only: trapezoid, error_table

    implicit none
    real(kind=8) :: a,b,int_true,k
    integer :: nvals(12), i, n

    a = 0.d0
    b = 2.d0
	k = 1000.d0
    int_true = (b-a) + (b**4 - a**4) / 4.d0 - &
			   ((cos(k*b) - cos(k*a)) / k)

!   values of n to test:
    do i=1,12
        nvals(i) = 5 * 2**(i-1)
    enddo

    call error_table(f, a, b, nvals, int_true)

contains

	! To return the the value of the function to integrate, f

    real(kind=8) function f(x)
        implicit none
        real(kind=8), intent(in) :: x 
        
        f = 1.d0 + x**3 + sin(1000*x)
    end function f

end program test2
