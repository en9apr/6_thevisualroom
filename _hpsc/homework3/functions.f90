! $MYHPSC/homework3/functions.f90

module functions
	
	! parameters for module
!	implicit none
!	real(kind=8) :: pi
!	save
	

contains

real(kind=8) function f_sqrt(x)
    implicit none
    real(kind=8), intent(in) :: x

    f_sqrt = x**2 - 4.d0

end function f_sqrt


real(kind=8) function fprime_sqrt(x)
    implicit none
    real(kind=8), intent(in) :: x
    
    fprime_sqrt = 2.d0 * x

end function fprime_sqrt

real(kind=8) function f_function(x)
    implicit none
    real(kind=8), intent(in) :: x
	real(kind=8) :: pi = 3.141592653589793

    f_function = (x*cos(pi*x))-(1.0d0-(0.6d0*(x**2)))

end function f_function

	
real(kind=8) function fprime_function(x)
    implicit none
    real(kind=8), intent(in) :: x
    real(kind=8) :: pi = 3.141592653589793

    fprime_function = ((-x*pi)*sin(pi*x))+cos(pi*x)+(1.2d0*x)

end function fprime_function


end module functions
