! ---
! This module returns the value of a function:
! sum(x**2) over all dimensions
! it also sets the number of evaulations of this function
! ---
module functions
! ---
! gevals is a module variable
! Example use: use functions, only gevals
! ---
    implicit none
    integer :: gevals
    save

contains
! ---
! Multiple dimension function
! Example use: f = g(xin,ndim)
! Inputs: x (vector of length ndim), ndim (number of dimensions)
! Output: g (the value of the function over all dimensions)	
! ---		
    function g(x,ndim)
	! ---
	! Declare parameters and vectors
	! ---
		implicit none
		integer, intent(in) :: ndim
		real(kind=8), intent(in) :: x(ndim)
		real(kind=8) :: g
		integer :: i
	! ---
	! We need to initialise the number of evaluations to zero
	! ---
		g = 0.d0
	! ---
	! Evaluate function in 20 dimensions
	! ---
		do i=1,ndim
		    g = g + x(i)**2
		enddo
	! ---
	! Update the number of function evaluations
	! ---
		gevals = gevals + 1

	end function g

end module functions
