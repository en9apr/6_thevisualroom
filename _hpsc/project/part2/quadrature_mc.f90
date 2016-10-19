! ---
! This module performs Monte Carlo Integration by:
! 1) Generating random input over a number of points
! 2) Evaluating the function at these points for 20 dimensions
! 3) Summing the function and dividing by the number of points for the average
! 4) Multiplying the average by the base for the integral
! ---
module quadrature_mc
! ---
! Uses random seed and functions modules
! ---
	use random_util, only: init_random_seed
	use functions, only: g

contains
! ---
! Monte Carlo Integration function
! Example use: int_mc = quad_mc(g,a,b,ndim,npoints)
! Inputs: g (from functions module), a and b (vectors of length ndim)
!		  ndim (number of dimensions) npoints (number of points)
! Output: quad_mc (the value of the integral)	
! ---		
	function quad_mc(g, a, b, ndim, npoints)
	! ---
	! Declare parameters and arrays.
	! ---
		implicit none
		real(kind=8), external :: g
		real(kind=8), dimension(ndim), intent(in) :: a, b
		integer, intent(in) :: ndim, npoints
		integer :: seed1, i, j
		real(kind=8), dimension(:), allocatable :: f, xin, x
		real(kind=8) :: quad_mc 
	! ---
	! Number of points & dimensions is unknown, so must allocate:
	! ---
		allocate(x(npoints)) ! random number array 0 to 1  
		allocate(xin(ndim))  ! random number array a+(b-a)x
		allocate(f(npoints)) ! result of sum(x**2) over all dimensions
							 ! for n points
	! ---
	! Generate random number array x once per call of quad_mc
	! Length of the array is the number of random points 
	! ---
		seed1 = 12345   ! seed1 = 12345 for repeatable seed
						! seed1 = 0 for random seed
 		call init_random_seed(seed1)
		call random_number(x)
	! ---
	! Generate random input based on a and b and 
	! compute sum(x**2) for all dimensions
	! ---
		do i=1,npoints
			do j=1,ndim
				xin(j)=a(j)+(b(j)-a(j))*x(i)
			enddo
			f(i) = g(xin,ndim)
		enddo
	! ---
	! Return the value of the Monte Carlo Integral by multiplying the
	! average value by the base
	! ---
		quad_mc = product(b-a)*(sum(f)/npoints)

	end function quad_mc

end module quadrature_mc
