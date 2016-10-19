! ---
! This program computes the Monte Carlo Integral by doubling the number
! of points every time
! ---
program test_quad_mc
! ---
! Uses functions, quadrature_mc and random_util modules
! ---
	use functions, only: g, gevals
    use quadrature_mc, only: quad_mc
    use random_util, only: init_random_seed
! ---
! Declare the parameters
! ---
    implicit none
    integer, parameter :: ndim = 20
    real(kind=8), dimension(ndim) :: a, b
    real(kind=8) :: volume, int_true, int_mc, int_mc_total, error
    integer :: i, seed1, n_total, npoints
! ---
! Need to initialize number of g evaulations to zero, because it's a module variable
! ---
    gevals = 0
! ---
! Open a file to write to
! ---
    open(unit=25, file='mc_quad_error.txt', status='unknown')
! ---
! a and b are the vectors for the integral limits in all dimensions
! and are the same in all directions
! ---
    do i=1,ndim
        a(i) = 2.d0
        b(i) = 4.d0
    enddo
! ---
! Compute the true integral for special case where
! g(x) = sum(x**2) over all dimensions
! True integral = base * height, 
! where height is the sum of the average value of the true integral in all dimensions
! ---
    volume = product(b-a)  ! =  product of b(i)-a(i) of ndim dimensions
    int_true = (volume) * sum((b**3 - a**3) / (3.d0*(b-a)))

    print '("Testing Monte Carlo quadrature in ",i2," dimensions")', ndim
    print '("True integral: ", es22.14)', int_true

! ---
! Start with Monte Carlo using only a few points.
! ---    
	npoints = 10
! ---
! Loop to successively and double the number of points used:
! ---
	do i=1,17
	! ---
    ! Compute integral using Monte Carlo method for a given set of points
	! compute relative error and write these to file with number of points
	! ---     
   		int_mc = quad_mc(g,a,b,ndim,npoints)
        error = abs(int_mc - int_true) / abs(int_true)
		write(25,'(i10,e23.15,e15.6)') npoints, int_mc, error
	! ---
	! Double the number of points
	! ---
        npoints = 2*npoints
    enddo

    print '("Final approximation to integral: ",es22.14)',int_mc
    print *, "Total g evaluations: ",gevals

end program test_quad_mc

