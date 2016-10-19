program test

	use quadrature_mc, only: quad_mc
	use functions, only: g

	implicit none
    integer :: seed1, npoints, i
	real(kind=8) :: int_mc
	integer, parameter :: ndim = 20
    real(kind=8), dimension(ndim) :: a, b

 	npoints = 10

 	do i=1,ndim
        a(i) = 2.d0
        b(i) = 4.d0
    enddo

    int_mc = quad_mc(g,a,b,ndim,npoints)
 	print *, "int_mc from random number generator:", int_mc 

end program test
