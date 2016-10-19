! $MYHPSC/homework3/problem7/test_quartic.f90

program test_quartic

    use newton, only: solve, tol
    use functions, only: f_quartic, fprime_quartic, alpha

    implicit none
    real(kind=8) :: x, x0, fx, dfx, xstar
    real(kind=8) :: tolvals(3), alphavals(3)
    integer :: iters, i, j
	logical :: debug         ! set to .true. or .false.

    ! values to test as x0:
    x0 = 4.0d0

    print 10, x0
10	format("Starting with initial guess x = ", es22.15)

    debug = .false.

	! values to use as tolerance:
	tolvals = (/1.0d-5, 1.0d-10, 1.0d-14/)

	! values to use as alpha:
	alphavals = (/1.0d-4, 1.0d-8, 1.0d-12/)

	print *, ' '  ! blank line
	print *, '   alpha        tol        iters          x              f(x)/f`(x)   x-xstar'

    do i=1,3
		do j =1,3
			tol = tolvals(j)
		    alpha = alphavals(i)
			xstar = 1+(alpha**0.25)
		    call solve(f_quartic, fprime_quartic, x0, x, iters, debug)
			fx = f_quartic(x)
			dfx = fprime_quartic(x)
		    print 11, alpha, tol, iters, x, fx/dfx, x-xstar
11 			format(2es13.3, i4, es24.15, 2es13.3)

        enddo
		print *, ' '  ! blank line
	enddo

end program test_quartic
