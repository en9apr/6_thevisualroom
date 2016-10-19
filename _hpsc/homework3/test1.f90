! $MYHPSC/homework3/test1.f90

program test1

    use newton, only: solve
    use functions, only: f_sqrt, fprime_sqrt

    implicit none
    real(kind=8) :: x, x0, fx
    real(kind=8) :: x0vals(3)
    integer :: iters, i
	logical :: debug         ! set to .true. or .false.

    print *, "Test routine for computing zero of f"
    debug = .true.

    ! values to test as x0:
    x0vals = (/1.d0, 2.d0, 100.d0 /)

    do i=1,3
        x0 = x0vals(i)
		print *, ' '  ! blank line
        call solve(f_sqrt, fprime_sqrt, x0, x, iters, debug)

        print 11, x, iters
11      format('solver returns x = ', es22.15, ' after', i3, ' iterations')

        fx = f_sqrt(x)
        print 12, fx
12      format('the value of f(x) is ', es22.15)

        if (abs(x-2.d0) > 1d-14) then
            print 13, x
13          format('*** Unexpected result: x = ', es22.15)
            endif
        enddo

end program test1
