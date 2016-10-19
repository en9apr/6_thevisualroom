! $MYHPSC/homework3/intersections.f90

program intersections

    use newton, only: solve
    use functions, only: f_function, fprime_function

    implicit none
    real(kind=8) :: x, x0, fx
    real(kind=8) :: x0vals(4)
    integer :: iters, i
	logical :: debug         ! set to .true. or .false.

    print *, "Test routine for computing zero of f"
    debug = .true.

    ! values to test as x0:
    x0vals = (/-2.2d0, -1.6d0, -0.8d0, 1.4d0 /)

    do i=1,4
        x0 = x0vals(i)
		print *, ' '  ! blank line
        call solve(f_function, fprime_function, x0, x, iters, debug)

        print 11, x, iters
11      format('solver returns x = ', es22.15, ' after', i3, ' iterations')

        fx = f_function(x)
        print 12, fx
12      format('the value of f(x) is ', es22.15)

!        if (abs(x-2.d0) > 1d-14) then
!            print 13, x
!13          format('*** Unexpected result: x = ', es22.15)
!            endif
        enddo

end program intersections
