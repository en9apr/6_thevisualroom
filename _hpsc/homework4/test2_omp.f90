! $MYHPSC/homework4/test2_omp.f90
! The purpose is to print a table for the effect of the number
! of integration points on the accuracy of the integration

! Example use: 
! $ gfortran -fopenmp quadrature_omp.f90 test2_omp.f90
! $ ./a.out

program test2_omp

	! To provide input values to print error_table to the screen
	use omp_lib
    use quadrature_omp, only: trapezoid, error_table

    implicit none
    real(kind=8) :: a,b,int_true,k
    integer :: nvals(20), i, n, nthreads
	real(kind=8) :: t1, t2, elapsed_time
    integer(kind=8) :: tclock1, tclock2, clock_rate
    a = 0.d0
    b = 2.d0
	k = 1000.d0
    int_true = (b-a) + (b**4 - a**4) / 4.d0 - &
			   ((cos(k*b) - cos(k*a)) / k)

	! Specify number of threads to use:
    nthreads = 1       ! need this value in serial mode
	!$ nthreads = 2    
	!$ call omp_set_num_threads(nthreads)
	!$ print "('Using OpenMP with ',i3,' threads')", nthreads

	! Specify number of threads to use:
	!$ call omp_set_num_threads(2)

!   values of n to test:
    do i=1,20
        nvals(i) = 5 * 2**(i-1)
    enddo

	call system_clock(tclock1)  ! start wall timer
    call cpu_time(t1)   ! start cpu timer

    call error_table(f, a, b, nvals, int_true)

	call cpu_time(t2)   ! end cpu timer

    print 10, t2-t1
10  format("CPU time = ",es13.3, " seconds")
    
    call system_clock(tclock2, clock_rate)
    elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
    print 11, elapsed_time
11  format("Elapsed time = ",es13.3, " seconds")

contains

	! To return the the value of the function to integrate, f

    real(kind=8) function f(x)
        implicit none
        real(kind=8), intent(in) :: x 
        
        f = 1.d0 + x**3 + sin(1000*x)
    end function f

end program test2_omp
