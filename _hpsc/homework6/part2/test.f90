
program test

    use mpi

    use quadrature, only: trapezoid
    use functions, only: f, fevals_proc, k

    implicit none
    real(kind=8) :: a,b,int_true, int_approx

    integer :: proc_num, num_procs, ierr, n, fevals_total
    integer, dimension(MPI_STATUS_SIZE) :: status

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)

    ! All processes set these values so we don't have to broadcast:
    k = 1.d3   ! functions module variable 
    a = 0.d0
    b = 2.d0
    int_true = (b-a) + (b**4 - a**4) / 4.d0 - (1.d0/k) * (cos(k*b) - cos(k*a))
    n = 1000

    ! Each process keeps track of number of fevals:
    fevals_proc = 0

    if (proc_num==0) then
        print '("Using ",i3," processes")', num_procs
        print '("true integral: ", es22.14)', int_true
        print *, " "  ! blank line
        endif

    call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! wait for process 0 to print

    ! Note: In this version all processes call trap and repeat the
    !       same work, so each should get the same answer.  
    int_approx = trapezoid(f,a,b,n)
    print '("Process ",i3," with n = ",i8," computes int_approx = ",es22.14)', &
            proc_num,n, int_approx

    call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! wait for all process to print

    ! print the number of function evaluations by each thread:
    print '("fevals by Process ",i2,": ",i13)',  proc_num, fevals_proc

    call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! wait for all process to print
	
	call MPI_REDUCE(fevals_proc, fevals_total, 1, &
			MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
		    MPI_COMM_WORLD, ierr)

    if (proc_num==0) then
        ! This is wrong -- part of homework is to fix this:
		
    !    fevals_total = 0   !! need to fix
        print '("Total number of fevals:  ",i10)', fevals_total
        endif

    call MPI_FINALIZE(ierr)

end program test
