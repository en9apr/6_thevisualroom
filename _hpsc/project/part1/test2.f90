
program test2

    use mpi

    use quadrature, only: trapezoid
    use functions, only: f, fevals_proc, k

    implicit none
    real(kind=8) :: a, b, int_true, int_approx, dx_sub, int_sub

    integer :: proc_num, num_procs, ierr, n, fevals_total, n_sub, j, jj
    integer, dimension(MPI_STATUS_SIZE) :: status
	real(kind=8), allocatable, dimension(:,:) :: ab_sub
	real(kind=8), allocatable, dimension(:) :: sub, int_array
 
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)

    ! All processes set these values so we don't have to broadcast:
    k = 1.d3   ! functions module variable 
    a = 0.d0
    b = 2.d0
    int_true = (b-a) + (b**4 - a**4) / 4.d0 - (1.d0/k) * (cos(k*b) - cos(k*a))
    n = 1000
	
	n_sub = num_procs-1 ! Set the number of Worker Processes

    ! Each process keeps track of number of fevals:
    fevals_proc = 0

    if (proc_num .GT. 0) then
        allocate(sub(2))   ! to hold a column vector sent from master
    endif 

    if (proc_num .EQ. 0) then
        print '("Using ",i3," processes")', num_procs
        print '("true integral: ", es22.14)', int_true
        print *, " "  ! blank line

		allocate(ab_sub(2,n_sub))    ! to hold a-b values for Workers
		allocate(int_array(n_sub))	 ! to hold integral from Workers

		ab_sub = 1.d0

		! Process 0 (Master) sends other processes (Workers) left and right edges
		! of the integral      		

		dx_sub = (b-a) / n_sub

    	do j=1,n_sub
      		ab_sub(1,j) = a + (j-1)*dx_sub
      		ab_sub(2,j) = a + j*dx_sub
      		call MPI_SEND(ab_sub(1,j), 2, MPI_DOUBLE_PRECISION, j, j, &
                    MPI_COMM_WORLD, ierr)
      	enddo

   		! Master recieves part of the integral from the Workers and stores it in
		! an array from any process

		do j=1,n_sub
			call MPI_RECV(int_sub, 1, MPI_DOUBLE_PRECISION, &
                    MPI_ANY_SOURCE, MPI_ANY_TAG, &
                    MPI_COMM_WORLD, status, ierr)
			jj = status(MPI_TAG)
			int_array(jj) = int_sub       	
		enddo

		int_approx = sum(abs(int_array))

    endif

	! Process 0 (Master) sends other processes (Workers) left and right edges
	! of the integral 

	if(proc_num .GT. 0) then

		call MPI_RECV(sub, 2, MPI_DOUBLE_PRECISION, &
		              0, MPI_ANY_TAG, &
		              MPI_COMM_WORLD, status, ierr)

		int_sub = trapezoid(f,sub(1),sub(2),n)

		j = status(MPI_TAG)   ! this is the row number
		                      ! (should agree with proc_num)

	 	call MPI_SEND(int_sub, 1, MPI_DOUBLE_PRECISION, &
	 	               0, j, MPI_COMM_WORLD, ierr)

		! print the number of function evaluations by each thread:
		print '("fevals by Process ",i2,": ",i13)', proc_num, fevals_proc
	
	endif

	! Master Process prints out the result:
	call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! wait for all processes to print

    if (proc_num .EQ. 0) then
   		print '("Trapezoid approximation with ",i8," total points: ",es22.14)',&
        n_sub*n, int_approx
    endif

    call MPI_FINALIZE(ierr)

end program test2
