! This program is designed to test communication using MPI, with a known
! solution to an integral. It must cope with:
! 1) No subdivisions specified - error (as master does nothing)
! 2) Subdivisions less than (the number of processes-1)
! 3) Subdivisions equal to (the number of processes-1)
! 4) Subdivisions greater than (the number of processes-1)

! Example Use: 
! $ make clean
! $ make test3

program test3

! ----
! MODULES AND FUNCTIONS
! ----
    use mpi
    use quadrature, only: trapezoid
    use functions, only: f, fevals_proc, k

! ----
! PARAMETERS, VECTORS AND ARRAYS
! ----
    implicit none
    real(kind=8) :: a, b, integral_true, integral_approx, dx_sub, integral_sub
    integer :: process_num, num_processes, ierr, n, fevals_total, num_sub, j, jj, num_sent, &
		next_sub, sender
    integer, dimension(MPI_STATUS_SIZE) :: status
	real(kind=8), allocatable, dimension(:,:) :: ab_master
	real(kind=8), allocatable, dimension(:) :: ab_sub, integral_master

! ---- 
! INITIALIZE MPI
! ----
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_processes, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, process_num, ierr)

! ----
! SET GLOBAL VARIABLES
! ----
    k = 1.d3   						! Functions module variable 
    a = 0.d0						! LHS of integral range
    b = 2.d0						! RHS of integral range
    integral_true = (b-a) + (b**4 - a**4) / 4.d0 - (1.d0/k) * (cos(k*b) - cos(k*a))   ! True integral
    n = 1000						! Each division uses 1000 points
 	fevals_proc = 0		   			! Each process keeps track of number of function evaluations

! ----
! INITIALIZATION IN THE MASTER PROCESS
! ----
	if (process_num .EQ. 0) then
	! ----
	! User must enter the number of subintervals and this cannot be zero:
	! ----
100		continue
	    print *, "Please enter the number of subintervals: "
        read *, num_sub

		if(num_sub .EQ. 0) then
			print *, "Error: number of subintervals cannot be zero"
			go to 100
		endif
	! ----
	! Now allocate size of the arrays using number of subintervals:
	! ----
        allocate(ab_master(2, num_sub))	! to hold a-b values to Workers
		ab_master = 1.d0
		allocate(integral_master(num_sub))	! to hold integral from Workers
	! ----
	! We must keep track of the number of times we send to the Workers:
	! ----
		num_sent = 0
	endif

! ----
! INITIALIZATION IN THE WORKER PROCESS(ES)
! ----
	! ----
	! Send Workers the number of subintervals:
	! ----
		call MPI_BCAST(num_sub, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
	! ----
	! Allocation of the a-b values vector:
	! ----
    	if (process_num .GT. 0) then
        	allocate(ab_sub(2))   ! to hold a column vector sent from master
    	endif

! ----
! MASTER PROCESS: MASTER-WORKER PARADIGM 
! ----
	if (process_num .EQ. 0) then
	! ----
	! Print the number of processes and the true integral
	! ----
        print '("Using ",i3," processes")', num_processes
        print '("true integral: ", es22.14)', integral_true
        print *, " "  ! blank line
	! ----
	! Master populates the matrix with values for all subintervals      		
	! ----
		dx_sub = (b-a) / num_sub

    	do j=1,num_sub
      		ab_master(1,j) = a + (j-1)*dx_sub
      		ab_master(2,j) = a + j*dx_sub
		enddo
	! ----
	! Send first set of a-b values to Workers - number sent can be less than the
	! number of Workers, so take the minimum, then record how many were sent 
	! ----
		do j=1,min(num_processes-1,num_sub)
      		call MPI_SEND(ab_master(1,j), 2, MPI_DOUBLE_PRECISION, &
				j, j, MPI_COMM_WORLD, ierr)
			num_sent = num_sent + 1
      	enddo
	! ----
	! Loop through the number of subintervals:	
	! ----
		do j=1,num_sub
		! ----
   		! Master recieves part of the integral from the Workers and stores it in
		! an array from any process:
		! ----
			call MPI_RECV(integral_sub, 1, MPI_DOUBLE_PRECISION, &
            	MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			sender = status(MPI_SOURCE)  		! Whichever process sent it is now free
			jj = status(MPI_TAG)		 		! The tag is j or next_sub
			integral_master(jj) = integral_sub	! Store the value in the matrix for summation later
		! ----
		! Test to see if there are more intervals to send -
		! i.e. if number sent is less than number of subintervals
		! ----
			if (num_sent .LT. num_sub) then
		        next_sub = num_sent + 1     ! The tag is the next subinterval
		        call MPI_SEND(ab_master(1,next_sub), 2, MPI_DOUBLE_PRECISION,&
		        	sender, next_sub, MPI_COMM_WORLD, ierr)
		        num_sent = num_sent + 1		! The number sent increases by 1
		! ----
		! If there are no more intervals to send, tag=0 and Workers return:
		! ----
          	else
		        call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, &
		        	sender, 0, MPI_COMM_WORLD, ierr)
          	endif
		enddo
	! ----
	! There are no more subintervals, so summate the parts of the integral and print:	
	! ----
		integral_approx = sum(abs(integral_master))
		print '("Trapezoid approximation with ",i8," total points: ",es22.14)',&
        num_sub*n, integral_approx
    endif

! ----
! WORKER PROCESS(ES): MASTER-WORKER PARADIGM 
! ----
	if(process_num .GT. 0) then
	! ----
	! If we have more processors than subintervals, some processors do no work:
	! ----
		if (process_num .GT. num_sub) go to 99
	! ----
	! Otherwise loop through until Master sends tag=0 (when we've sent all the intervals):
	! ----
		do while (.true.)
        ! ----
		! Workers recieve the a-b matrix as a vector, unpack it and send it to
		! the trapezoid function to evaluate the integral with 1000 divisions:
		! ----
			call MPI_RECV(ab_sub, 2, MPI_DOUBLE_PRECISION, &
				0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			integral_sub = trapezoid(f,ab_sub(1),ab_sub(2),n)
			j = status(MPI_TAG)   ! this is the interval number or the number of Workers
								  ! (whichever is smaller)
		! ----
		! If the Master sends tag=0, then the Workers return because there are no
		! more intervals to compute:
		! ----
			if (j .EQ. 0) go to 99
		! ----
		! Otherwise we must send the Master the computed integral and print it:
		! ----
		 	call MPI_SEND(integral_sub, 1, MPI_DOUBLE_PRECISION, &
		 		0, j, MPI_COMM_WORLD, ierr)
			print '("fevals by Process ",i2,": ",i13)', process_num, fevals_proc
		enddo
	endif

! ----
! COMPLETE MPI
! ----
99  continue   ! Might jump to here if number of Workers > number of subintervals
    call MPI_FINALIZE(ierr)

end program test3
