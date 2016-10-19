! ---
! This module performs a random walk and also many random walks
! ---
module mc_walk
! ---
! ALL PROCESSES: Use random seed module & mpi
! ---
 	use mpi
	use problem_description, only: uboundary, ax, ay, dx, dy, nx, ny
! ---
! ALL PROCESSES: nwalks & n_success_proc are module variables
! Example use: use mc_walk, only nwalks, n_success_proc
! ---
    implicit none
    integer :: n_walks, n_success_proc
	real(kind=8), dimension(:), allocatable :: x
    save

contains
! ---
! ALL PROCESSES: A single random walk
! Inputs: i0 (step in x), j0 (step in y), max_steps (maximum number of steps)
! Outputs: ub (value at boundary), iabort (0 for success, 1 for failure)
! Example use: random_walk(0, 0, 100, ub, iabort) 
!     
!    i:  0,  1,.. nx, nx+1
! j:
! 0,   (BC) (BC) (BC) (BC)
! 1,   (BC)           (BC)
! ...  (BC)           (BC)
! nj   (BC)           (BC)
! nj+1 (BC) (BC) (BC) (BC)
! ---
	subroutine random_walk(i0, j0, max_steps, ub, iabort, debug_random)
	! ---
	! ALL PROCESSES: Declare parameters
	! ---
		implicit none
		integer, intent(in) :: i0, j0, max_steps
		logical, intent(in) :: debug_random
		real(kind=8), intent(out) :: ub
		integer, intent(out) :: iabort
		integer :: i, j, istep, seed1, start
		real(kind=8) :: xb, yb
	! ---
	! ALL PROCESSES: Initialize the walk from it's starting point
	! ---
		i = i0
		j = j0
		start = n_walks
	! ---
	! ALL PROCESSES: Go through the vector of random numbers and
	! make moves depending on the value of x
	! ---
		do istep=start,max_steps
			if (x(istep) .LT. 0.25) then
				i = i-1    ! go left
				if (debug_random) then
					print *, "go left"
					print *, "i=", i
					print *, "istep=", istep
					print *, "x=", x(istep)
				endif
				n_walks = n_walks + 1
			else if (x(istep) .LT. 0.5) then
				i = i+1    ! go right
				if (debug_random) then
					print *, "go right"
					print *, "i=", i
					print *, "istep=", istep
					print *, "x=", x(istep)
				endif
				n_walks = n_walks + 1
			else if (x(istep) .LT. 0.75) then
				j = j-1    ! go down
				if (debug_random) then
					print *, "go up"
					print *, "j=", j
					print *, "istep=", istep
					print *, "x=", x(istep)
				endif
				n_walks = n_walks + 1
			else
				j = j+1    ! go up
				if (debug_random) then
					print *, "go down"
					print *, "j=", j
					print *, "istep=", istep
					print *, "x=", x(istep)
				endif
				n_walks = n_walks + 1
			endif
		! ---
		! ALL PROCESSES: What if we are at the boundary?
		! ---			
			if (i*j*(nx+1 - i)*(ny+1 - j) .EQ. 0) then
				xb = ax + i*dx 			! x is ax or ax+i*dx
				yb = ay + j*dy			! y is ay or ay+j*dy
				ub = uboundary(xb, yb)	! the boundary value is then computed
				iabort = 0				! success
				n_success_proc = n_success_proc + 1
				if (debug_random) then
					print *, "Reached the boundary"
				endif				
				go to 99				! end do loop
			endif
		! ---
		! ALL PROCESSES: What if we never reach the boundary within max steps?		
		! ---
			if (istep .EQ. max_steps) then 
				ub = 0					! set ub to zero
				iabort = 1				! failure
				n_success_proc = n_success_proc
				if (debug_random) then
					print *, "Failed to reach boundary"
				endif	
			endif
		enddo

99  continue

	end subroutine random_walk 
! ---
! ALL PROCESSES: A single random walk
! Inputs: i0 (step in x), j0 (step in y), max_steps (maximum number of steps)
! Outputs: ub (value at boundary), iabort (0 for success, 1 for failure)
! Example use: many_walks(0, 0, 100, ub, iabort)
! ---
	subroutine many_walks(i0, j0, max_steps, n_mc, u_mc, n_success, debug_many)
	! ---
	! ALL PROCESSES: Declare parameters
	! ---
		implicit none
		integer, intent(in) :: i0, j0, max_steps, n_mc
		real(kind=8), intent(out) :: u_mc
		integer, intent(out) :: n_success
		logical, intent(in) :: debug_many
		real(kind=8) :: ub_sum, ub, answer_r
		real(kind=8), allocatable, dimension(:) :: int_array
		integer :: i, j, k, iabort, num_processes, process_num, &
		ierr, m, r, sender, iabort_r, answer, num_sent, next_walk, jj
		integer, dimension(MPI_STATUS_SIZE) :: status
		logical :: debug_random
    	call MPI_COMM_SIZE(MPI_COMM_WORLD, num_processes, ierr)
    	call MPI_COMM_RANK(MPI_COMM_WORLD, process_num, ierr)

		ub_sum = 0
		n_success = 0
		num_sent = 0
		n_success_proc = 0
	! ---
	! MASTER PROCESS
	! ---
		if(process_num .EQ. 0) then
			allocate(int_array(n_mc))	 ! to hold answer from Workers
 			! ---
			! INITIALIZE THE PROCESS
			! Send: 0, Tag: 1, To: Processes 1 to 3
			! Increment number sent by 1
 			! ---
				do m=1,num_processes-1
					call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, &
				    	m, 1, MPI_COMM_WORLD, ierr)
					num_sent = num_sent + 1
					if (debug_many) then
						print '("Process ",i1," sent tag = ",i1, " to process ",i1)',process_num,1,m
					endif				
				enddo
				if (debug_many) then
					print '("Number of values sent initially ",i1)',num_sent
				endif			
			! ---
			! RECIEVE ALL POSSIBLE CALLS TO RANDOM WALK
			! Receive: answer_r, Tag: iabort (0 or 1), From: Any Process (1 to 3)
			! Sender is now free
			! ---
				do m=1,n_mc
					call MPI_RECV(answer_r, 1, MPI_DOUBLE_PRECISION, &
				        MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
					iabort_r = status(MPI_TAG)
					sender = status(MPI_SOURCE)    ! Whichever process sent it is now free
				! ---
				! HAVE WE REACHED THE BOUNDARY?
				! If the walk reached the boundary, iabort = 0, so add it to the array
				! ---
					if(iabort_r .EQ. 0) then
						int_array(m) = answer_r
						n_success = n_success + 1
					else
						int_array(m) = 0
						n_success = n_success
					endif
					if (debug_many) then
						print '("Process ",i1," recieved answer = ",es22.14," tag = ",i3, " from process ",i1)', &
						process_num,int_array(m),iabort_r,sender
					endif
				! ---
				! SEND ADDITIONAL CALLS TO RANDOM WALK 
				! If the number sent is less than the number of Monte Carlo points
				! Send: 0, Tag: 1, To: Processes 1 to 3 (whichever is free)
				! Increment number sent by 1
				! ---
					if (num_sent .LT. n_mc) then
						call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, &
				    	sender, 1, MPI_COMM_WORLD, ierr)
						num_sent = num_sent + 1
						if (debug_many) then
							print '("Process ",i1," sent tag = ",i1, " to process ",i1)', &
							process_num,1,sender
						endif
				! ---
				! END WORKER PROCESSES
				! Send: 0, Tag: 0, To: Processes 1 to 3 (whichever is free)
				! ---
				  	else
						call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, &
						sender, 0, MPI_COMM_WORLD, ierr)
				  	endif
				enddo
			! ---
			! SUM OVER THE ARRAY TO RETURN U_MC TO MANY_WALKS
			! ---
				ub_sum=sum(int_array) 		! Summation of the values from each process
				u_mc = ub_sum/n_success		! Returned to many_walks
				if (debug_many) then
			 		print '("Returned answer ",es22.14)', u_mc
		 			print '("Number of values sent finally ",i10)', num_sent
				endif
 		endif

	! ---
	! WORKER PROCESS
	! ---
 		if(process_num .GT. 0) then

			do while (.true.)
				! ---
				! RECIEVE CALLS FROM MASTER PROCESS
				! Receive: 0, Tag: 1, From: Master (0)
				! Sender is now free
				! ---
					call MPI_RECV(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, &
				        MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
					r = status(MPI_TAG)
					if (debug_many) then
				 		print '("Process ",i1," recieved tag = ",i1, " from process ",i1)',process_num,r,0
					endif
				! ---
				! THERE IS ANOTHER WALK TO DO
				! Call Random Walk again
				! ---
					if(r .EQ. 1) then
						i = i0
						j = j0
						debug_random = .False.	! Change n_mc to 3 if debugging this
						call random_walk(i, j, max_steps, ub, iabort, debug_random)
					! ---
					! SEND THE ANSWER TO THE MASTER
					! Send: ub, Tag: iabort(0 or 1), To: Master (Process 0)
					! ---
						call MPI_SEND(ub, 1, MPI_DOUBLE_PRECISION, &
							0, iabort, MPI_COMM_WORLD, ierr)
						if (debug_many) then
							print '("Process ",i1," sent answer = ",es22.14, " & tag = ",i3, " to process ",i1)', &
						 	process_num,ub,iabort,0
						endif
				! ---
				! THERE ISN'T ANOTHER WALK TO DO
				! Worker finishes
				! ---
					else if(r .EQ. 0) then
						if (debug_many) then
							print '("Process ",i1," ended ")', process_num
						endif
						go to 100
					endif
			enddo
 		endif
! ----
! COMPLETE MPI
! ----
100 	continue ! Jumps here if the Worker finishes
 
	end subroutine many_walks

end module mc_walk

