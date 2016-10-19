! ---
! This program loops through the Monte Carlo solution of Laplace's equation,
! doubling the number of sample points every loop using MPI
! ---
program laplace_mc
! ---
! ALL PROCESSES: Modules and functions
! ---
    use mpi
	use random_util, only: init_random_seed
	use mc_walk, only: random_walk, many_walks, n_walks, x, n_success_proc
	use problem_description, only: utrue, ax, ay, bx, by, nx, ny, dx, dy
! ---
! ALL PROCESSES: Parameters, vectors and arrays
! ---	
	implicit none
	real(kind=8) :: u_mc, x0, y0, rel_error, u_true
	integer :: n_mc, n_success, i0, j0, max_steps, i, &
	process_num, num_processes, ierr, &
	istep, seed1, n_success_total
	logical :: debug_many
! ---
! ALL PROCESSES: Number of steps is unknown, so we must allocate -
! 200000 is the maximum number of moves
! ---
	n_mc = 10		! This must be at least 3 as there are 3 processes
					! Use n_mc = 10 if not debugging random_walk
	max_steps = 200000*max(nx, ny)
	allocate(x(max_steps)) ! random number array 0 to 1 
! ---- 
! ALL PROCESSES: Initialize MPI
! ----
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_processes, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, process_num, ierr) 
! ---
! ALL PROCESSES: Generate random number array x (0 to 1)
! ---
	seed1 = 12345   ! seed1 = 12345 for repeatable seed
					! seed1 = 0 for random seed
	seed1 = seed1 + 97*process_num 	! unique for each process
	print '(" Seed for process",i2, " is ",i5)', process_num, seed1
	call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! wait for processes to print
	call init_random_seed(seed1)
	call random_number(x)
! ----
! ALL PROCESSES: Try it out from a specific (x0,y0):
! ----
	x0 = 0.9
	y0 = 0.6
	i0 = nint((x0-ax)/dx)
	j0 = nint((y0-ay)/dy)
! ----
! ALL PROCESSES: Shift (x0,y0) to a grid point if it wasn't already:
! ----
	x0 = ax + i0*dx
	y0 = ay + j0*dy
	u_true = utrue(x0,y0)
	n_walks = 0
! ---
! MASTER PROCESS: Open a file to write to, provide true solution and print output headers
! ---
	if (process_num .EQ. 0) then
		open(unit=25, file='laplace_mc.txt', status='unknown')
		print '(" True solution of PDE at x =",f10.4," y =",f10.4)', x0, y0
		print '(" u_true =",e23.15)', u_true
		print *, ""
		print *, "  Successes    Monte Carlo Solution   Relative Error"
	endif
! ---
! ALL PROCESSES: Loop through and call many_walks, doubling the number
! of Monte Carlo points by two each time
! ---
! MASTER PROCESS: Compute relative error and print output of many_walks  
! ---
	debug_many = .False.		! Best to remove the do loop if you're trying this
	do i=1,13
		call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success, debug_many)
		n_mc = 2.0d0 * n_mc
		if (process_num .EQ. 0) then
			rel_error = abs(u_true-u_mc)/u_true
			print '("  ",i10," ",e23.15," ",e23.15)', n_success, u_mc, rel_error
			write(25,'(i10,e23.15,e15.6)') n_success, u_mc, rel_error
		endif
	enddo
! ---
! ALL PROCESSES: Wait for process 0 to print
! ---
	call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
! ---
! MASTER PROCESS: Reduce n_success_proc to n_success_total  
! ---
	call MPI_REDUCE(n_success_proc, n_success_total, 1, &
		MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
! ---
! MASTER PROCESS: Reduce n_success_proc to n_success_total  
! ---
	if (process_num .EQ. 0) then
		print '("Final approximation to u(x0,y0):",e23.15)', u_mc
		print '("Total number of walks by all processes:",i10)', n_success_total
	endif
! ---
! ALL PROCESSES: Wait for process 0 to print
! ---
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! ---
! ALL PROCESSES: Print out process number and number of successes
! ---
	print '("Walks by process ",i1,":",i8)', process_num, n_success_proc
! ----
! ALL PROCESSES: COMPLETE MPI
! ----
    call MPI_FINALIZE(ierr)

end program laplace_mc
