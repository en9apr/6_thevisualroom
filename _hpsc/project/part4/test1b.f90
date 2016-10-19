program test1b
! ---
! Uses random seed module
! ---
	use mc_walk, only: random_walk, n_walks, x
	use random_util, only: init_random_seed
	use problem_description, only: utrue, ax, ay, bx, by, nx, ny, dx, dy
	real(kind=8) :: ub, rel_error, x0, y0
	integer :: iabort, i, j, max_steps, seed1, n_success, n_mc
!	integer, parameter :: nx = 19
!	integer, parameter :: ny = 11
! ---
! Set the number of samples for the Monte Carlo Solution
! ---
	n_mc=10
! ---
! Number of steps is unknown, so we must allocate:
! ---
	max_steps = 20*max(nx, ny)*n_mc
!	max_steps = 10
	allocate(x(max_steps)) ! random number array 0 to 1 
! ---
! Generate random number array x (0 to 1)
! ---
	seed1 = 12345   	! seed1 = 12345 for repeatable seed
						! seed1 = 0 for random seed
	call init_random_seed(seed1)
	call random_number(x)
	! Try it out from a specific (x0,y0):
		x0 = 0.9
		y0 = 0.6
		
		i0 = nint((x0-ax)/dx)
		j0 = nint((y0-ay)/dy)

		! shift (x0,y0) to a grid point if it wasn't already:
		x0 = ax + i0*dx
		y0 = ay + j0*dy

		u_true = utrue(x0,y0)

		i = i0
		j = j0
! ---
! Record the number of times we successfully hit the boundary
! ---
	n_success = 0
! ---
! Record the number of random walks we make
! ---
	n_walks = 0
! ---
! Loop through and call a single random walk each time, indexing
! You must start the walk using the last successful walk point
! ---
	do k=1,n_mc 0
		call random_walk(i, j, max_steps, ub, iabort) 
		if(iabort .NE. 1) then
			ub_sum = ub_sum + ub
			n_success = n_success + 1
		endif
		u_mc = ub_sum/n_success
		print *, ""
		print *, "ub=", ub
	!	print *, "x=",x(k)
		print *, "iabort=", iabort
		print *, "nwalks=", n_walks
		print *, "n_success=", n_success
		print *, "u_mc=", u_mc
	enddo
end program test1b
