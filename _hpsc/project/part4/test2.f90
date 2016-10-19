program test2
! ---
! Uses random seed module
! ---
	use mc_walk, only: random_walk, many_walks, n_walks, x
	use random_util, only: init_random_seed
	use problem_description, only: utrue, ax, ay, bx, by, nx, ny, dx, dy

		real(kind=8) :: u_mc, rel_error, x0, y0
		integer :: n_mc, n_success, i, j, max_steps, seed1, k, i0, j0
	! ---
	! Number of steps is unknown, so we must allocate
	! 200000 is the maximum number of moves
	! ---
		n_mc = 10
		max_steps = 200000*max(nx, ny)
		!max_steps = 25
		allocate(x(max_steps)) ! random number array 0 to 1 
	! Try it out from a specific (x0,y0):
		x0 = 0.9
		y0 = 0.6
		
		i0 = nint((x0-ax)/dx)
		j0 = nint((y0-ay)/dy)

		! shift (x0,y0) to a grid point if it wasn't already:
		x0 = ax + i0*dx
		y0 = ay + j0*dy

		u_true = utrue(x0,y0)

	! ---
	! Generate random number array x (0 to 1)
	! ---
		seed1 = 12345   	! seed1 = 12345 for repeatable seed
							! seed1 = 0 for random seed
		call init_random_seed(seed1)
		call random_number(x)
	! ---
	! Set the point of interest
	! ---
		i = i0
		j = j0
	! ---
	! Record the number of random walks we make
	! ---
		n_walks = 0
		print *, ""
		print '(" True solution of PDE at x =",f10.4," y =",f10.4)', x0, y0
		print '(" u_true =   ",es13.3)', u_true
		print *, ""
		print *, "    Successes  Walks    Monte Carlo Solution    Relative Error"
		do k=1,13
			call many_walks(i, j, max_steps, n_mc, u_mc, n_success)
			n_mc = 2 * n_mc
			rel_error = abs(u_true-u_mc)/u_true
			print '("  ",i8," ",i8,"      ",es13.3,"        ",es13.3)', &
			n_success, n_walks, u_mc, rel_error
		enddo

end program test2
