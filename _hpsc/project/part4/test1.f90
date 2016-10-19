program test1
! ---
! Uses random seed module
! ---
	use mc_walk, only: random_walk, n_walks, x
	use random_util, only: init_random_seed
	use problem_description, only: utrue, ax, ay, bx, by, nx, ny, dx, dy
		real(kind=8) :: ub, rel_error, x0, y0
		integer :: iabort, i, j, max_steps, seed1, n_success, i0, j0

	! ---
	! Number of steps is unknown, so we must allocate:
	! ---
		max_steps = 100*max(nx, ny)
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

		i = i0
		j = j0

    	n_success = 0
		print *, ""
		print '(" True solution of PDE at x =",f10.4," y =",f10.4)', x0, y0
		print '(" i =",i8," j =",i8)', i, j
		print '(" u_true =   ",es13.3)', u_true
		call random_walk(i, j, max_steps, ub, iabort) 
		print *, "ub=", ub
		print *, "x=",x(1)
		print *, "iabort=", iabort
		print *, "nwalks=", n_walks

end program test1
