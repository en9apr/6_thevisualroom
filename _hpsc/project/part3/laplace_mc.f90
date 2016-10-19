program laplace_mc
! ---
! Uses random seed module
! ---
	use mc_walk, only: random_walk, many_walks, nwalks
	use problem_description, only: utrue, ax, ay, bx, by, nx, ny, dx, dy

		real(kind=8) :: u_mc, x0, y0, rel_error
		integer :: n_mc, n_success, i0, j0, max_steps

		! Try it out from a specific (x0,y0):
		x0 = 0.9
		y0 = 0.6
		
		i0 = nint((x0-ax)/dx)
		j0 = nint((y0-ay)/dy)

		! shift (x0,y0) to a grid point if it wasn't already:
		x0 = ax + i0*dx
		y0 = ay + j0*dy

		u_true = utrue(x0,y0)
		!nwalks = 0
		n_mc = 10
    	max_steps = 100*max(nx, ny)
		! do i=1,13
		! call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)

		! ---
		! Open a file to write to
		! ---
    		open(unit=25, file='laplace_mc.txt', status='unknown')
		!rel_error = abs(u_true-u_mc)/u_true
		print *, ""
		print '(" True solution of PDE at x =",f10.4," y =",f10.4)', x0, y0
		print '(" u_true =   ",es13.3)', u_true
		print *, ""
		print *, "    Walks    Monte Carlo Solution    Relative Error"
		do i=1,13
			nwalks=0
			call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)
			n_mc = 2.0d0 * n_mc
			rel_error = abs(u_true-u_mc)/u_true
			print '("  ",i8,"      ",es13.3,"        ",es13.3)', nwalks, u_mc, rel_error
			write(25,'(i10,e23.15,e15.6)') nwalks, u_mc, rel_error
		enddo
end program laplace_mc
