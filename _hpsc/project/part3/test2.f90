program test2
! ---
! Uses random seed module
! ---
	use mc_walk, only: random_walk, many_walks, nwalks
		real(kind=8) :: u_mc
		integer :: n_mc, n_success, i, j, max_steps
 		integer, parameter :: nx = 19
    	integer, parameter :: ny = 11
		! range of i is 1 to nx
		! range of j is 1 to nj
		
		i = 1
		j = 1
    	max_steps = 100*max(nx, ny)
		call many_walks(i, j, max_steps, n_mc, u_mc, n_success)
		print *, "u_mc=", u_mc
		print *, "n_success=", n_success

end program test2
