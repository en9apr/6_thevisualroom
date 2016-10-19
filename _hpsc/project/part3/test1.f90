program test1
! ---
! Uses random seed module
! ---
	use mc_walk, only: random_walk, nwalks
		real(kind=8) :: ub
		integer :: iabort, i, j, max_steps
 		integer, parameter :: nx = 19
    	integer, parameter :: ny = 11
		! range of i is 1 to nx
		! range of j is 1 to nj
		
		i = 1
		j = 1
    	max_steps = 100*max(nx, ny)
		call random_walk(i, j, max_steps, ub, iabort) 
		print *, "ub=", ub
		print *, "iabort=", iabort
		print *, "nwalks=", nwalks

end program test1
