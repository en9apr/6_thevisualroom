! ---
! This has two functions random_walk and many_walks
! ---
module mc_walk
! ---
! Uses random seed module
! ---
	use random_util, only: init_random_seed
	use problem_description, only: uboundary, ax, ay, dx, dy, nx, ny
! ---
! nwalks is a module variable
! Example use: use mc_walk, only nwalks
! ---
    implicit none
    integer :: nwalks
	!integer :: iabort 
    save

contains

! ---
! A single random walk
! Inputs: i (step in x), j (step in y), max_steps (maximum number of steps)
! Outputs: ub (value at boundary), iabort (0 for success, 1 for failure)
! Example use: ub, iabort = random_walk(0, 0, 100, ub, iabort) 
!     
!    i:  0,  1,.. nx, nx+1
! j:
! 0,   (BC) (BC) (BC) (BC)
! 1,   (BC)           (BC)
! ...  (BC)           (BC)
! nj   (BC)           (BC)
! nj+1 (BC) (BC) (BC) (BC)
! ---
	subroutine random_walk(i0, j0, max_steps, ub, iabort)
	! ---
	! Declare parameters
	! ---
		implicit none
		integer, intent(in) :: i0, j0, max_steps 	! i0 and j0 are inputs (static)
		integer :: i, j 							! i and j are dummy variables
		real(kind=8), intent(out) :: ub
		integer, intent(out) :: iabort
		real(kind=8), dimension(:), allocatable :: x
		integer :: istep, seed1
		real(kind=8) :: xb, yb
	! ---
	! Number of steps is unknown, so we must allocate:
	! ---
		allocate(x(max_steps)) ! random number array 0 to 1 
	! ---
	! Generate random number array x (0 to 1)
	! ---
		seed1 = 0   ! seed1 = 12345 for repeatable seed
						! seed1 = 0 for random seed
		call init_random_seed(seed1)
		call random_number(x)
		i = i0
		j = j0
		!nwalks = 0
	! ---
	! Go through the vector of random numbers and
	! make moves depending on the value of x
	! ---
		do istep=1,max_steps
			if (x(istep) .LT. 0.25) then
				i = i-1    ! go left
				! print *, "go left"
				! print *, "i=", i
				! print *, "istep=", istep
				! print *, "x=", x(istep)
			else if (x(istep) .LT. 0.5) then
				i = i+1    ! go right
				! print *, "go right"
				! print *, "i=", i
				! print *, "istep=", istep
				! print *, "x=", x(istep)
			else if (x(istep) .LT. 0.75) then
				j = j-1    ! go down
				! print *, "go up"
				! print *, "j=", j
				! print *, "istep=", istep
				! print *, "x=", x(istep)
			else
				j = j+1    ! go up
				! print *, "go down"
				! print *, "j=", j
				! print *, "istep=", istep
				! print *, "x=", x(istep)
			endif
		
		! ---
		! What if we are at the boundary?
		! ---			
			if (i*j*(nx+1 - i)*(ny+1 - j) .EQ. 0) then
				! print *, "we are at the boundary"
				! print *, " "
				xb = ax + i*dx 			! x is ax or ax+i*dx
				yb = ay + j*dy			! y is ay or ay+j*dy
				! print *, "xb=",xb
				! print *, "yb=",yb
				ub = uboundary(xb, yb)	! the boundary value is then computed
				iabort = 0				! success
				go to 99				! end do loop
			endif
		! ---
		! What if we never reach the boundary within max steps?		
		! ---
			if (istep .EQ. max_steps) then 
			! 	print *, "we aborted"
			! 	print *, " "
				ub = 0					! set ub to zero
				iabort = 1				! failure
			endif
		enddo

99  continue


	end subroutine random_walk 

	subroutine many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)
	! ---
	! Declare parameters
	! ---
		implicit none
		integer, intent(in) :: i0, j0, max_steps 		! i0 and j0 are inputs (static)
		integer, intent(in) :: n_mc
		integer :: i, j, k, iabort !, nwalks		! i and j are dummy variables
		real(kind=8), intent(out) :: u_mc
		integer, intent(out) :: n_success
		real(kind=8) :: ub_sum, ub

		ub_sum = 0
		n_success = 0
		!nwalks = 0
		do k=1,n_mc
			!nwalks = 0
			i = i0
			j = j0
			call random_walk(i, j, max_steps, ub, iabort)
			if(iabort .NE. 1) then
				ub_sum = ub_sum + ub
				n_success = n_success + 1
				nwalks = nwalks +1	
			endif
			u_mc = ub_sum/n_success
	
		enddo

	end subroutine many_walks

end module mc_walk

