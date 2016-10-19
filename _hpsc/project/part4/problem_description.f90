! ---
! This module returns the true value and the value at the boundary
! ---
module problem_description
! ---
! Initilize parameters
! ---
    implicit none
    real(kind=8), parameter :: ax = 0.0d0
    real(kind=8), parameter :: bx = 1.0d0
    real(kind=8), parameter :: ay = 0.4d0
    real(kind=8), parameter :: by = 1.0d0
    integer, parameter :: nx = 19
    integer, parameter :: ny = 11
    real(kind=8), parameter :: dx = (bx - ax) / (nx+1)
    real(kind=8), parameter :: dy = (by - ay) / (ny+1)
	save

contains
! ---
! Function for the true solution (if known)
! Inputs: x and y
! Output: z
! ---
	function utrue(x, y)

		implicit none
		real(kind=8), intent(in) :: x,y
		real(kind=8) :: utrue

		utrue = x**2 - y**2

	end function utrue
! ---
! Function for the value at the boundary
! Inputs: x and y
! Output: z
! ---
	function uboundary(x, y)

		implicit none
		real(kind=8), intent(in) :: x,y
		real(kind=8) :: uboundary
	! ---
	! If we are not at the boundary we print an error and stop
	! ---
		if ((x-ax)*(x-bx)*(y-ay)*(y-by) .ne. 0.d0) then
		    print *, "*** Error -- called uboundary at non-boundary point"
		    stop
		endif
	! ---
	! If we are at the boundary we compute the solution
	! ---
		uboundary = utrue(x,y)   ! assuming we know this

	end function uboundary
    
end module problem_description

