! ---
! This module returns seed for the random_number function
! ---
module random_util

contains

	subroutine init_random_seed(seed1)
	! ---
	! Declare the parameters
	! ---
		integer, intent(inout) :: seed1
		integer :: clock 
		integer, dimension(:), allocatable :: seed
	! ---
	! Determine how many numbers needed to seed and allocate
	! ---
		call random_seed(size = nseed)  
		allocate(seed(nseed))
	! ---  
	! If seed1 = 0 then set seed1 randomly using the system clock.
	! This will give different sequences of random numbers
	! from different runs, so results are not reproducible.
	! ---
		if (seed1 == 0) then
		    call system_clock(count = clock)
		    seed1 = clock
		endif
	! ---
	! If seed1 > 0 then results will be reproducible since calling this
	! twice with the same seed will initialize the random
	! number generator so the same sequence is generated.
	! Once seed1 is set, set the other elements of the seed array by adding
	! multiples of 37 as suggested in the documentation. 
	! --- 
		do i=1,nseed
		    seed(i) = seed1 + 37*(i-1)  
		enddo
	! ---
	! Seed the generator
	! ---
		call random_seed(put = seed)
	! ---
	! Deallocate the seed
	! ---
		deallocate(seed)

	end subroutine init_random_seed

end module random_util

