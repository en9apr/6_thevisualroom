
module functions

    use omp_lib
    implicit none
    integer :: fevals(0:7), gevals(0:7)
    real(kind=8) :: k
    save

contains

    real(kind=8) function f(g, x, c, d)
        implicit none
		real(kind=8), intent(in) :: c,d
        real(kind=8), intent(in) :: x
		real(kind=8), external :: g

        integer thread_num, n, j
		real(kind=8) :: h, trap_sum, y

        ! keep track of number of function evaluations by 
        ! each thread: 
        thread_num = 0   ! serial mode
        !$ thread_num = omp_get_thread_num()
        fevals(thread_num) = fevals(thread_num) + 1
        n = 1000
		h = (d-c)/(n-1)
    	trap_sum = 0.5d0*(g(x,c) + g(x,d))  ! endpoint contributions
    
!    	!$omp parallel do private(y) reduction(+ : trap_sum) 
    	do j=2,n-1
       		y = c + (j-1)*h
        	trap_sum = trap_sum + g(x,y)
        	enddo

    	f = h * trap_sum

    end function f

	real(kind=8) function g(x,y)
		implicit none
		real(kind=8), intent(in) :: x,y
		integer thread_num_g

		thread_num_g = 0

		!$ thread_num_g = omp_get_thread_num()
		gevals(thread_num_g) = gevals(thread_num_g) + 1

		g = sin(x+y)

	end function g

end module functions
