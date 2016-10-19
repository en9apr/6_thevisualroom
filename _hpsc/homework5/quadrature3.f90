
module quadrature3

    use omp_lib

contains

real(kind=8) function trapezoid(f, a, b, n)

    ! Estimate the integral of f(x) from a to b using the
    ! Trapezoid Rule with n points.

    ! Input:
    !   f:  the function to integrate
    !   a:  left endpoint
    !   b:  right endpoint
    !   n:  number of points to use
    ! Returns:
    !   the estimate of the integral
     
    implicit none
    real(kind=8), intent(in) :: a,b
    real(kind=8), external :: f
    integer, intent(in) :: n

    ! Local variables:
    integer :: j
    real(kind=8) :: h, trap_sum, xj

    h = (b-a)/(n-1)
    trap_sum = 0.5d0*(f(a) + f(b))  ! endpoint contributions
    
    do j=2,n-1
        xj = a + (j-1)*h
        trap_sum = trap_sum + f(xj)
        enddo

    trapezoid = h * trap_sum

end function trapezoid


subroutine error_table(f,a,b,nvals,int_true,method)

    ! Compute and print out a table of errors when the quadrature
    ! rule specified by the input function method is applied for
    ! each value of n in the array nvals.

    implicit none
    real(kind=8), intent(in) :: a,b, int_true
    real(kind=8), external :: f, method
    integer, dimension(:), intent(in) :: nvals

    ! Local variables:
    integer :: j, n
    real(kind=8) :: ratio, last_error, error, int_approx
	integer thread_num

    print *, "      n         approximation        error       ratio     thread no."
    last_error = 0.d0   
	!$omp parallel do firstprivate(last_error) private(n, int_approx, error, ratio) & 
    !$omp schedule(dynamic)
	do j=size(nvals),1,-1
        n = nvals(j)
        int_approx = method(f,a,b,n)
        error = abs(int_approx - int_true)
        ratio = last_error / error
        last_error = error  ! for next n -- last_error is firstprivate variable - 
! each thead must start with last_error as zero - so ratio will be wrong.
	 	thread_num = 0   ! serial mode
        !$ thread_num = omp_get_thread_num()
        print 11, n, int_approx, error, ratio, thread_num
 11     format(i8, es22.14, es13.3, es13.3, i8)
        enddo

end subroutine error_table


end module quadrature3

