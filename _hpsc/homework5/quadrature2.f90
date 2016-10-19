
module quadrature2

    use omp_lib

contains

real(kind=8) function simpson(f, a, b, n)

    ! Estimate the integral of f(x) from a to b using the
    ! Simpson Rule with n points.

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
    real(kind=8) :: h, simp_sum, xj, xj2

    h = (b-a)/(n-1)
    simp_sum = (1.d0/6.d0)*(f(a) + f(b))  ! endpoint contributions
    simp_sum = simp_sum + (2.d0/3.d0)*f(a+0.5d0*h) !for the half before the loop
    !$omp parallel do private(xj,xj2) reduction(+ : simp_sum) 
    do j=2,n-1
        xj = a + (j-1.d0)*h
		xj2 = a + (j-0.5d0)*h
        simp_sum = simp_sum + (1.d0/3.d0)*f(xj) + (2.d0/3.d0)*f(xj2)
        enddo

    simpson = h * simp_sum

end function simpson


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

    print *, "      n         approximation        error       ratio"
    last_error = 0.d0   
    do j=1,size(nvals)
        n = nvals(j)
        int_approx = method(f,a,b,n)
        error = abs(int_approx - int_true)
        ratio = last_error / error
        last_error = error  ! for next n

        print 11, n, int_approx, error, ratio
 11     format(i8, es22.14, es13.3, es13.3)
        enddo

end subroutine error_table


end module quadrature2

