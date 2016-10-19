module quadrature_omp

! Performs integration using quadrature integration 
! and creates a table of the error between this and
! the known solution. 

contains

real(kind=8) function trapezoid(f,a,b,n)

	use omp_lib
	implicit none
	real(kind=8), external :: f
    real(kind=8), intent(in) :: a, b
	integer, intent(in) :: n
	real(kind=8) :: h
	integer :: i
	real(kind=8), dimension(n) :: xpoints, ypoints


	h=(b-a)/(n-1)
   
    !$omp parallel private(i)

	!$omp do reduction(+ : trapezoid)

	do i=1,n
		xpoints(i)=(i-1)*h+a
		ypoints(i)=f(xpoints(i))
	enddo

	trapezoid = h*sum(ypoints)-0.5*h*(ypoints(1)+ypoints(n))

	!$omp end parallel

end function trapezoid

subroutine error_table(f,a,b,nvals,int_true)
	implicit none
	real(kind=8), external :: f
	real(kind=8), intent(in) :: a, b, int_true
	integer, dimension(:), intent(in) :: nvals
	integer :: n
	real(kind=8) :: last_error, int_trap, error, ratio

	print *, "      n  trapezoid               error        ratio"
	last_error = 0.d0
	do n=1,size(nvals)
		int_trap = trapezoid(f, a, b, nvals(n))
		error = abs(int_trap - int_true)
		ratio = last_error / error
		last_error = error
		print 11, nvals(n), int_trap, error, ratio
11		format(i8, es22.14, es13.3, es13.3)
	enddo

end subroutine error_table

end module quadrature_omp
