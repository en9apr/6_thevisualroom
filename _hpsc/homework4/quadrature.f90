module quadrature

! Performs integration using quadrature integration 
! and creates a table of the error between this and
! the known solution. 

contains

real(kind=8) function trapezoid(f,a,b,n)

	implicit none
	real(kind=8), external :: f
    real(kind=8), intent(in) :: a, b
	integer, intent(in) :: n
	real(kind=8) :: h
	integer :: i
	real(kind=8), dimension(n) :: xpoints, ypoints

	h=(b-a)/(n-1)
!	print *, " "  ! blank line
!	print 12, h

12  format("h: ", es22.14)
	xpoints=(/((i*h+a),i=0,n-1)/)
	ypoints=(/(f(xpoints(i)),i=1,n)/)
	
	trapezoid = h*sum(ypoints)-0.5*h*(ypoints(1)+ypoints(n))
!	do i=1,5	
!		print 13, xpoints(i), ypoints(i)
!13		format("x,y = " es22.14, es22.14)
!	enddo

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

end module quadrature
