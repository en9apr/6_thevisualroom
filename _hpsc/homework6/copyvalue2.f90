! $UWHPSC/codes/mpi/copyvalue2.f90
!
! Set value in Process numprocs-1 and copy this through a chain of processes
! and finally print result from Process 0.
!

program copyvalue2

    use mpi

    implicit none

    integer :: i, proc_num, num_procs,ierr
    integer, dimension(MPI_STATUS_SIZE) :: status

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)

    if (num_procs .EQ. 1) then
        print *, "Only one process, cannot do anything"
        call MPI_FINALIZE(ierr)
        stop
        endif


    if (proc_num .EQ. num_procs-1) then
        i = 55
        print '("Process ",i3," setting      i = ",i3)', proc_num, i

        call MPI_SEND(i, 1, MPI_INTEGER, num_procs-2, 21, &
                      MPI_COMM_WORLD, ierr)

      else if (proc_num .GT. 0) then

        call MPI_RECV(i, 1, MPI_INTEGER, proc_num+1, 21, &
                      MPI_COMM_WORLD, status, ierr)

        print '("Process ",i3," receives     i = ",i3)', proc_num, i
        print '("Process ",i3," sends        i = ",i3)', proc_num, i-1

        call MPI_SEND(i-1, 1, MPI_INTEGER, proc_num-1, 21, &
                      MPI_COMM_WORLD, ierr)


      else if (proc_num .EQ. 0) then

        call MPI_RECV(i, 1, MPI_INTEGER, 1, 21, &
                      MPI_COMM_WORLD, status, ierr)

        print '("Process ",i3," ends up with i = ",i3)', proc_num, i
      endif

    call MPI_FINALIZE(ierr)

end program copyvalue2
