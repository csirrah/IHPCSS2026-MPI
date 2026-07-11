! start of program
program unblock

    ! prevent typographical errors from declaring new variables
    implicit none

    ! include the MPI Fortran headers 
    include "mpif.h"

    ! declare variables 
    integer :: rank                              ! value corresponding to this MPI process
    integer :: total                             ! total number of MPI processes
    integer :: err                               ! error code returned by MPI calls (not checked) 
    integer, parameter :: tag=0                  ! tag for MPI datas (set to 0 and not used)
    integer, parameter :: data_size=16384        ! size of the data
    real*8, dimension(data_size) :: data_send    ! the data the rank sends
    real*8, dimension(data_size) :: data_recv    ! the data the rank recieves
    integer :: partner                           ! the rank of the other process

    ! initialise the MPI implementation 
    call MPI_INIT(err)

    ! determine the size of the MPI_COMM_WORLD communicator 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, total, err)

    ! check there are two MPI processes
    if (total .ne. 2) then
        print*, "This program must be launched with 2 MPI processes, instead of ", total
        error stop
    end if

    ! determine the rank of this MPI process in the MPI_COMM_WORLD communicator
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)

    ! set initial value to rank
    data_send(:) = rank
    partner = mod(rank + 1, 2)

    ! branch based on rank
    if (rank .eq. 0) then

        ! send data and print that it has been sent
        call MPI_SEND(data_send, data_size, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD, err) 
        print*, "Rank ", rank, " sent ", data_send(1)

        ! recieve data and print that it has been received
        call MPI_RECV(data_recv, data_size, MPI_DOUBLE, partner, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
        print*, "Rank ", rank, " received ", data_recv(1)

    else ! (rank is 1)

        ! recieve data and print that it has been received
        call MPI_RECV(data_recv, data_size, MPI_DOUBLE, partner, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err)
        print*, "Rank ", rank, " received ", data_recv(1)

        ! send data and print that it has been sent
        call MPI_SEND(data_send, data_size, MPI_DOUBLE, partner, tag, MPI_COMM_WORLD, err) 
        print*, "Rank ", rank, " sent ", data_send(1)

    endif

    ! finalise the MPI implementation 
    call MPI_FINALIZE(err)
   
! end of program 
end
