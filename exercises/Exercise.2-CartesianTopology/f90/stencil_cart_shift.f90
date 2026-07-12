program stencil_cart_shift
    use mpi
    implicit none

    integer :: rank, p, ierr
    integer :: n, energy, niters, px, py
    integer :: rx, ry
    integer :: north, south, east, west
    integer :: bx, by, offx, offy
    integer :: out_freq
    integer :: comm_cart
    integer :: dims(2), coords(2)
    logical :: periods(2)
    integer :: reqs(8)
    integer :: statuses(MPI_STATUS_SIZE, 8)

    integer, parameter :: nsources = 3
    integer :: sources(nsources, 2), locsources(nsources, 2)
    integer :: locnsources, locx, locy

    double precision, pointer :: aold(:,:), anew(:,:), tmp_ptr(:,:)
    double precision, allocatable :: sbufnorth(:), sbufsouth(:)
    double precision, allocatable :: sbufeast(:), sbufwest(:)
    double precision, allocatable :: rbufnorth(:), rbufsouth(:)
    double precision, allocatable :: rbufeast(:), rbufwest(:)

    double precision :: heat, rheat, t1, t2
    integer :: iter, i, j, final_flag

    character(len=32)  :: carg

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)

    ! Argument parsing and validation
    final_flag = 0
    if (command_argument_count() < 5) then
        if (rank == 0) print *, "usage: stencil_mpi <n> <energy> <niters> <px> <py>"
        final_flag = 1
    else
        call get_command_argument(1, carg); read(carg,*) n
        call get_command_argument(2, carg); read(carg,*) energy
        call get_command_argument(3, carg); read(carg,*) niters
        call get_command_argument(4, carg); read(carg,*) px
        call get_command_argument(5, carg); read(carg,*) py
        if (px * py /= p) then
            if (rank == 0) print *, "px * py must equal to the number of processes."
            call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if (mod(n, px) /= 0) then
            if (rank == 0) print *, "grid size n must be divisible by px."
            call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
        end if
        if (mod(n, py) /= 0) then
            if (rank == 0) print *, "grid size n must be divisible by py."
            call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
        end if
    end if

    if (final_flag == 1) then
        call MPI_Finalize(ierr)
        stop
    end if

    out_freq = 1000

    ! =========================================================================
    ! Step 1. Create Cartesian communicator
    !
    !   Summary:
    !     Creates a new communicator with Cartesian topology information
    !     attached.
    !
    !   Function Call:
    !     call MPI_Cart_create(comm_old, ndims, dims, periods, reorder, &
    !                          comm_cart, ierr)
    !
    !   Input Parameters:
    !     comm_old
    !       Input communicator (MPI communicator).
    !     ndims
    !       Number of dimensions of the Cartesian grid (integer).
    !     dims
    !       Integer array of size ndims specifying the number of
    !       processes in each dimension.
    !     periods
    !       Logical array of size ndims specifying whether each
    !       dimension is periodic (.true.) or non-periodic (.false.).
    !     reorder
    !       Logical flag indicating whether MPI may reorder process
    !       ranks (.true.) or not (.false.).
    !
    !   Output Parameters:
    !     comm_cart
    !       Communicator with the new Cartesian topology (MPI communicator).
    !     ierr
    !       Error status returned by the MPI routine (integer).
    !
    periods = ! TODO
    dims(1) = ! TODO
    dims(2) = ! TODO
    ! TODO create comm_cart communicator

    ! =========================================================================

    ! =========================================================================
    ! Step 2. Get neighbor ranks via Cart_shift
    !
    !   Summary:
    !     Returns the source and destination ranks corresponding to a shift
    !     along a specified Cartesian coordinate direction.
    !
    !   Function Call:
    !     call MPI_Cart_shift(comm, direction, disp, rank_source, &
    !                         rank_dest, ierr)
    !
    !   Input Parameters:
    !     comm
    !       Communicator with Cartesian topology (MPI communicator).
    !     direction
    !       Coordinate dimension along which to perform the shift (integer).
    !     disp
    !       Displacement of the shift:
    !         > 0 : positive direction
    !         < 0 : negative direction
    !       (integer).
    !
    !   Output Parameters:
    !     rank_source
    !       Rank of the source process (integer).
    !     rank_dest
    !       Rank of the destination process (integer).
    !     ierr
    !       Error status returned by the MPI routine (integer).
    !
    !   dim 0 (x): west and east neighbors
    !   dim 1 (y): north and south neighbors
    !
    ! TODO: set west,east,south,north using MPI_Cart_shift
   
  
    ! =========================================================================

    ! =========================================================================
    ! Step 3. Get my (rx, ry) coordinates
    !
    !   Summary:
    !     Returns the Cartesian coordinates of a specified process within
    !     a communicator with Cartesian topology.
    !
    !   Function Call:
    !     call MPI_Cart_coords(comm, rank, maxdims, coords, ierr)
    !
    !   Input Parameters:
    !     comm
    !       Communicator with Cartesian topology (MPI communicator).
    !     rank
    !       Rank of the process within the communicator (integer).
    !     maxdims
    !       Length of the coords array in the calling program (integer).
    !
    !   Output Parameters:
    !     coords
    !       Integer array of size ndims containing the Cartesian
    !       coordinates of the specified process.
    !     ierr
    !       Error status returned by the MPI routine (integer).
    !
    ! TODO: get coords
    
    ! TODO: set rx and ry from coords
    rx = ! TODO
    ry = ! TODO
    ! =========================================================================

    ! Decompose domain
    bx   = n / px
    by   = n / py
    offx = rx * bx
    offy = ry * by

    ! Print neighbor info matching C format
    write(*,'(A,I2,A,I0,A,I0,A)', advance='no') &
        'Rank[', rank, ']' // achar(9) // ' (', ry, ',', rx, ') - '
    if (west  < 0) then
        write(*,'(A,I0,A)', advance='no') '[w]: ', west,  ', '
    else
        write(*,'(A,I2,A)', advance='no') '[w]: ', west,  ', '
    end if
    if (east  < 0) then
        write(*,'(A,I0,A)', advance='no') '[e]: ', east,  ', '
    else
        write(*,'(A,I2,A)', advance='no') '[e]: ', east,  ', '
    end if
    if (north < 0) then
        write(*,'(A,I0,A)', advance='no') '[n]: ', north, ', '
    else
        write(*,'(A,I2,A)', advance='no') '[n]: ', north, ', '
    end if
    if (south < 0) then
        write(*,'(A,I0)') '[s]: ', south
    else
        write(*,'(A,I2)') '[s]: ', south
    end if

    ! Initialize heat sources
    sources(1,:) = [n/2,   n/2  ]
    sources(2,:) = [n/3,   n/3  ]
    sources(3,:) = [n*4/5, n*8/9]
    locnsources = 0
    do i = 1, nsources
        locx = sources(i,1) - offx
        locy = sources(i,2) - offy
        if (locx >= 0 .and. locx < bx .and. locy >= 0 .and. locy < by) then
            locnsources = locnsources + 1
            locsources(locnsources, 1) = locx + 1
            locsources(locnsources, 2) = locy + 1
        end if
    end do

    ! Allocate working arrays (with 1-wide halo)
    allocate(aold(0:bx+1, 0:by+1))
    allocate(anew(0:bx+1, 0:by+1))
    aold = 0.0d0
    anew = 0.0d0

    ! Allocate communication buffers
    allocate(sbufnorth(bx), sbufsouth(bx), sbufeast(by), sbufwest(by))
    allocate(rbufnorth(bx), rbufsouth(bx), rbufeast(by), rbufwest(by))
    sbufnorth = 0.0d0; sbufsouth = 0.0d0
    sbufeast  = 0.0d0; sbufwest  = 0.0d0
    rbufnorth = 0.0d0; rbufsouth = 0.0d0
    rbufeast  = 0.0d0; rbufwest  = 0.0d0

    t1 = MPI_Wtime()

    do iter = 0, niters - 1

        ! Refresh heat sources (set value, not accumulate)
        do i = 1, locnsources
            aold(locsources(i,1), locsources(i,2)) = energy
        end do

        ! Pack data into send buffers
        sbufnorth(1:bx) = aold(1:bx, 1)
        sbufsouth(1:bx) = aold(1:bx, by)
        sbufeast(1:by)  = aold(bx, 1:by)
        sbufwest(1:by)  = aold(1, 1:by)

        ! Nonblocking exchange with all four neighbors
        call MPI_Irecv(rbufnorth, bx, MPI_DOUBLE_PRECISION, north, 9, MPI_COMM_WORLD, reqs(1), ierr)
        call MPI_Irecv(rbufsouth, bx, MPI_DOUBLE_PRECISION, south, 9, MPI_COMM_WORLD, reqs(2), ierr)
        call MPI_Irecv(rbufeast,  by, MPI_DOUBLE_PRECISION, east,  9, MPI_COMM_WORLD, reqs(3), ierr)
        call MPI_Irecv(rbufwest,  by, MPI_DOUBLE_PRECISION, west,  9, MPI_COMM_WORLD, reqs(4), ierr)

        call MPI_Isend(sbufnorth, bx, MPI_DOUBLE_PRECISION, north, 9, MPI_COMM_WORLD, reqs(5), ierr)
        call MPI_Isend(sbufsouth, bx, MPI_DOUBLE_PRECISION, south, 9, MPI_COMM_WORLD, reqs(6), ierr)
        call MPI_Isend(sbufeast,  by, MPI_DOUBLE_PRECISION, east,  9, MPI_COMM_WORLD, reqs(7), ierr)
        call MPI_Isend(sbufwest,  by, MPI_DOUBLE_PRECISION, west,  9, MPI_COMM_WORLD, reqs(8), ierr)

        call MPI_Waitall(8, reqs, statuses, ierr)

        ! Unpack received data into halo zones
        aold(1:bx, 0)    = rbufnorth(1:bx)
        aold(1:bx, by+1) = rbufsouth(1:bx)
        aold(bx+1, 1:by) = rbufeast(1:by)
        aold(0, 1:by)    = rbufwest(1:by)

        ! Update grid points
        heat = 0.0d0
        do j = 1, by
            do i = 1, bx
                anew(i,j) = anew(i,j) / 2.0d0 + &
                            (aold(i-1,j) + aold(i+1,j) + &
                             aold(i,j-1) + aold(i,j+1)) / 4.0d0 / 2.0d0
                heat = heat + anew(i,j)
            end do
        end do

        ! Swap working arrays (pointer swap — no data copy)
        tmp_ptr => anew
        anew    => aold
        aold    => tmp_ptr

        ! Output BMP at requested intervals
        if (mod(iter, out_freq) == 0 .or. iter == niters - 1) then
            if (rank == 0) write(*,'(A,I0)') 'Outputting state: ', iter
            call printarr_par(iter, anew, n, px, py, rx, ry, bx, by, offx, offy, MPI_COMM_WORLD, p)
        end if

    end do

    t2 = MPI_Wtime()

    ! Free working arrays and buffers
    deallocate(aold, anew)
    deallocate(sbufnorth, sbufsouth, sbufeast, sbufwest)
    deallocate(rbufnorth, rbufsouth, rbufeast, rbufwest)

    ! Sum heat across all processes
    call MPI_Allreduce(heat, rheat, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (rank == 0) write(*,'(A,I0,A,F18.6,A,F10.6)') &
        '[', rank, '] last heat: ', rheat, ' time: ', t2 - t1

    call MPI_Finalize(ierr)

contains

    subroutine printarr_par(iter, array, size_n, px_in, py_in, rx_in, ry_in, &
                             bx_in, by_in, offx_in, offy_in, comm, nproc)
        integer, intent(in) :: iter, size_n, px_in, py_in, rx_in, ry_in
        integer, intent(in) :: bx_in, by_in, offx_in, offy_in, comm, nproc
        double precision, intent(in) :: array(0:bx_in+1, 0:by_in+1)

        integer :: myrank, fh, ierr2
        integer :: linesize, padding, col, rgb
        integer :: xcnt, ycnt, my_xcnt, my_ycnt
        integer :: pos
        integer(kind=1) :: hdr(54), b4(4), b2(2)
        integer(kind=4) :: i4
        integer(kind=2) :: i2
        integer(kind=1), allocatable :: myline(:)
        character(len=128) :: fname
        character(len=8)   :: nprocstr

        call MPI_Comm_rank(comm, myrank, ierr2)

        ! Create output directory
        write(nprocstr, '(I0)') nproc
        if (myrank == 0) then
            call execute_command_line('mkdir -p WRK/' // trim(nprocstr), wait=.true.)
        end if
        call MPI_Barrier(comm, ierr2)

        write(fname, '(A,A,A,I0,A)') './WRK/', trim(nprocstr), '/output-', iter, '.bmp'
        call MPI_File_open(comm, trim(fname), &
                           MPI_MODE_SEQUENTIAL + MPI_MODE_CREATE + MPI_MODE_WRONLY, &
                           MPI_INFO_NULL, fh, ierr2)

        if (myrank == 0) then
            hdr = 0_1
            pos = 1

            ! Magic: 'B','M'
            hdr(1) = int(66, kind=1)
            hdr(2) = int(77, kind=1)
            pos = 3

            ! filesz
            i4 = int(54 + size_n*(size_n*3 + mod(size_n*3, 4)), kind=4)
            b4 = transfer(i4, b4); hdr(pos:pos+3) = b4; pos = pos + 4
            ! creator1
            i2 = int(z'FE', kind=2)
            b2 = transfer(i2, b2); hdr(pos:pos+1) = b2; pos = pos + 2
            ! creator2
            i2 = int(z'FE', kind=2)
            b2 = transfer(i2, b2); hdr(pos:pos+1) = b2; pos = pos + 2
            ! bmp_offset = 54
            i4 = 54_4
            b4 = transfer(i4, b4); hdr(pos:pos+3) = b4; pos = pos + 4

            ! Info header
            i4 = 40_4
            b4 = transfer(i4, b4); hdr(pos:pos+3) = b4; pos = pos + 4
            i4 = int(size_n, kind=4)                        ! width
            b4 = transfer(i4, b4); hdr(pos:pos+3) = b4; pos = pos + 4
            i4 = int(size_n, kind=4)                        ! height
            b4 = transfer(i4, b4); hdr(pos:pos+3) = b4; pos = pos + 4
            i2 = 1_2                                        ! nplanes
            b2 = transfer(i2, b2); hdr(pos:pos+1) = b2; pos = pos + 2
            i2 = 24_2                                       ! bitspp
            b2 = transfer(i2, b2); hdr(pos:pos+1) = b2; pos = pos + 2
            i4 = 0_4                                        ! compress_type
            b4 = transfer(i4, b4); hdr(pos:pos+3) = b4; pos = pos + 4
            i4 = int(size_n*(size_n*3 + mod(size_n*3,4)), kind=4) ! bmp_bytesz
            b4 = transfer(i4, b4); hdr(pos:pos+3) = b4; pos = pos + 4
            i4 = int(size_n, kind=4)                        ! hres
            b4 = transfer(i4, b4); hdr(pos:pos+3) = b4; pos = pos + 4
            i4 = int(size_n, kind=4)                        ! vres
            b4 = transfer(i4, b4); hdr(pos:pos+3) = b4; pos = pos + 4
            i4 = 0_4                                        ! ncolors
            b4 = transfer(i4, b4); hdr(pos:pos+3) = b4; pos = pos + 4
            i4 = 0_4                                        ! nimpcolors
            b4 = transfer(i4, b4); hdr(pos:pos+3) = b4; pos = pos + 4

            call MPI_File_write_shared(fh, hdr, 54, MPI_BYTE, MPI_STATUS_IGNORE, ierr2)
        end if

        linesize = bx_in * 3
        padding  = 0
        if (mod(rx_in + 1, px_in) == 0) padding = mod(size_n * 3, 4)
        allocate(myline(linesize + padding))
        myline = 0_1

        my_xcnt = 0
        my_ycnt = 0
        xcnt    = 0
        ycnt    = size_n

        do while (ycnt >= 0)
            call MPI_Barrier(comm, ierr2)
            if (xcnt == offx_in .and. ycnt >= offy_in .and. ycnt < offy_in + by_in) then
                do col = 0, bx_in - 1
                    rgb = nint(255.0d0 * array(col + 1, by_in - my_ycnt))
                    if (col == 0 .or. col == bx_in - 1 .or. &
                        my_ycnt == 0 .or. my_ycnt == by_in - 1) rgb = 255
                    if (rgb > 255) rgb = 255
                    myline(col*3 + 1) = int(-128, kind=1)   ! blue = 128
                    myline(col*3 + 2) = 0_1                  ! green = 0
                    if (rgb > 127) then
                        myline(col*3 + 3) = int(rgb - 256, kind=1)
                    else
                        myline(col*3 + 3) = int(rgb, kind=1)
                    end if
                end do
                my_xcnt = my_xcnt + bx_in
                my_ycnt = my_ycnt + 1
                call MPI_File_write_shared(fh, myline, linesize + padding, &
                                            MPI_BYTE, MPI_STATUS_IGNORE, ierr2)
            end if
            xcnt = xcnt + bx_in
            if (xcnt >= size_n) then
                xcnt = 0
                ycnt = ycnt - 1
            end if
        end do

        deallocate(myline)
        call MPI_File_close(fh, ierr2)

    end subroutine printarr_par

end program stencil_cart_shift
