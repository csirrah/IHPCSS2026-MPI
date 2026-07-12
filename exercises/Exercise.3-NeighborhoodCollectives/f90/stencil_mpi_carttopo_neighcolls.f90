program stencil_mpi_carttopo_neighcolls
    use mpi
    implicit none

    integer :: r, p, ierr
    integer :: comm, topocomm
    integer :: n, energy, niters
    integer :: args(3)
    integer :: pdims(2), coords(2)
    logical :: periods(2)
    integer :: px, py, rx, ry
    integer :: bx, by, offx, offy
    integer :: west, east, north, south
    integer :: iter, i, j
    integer :: locnsources, locx, locy
    integer :: counts(4), displs(4)
    integer :: req
    integer :: status(MPI_STATUS_SIZE)
    double precision :: heat, rheat

    integer, parameter :: nsources = 3
    integer :: sources(nsources, 2), locsources(nsources, 2)

    double precision, pointer :: aold(:,:), anew(:,:), tmp_ptr(:,:)
    double precision, allocatable :: sbuf(:), rbuf(:)

    character(len=32)  :: carg
    character(len=128) :: fname

    call MPI_Init(ierr)
    comm = MPI_COMM_WORLD
    call MPI_Comm_rank(comm, r,   ierr)
    call MPI_Comm_size(comm, p,   ierr)

    ! Parse arguments on rank 0 and broadcast
    if (r == 0) then
        if (command_argument_count() < 3) then
            print *, "usage: stencil_mpi <n> <energy> <niters>"
            call MPI_Abort(comm, 1, ierr)
        end if
        call get_command_argument(1, carg); read(carg,*) n
        call get_command_argument(2, carg); read(carg,*) energy
        call get_command_argument(3, carg); read(carg,*) niters
        args(1) = n; args(2) = energy; args(3) = niters
        call MPI_Bcast(args, 3, MPI_INTEGER, 0, comm, ierr)
    else
        call MPI_Bcast(args, 3, MPI_INTEGER, 0, comm, ierr)
        n = args(1); energy = args(2); niters = args(3)
    end if

    ! Compute a good 2D domain decomposition
    pdims = 0
    call MPI_Dims_create(p, 2, pdims, ierr)
    px = pdims(1)
    py = pdims(2)

    ! Create Cartesian topology
    periods = .false.
    call MPI_Cart_create(comm, 2, pdims, periods, .false., topocomm, ierr)

    ! Get my (rx, ry) coordinates
    call MPI_Cart_coords(topocomm, r, 2, coords, ierr)
    rx = coords(1)
    ry = coords(2)

    ! Get neighbor ranks
    call MPI_Cart_shift(topocomm, 0, 1, west,  east,  ierr)
    call MPI_Cart_shift(topocomm, 1, 1, north, south, ierr)

    ! Decompose domain
    bx   = n / px
    by   = n / py
    offx = rx * bx
    offy = ry * by

    ! Allocate working arrays (with 1-wide halo in each direction)
    allocate(aold(0:bx+1, 0:by+1))
    allocate(anew(0:bx+1, 0:by+1))
    aold = 0.0d0
    anew = 0.0d0

    ! Determine which of the 3 global heat sources fall in my patch
    sources(1,:) = [n/2,     n/2    ]
    sources(2,:) = [n/3,     n/3    ]
    sources(3,:) = [n*4/5,   n*8/9  ]
    locnsources = 0
    do i = 1, nsources
        locx = sources(i,1) - offx
        locy = sources(i,2) - offy
        if (locx >= 0 .and. locx < bx .and. locy >= 0 .and. locy < by) then
            locnsources = locnsources + 1
            locsources(locnsources, 1) = locx + 1   ! shift into halo-zone offset
            locsources(locnsources, 2) = locy + 1
        end if
    end do

    ! Communication buffers (west, east, north, south interleaved)
    allocate(sbuf(2*by + 2*bx))
    allocate(rbuf(2*by + 2*bx))

    ! ===================================================================
    ! Step 1: Neighborhood collective parameters:
    !   [west: by] [east: by] [north: bx] [south: bx]
    ! TODO: fillout 'counts' array: how much elements of data are we sending to each neighbor (hint: look at the packing of sbuf above)
    counts = [!TODO, TODO, TODO, TODO]
    ! TODO: fillout 'displs' array: What are the starting indexes of each "chunk" (hint: look at the packing of sbuf above)
    displs = [!TODO, TODO, TODO, TODO]
    ! ===================================================================

    do iter = 0, niters-1

        ! Refresh heat sources
        do i = 1, locnsources
            aold(locsources(i,1), locsources(i,2)) = &
                aold(locsources(i,1), locsources(i,2)) + energy
        end do

        ! Pack send buffer: west col, east col, north row, south row
        do i = 1, by
            sbuf(i)      = aold(1,  i)
            sbuf(by+i)   = aold(bx, i)
        end do
        do i = 1, bx
            sbuf(2*by+i)    = aold(i, 1)
            sbuf(2*by+bx+i) = aold(i, by)
        end do

        ! ===================================================================
        ! Step 2. Nonblocking neighborhood collective
        !
        !   MPI_INEIGHBOR_ALLTOALLV(SENDBUF, SENDCOUNTS, SDISPLS, SENDTYPE,
        !                           RECVBUF, RECVCOUNTS, RDISPLS, RECVTYPE,
        !                           COMM, REQUEST, IERROR)
        !
        !   sendcounts / recvcounts : elements to send/recv per neighbor
        !   sdispls    / rdispls    : offsets into send/recv buffer per neighbor
        !   Neighbor order for 2D Cart: [dim0-, dim0+, dim1-, dim1+]
        !                              = [west,  east,  north, south]
        ! TODO: perform nonblocking neighborhood collective call
        !       (counts serves for both sendcounts and recvcounts;
        !        displs serves for both sdispls and rdispls)
                                     
        !
        !   MPI_WAIT( REQUEST, STATUS, IERROR)
        ! TODO: wait for request to complete
        
        ! ===================================================================

        ! Unpack receive buffer into halo zones
        do i = 1, by
            aold(0,    i) = rbuf(i)
            aold(bx+1, i) = rbuf(by+i)
        end do
        do i = 1, bx
            aold(i, 0)    = rbuf(2*by+i)
            aold(i, by+1) = rbuf(2*by+bx+i)
        end do

        ! Update interior grid points
        heat = 0.0d0
        do j = 1, by
            do i = 1, bx
                anew(i,j) = aold(i,j) / 2.0d0 + &
                            (aold(i-1,j) + aold(i+1,j) + &
                             aold(i,j-1) + aold(i,j+1)) / 4.0d0 / 2.0d0
                heat = heat + anew(i,j)
            end do
        end do

        ! Swap working arrays (pointer swap — no data copy)
        tmp_ptr => anew
        anew    => aold
        aold    => tmp_ptr

        ! Write BMP output on final iteration
        if (iter == niters-1) then
            call printarr_par(iter, anew, n, px, py, rx, ry, bx, by, offx, offy, topocomm)
        end if

    end do

    ! Sum heat across all processes
    call MPI_Allreduce(heat, rheat, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    if (r == 0) write(*,'(A,I0,A,F18.6)') "[", r, "] last heat: ", rheat

    deallocate(sbuf, rbuf)
    deallocate(aold, anew)
    call MPI_Finalize(ierr)

contains

    subroutine printarr_par(iter, array, size_n, px, py, rx, ry, bx, by, offx, offy, comm)
        integer, intent(in) :: iter, size_n, px, py, rx, ry, bx, by, offx, offy, comm
        double precision, intent(in) :: array(0:bx+1, 0:by+1)

        integer :: myrank, fh, ierr2
        integer :: linesize, padding, i, rgb
        integer :: xcnt, ycnt, my_xcnt, my_ycnt
        integer :: pos
        integer(kind=1) :: hdr(54), b4(4), b2(2)
        integer(kind=4) :: i4
        integer(kind=2) :: i2
        integer(kind=1), allocatable :: myline(:)
        integer :: wstatus(MPI_STATUS_SIZE)
        character(len=128) :: fname

        call MPI_Comm_rank(comm, myrank, ierr2)

        write(fname, '(A,I0,A)') './output-', iter, '.bmp'
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
            i4 = 40_4                                       ! header_sz
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

        linesize = bx * 3
        padding  = 0
        if (mod(rx+1, px) == 0) padding = mod(size_n*3, 4)
        allocate(myline(linesize + padding))
        myline = 0_1

        my_xcnt = 0
        my_ycnt = 0
        xcnt    = 0
        ycnt    = size_n

        do while (ycnt >= 0)
            call MPI_Barrier(comm, ierr2)
            if (xcnt == offx .and. ycnt >= offy .and. ycnt < offy+by) then
                do i = 0, linesize-1, 3
                    if (i/3 > bx) then
                        rgb = 0
                    else
                        rgb = nint(255.0d0 * array(i/3, by-my_ycnt))
                    end if
                    if (i == 0 .or. i/3 == bx-1 .or. &
                        my_ycnt == 0 .or. my_ycnt == by-1) rgb = 255
                    if (rgb > 255) rgb = 255
                    myline(i+1) = 0_1
                    myline(i+2) = 0_1
                    if (rgb > 127) then
                        myline(i+3) = int(rgb - 256, kind=1)
                    else
                        myline(i+3) = int(rgb, kind=1)
                    end if
                end do
                my_xcnt = my_xcnt + bx
                my_ycnt = my_ycnt + 1
                call MPI_File_write_shared(fh, myline, linesize+padding, &
                                            MPI_BYTE, MPI_STATUS_IGNORE, ierr2)
            end if
            xcnt = xcnt + bx
            if (xcnt >= size_n) then
                xcnt = 0
                ycnt = ycnt - 1
            end if
        end do

        deallocate(myline)
        call MPI_File_close(fh, ierr2)

    end subroutine printarr_par

end program stencil_mpi_carttopo_neighcolls
