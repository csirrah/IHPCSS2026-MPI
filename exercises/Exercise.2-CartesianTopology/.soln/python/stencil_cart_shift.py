from mpi4py import MPI
import numpy as np
import sys
import os
import struct


def printarr_par(iteration, array, n, px, py, rx, ry, bx, by, offx, offy, comm, nproc):
    myrank = comm.Get_rank()

    nprocstr = str(nproc)
    if myrank == 0:
        os.makedirs(f'WRK/{nprocstr}', exist_ok=True)
    comm.Barrier()

    fname = f'./WRK/{nprocstr}/output-{iteration}.bmp'
    fh = MPI.File.Open(comm, fname,
                       MPI.MODE_SEQUENTIAL | MPI.MODE_CREATE | MPI.MODE_WRONLY)

    if myrank == 0:
        linesize_full = n * 3
        padding_full  = (4 - linesize_full % 4) % 4
        bmp_data_size = n * (linesize_full + padding_full)
        hdr = struct.pack('<2sIHHI',
                          b'BM',
                          54 + bmp_data_size,
                          0xFE, 0xFE,
                          54)
        hdr += struct.pack('<IiiHHIIiiII',
                           40,
                           n, n,
                           1, 24,
                           0,
                           bmp_data_size,
                           n, n,
                           0, 0)
        fh.Write_shared(np.frombuffer(hdr, dtype=np.uint8))

    linesize = bx * 3
    padding  = (n * 3) % 4 if (rx + 1 == px) else 0
    myline   = np.zeros(linesize + padding, dtype=np.uint8)

    xcnt    = 0
    ycnt    = n
    my_ycnt = 0

    while ycnt >= 0:
        comm.Barrier()
        if xcnt == offx and offy <= ycnt < offy + by:
            row = by - my_ycnt
            for col in range(bx):
                val = array[col + 1, row]
                rgb = min(255, max(0, round(255.0 * val)))
                if col == 0 or col == bx - 1 or my_ycnt == 0 or my_ycnt == by - 1:
                    rgb = 255
                myline[col * 3 + 0] = 128   # blue
                myline[col * 3 + 1] = 0     # green
                myline[col * 3 + 2] = rgb   # red
            my_ycnt += 1
            fh.Write_shared(myline)
        xcnt += bx
        if xcnt >= n:
            xcnt = 0
            ycnt -= 1

    fh.Close()


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    p    = comm.Get_size()

    if len(sys.argv) < 6:
        if rank == 0:
            print("usage: stencil_mpi <n> <energy> <niters> <px> <py>")
        MPI.Finalize()
        sys.exit(0)

    n      = int(sys.argv[1])
    energy = int(sys.argv[2])
    niters = int(sys.argv[3])
    px     = int(sys.argv[4])
    py     = int(sys.argv[5])

    if px * py != p:
        if rank == 0:
            print("px * py must equal to the number of processes.")
        comm.Abort(1)
    if n % px != 0:
        if rank == 0:
            print("grid size n must be divisible by px.")
        comm.Abort(2)
    if n % py != 0:
        if rank == 0:
            print("grid size n must be divisible by py.")
        comm.Abort(3)

    out_freq = 1000

    # =========================================================================
    # Step 1. Create Cartesian communicator
    topocomm = comm.Create_cart(dims=[px, py], periods=[False, False], reorder=False)
    # =========================================================================

    # =========================================================================
    # Step 2. Get neighbor ranks via Shift
    west,  east  = topocomm.Shift(0, 1)
    north, south = topocomm.Shift(1, 1)
    # =========================================================================

    # =========================================================================
    # Step 3. Get my (rx, ry) coordinates
    coords = topocomm.Get_coords(rank)
    rx, ry = coords[0], coords[1]
    # =========================================================================

    bx   = n // px
    by   = n // py
    offx = rx * bx
    offy = ry * by

    # Print neighbor info matching C format
    def fmt(label, val, sep):
        if val < 0:
            return f'[{label}]: {val},{sep}'
        else:
            return f'[{label}]: {val:2d},{sep}'

    line = (f'Rank[{rank:2d}]\t ({ry},{rx}) - '
            + fmt('w', west,  ' ')
            + fmt('e', east,  ' ')
            + fmt('n', north, ' ')
            + (f'[s]: {south}' if south < 0 else f'[s]: {south:2d}') + '  ')
    print(line)

    # Heat sources
    sources    = [(n // 2, n // 2), (n // 3, n // 3), (n * 4 // 5, n * 8 // 9)]
    locsources = []
    for sx, sy in sources:
        locx = sx - offx
        locy = sy - offy
        if 0 <= locx < bx and 0 <= locy < by:
            locsources.append((locx + 1, locy + 1))

    # Working arrays with 1-wide halo
    aold = np.zeros((bx + 2, by + 2), dtype=np.float64)
    anew = np.zeros((bx + 2, by + 2), dtype=np.float64)

    # Communication buffers
    sbufnorth = np.zeros(bx, dtype=np.float64)
    sbufsouth = np.zeros(bx, dtype=np.float64)
    sbufeast  = np.zeros(by, dtype=np.float64)
    sbufwest  = np.zeros(by, dtype=np.float64)
    rbufnorth = np.zeros(bx, dtype=np.float64)
    rbufsouth = np.zeros(bx, dtype=np.float64)
    rbufeast  = np.zeros(by, dtype=np.float64)
    rbufwest  = np.zeros(by, dtype=np.float64)

    t1 = MPI.Wtime()

    for iteration in range(niters):

        # Refresh heat sources (set value, not accumulate)
        for lx, ly in locsources:
            aold[lx, ly] = energy

        # Pack
        sbufnorth[:] = aold[1:bx+1, 1]
        sbufsouth[:] = aold[1:bx+1, by]
        sbufeast[:]  = aold[bx, 1:by+1]
        sbufwest[:]  = aold[1, 1:by+1]

        # Nonblocking exchange
        reqs = []
        reqs.append(comm.Irecv([rbufnorth, MPI.DOUBLE], source=north, tag=9))
        reqs.append(comm.Irecv([rbufsouth, MPI.DOUBLE], source=south, tag=9))
        reqs.append(comm.Irecv([rbufeast,  MPI.DOUBLE], source=east,  tag=9))
        reqs.append(comm.Irecv([rbufwest,  MPI.DOUBLE], source=west,  tag=9))
        reqs.append(comm.Isend([sbufnorth, MPI.DOUBLE], dest=north, tag=9))
        reqs.append(comm.Isend([sbufsouth, MPI.DOUBLE], dest=south, tag=9))
        reqs.append(comm.Isend([sbufeast,  MPI.DOUBLE], dest=east,  tag=9))
        reqs.append(comm.Isend([sbufwest,  MPI.DOUBLE], dest=west,  tag=9))
        MPI.Request.Waitall(reqs)

        # Unpack into halo zones
        aold[1:bx+1, 0]    = rbufnorth[:]
        aold[1:bx+1, by+1] = rbufsouth[:]
        aold[bx+1, 1:by+1] = rbufeast[:]
        aold[0, 1:by+1]    = rbufwest[:]

        # Update grid
        anew[1:bx+1, 1:by+1] = (
            anew[1:bx+1, 1:by+1] / 2.0 +
            (aold[0:bx,   1:by+1] + aold[2:bx+2, 1:by+1] +
             aold[1:bx+1, 0:by]   + aold[1:bx+1, 2:by+2]) / 4.0 / 2.0
        )
        heat = np.sum(anew[1:bx+1, 1:by+1])

        aold, anew = anew, aold

        if iteration % out_freq == 0 or iteration == niters - 1:
            if rank == 0:
                print(f'Outputting state: {iteration}')
            printarr_par(iteration, anew, n, px, py, rx, ry, bx, by, offx, offy, topocomm, p)

    t2 = MPI.Wtime()

    rheat = comm.allreduce(heat, op=MPI.SUM)
    if rank == 0:
        print(f'[{rank}] last heat: {rheat:.6f} time: {t2 - t1:.6f}')


if __name__ == '__main__':
    main()
