#!/bin/bash
mpirun -np 4 ./stencil_mpi_carttopo_neighcolls 512 10 10000 2 2
