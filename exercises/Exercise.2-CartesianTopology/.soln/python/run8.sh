#!/bin/bash
mpirun -np 8 python3 stencil_cart_shift.py 512 10 10000 4 2
