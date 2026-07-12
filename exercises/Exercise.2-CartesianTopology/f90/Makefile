MPIFC  = mpif90
FFLAGS =
BINS   = mpi_stencil_cart

all: $(BINS)

mpi_stencil_cart: stencil_cart_shift.f90 Makefile
	$(MPIFC) $(FFLAGS) -o $@ $<

clean:
	rm -f $(BINS) *.o *.mod
	rm -rf WRK
