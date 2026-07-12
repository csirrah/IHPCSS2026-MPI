MPIFC=mpif90
FFLAGS=
BINS=stencil_mpi_carttopo_neighcolls

all: $(BINS)

stencil_mpi_carttopo_neighcolls: stencil_mpi_carttopo_neighcolls.f90 Makefile
	$(MPIFC) $(FFLAGS) -o $@ $<

clean:
	rm -f $(BINS) *.o *.mod
	rm -f output*.bmp
