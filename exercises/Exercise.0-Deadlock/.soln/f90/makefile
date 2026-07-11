FC = mpif90
CFLAGS = -g
EXEC = deadlock

all:
	$(FC) $(FFLAGS) deadlock.f90 -o $(EXEC)
clean:
	@echo "*** Clean : removing object files"
	@( $(RM) *.o *.mod *~)
	@( $(RM) $(EXEC))
