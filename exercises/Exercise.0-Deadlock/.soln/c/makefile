CXX = mpicxx
CFLAGS = -g
EXEC = deadlock

all:
	$(CXX) $(CFLAGS) deadlock.c -o $(EXEC)
clean:
	@echo "*** Clean : removing object files"
	@( $(RM) *.o *.mod *~)
	@( $(RM) $(EXEC))
