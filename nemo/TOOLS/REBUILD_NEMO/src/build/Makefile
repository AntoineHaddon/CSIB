# define compiler options
#   - cray-netcdf must be loaded to handle netcdf linking
#	- openmp flag will need to be changed for different compilers (its using the intel settings here)
FC		:= ftn
FFLAGS	:= -static -qopenmp

# define source files which will be brought in from $SRCDIR
SRCDIR			:= $(realpath $(CURDIR)/..)
REBLD_SRC		:= rebuild_nemo.f90
REBLD_MPI_SRC	:= rebuild_nemo_mpi.f90
REBLD_OBJ		:= $(REBLD_SRC:.f90=.o)
REBLD_MPI_OBJ	:= $(REBLD_MPI_SRC:.f90=.o)

# make all targets by default
.DEFAULT_GOAL := all

%.o : $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $<

all : rebuild_nemo.exe rebuild_nemo_mpi.exe

# non-mpi
rebuild_nemo.exe: $(REBLD_OBJ)
	$(FC) $(FFLAGS) -o $@ $(REBLD_OBJ)

# mpi version
rebuild_nemo_mpi.exe: $(REBLD_MPI_OBJ)
	$(FC) $(FFLAGS) -o $@ $(REBLD_MPI_OBJ)

clean: 
	rm -f *.o *.mod rebuild_nemo.exe rebuild_nemo_mpi.exe
