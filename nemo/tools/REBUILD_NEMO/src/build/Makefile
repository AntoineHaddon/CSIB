# define default compiler template - should be provided at call line
COMPILER_TEMPLATE ?= 

# warn if no compiler template was given
ifndef COMPILER_TEMPLATE
$(warning For all targets but 'clean' COMPILER_TEMPLATE must be given to the make call)
else
# bring in compiler flags
include $(COMPILER_TEMPLATE)
endif

# define flags to be used 
FFLAGS = $(FFLAGS_STATIC) $(FFLAGS_OMP)
LDFLAGS = $(INCLUDEFLAGS_USR) $(LDFLAGS_USR) $(LDFLAGS_NETCDF)

# define source files which will be brought in from $SRCDIR
SRCDIR          := $(realpath $(CURDIR)/..)
REBLD_SRC       := rebuild_nemo.f90
REBLD_MPI_SRC   := rebuild_nemo_mpi.f90
REBLD_OBJ       := $(REBLD_SRC:.f90=.o)
REBLD_MPI_OBJ   := $(REBLD_MPI_SRC:.f90=.o)

# create implicit rule for `.o` files so it knows to look in SRCDIR
%.o : $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< $(LDFLAGS)
	
# make all targets by default
.DEFAULT_GOAL := all

all : rebuild_nemo.exe rebuild_nemo_mpi.exe

# non-mpi
rebuild_nemo.exe: $(REBLD_OBJ)
	$(FC) $(FFLAGS) -o $@ $(REBLD_OBJ) $(LDFLAGS)

# mpi version
rebuild_nemo_mpi.exe: $(REBLD_MPI_OBJ)
	$(FC) $(FFLAGS) -o $@ $(REBLD_MPI_OBJ) $(LDFLAGS)

clean: 
	rm -f *.o *.mod rebuild_nemo.exe rebuild_nemo_mpi.exe
