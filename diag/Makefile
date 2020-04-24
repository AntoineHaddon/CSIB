FC      := ftn
SRCDIR  := ".."
FFLAGS	:= -static -I${NETCDF_DIR}/include
LDFLAGS := -L${NETCDF_DIR}/lib -lnetcdff -lnetcdf

# define source files
NETCDF_SRC		:= uvic_netcdf.f
PHYS_DIAG_FILES := nemo_diag_glovars.F90 nemo_diag_cal.F90 nemo_diag.F90

# make object substitution on source file
NETCDF_OBJ		:= $(NETCDF_SRC:.f=.o)
PHYS_DIAG_OBJS	:= $(PHYS_DIAG_FILES:.F90=.o)

all: nemo_diag.exe

%.o : %.F90
	$(FC) $(FFLAGS) -c $<

# physical diagnostics
nemo_diag.exe: $(PHYS_DIAG_OBJS) $(NETCDF_OBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(PHYS_DIAG_OBJS) $(NETCDF_OBJ)

# clean
clean :
	rm -f *.o *.mod nemo_diag.exe
