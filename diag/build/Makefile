# define compiler options
#	- NETCDF_DIR gets inherited from the cray-netcdf module
FC      := ftn
FFLAGS	:= -static -I${NETCDF_DIR}/include
LDFLAGS := -L${NETCDF_DIR}/lib -lnetcdff -lnetcdf

# define source files which will be brought in from $SRCDIR
SRCDIR				:= $(realpath $(CURDIR)/..)
NETCDF_SRC			:= uvic_netcdf.f
PHYS_DIAG_FILES		:= nemo_diag_glovars.F90 nemo_diag_cal.F90 nemo_diag.F90
CMOC_DIAG_FILES		:= nemo_diag_glovars_cmoc.F90 nemo_diag_cal_cmoc.F90 nemo_diag_cmoc.F90
CANOE_DIAG_FILES	:= nemo_diag_glovars_canoe.F90 nemo_diag_cal_canoe.F90 nemo_diag_canoe.F90

# make object substitution on source file
NETCDF_OBJ		:= $(NETCDF_SRC:.f=.o)
PHYS_DIAG_OBJS	:= $(PHYS_DIAG_FILES:.F90=.o)
CMOC_DIAG_OBJS	:= $(CMOC_DIAG_FILES:.F90=.o)
CANOE_DIAG_OBJS	:= $(CANOE_DIAG_FILES:.F90=.o)

# make all targets by default
.DEFAULT_GOAL := all

# define implicit rules such that make knows to look in SRCDIR
%.o : $(SRCDIR)/%.F90
	$(FC) $(FFLAGS) -c $<
%.o : $(SRCDIR)/%.f
	$(FC) $(FFLAGS) -c $<

# physical diagnostics
nemo_diag.exe: $(PHYS_DIAG_OBJS) $(NETCDF_OBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(PHYS_DIAG_OBJS) $(NETCDF_OBJ)

# cmoc diagnostics
nemo_diag_cmoc.exe: $(CMOC_DIAG_OBJS) $(NETCDF_OBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(CMOC_DIAG_OBJS) $(NETCDF_OBJ)

# canoe diagnostics
nemo_diag_canoe.exe: $(CANOE_DIAG_OBJS) $(NETCDF_OBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(CANOE_DIAG_OBJS) $(NETCDF_OBJ)

all : nemo_diag.exe nemo_diag_canoe.exe nemo_diag_cmoc.exe

# clean
clean :
	rm -f *.o *.mod nemo_diag.exe nemo_diag_canoe.exe nemo_diag_cmoc.exe
