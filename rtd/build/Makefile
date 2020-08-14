# define compiler options
#	- cray-netcdf must be loaded to handle netcdf linking
FC      := ftn
FFLAGS  := -static

# define source files which will be brought in from $SRCDIR
SRCDIR				:= $(realpath $(CURDIR)/..)
NETCDF_SRC          := uvic_netcdf.f
PHYS_RTD_FILES      := ccc_nemo_rtd_utils.F90 ccc_nemo_physical_rtd.F90
ICE_RTD_FILES       := ccc_nemo_rtd_utils.F90 ccc_nemo_ice_rtd.F90
CMOC_RTD_FILES      := ccc_nemo_rtd_utils.F90 ccc_nemo_carbon-cmoc_rtd.F90
CANOE_RTD_FILES     := ccc_nemo_rtd_utils.F90 ccc_nemo_carbon-canoe_rtd.F90

# make object substitution on source file
NETCDF_OBJ      := $(NETCDF_SRC:.f=.o)
PHYS_RTD_OBJS   := $(PHYS_RTD_FILES:.F90=.o)
CMOC_RTD_OBJS   := $(CMOC_RTD_FILES:.F90=.o)
CANOE_RTD_OBJS  := $(CANOE_RTD_FILES:.F90=.o)
ICE_RTD_OBJS    := $(ICE_RTD_FILES:.F90=.o)

# make all targets by default
.DEFAULT_GOAL := all

# define implicit rules for handling automatic compiling of the object files
%.o : $(SRCDIR)/%.F90
	$(FC) $(FFLAGS) -c $<
%.o : $(SRCDIR)/%.f
	$(FC) $(FFLAGS) -c $<

# physical rtds
nemo_physical_rtd.exe: $(PHYS_RTD_OBJS) $(NETCDF_OBJ)
	$(FC) $(FFLAGS) -o $@ $(PHYS_RTD_OBJS) $(NETCDF_OBJ)

# cmoc rtds
nemo_carbon-cmoc_rtd.exe: $(CMOC_RTD_OBJS) $(NETCDF_OBJ)
	$(FC) $(FFLAGS) -o $@ $(CMOC_RTD_OBJS) $(NETCDF_OBJ)

# canoe rtds
nemo_carbon-canoe_rtd.exe: $(CANOE_RTD_OBJS) $(NETCDF_OBJ)
	$(FC) $(FFLAGS) -o $@ $(CANOE_RTD_OBJS) $(NETCDF_OBJ)

# ice rtds
nemo_ice_rtd.exe: $(ICE_RTD_OBJS) $(NETCDF_OBJ)
	$(FC) $(FFLAGS) -o $@ $(ICE_RTD_OBJS) $(NETCDF_OBJ)

all : nemo_physical_rtd.exe nemo_carbon-cmoc_rtd.exe nemo_carbon-canoe_rtd.exe nemo_ice_rtd.exe

clean :
	rm -f *.o *.mod nemo_physical_rtd.exe nemo_carbon-cmoc_rtd.exe nemo_carbon-canoe_rtd.exe nemo_ice_rtd.exe
