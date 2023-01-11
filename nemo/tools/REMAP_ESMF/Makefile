FC=mpif90
FFLAGS=-traceback -O2 -g -I/home/scrd102/canesm_libs/esmf/ESMF_8_1_1/include -I/home/scrd102/canesm_libs/esmf/ESMF_8_1_1/mod -I/fs/ssm/main/opt/hdf5-netcdf4/master/netcdf-fortran_4.4.5-inteloneapi-2022.1.2-serial-shared_rhel-8-amd64-64/include
LDFLAGS=-L/home/scrd102/canesm_libs/esmf/ESMF_8_1_1/lib -lesmf -L/fs/ssm/main/opt/hdf5-netcdf4/master/netcdf-fortran_4.4.5-inteloneapi-2022.1.2-serial-shared_rhel-8-amd64-64/lib -lnetcdff  -lnetcdf -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz -lcurl -lm
OBJ=remap.exe

all : $(OBJ)

%.mod : %.f90
	$(FC) $(FFLAGS) -c $^

%.o : %.f90 ncdf.mod esmf_utils.mod file_readers.mod river_spreading.mod
	$(FC) $(FFLAGS) -c $< 

esmf_utils.o : esmf_utils.f90 file_readers.mod
	$(FC) $(FFLAGS) -c $< 

file_readers.o : file_readers.f90 ncdf.mod
	$(FC) $(FFLAGS) -c $< 

river_spreading.o : river_spreading.f90
	$(FC) $(FFLAGS) -c $< 

%.exe : %.o ncdf.o xit.o esmf_utils.o file_readers.o river_spreading.o
	$(FC) $(FFLAGS)  -o $@ $^ $(LDFLAGS)

clean: #cleans all the old compilation files
	@rm -f *.mod *.o *.exe *.dbg
