#!/bin/bash

# An example for compilation and execution of nemo offline diagnostics 
# including time mean (1d->1m) and time series generation
# Sep/2018, D. Yang

# netcdf? Does this do anything?
. ssmuse-sh -d /fs/ssm/hpco/tmp/eccc/201402/04/intel-2016.1.150

# ifort
. ssmuse-sh -d /fs/ssm/hpco/tmp/eccc/201402/03/base  -d main/opt/intelcomp/intelcomp-2016.1.156

rm -f nemo_diag_cmoc.exe

ifort -o nemo_diag_cmoc.exe nemo_diag_glovars_cmoc.F90 nemo_diag_cal_cmoc.F90 nemo_diag_cmoc.F90 uvic_netcdf.f -I /fs/ssm/hpco/tmp/eccc/201402/04/intel-2016.1.150/ubuntu-14.04-amd64-64/include/ -L /fs/ssm/hpco/tmp/eccc/201402/04/intel-2016.1.150/ubuntu-14.04-amd64-64/lib/ -lnetcdf -lnetcdff

rm -f *.mod

#rm -f solproda.nc solprodc.nc o2sol.nc

./nemo_diag_cmoc.exe

exit 0
