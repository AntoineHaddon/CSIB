#!/bin/bash
# O Riche July 19th 2022
# remove files from the database
# O Riche July 7th 2022
# moving onto NO3 and other external sources
# O Riche May 26th 2022
# interpolate vertically to 75 levels in a batch
# set -e
# set -x

source /home/ror001/bin/nemo_env_setup_file

parname=(dic alkalini o2 no3 si po4 doc fer epsdb bathy)

for ((pp=0; pp<${#parname[@]}; pp++))
do
  echo \# nco/cdo command lines for ${parname[$pp]}
	echo release uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  echo fdb delete uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  echo rm -f ../../echo fdb delete uror_nemo4_eorca1_${parname[$pp]}_nomask.nc.001\;
  echo \#
done

# set +e
# set +x
