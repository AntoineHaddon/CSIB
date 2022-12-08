#!/bin/bash
# O Riche July 28th 2022
# Extract the surface map
set -e
parname=(no3)

# for ((pp=0; pp<1 pp++))

for ((pp=0; pp<${#parname[@]}; pp++))
do
  echo \# extract surface ${parname[$pp]}
  echo ncks -F -O -d deptht,1,1 uror_nemo4_eorca1_${parname[$pp]}_nomask.nc -o surf_${parname[$pp]}_orca1.nc\;
  echo ncwa -O -a deptht surf_${parname[$pp]}_orca1.nc surf_${parname[$pp]}_orca1.nc\;
  echo ncks -C -O -x -v deptht surf_${parname[$pp]}_orca1.nc surf_${parname[$pp]}_orca1.nc\;
  echo ncpdq -O surf_${parname[$pp]}_orca1.nc surf_${parname[$pp]}_orca1.nc\;
  echo cp surf_${parname[$pp]}_orca1.nc uror_nemo4_eorca1_surf_${parname[$pp]}_nomask.nc\;
  echo save uror_nemo4_eorca1_surf_${parname[$pp]}_nomask.nc uror_nemo4_eorca1_surf_${parname[$pp]}_nomask.nc\;
  echo \#
done

set +e