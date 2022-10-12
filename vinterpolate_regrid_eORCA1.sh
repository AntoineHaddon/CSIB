#!/bin/bash
# O Riche July 7th 2022
# moving onto NO3 and other external sources
# O Riche May 26th 2022
# interpolate vertically to 75 levels in a batch
set -e
# set -x

gridtarget=urdy_eorca1_en4_v1.1.1995_2014.monthlymean_eorca1t_nemo_l75_jstart.nc.001

parname=(dic alkalini o2 no3 si po4 doc fer epsdb bathy)
# parname=(o2)
for ((pp=0; pp<${#parname[@]}; pp++))
do
  echo \# nco/cdo command lines for ${parname[$pp]}
	echo cdo intlevelx,0.50576,1.555855,2.667682,3.85628,5.140361,6.543034,8.092519,9.82275,11.77368,13.99104,16.52532,19.4298,22.75762,26.5583,30.87456,35.7402,41.18002,47.21189,53.85064,61.11284,69.02168,77.61116,86.92943,97.04131,108.0303,120,133.0758,147.4062,163.1645,180.5499,199.79,221.1412,244.8906,271.3564,300.8875,333.8628,370.6885,411.7939,457.6256,508.6399,565.2923,628.026,697.2587,773.3683,856.679,947.4479,1045.854,1151.991,1265.861,1387.377,1516.364,1652.568,1795.671,1945.296,2101.027,2262.422,2429.025,2600.38,2776.039,2955.57,3138.565,3324.641,3513.446,3704.657,3897.982,4093.159,4289.953,4488.155,4687.581,4888.07,5089.479,5291.683,5494.575,5698.061,5902.058 -remapnn,$gridtarget uncs_orca1_data_${parname[$pp]}_nomask.nc.001 uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
	echo ncatted -O -a units,deptht,d,c,""     uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  echo ncatted -O -a axis,deptht,d,c,""      uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  echo ncatted -O -a long_name,deptht,d,c,"" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
	echo ncatted -O -a open_ocean_jstart,global,a,d,41 uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
  echo ncrename -d time_counter,time uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  echo ncrename -v time_counter,time uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  echo ncatted -O -a standard_name,time,o,c,\"time\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  echo ncatted -O -a calendar,time,o,c,\"360_day\"   uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  echo ncatted -O -a axis,time,o,c,\"T\"             uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  # if the variable is a monthly mean
  case ${parname[$pp]} in
    o2|si|no3|po4)
      echo ncatted -O -a units,time,o,c,"days" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
      echo cp uror_nemo4_eorca1_${parname[$pp]}_nomask.nc tmp_orca1_00.nc\;
      echo rm uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
      echo "ncap2 -s 'time(:)={15,45,75,105,135,165,195,225,255,285,315,345}' tmp_orca1_00.nc uror_nemo4_eorca1_${parname[$pp]}_nomask.nc;"
      echo rm tmp_orca1_00.nc\;
      ;;
    *)  # anything else do nothing
      : # null command is a colon in bash
      ;;
  esac
  echo save uror_nemo4_eorca1_${parname[$pp]}_nomask.nc uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  echo \#

done

set +e
# set +x


# O Riche Aug 02nd 2022
#### Note that I had to prep hydrofe_orca2.nc aka uncs*epsdb*.nc
#### to be process-able, if that's recuring I will have to add more case statements:

# ncrename -d z,deptht hydrofe_orca2.nc ;
# ncatted -O -a standard_name,nav_lon,o,c,"longitude" hydrofe_orca2.nc ;
# ncatted -O -a long_name,nav_lon,o,c,"Longitude"     hydrofe_orca2.nc ;
# ncatted -O -a units,nav_lon,o,c,"degrees_east"      hydrofe_orca2.nc ;
# ncatted -O -a _CoordinateAxisType,nav_lon,o,c,"Lon" hydrofe_orca2.nc ;
# ncatted -O -a standard_name,nav_lat,o,c,"latitude"  hydrofe_orca2.nc ;
# ncatted -O -a long_name,nav_lat,o,c,"Latitude"      hydrofe_orca2.nc ;
# ncatted -O -a units,nav_lat,o,c,"degrees_north"     hydrofe_orca2.nc ;
# ncatted -O -a _CoordinateAxisType,nav_lat,o,c,"lat" hydrofe_orca2.nc ;
# ncatted -O -a coordinates,epsdb,o,c,"nav_lon nav_lat" hydrofe_orca2.nc ;
