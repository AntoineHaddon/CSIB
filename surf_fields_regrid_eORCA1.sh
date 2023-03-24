#!/bin/bash
# O Riche July 7th 2022
# moving onto NO3 and other external sources
# O Riche May 26th 2022
# interpolate vertically to 75 levels in a batch
# set -x
set -e

coordorca2=surf_river_orca2.nc #### used as a reference file to copy nav_lat nav_lon if necessary

gridtarget=urdy_eorca1_en4_v1.1.1995_2014.monthlymean_eorca1t_nemo_l75_jstart.nc.001

parname=( surf_chla                # chl-a
          surf_dust                # dust
          surf_par                 # PAR
          surf_femask              # iron mask
          surf_ndepo               # N atm depo
          surf_rivercarbon         # river carbon
          surf_riverpisces         # river PISCES source
          surf_fesolub)            # iron solub.
          
filename=(surf_chlorophyll_orca1    # chl-a
          surf_dust_orca1           # dust
          surf_par_orca2            # PAR
          surf_femask_orca1         # iron mask
          surf_ndepo_orca1          # N atm depo
          surf_river_orca2          # river carbon
          surf_river_orca2_pisces   # river PISCES source
          surf_fesolub_orca2)       # iron solub.
          
varname=( CHLA                    # chl-a
          dust                    # dust
          fr_par                  # PAR
          femask                  # iron mask
          ndep                    # N atm depo
          river                   # river carbon
          river2                  # river PISCES source
          solub)                  # iron solub.
          
varunts=( "mg Chl-a m^-3"                # chl-a
          "-"                            # dust
          "-"                            # PAR
          "-"                            # iron mask
          "mmol N m^-2 y^-1?"            # N atm depo
          "mmol ? m^-3"                  # river carbon
          "mmol ? m^-3"                  # river PISCES source
          "nmol Fe m^-3 uatm^-1?")       # iron solub.

for ((pp=0; pp<${#parname[@]}; pp++))
do
  echo \# nco/cdo command lines for ${parname[$pp]}
  echo cp ${filename[$pp]}.nc uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\; 
	echo ncatted -O -a open_ocean_jstart,global,a,d,41       uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
  case "${varname[$pp]}" in
    river2)
      echo ncks -A -v nav_lat,nav_lon ${coordorca2}           uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
      ;;
    *)
      :
      ;;
  esac
  echo ncatted -O -a standard_name,nav_lon,o,c,"longitude" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
	echo ncatted -O -a long_name,nav_lon,o,c,"Longitude"     uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
	echo ncatted -O -a _CoordinateAxisType,nav_lon,o,c,"Lon" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
  echo ncatted -O -a standard_name,nav_lat,o,c,"latitude"  uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
	echo ncatted -O -a long_name,nav_lat,o,c,"Latitude"      uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
	echo ncatted -O -a _CoordinateAxisType,nav_lat,o,c,"lat" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
  case "${varname[$pp]}" in
  solub)
    for varname0 in solubility1 solubility2
    do
      echo ncatted -O -a coordinates,${varname0},o,c,\"nav_lon nav_lat\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
      echo ncatted -O -a units,${varname0},o,c,\"nmol-Fe/L/uamt?\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   
    done
    ;;
  river)
    for varname0 in riverdic riverdoc riverpoc
    do
      echo ncatted -O -a coordinates,${varname0},o,c,\"nav_lon nav_lat\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
      echo ncatted -O -a units,${varname0},o,c,\"mmol-C/L\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   
    done
    ;;
  river2)
    vars_list=(riverdin riverdip riverdon riverdop riverdic riverdoc riverdsi)
    crnc_list=(N P N P C C Si)
    for ((ii = 0 ; ii < ${#vars_list[@]} ; ii++ ))
    do
      echo ncatted -O -a coordinates,${vars_list[$ii]},o,c,\"nav_lon nav_lat\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
      echo ncatted -O -a units,${vars_list[$ii]},o,c,\"mmol-${crnc_list[$ii]}/L\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   
    done
    ;;
  *)
    echo ncatted -O -a coordinates,${varname[$pp]},o,c,\"nav_lon nav_lat\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
    echo ncatted -O -a units,${varname[$pp]},o,c,\"${varunts[$pp]}\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   
    ;;
  esac
#### time_counter to time
  echo ncrename -d time_counter,time uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  echo ncrename -v time_counter,time uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;

  echo ncatted -O -a standard_name,time,o,c,\"time\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
  echo ncatted -O -a units,time,o,c,\"days\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
  echo ncatted -O -a calendar,time,o,c,\"360_day\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
  echo ncatted -O -a axis,time,o,c,\"T\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
  echo ncatted -O -a long_name,time,o,c,\"Time Axis\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly

#### keep time_counter
 
  # echo ncatted -O -a standard_name,time_counter,o,c,\"time\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
  # echo ncatted -O -a units,time_counter,o,c,\"days\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
  # echo ncatted -O -a calendar,time_counter,o,c,\"360_day\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
  # echo ncatted -O -a axis,time_counter,o,c,\"T\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly
  # echo ncatted -O -a long_name,time_counter,o,c,\"Time Axis\" uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;   #### need this attribute to work properly

  echo mv uror_nemo4_eorca1_${parname[$pp]}_nomask.nc surf_tmpry_orca1_00.nc\;
  echo cdo -remapnn,$gridtarget surf_tmpry_orca1_00.nc uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  # echo "ncap2 -O -s 'defdim(\"deptht\",0.50576); "${varname[$pp]}"[time,y,x,deptht]="${varname[$pp]}"' surf_tmpry_orca1_00.nc uror_nemo4_eorca1_${parname[$pp]}_nomask.nc;"

  case "${parname[$pp]}" in
    surf_par)
      echo "ncap2 -O -s 'time(:)={ " `seq -s "," 1 365` " }' uror_nemo4_eorca1_${parname[$pp]}_nomask.nc surf_tmpry_orca1_00.nc;"
      # echo "ncap2 -O -s 'time_counter(:)={ " `seq -s "," 1 365` " }' uror_nemo4_eorca1_${parname[$pp]}_nomask.nc surf_tmpry_orca1_00.nc;"
      ;;
    surf_chla|surf_dust)  
      echo "ncap2 -O -s 'time(:)={15,45,75,105,135,165,195,225,255,285,315,345}' uror_nemo4_eorca1_${parname[$pp]}_nomask.nc surf_tmpry_orca1_00.nc;"
      # echo "ncap2 -O -s 'time_counter(:)={15,45,75,105,135,165,195,225,255,285,315,345}' uror_nemo4_eorca1_${parname[$pp]}_nomask.nc surf_tmpry_orca1_00.nc;"
      ;;
    *)
      echo mv uror_nemo4_eorca1_${parname[$pp]}_nomask.nc surf_tmpry_orca1_00.nc\;
      ;;
  esac

  echo mv surf_tmpry_orca1_00.nc uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;
  # echo rm surf_tmpry_orca1_00.nc
  
  echo save uror_nemo4_eorca1_${parname[$pp]}_nomask.nc uror_nemo4_eorca1_${parname[$pp]}_nomask.nc\;

  echo \#
 
done

set +e
# set +x