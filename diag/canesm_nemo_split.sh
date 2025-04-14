#########################################################
# nemo & time series
# N. Lambert, Jan 2025
#
#########################################################

set -e

# Note that nemo_rtd_mons used below is first month of the time chunk. 
# nemo_rtd_mons=1 for a run starting from January in a single 12-month chunk;
# nemo_rtd_mons=6 for a run starting from June in a single 12-month chunk;
# nemo_rtd_mons='1 7' for a run starting from January in two 6-month chunks;
# nemo_rtd_mons='6 12' for a run starting from June in two 6-month chunks.

# First and last month/year of 12-month period
  fmon=`echo $nemo_rtd_mons | cut -f1 -d' '`
  nmon=`echo $nemo_rtd_mons | wc -w` # number of chunks in 12-month period
  if [ $fmon -eq 1 ] ; then
    fyear=$year
    lyear=$year
    lmon=12
  else
    fyear=`echo $year | awk '{printf "%04d", $1 - 1}'`
    lyear=$year
    lmon=`echo $fmon | awk '{printf "%02d", $1 - 1}'`
  fi

# Access file containing grid information
  mask_mon=$(echo $nemo_rtd_mons | awk '{printf "%02d", $1}')  # get last element of nemo_rtd_mons, printed as 2 digit number
  orca_grid_info=mc_${runid}_${fyear}_m${mask_mon}_mesh_mask.nc
  [ -s orca_mesh_mask ] || access orca_mesh_mask $orca_grid_info || bail "Failed to access $orca_grid_info"

##################################################
# Split historical files to time series and save #
##################################################
  for sfx in $nemo_diag_file_suffix_list $nemo_diag_file_1y_suffix_list ; do
    yr=$fyear
    mp=0
    for mm in $nemo_rtd_mons ; do
      if [ $mm -lt $mp ] ; then
        # increment year by 1 if the current month is smaller than the previous month
        yr=`echo $yr | awk '{printf "%04d", $1 + 1}'`;
      fi
      diag_hist="mc_${runid}_${yr}_m${mm}_${sfx}.nc"
      # Access the history files
      access ${sfx}_${mm} $diag_hist 
      mp=$mm

      # split to time series
      # Get the list of variable (use "coordinates" as the key word because it won't apply to the coordinates variables, e.g. nav_lat)
      vars_list=$( (ncdump -h ${sfx}_${fmon} | grep coordinates | awk -F: '{print $1}') )
      [ -z "$vars_list" ] && sleep 30 && vars_list=$( (ncdump -h ${sfx}_${fmon} | grep coordinates | awk -F: '{print $1}') )
      [ -z "$vars_list" ] && bail 'error: fail to make var_lists'
      for var in $vars_list; do
          # produce one file per variables (3 times because the file system can failed)
        ncks -v $var ${sfx}_${mm} xxx-${sfx}_${var}.nc || \
            (sleep 5 && ncks -v $var ${sfx}_${mm} xxx-${sfx}_${var}.nc) || \
            (sleep 5 && ncks -v $var ${sfx}_${mm} xxx-${sfx}_${var}.nc)
      done
      # release the history files
      release  ${sfx}_${mm}
    done
  done

  # Save time series (now, at the end when all the time series are produced).
  tslist=`ls -1 xxx-*`
  for ts in $tslist ; do
     tssfx=`echo $ts |cut -f 2 -d '-' |sed 's/_ar6//'`
     save $ts sc_${runid}_${fyear}${fmon}_${year}${lmon}_${tssfx}
     release $ts
  done
  #

  # Save orca grid mask with consistent name as TS files
  save orca_mesh_mask sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_mesh_mask.nc

