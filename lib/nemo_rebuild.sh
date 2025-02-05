#!/bin/bash
set -e
#~~~~~~~~~~~~~~~
# Function Defs
#~~~~~~~~~~~~~~~
source $CANESM_SRC_ROOT/CCCma_tools/tools/CanESM_shell_functions.sh
function rebuild_nemo_tiles(){
   # Function to rebuild NEMO tiles, for files matching
   # the pattern in $pfx. The number of tiles is specified
   # as $ndomain. It is assumed the rebuild_nemo.exe is locally
   # available.

   fnpatt=${pfx}_0000.nc
   ndomain=`ls -l ${pfx}_[0-9]*.nc |wc -l`
   if [ -s "$fnpatt" ]; then
      # Create a namelist file for use with this rebuild
      rbnl_file=nam_rebuild
      echo "&nam_rebuild"           > $rbnl_file
      echo "    filebase='${pfx}'" >> $rbnl_file
      echo "    ndomain=$ndomain"  >> $rbnl_file
      echo "/"                     >> $rbnl_file

      # Do the rbld
      [ -s rebuild_nemo.exe ] && ./rebuild_nemo.exe || bail "rebuild_nemo.exe not available"
      [ -s ${pfx}.nc ] || bail "Failed to rebuild file ${pfx}.nc"
  
      # preserve time stamp from tiles
      touch -r $fnpatt ${pfx}.nc

      # remove individual tiles after rbld
      rm -f ${pfx}_[0-9][0-9][0-9][0-9].nc
    else
        bail "No files found for pattern $pfx"
    fi
}

#~~~~~~~~~~~~~
# Main Script
#~~~~~~~~~~~~~
# Figure out NEMO dates. We are operating 1 month at a time in coupled mode.
# Determine run start and stop dates in YYYYMMDD format

# Currently start_day is always 1
start_day=1

# Determine the number of days in the last month of the run
monlast=$(echo $mon $months_gcm | awk '{printf "%02d",$1+$2-1}')
if [ $monlast -gt 12 ] ; then
  monlast=$(echo $monlast | awk '{printf "%02d",$1-12}')
  yearlast=$(echo $year | awk '{printf "%04d",$1+1}')
else
  yearlast=$year
fi
stop_day=$( days_in_month $monlast )

# First and last calendar date YYYYMMDD determined from user input start/stop
start_date=$(echo $year $mon $start_day|awk '{printf "%04d%02d%02d",$1,$2,$3}' -)
stop_date=$(echo $yearlast $monlast $stop_day|awk '{printf "%04d%02d%02d",$1,$2,$3}' -)

# History file suffix and frequency lists
# Get bash arrays with the freq and suffix list
nemo_hist_file_suffix_list_array=($nemo_rbld_hist_file_suffix_list)
n_suffix=${#nemo_hist_file_suffix_list_array[@]}

nemo_file_suffixes_array=()
nemo_file_freqs_array=()
for i in $(seq 0 $((n_suffix-1))); do
    fs=${nemo_hist_file_suffix_list_array[$i]}
    IFS='_' read -r freq param <<< $fs
    if [[ $freq =~ ^[0-9]+[hdmyt] ]]; then
        nemo_file_suffixes_array+=("$param")
        nemo_file_freqs_array+=("$freq")
    else
        echo "file list entry does not match required format {\$freq_\$suffix}; skipping: $a"
    fi
done

# The rebuild executable must be accessable at run time and namelist files present in cwd
cp ${EXEC_STORAGE_DIR}/rebuild_nemo.exe . || bail "Unable to get rebuild_nemo.exe"

# Can use Open MP. But probably only running on one processor.
export OMP_NUM_THREADS=2

# the tmpdir we are working in
wrkdir=$(pwd)

#access the coordinates files (used for lat/lon later)
access coor.nc $nemo_coordinates  nocp=no  #force copy because we make temporary changes
ncrename -h -O -d t,time_counter coor.nc coor.nc || true #no error if already done
ncks -h -O -v e1.,e2.,nav_lon,nav_lat,glam.,gphi. coor.nc coor.nc
ncwa -h -O -a time_counter coor.nc coor.nc && ncks -h -O -x -v  time_counter coor.nc coor.nc


# A list of directories to delete from RUNPATH at the end
if (( with_rbld_nemo == 1 )) ; then
   # Loop over the list of history files/freqs to rebuild
   for i in $(seq 0 $(($n_suffix-1))); do
      cd $wrkdir
      # Get the directory with the tiles, and extract them
      sfx=${nemo_file_suffixes_array[$i]}
      freq=${nemo_file_freqs_array[$i]}
      lsfx=$(echo "$sfx" | tr '[:upper:]' '[:lower:]')
      indir=${model1}_${freq}_${lsfx}
      access $indir $indir nocp=off na
      if [ -d "$indir" ] ; then 
         # if don't exist, re-tile probably done by NEMO
         cd $indir

         # Define the pattern, get the exe, do the rbld, and save.
         pfx=${runid}_${freq}_${start_date}_${stop_date}_$sfx
         ln -s ../rebuild_nemo.exe .
         rebuild_nemo_tiles
         ncsave=${model1}_${freq}_${sfx}.nc
         save ${pfx}.nc $ncsave

         # Move back up and cleanup
         cd $wrkdir
         rm -rf $indir
      fi
               # Replace the lat/lon to remove the hold made by the land processors elimination
      ncsave=${freq}_${lsfx}
      access  $ncsave.nc $indir.nc na 
      if [ -e "$ncsave.nc" ] ; then
        chmod u+w $(readlink -f "$ncsave.nc")
        # detect the grid (U/V/F/T) with the suffix
        if [[ ${sfx,,} == *"grid_u"*  ]];then
                  ( ncks -A -h -v glamu,gphiu coor.nc $ncsave.nc && 
                    ncap2 -h -O -s "nav_lon=glamu;nav_lat=gphiu"  $ncsave.nc  $ncsave.nc )
        elif [[ ${sfx,,} == *"grid_v"*  ]];then
                  ( ncks -A -h -v glamv,gphiv coor.nc $ncsave.nc && 
                    ncap2 -h -O -s "nav_lon=glamv;nav_lat=gphiv"  $ncsave.nc  $ncsave.nc )
        elif [[ ${sfx,,} == *"grid_f"*  ]];then
                  ( ncks -A -h -v glamf,gphif coor.nc $ncsave.nc && 
                    ncap2 -h -O -s "nav_lon=glamf;nav_lat=gphif"  $ncsave.nc  $ncsave.nc )
        elif [[ ${sfx,,} == *"diaptr"*  ]];then
                  (  release $ncsave.nc &&
                   continue )
        else # grid T is the default 
                  ( ncks -A -h -v glamt,gphit coor.nc $ncsave.nc && 
                    ncap2 -h -O -s "nav_lon=glamt;nav_lat=gphit"  $ncsave.nc  $ncsave.nc )
        fi
        ncks -h -O -x -v gphi.,glam.  $ncsave.nc  $ncsave.nc
        chmod u-w $(readlink -f "$ncsave.nc")
        release $ncsave.nc
      fi

   done
fi

# Rebuild the mesh_mask file created by the model, if present
# Access the directory, if it is successfull, cd into it
indir=${model1}_mesh_mask
access $indir $indir nocp=off na
if [ -s "$indir" ] ; then
   cd $indir

   # Define the pattern and do the rebld
   pfx=mesh_mask
   ln -s ../rebuild_nemo.exe .
   rebuild_nemo_tiles
   ncsave=${model1}_${pfx}.nc
        # Replace the global lat/lon to remove the hold made by the land processors elimination
   ncks -x -h -O -v  nav_lon,nav_lat,glamf,gphif,glamv,gphiv,glamu,gphiu,glamt,gphit ${pfx}.nc ${pfx}.nc 
   ncks -A -h -v nav_lon,nav_lat,glamf,gphif,glamv,gphiv,glamu,gphiu,glamt,gphit  ${wrkdir}/coor.nc ${pfx}.nc  
   ncks -x -h -O -v  e1f,e2f,e1v,e2v,e1u,e2u,e1t,e2t ${pfx}.nc ${pfx}.nc 
   ncks -A -h -v e1f,e2f,e1v,e2v,e1u,e2u,e1t,e2t  ${wrkdir}/coor.nc ${pfx}.nc  
   save ${pfx}.nc $ncsave

   # cleanup
   cd $wrkdir
   rm -rf $indir
fi

# Rebuild the restart files. These will be tarred below and saved alltogether, as
# is custom for NEMO rs' historically.

# Determine the restart to use
#   - note that we expect to use tiled_nemors for most executions,
#       but at the beginning of a run, we may have a restart that has tiled files
#       but is still using the nemors suffix
modellast="mc_${runid}_${yearlast}_m${monlast}";
if fdb exists ${modellast}_tiled_nemors; then
    inrs=${modellast}_tiled_nemors
else
    inrs=${modellast}_nemors
fi
outrs=${modellast}_nemors

# Access the restart directory, and cd into it.
#   - we add "in_*" to the input directory to differentiate
#       between input/output restarts incase we have a tiled restart
#       that doesn't have the new "tiled_nemors" suffix,
#       which would result in $inrs=$outrs
access in_${inrs} $inrs nocp=off
cd in_${inrs}
ln -sf ../rebuild_nemo.exe .
# Figure out the last time step, which is needed for the rs tile names.
nn_itend=$(cat rs_time.step)
start_step=$(grep -m 1 -w nn_it000 rs_namelist_cfg | awk '{printf "%8.8d",$3 - 1}')
end_step=$(echo $nn_itend | awk '{printf "%8.8d",$1}')

# The initial ice state files
pfx=output.init_ice
# Check if the RS is already rebuilt, in which case do nothing.
if [ -s "${pfx}_0000.nc" ]; then
   rebuild_nemo_tiles
   # Replace the global lat/lon to remove the hold made by the land processors elimination
   ncks -x -h -O -v  nav_lon,nav_lat $pfx.nc $pfx.nc
   ncks -A -h -v nav_lon,nav_lat ${wrkdir}/coor.nc $pfx.nc
   ncsave=${runid}_${start_step}_initial_ice.nc
   mv  $pfx.nc $ncsave
fi

# The physics rs file
pfx=${runid}_${end_step}_restart

if [ ! -s "${pfx}_0000.nc" ]; then
   # Look for files, which might not have the same name as the run
   found_rs=`(ls -1 *_restart_[0-9][0-9][0-9][0-9].nc || : ) 2>/dev/null`
   # Rename them
   for rsfile in $found_rs; do
      sfx=`echo $rsfile|sed 's/^.*\(_[0-9][0-9][0-9][0-9].nc\).*$/\1/'`
      [ "x$sfx" = "x$rsfile" ] && sfx=`echo $rsfile|sed 's/^.*\(\.nc\).*$/\1/'`
      mv $rsfile $pfx$sfx
    done
   # an already rebuilt rs with a different name
   found_rs=`(ls -1 *_restart.nc || : ) 2>/dev/null`
   [ -z "$found_rs" ] || mv -n $found_rs $pfx.nc
fi

# Check if the RS is already rebuilt, in which case do nothing.
fnpatt=${pfx}_0000.nc
if [ -s "$fnpatt" ]; then
   rebuild_nemo_tiles
fi

# The ice rs file
pfx=${runid}_${end_step}_restart_ice
if [ ! -s "${pfx}_0000.nc" ]; then
   # Look for files, which might not have the same name as the run
   found_rs=`(ls -1 *_restart_ice_[0-9][0-9][0-9][0-9].nc || : ) 2>/dev/null`
   # Rename them
   for rsfile in $found_rs; do
      sfx=`echo $rsfile|sed 's/^.*\(_[0-9][0-9][0-9][0-9].nc\).*$/\1/'`
      [ "x$sfx" = "x$rsfile" ] && sfx=`echo $rsfile|sed 's/^.*\(\.nc\).*$/\1/'`
      mv $rsfile $pfx$sfx
    done
   # an already rebuilt rs with a different name
   found_rs=`(ls -1 *_restart_ice.nc || : ) 2>/dev/null`
   [ -z "$found_rs" ] || mv -n $found_rs $pfx.nc
fi

fnpatt=${pfx}_0000.nc
if [ -s "$fnpatt" ]; then
   rebuild_nemo_tiles
fi

# The trc rs file
pfx=${runid}_${end_step}_restart_trc
if [ ! -s "${pfx}_0000.nc" ]; then
   # Look for files, which might not have the same name as the run
   found_rs=`(ls -1 *_restart_trc_[0-9][0-9][0-9][0-9].nc || : ) 2>/dev/null`
   # Rename them
   for rsfile in $found_rs; do
      sfx=`echo $rsfile|sed 's/^.*\(_[0-9][0-9][0-9][0-9].nc\).*$/\1/'`
      [ "x$sfx" = "x$rsfile" ] && sfx=`echo $rsfile|sed 's/^.*\(\.nc\).*$/\1/'`
      mv $rsfile $pfx$sfx
    done
   # an already rebuilt rs with a different name
   found_rs=`(ls -1 *_restart_trc.nc || : ) 2>/dev/null`
   [ -z "$found_rs" ] || mv -n $found_rs $pfx.nc
fi

fnpatt=${pfx}_0000.nc
if [ -s "$fnpatt" ]; then
   rebuild_nemo_tiles
fi

# The physics init file
pfx=output.init
# Check if the init files are present 
fnpatt=${pfx}_0000.nc
if [ -s "$fnpatt" ]; then
   rebuild_nemo_tiles
   mv $pfx.nc ${runid}_${start_step}_initial.nc
fi

# The trc init file
pfx=output_trc.init
# Check if the init files are present 
fnpatt=${pfx}_0000.nc
if [ -s "$fnpatt" ]; then
   rebuild_nemo_tiles
   mv $pfx.nc ${runid}_${start_step}_initial_trc.nc
fi

release rebuild_nemo.exe $rbnl_file
cd $wrkdir

# since rebuild has gone successfully, cleanup tile directories from RUNPATH,
#   removing the input restart (inrs) if inrs==outrs (which should only happen
#   for the initial restart)
if [[ ${inrs} == ${outrs} ]]; then
   fdb mdelete $inrs
fi

# Finally, save new directory with the rebuilt files
mkdir out_${outrs}
mv in_${inrs}/* out_${outrs}/
save out_${outrs} ${outrs}
rm -rf in_${inrs} out_${outrs}
