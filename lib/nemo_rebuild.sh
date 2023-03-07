#!/bin/bash

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
        # not an error anymore because most of the files are re-tiled by NEMO+XIOS
        echo "No tiles files found for pattern $pfx"
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

# History file suffix and frequency lists. Must be the same lengths. I.E. repeat filenames for multiple freq.
# Get bash arrays with the freq and suffix list
nemo_hist_file_suffix_list_array=($nemo_rbld_hist_file_suffix_list)
nemo_hist_file_freq_list_array=($nemo_rbld_hist_file_freq_list)

# Check that the lists are the same length. Use this number to loop below..
n_suffix=${#nemo_hist_file_suffix_list_array[@]}
n_freq=${#nemo_hist_file_freq_list_array[@]}
if [ "$n_suffix" -ne "$n_freq" ]; then
   bail "ERROR: nemo_hist_file_suffix_list_array and nemo_hist_file_freq_list_array MUST have the same number of elements."
fi

# The rebuild executable must be accessable at run time and namelist files present in cwd
if [[ ! -f rebuild_nemo.exe ]]; then
  cp ${EXEC_STORAGE_DIR}/rebuild_nemo.exe . || bail "Unable to get rebuild_nemo.exe"
fi

# Can use Open MP. But probably only running on one processor.
export OMP_NUM_THREADS=2

# the tmpdir we are working in
wrkdir=$(pwd)
dir_del_list=""

# Access the restart directory 
# done ealier to have access to the namlist_cdf stored in $inrs
# Restart will be rebuilt later.
modellast="mc_${runid}_${yearlast}_m${monlast}";
inrs=${modellast}_nemors
access ${inrs}.tar ${inrs}.tar nocp=off na
if [ -s "${inrs}.tar" ]; then
  mkdir  ${inrs} 
  tar -xf ${inrs}.tar -C ${inrs}
  dir_del_list+=" ${inrs}.tar" #delete the tar and save the directory
else
  access $inrs $inrs nocp=on na #make a link to update the files in the directory 
  #dir_del_list+=" $inrs" # DON'T delete the directory, files inside updated
fi
[ -s "${inrs}" ]|| bail "Could not find ${inrs}"

#access the coordinates files (used for lat/lon later)
access coor.nc $nemo_coordinates  nocp=no  #force copy because we make temporary changes
ncrename -h -O -d t,time_counter coor.nc coor.nc || true #no error if already done
ncks -h -O -v e1.,e2.,nav_lon,nav_lat,glam.,gphi. coor.nc coor.nc
ncwa -h -O -a time_counter coor.nc coor.nc && ncks -h -O -x -v  time_counter coor.nc coor.nc

# remove the jstart if  ln_use_jatt is true in the namelist
if [ $( get_namelist_var ln_use_jattr $inrs/rs_namelist_cfg ) == ".true." ];then
  jstart=$( ncdump -h coor.nc | grep  --color=never "open_ocean_jstart\s*=" )
  jstart=${jstart##*=}
  jstart=$( trim_whitespace $jstart )
  jstart=$( echo $jstart | sed 's/[,;.]$//g' )
  # cut coor.nc according to open_ocean_jstart
  # remember that coor.nc is a temporary file
  ncks -h -O -d y,$(expr $jstart - 1), coor.nc coor.nc 
fi

# A list of directories to delete from RUNPATH at the end
if [ $nemo_save_hist == "on" ] ; then
   # Loop over the list of history files/freqs to rebuild
   for i in $(seq 0 $(($n_suffix-1))); do
      cd $wrkdir
      # Get the directory with the tiles, and extract them
      sfx=${nemo_hist_file_suffix_list_array[$i]}
      freq=${nemo_hist_file_freq_list_array[$i]}
      lsfx=$(echo "$sfx" | tr '[:upper:]' '[:lower:]')
      indir=${model1}_${freq}_${lsfx}
      access $indir $indir nocp=off na 
      if [ -d "$indir" ] ; then 
        # if don't exist, re-tile probably done by NEMO
        dir_del_list+=" $indir"
        cd $indir

        # Define the pattern, get the exe, do the rbld, and save.
        pfx=${runid}_${freq}_${start_date}_${stop_date}_$sfx
        ln -s ../rebuild_nemo.exe .
        rebuild_nemo_tiles
        save ${pfx}.nc $indir.nc

        # Move back up and cleanup
        cd $wrkdir
        rm -rf $indir
      fi
         # Replace the lat/lon to remove the hold made by the land processors elimination
        ncsave=${freq}_${lsfx}
        access  $ncsave.nc $indir.nc nocp=no na #force copy because we make temporary changes
      if [ -e "$ncsave.nc" ] ; then
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
        access $indir.nc $indir.nc na 
        delete  $indir.nc #delete and re-save to avoid new version
        save $ncsave.nc $indir.nc
        release $ncsave.nc
      fi
   done
fi

# Rebuild the mesh_mask file created by the model, if present
# Access the directory, if it is successfull, cd into it
indir=${model1}_mesh_mask
pfx=mesh_mask
ncsave=${model1}_${pfx}.nc
access $indir $indir nocp=off na
if [ -s "$indir" ] ; then
   cd $indir
   dir_del_list+=" $indir"

   # Define the pattern and do the rebld
   ln -s ../rebuild_nemo.exe .
   rebuild_nemo_tiles
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

# Rebuild the restart files. These will be saved alltogether, as
# is custom for NEMO rs' historically.

cp rebuild_nemo.exe ${inrs}/
cd $inrs
# Figure out the last time step, which is needed for the rs tile names.
ls -l  rs_time.step
nn_itend=$(cat rs_time.step)
end_step=$(echo $nn_itend | awk '{printf "%8.8d",$1}')

# The initial state files
pfx=output.init
# Check if the RS is already rebuilt, in which case do nothing.
if [ -s "${pfx}_0000.nc" ]; then
   rebuild_nemo_tiles
   # Replace the global lat/lon to remove the hold made by the land processors elimination
   ncks -x -h -O -v  nav_lon,nav_lat $pfx.nc $pfx.nc
   ncks -A -h -v nav_lon,nav_lat ${wrkdir}/coor.nc $pfx.nc
   ncsave=${model1}_istate_$start_date.nc
   mv  $pfx.nc $ncsave
fi

# The physics rs file
pfx=${runid}_${end_step}_restart

if [ ! -s "${pfx}_0000.nc" -a ! -e "${pfx}.nc" ]; then
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
   [ -z "$found_rs" ] || mv $found_rs $pfx.nc
fi

# Check if the RS is already rebuilt, in which case do nothing.
if [ -s "${pfx}_0000.nc" ]; then
   rebuild_nemo_tiles
   # Replace the global lat/lon to remove the hold made by the land processors elimination
   ncks -x -h -O -v  nav_lon,nav_lat $pfx.nc $pfx.nc
   ncks -A -h -v nav_lon,nav_lat ${wrkdir}/coor.nc $pfx.nc
fi

# The ice rs file
pfx=${runid}_${end_step}_restart_ice
if [ ! -s "${pfx}_0000.nc"  -a ! -e "${pfx}.nc" ]; then
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
   [ -z "$found_rs"  ] || mv $found_rs $pfx.nc
fi

if [ -s "${pfx}_0000.nc" ]; then
   rebuild_nemo_tiles
        # Replace the global lat/lon to remove the hold made by the land processors elimination
   ncks -x -h -O -v  nav_lon,nav_lat $pfx.nc $pfx.nc
   ncks -A -h -v nav_lon,nav_lat ${wrkdir}/coor.nc $pfx.nc
fi


# The trc rs file
pfx=${runid}_${end_step}_restart_trc
if [ ! -s "${pfx}_0000.nc"  -a ! -e "${pfx}.nc" ]; then
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
   [ -z "$found_rs" ] || mv $found_rs $pfx.nc
fi

if [ -s "${pfx}_0000.nc" ]; then
   rebuild_nemo_tiles
        # Replace the global lat/lon to remove the hold made by the land processors elimination
   ncks -x -h -O -v  nav_lon,nav_lat $pfx.nc $pfx.nc
   ncks -A -h -v nav_lon,nav_lat ${wrkdir}/coor.nc $pfx.nc
fi
rm ${wrkdir}/coor.nc 

cd $wrkdir

#save the new untar ${inrs} (if already untar, file in direcotory are just kept)
 [ -e "${inrs}.tar" ] && ( save ${inrs} ${inrs} || bail "Could not save ${inrs}" )

# since everything has gone successfully, cleanup tile directories from RUNPATH
mkdir cleanup
cd cleanup
for fil in $dir_del_list; do
   access xxx-to-del $fil
   delete xxx-to-del
done
cd $wrkdir
rm -rf cleanup
