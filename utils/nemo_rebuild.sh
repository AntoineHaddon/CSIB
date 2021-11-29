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

# A list of directories to delete from RUNPATH at the end
dir_del_list=""
if (( canesm_nemo_rbld_save_hist == 1 )) ; then
   # Loop over the list of history files/freqs to rebuild
   for i in $(seq 0 $(($n_suffix-1))); do
      cd $wrkdir
      # Get the directory with the tiles, and extract them
      sfx=${nemo_hist_file_suffix_list_array[$i]}
      freq=${nemo_hist_file_freq_list_array[$i]}
      lsfx=$(echo "$sfx" | tr '[:upper:]' '[:lower:]')
      indir=${model1}_${freq}_${lsfx}
      access $indir $indir nocp=off
      dir_del_list+=" $indir"
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
   done
fi

# Rebuild the mesh_mask file created by the model, if present
# Access the directory, if it is successfull, cd into it
indir=${model1}_mesh_mask
access $indir $indir nocp=off na
if [ -s "$indir" ] ; then
   cd $indir
   dir_del_list+=" $indir"

   # Define the pattern and do the rebld
   pfx=mesh_mask
   ln -s ../rebuild_nemo.exe .
   rebuild_nemo_tiles
   ncsave=${model1}_${pfx}.nc
   save ${pfx}.nc $ncsave

   # cleanup
   cd $wrkdir
   rm -rf $indir
fi

# Rebuild the restart files. These will be tarred below and saved alltogether, as
# is custom for NEMO rs' historically.

# Access the restart directory, and cd into it.
modellast="mc_${runid}_${yearlast}_m${monlast}";
inrs=${modellast}_nemors
access $inrs $inrs nocp=off
dir_del_list+=" $inrs"
cd $inrs
ln -s ../rebuild_nemo.exe .
# Figure out the last time step, which is needed for the rs tile names.
nn_itend=$(cat rs_time.step)
end_step=$(echo $nn_itend | awk '{printf "%8.8d",$1}')

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
   [ -z "$found_rs" ] || mv $found_rs $pfx.nc
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
   [ -z "$found_rs" ] || mv $found_rs $pfx.nc
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
   [ -z "$found_rs" ] || mv $found_rs $pfx.nc
fi

fnpatt=${pfx}_0000.nc
if [ -s "$fnpatt" ]; then
   rebuild_nemo_tiles
fi

# Create the tar archive for the nemors and save it.
release rebuild_nemo.exe $rbnl_file
tar -cf ${inrs}.tar *
# preserve time stamp from restart.nc
touch -r ${runid}_${end_step}_restart.nc $fnpatt ${inrs}.tar
save ${inrs}.tar ${inrs}.tar || bail "Could not save ${indir}.tar"
cd $wrkdir
rm -rf $inrs

# since everything has gone successfully, cleanup tile directories from RUNPATH
mkdir cleanup
cd cleanup
for fil in $dir_del_list; do
   access xxx-to-del $fil
   delete xxx-to-del
done
cd $wrkdir
rm -rf cleanup
