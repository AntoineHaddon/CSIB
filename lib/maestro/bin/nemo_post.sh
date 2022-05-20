#!/bin/bash
#=========================================================
# This script contains two functions:
# * Rebuild the tiled history and restart files.
# * Create the restart tar archive.
# 
# This script is peeled off from nemo_post.tsk and adapted 
# to facilitate migrations (D. Yang, Sep 2020).
# 
# This script is adapted for NEMO4.0.3 (D. Yang, OCT 2020)
#=========================================================
 
# Get input files and executables
ln -s ${task_input_modfiles_dir}/* .
cp ${storage_dir}/executables/rebuild_nemo_mpi.exe .
acc_cp mpi_rebuild ${seq_exp_home_dir}/bin/mpi_rebuild

# Establish a hist file list for moving to the output directory 

  # Get bash arrays with the freq and suffix list
  nemo_hist_file_suffix_list_array=($nemo_hist_file_suffix_list)
  nemo_hist_file_freq_list_array=($nemo_hist_file_freq_list)

  # Check that the lists are the same length
  n_suffix=${#nemo_hist_file_suffix_list_array[@]}
  n_freq=${#nemo_hist_file_freq_list_array[@]}
  if [ "$n_suffix" -ne "$n_freq" ]; then
      bail "ERROR: nemo_hist_file_suffix_list_array and nemo_hist_file_freq_list_array MUST have the same number of elements."
  fi

  # Loop through each element, and setup a rebuild namelist.
  histfile_list=""
  if [ $nemo_save_hist -eq 1 ] || [ $nemo_rtd -eq 1 ]; then
      for i in $(seq 0 $(($n_suffix-1))); do
          echo "${nemo_hist_file_suffix_list_array[$i]}"
          echo "${nemo_hist_file_freq_list_array[$i]}"
          histfile_list="$histfile_list ${cn_exp}_${nemo_hist_file_freq_list_array[$i]}_${chunk_start_date}_${chunk_end_date}_${nemo_hist_file_suffix_list_array[$i]}"
      done
  fi

# If used, the below section does the parallel rebuild of tiles
if [ "0" -eq "0" ] ; then
  rebuild_list=""
  n_rebuild=0
  # Depending on parameter choices, add relevant restarts to rebuild_list
  end_step=${chunk_nn_itend}
  if [ $pisces_offline -eq 0 ]; then
      rebuild_list="$rebuild_list ${cn_exp}_${end_step}_restart"
      n_rebuild=$(($n_rebuild + 1))
      [ -s ${cn_exp}_icebergs_${end_step}_restart ] && rebuild_list="$rebuild_list ${cn_exp}_icebergs_${end_step}_restart"
      [ -s ${cn_exp}_icebergs_${end_step}_restart ] && n_rebuild=$(($n_rebuild + 1))
      if [ $nn_ice -eq 2 ]; then
          rebuild_list="$rebuild_list ${cn_exp}_${end_step}_restart_ice"
          n_rebuild=$(($n_rebuild + 1))
      fi
  fi
  if [ $nemo_trc -eq 1 ] || [ $nemo_carbon -eq 1 ]; then
      rebuild_list="$rebuild_list ${cn_exp}_${end_step}_restart_trc"
      n_rebuild=$(($n_rebuild + 1))
  fi
  # setup the namelists. The script mpi_rebuild makes the namelists, based on the input of "rebuild_list". This script is adjusted from the
  # default nemo version of rebuild, and produces mulitple namelist files (one for each entry in "rebuild_list" - which corresponds to
  # each mpi task in the parallel rebuild, which gets its own input. (There is one file-type/task [e.g. grid_T] and one namelist/task).
  ${task_work_dir}/mpi_rebuild "${rebuild_list}" $jpnij
  
  # run the parallel rebuild
  export OMP_NUM_THREADS=5
  #${OMP_NUM_THREADS:-3}
  # The total number for "n" must be <=36 to fit on one node of the XC40
  # This could be tuned, and defitely changed for other machines.
  while [ $(($OMP_NUM_THREADS * $n_rebuild)) -gt 36 ]; do
      OMP_NUM_THREADS=$(($OMP_NUM_THREADS -1))
  done

  # Parallel rebuild execution
  mpirun -n $n_rebuild -env OMP_NUM_THREADS ${OMP_NUM_THREADS} ${task_work_dir}/rebuild_nemo_mpi.exe
  wait
fi

if [ $with_nemo_diag -eq 1 ]; then
   # rebuild 1d diaptr files
   rebuild_list2=""
   # number of tiles 
   n_diaptr=`ls -l *diaptr* |wc -l`
   rebuild_list2="$rebuild_list2 ${cn_exp}_1d_${chunk_start_date}_${chunk_end_date}_diaptr"
   # setup the rebuild namelist
   ${task_work_dir}/mpi_rebuild "${rebuild_list2}" $n_diaptr
   # rebuild diaptr file
   mpirun -n 1 -env ${OMP_NUM_THREADS} ${task_work_dir}/rebuild_nemo_mpi.exe
   # add scalar and diaptr files to the list being copied to sitestore. 
   histfile_list="$histfile_list ${cn_exp}_1d_${chunk_start_date}_${chunk_end_date}_diaptr"
   if [ ${output_level} -ge 1 ] ; then
     histfile_list="$histfile_list ${cn_exp}_1m_${chunk_start_date}_${chunk_end_date}_scalar_ar6"
   fi
fi

# Since the rebuild is done we can remove links to the tiled files (they will still exist in the run tmp dir)
rm *[0-9][0-9].nc

#=================================
# Create the restart tar archive
#=================================

# tarlist contains a list of additional files to add to the restart archive
tarlist=''

cp -f ${runwrk}/namelist_cfg     rs_namelist_cfg        || :
cp -f ${runwrk}/namelist_ref     rs_namelist_ref        || :
cp -f ${runwrk}/namelist_ice_cfg rs_namelist_ice_cfg    || :
cp -f ${runwrk}/namelist_ice_ref rs_namelist_ice_ref    || :

[ -s rs_namelist_cfg ]     && tarlist="$tarlist rs_namelist_cfg"
[ -s rs_namelist_ref ]     && tarlist="$tarlist rs_namelist_ref"
[ -s rs_namelist_ice_cfg ] && tarlist="$tarlist rs_namelist_ice_cfg"
[ -s rs_namelist_ice_ref ] && tarlist="$tarlist rs_namelist_ice_ref"

if [ $pisces_offline -eq 0 ]; then
  [ -s ${cn_exp}_${chunk_nn_itend}_restart.nc ] && 
     tarlist="$tarlist ${cn_exp}_${chunk_nn_itend}_restart.nc"
  [ -s ${cn_exp}_icebergs_${chunk_nn_itend}_restart.nc ] &&
     tarlist="$tarlist ${cn_exp}_icebergs_${chunk_nn_itend}_restart.nc"
  [ -s ${cn_exp}_${chunk_nn_itend}_restart_ice.nc ]  && 
     tarlist="$tarlist ${cn_exp}_${chunk_nn_itend}_restart_ice.nc"
  [ -s nemo_physical_rtd.nc ]  && tarlist="$tarlist nemo_physical_rtd.nc"
  [ -s nemo_ice_rtd.nc ]  && tarlist="$tarlist nemo_ice_rtd.nc"
fi

if [ $nemo_trc -eq 1 ] || [ $nemo_carbon -eq 1 ]; then
  cp -f ${runwrk}/namelist_top_cfg rs_namelist_top_cfg    || :
  cp -f ${runwrk}/namelist_top_ref rs_namelist_top_ref    || :
  [ -s rs_namelist_top_cfg ] && tarlist="$tarlist rs_namelist_top_cfg"
  [ -s rs_namelist_top_ref ] && tarlist="$tarlist rs_namelist_top_ref"
  tarlist="$tarlist ${cn_exp}_${chunk_nn_itend}_restart_trc.nc"
fi

if [ $nemo_carbon -eq 1 ]; then
  cp -f ${runwrk}/namelist_pisces_cfg rs_namelist_pisces_cfg || :
  cp -f ${runwrk}/namelist_pisces_ref rs_namelist_pisces_ref || :
  [ -s rs_namelist_pisces_cfg ] && tarlist="$tarlist rs_namelist_pisces_cfg"
  [ -s rs_namelist_pisces_ref ] && tarlist="$tarlist rs_namelist_pisces_ref"
  [ -s nemo_carbon_rtd.nc ]  && tarlist="$tarlist nemo_carbon_rtd.nc"

  if [ $nemo_cmoc -eq 1 ]; then
    cp -f ${runwrk}/namelist_cmoc_cfg rs_namelist_cmoc_cfg || :
    cp -f ${runwrk}/namelist_cmoc_ref rs_namelist_cmoc_ref || :
    [ -s rs_namelist_cmoc_cfg ]  && tarlist="$tarlist rs_namelist_cmoc_cfg"
    [ -s rs_namelist_cmoc_ref ]  && tarlist="$tarlist rs_namelist_cmoc_ref"
  fi
fi

# Save a copy of the XIOS executable that was just run to the restart archive
cp -f ${runwrk}/xios_server.exe rs_xios_exec || :
[ -s rs_xios_exec ] && tarlist="$tarlist rs_xios_exec"

# Save a copy of the NEMO executable that was just run to the restart archive
cp -f ${runwrk}/nemo.exe rs_nemo_exec || :
[ -s rs_nemo_exec ] && tarlist="$tarlist rs_nemo_exec"

# Add a file containing the time step
cp -f ${runwrk}/time.step rs_time.step || :
[ -s rs_time.step ] && tarlist="$tarlist rs_time.step"

# Add ocean output and timing info
cp -f ${runwrk}/ocean.output rs_ocean.output || :
[ -s rs_ocean.output ] && tarlist="$tarlist rs_ocean.output"

cp -f ${runwrk}/timing.output rs_timing.output || :
[ -s rs_timing.output ] && tarlist="$tarlist rs_timing.output"

#cp -f rebuild_list rs_saved_history_file_names || :
#[ -s rs_saved_history_file_names ] && tarlist="$tarlist rs_saved_history_file_names"

cp -f ${runwrk}/xnemo_access_log rs_xnemo_access_log || :
[ -s rs_xnemo_access_log ] && tarlist="$tarlist rs_xnemo_access_log"

rm -f NEXT_RESTART_ARC
tar cf NEXT_RESTART_ARC $tarlist ||
  bail "Unable to create restart archive"
rs_out=${cn_exp}_${chunk_end_date}_restart.tar

mv NEXT_RESTART_ARC ${task_output_dir}/$rs_out

for i in $histfile_list; do
  mv ${i}.nc ${task_output_dir}/
done
