#########################################################
# nemo merge 1y script
# N. Lambert, Sept 2022
#
# Merge sub-yearly files inot a full 1y file
#
#########################################################

set -x

  # ======================================================
  # Make sure we have necessary variables from parent code
  # ======================================================
  is_defined $nemo_diag_file_suffix_list   || bail "nemo_merge_1h.sh: nemo_diag_file_suffix_list must be defined!"

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

 # Access the history files (and add the mesh_grid)
   nemo_diag_file_suffix_list_array_save=($nemo_diag_file_suffix_list mesh_mask)
    n_suffix=${#nemo_diag_file_suffix_list_array_save[@]}
   for ifile in $(seq 0 $(($n_suffix-1))); do
     sfx=${nemo_diag_file_suffix_list_array_save[$ifile]}
     [ $sfx == "1ts_cfg" ] || [ $sfx == "mesh_mask" ] && continue # skip the files (without temporal records)
     yr=$fyear
     mp=0
     for mm in $nemo_rtd_mons ; do
       if [ $mm -lt $mp ] ; then
         # increment year by 1 if the current month is smaller than the previous month
         yr=`echo $yr | awk '{printf "%04d", $1 + 1}'`;
       fi
       diag_hist="mc_${runid}_${yr}_m${mm}_${sfx}.nc"
       access ${sfx}_${mm} $diag_hist na
       mp=$mm
     done
 # Merge sub-yearly files and save it (delete the sub-year files)
     if [ $nmon -gt 1 -a -e "${sfx}_$fmon" ] ; then
       diag_hist="mc_${runid}_${yr}_m${fmon}_${sfx}.nc"
       cdo mergetime  ${sfx}_?? ${sfx}_merged
       for dfile in $(ls  ${sfx}_??)
       do
          delete ${dfile}
       done
       save ${sfx}_merged $diag_hist
     fi
   done

