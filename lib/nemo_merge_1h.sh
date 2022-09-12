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
  is_defined $nemo_hist_file_suffix_list_save   || bail "nemo_merge_1h.sh: nemo_hist_file_suffix_list_save must be defined!"
  is_defined $nemo_hist_file_freq_list_save   || bail "nemo_merge_1h.sh: nemo_hist_file_freq_list_save must be defined!"

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

 # Access the history files
   nemo_hist_file_suffix_list_array_save=($nemo_hist_file_suffix_list_save)
   nemo_hist_file_freq_list_array_save=($nemo_hist_file_freq_list_save)
    n_suffix=${#nemo_hist_file_suffix_list_array_save[@]}
   for ifile in $(seq 0 $(($n_suffix-1))); do
     sfx=${nemo_hist_file_freq_list_array_save[$ifile]}_${nemo_hist_file_suffix_list_array_save[$ifile]}
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
 # Merge sub-yearly files
     if [ $nmon -gt 1 -a -e ${sfx}_m$fmon ] ; then
       cdo mergetime ${sfx}_?? ${sfx}_m$fmon
       rm -f ${sfx}_??
       mv ${sfx}_m$fmon ${sfx}_$fmon
     fi
   done

