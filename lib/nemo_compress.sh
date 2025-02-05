#!/bin/bash

set -e

# Compress stored history files

# only proceed if files will actually be kept!
if (( with_delhist==0 )) && (( with_nemo_compress==1)) ; then
  echo "Compressing NEMO history files."
  # keep track of files that have been saved with compression and files that will be deleted if successful
  saved_list=
  compress_failed=0 # successful if compress_failed==0
  # loop through each saved file
  nemo_diag_file_suffix_list_array_save=($nemo_diag_file_suffix_list)
  #for sfx in $nemo_diag_file_suffix_list_array_save ; do
  for sfx in $nemo_diag_file_suffix_list ; do
    # loop through each stored file and attemp to compress using CDO
    fname=mc_${runid}_${job_stop_date_array}_m01_${sfx}
    # overwrite original file if compression successful
    # otherwise, leave file as-is
    access $fname.lnk.nc $fname.nc >> ${WRK_DIR}/cmpr.out
    if [[ -e $fname.lnk.nc ]] ; then
      # compress history file
      cdo -f nc4c -z zip_2 copy  $(readlink -f "${fname}.lnk.nc") ${fname}.zip.nc
      if [[ ! -e  "${fname}.zip.nc" ]] ; then
        echo "Compression failed for $fname"
        compress_failed=1
        break
      else
        # save compressed file with new number
        save $fname.zip.nc ${fname}.nc || sabrt=1
        if (( sabrt!=0 )) ; then 
          # if saving failed, revert
          echo "Saving compressed files failed on ${fname}! Reverting."
          compress_failed=1
          break
        else
          # was successful; remove the local compressed file
          saved_list="${saved_list} ${fname}"
          release ${fname}.zip.nc
        fi
      fi
    fi
  done
 
  # revert to original files if unsuccessful in compression stage
  if (( compress_failed==1 )) ; then
    # remove compressed files
    rm -f *.zip.nc
    # remove any saved files
    for sname in $saved_list ; do
      access lnk.nc ${sname}.nc
      delete lnk.nc
      release lnk.nc
    done
    echo "Compression failed."
    exit 10
  else
    # compression was successful and everything was saved
    # delete the original files
    echo "Saved list: $saved_list"
    for sname in $saved_list ; do
      # remove existing file
      delete $sname.lnk.nc
      release $sname.lnk.nc
    done
    echo "Compression successful."
  fi
fi
