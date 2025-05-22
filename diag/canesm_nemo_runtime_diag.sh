#!/bin/bash
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Local function/variable defs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  stamp=$(date "+%j%H%M%S"$$)


  function rtd_staging() {
    # Copy a single rtd file to a local staging directory
    # Usage: rtd_staging local_rtd_file remote_rtd_file
    [ -z "$1" ] && bail "rtd_staging requires a file name as the first arg."
    [ -z "$2" ] && bail "rtd_staging requires a file name as the second arg."
    [ -e $1 ] || bail "rtd_staging: File --> $1 <-- is missing."


    rtd_staging_dir=$cannemo_rtd_staging_dir

    echo "staging" $cannemo_rtd_staging_dir
    [ -z $cannemo_rtd_staging_dir ] && bail "rtd_staging_dir not defined"
    cp $1 $cannemo_rtd_staging_dir/$2
  }

#~~~~~~~~~~~~~~~~~~~
#   Main Script
#~~~~~~~~~~~~~~~~~~~

# NCS, Feb 2018
# We are reverting to the AGCM RTD approach, since this is being run in the
# coupled model, and it is easiest to be consistent across models
# What this means is that we are going to run RTD in one year chunks. We
# are going to make the hard assumption that the coupled model was run
# in monthly chunks, with monthly output. Hence, this RTD job will loop
# over the 12 months for the current year, computing the RTD for each
# month progressively.
# We are also going to use the rtdiag_start_year paradigm. That is, we are
# going to save files with names including the dates, which for the current
# segment will be from rtdiag_start_year to year.

# first and last months
  mon1=`echo ${month_rtdiag_start:-1} | awk '{printf "%02d",$1}'`
  monl=`expr $mon1 + 11`
  if [ $monl -gt 12 ] ; then
    mon2=`echo $monl | awk '{printf "%02d",$1-12}'`
    yr1=`echo $year | awk '{printf "%04d", $1-1}'`;
  else
    mon2=`echo $monl | awk '{printf "%02d",$1}'`
    yr1=`echo $year | awk '{printf "%04d", $1}'`;
  fi
  nmon=`expr $monl - $mon1 + 1`

# Timing: figure out previous year, and previous RTD file
  keep_old_rtdiag_number=${keep_old_rtdiag_number:=3}
  echo keep_old_rtdiag_number=$keep_old_rtdiag_number

  yearm1=`echo $year | awk '{printf "%04d", $1 - 1}'`
  yearmo=`echo $year $keep_old_rtdiag_number | awk '{printf "%04d", $1 - $2}'`

  # current year rtd file names
  physical_rtdfile="sc_${runid}_${year_rtdiag_start}${mon1}_${year}${mon2}_nemo_physical_rtd.nc"
  ice_rtdfile="sc_${runid}_${year_rtdiag_start}${mon1}_${year}${mon2}_nemo_ice_rtd.nc"
  carbon_rtdfile="sc_${runid}_${year_rtdiag_start}${mon1}_${year}${mon2}_nemo_carbon_rtd.nc"

  # previous year rtd file names
  physical_rtdfile1="sc_${runid}_${year_rtdiag_start}${mon1}_${yearm1}${mon2}_nemo_physical_rtd.nc"
  ice_rtdfile1="sc_${runid}_${year_rtdiag_start}${mon1}_${yearm1}${mon2}_nemo_ice_rtd.nc"
  carbon_rtdfile1="sc_${runid}_${year_rtdiag_start}${mon1}_${yearm1}${mon2}_nemo_carbon_rtd.nc"

  # older rtd file to be deleted (depends on keep_old_rtdiag_number)
  physical_rtdfileo="sc_${runid}_${year_rtdiag_start}${mon1}_${yearmo}${mon2}_nemo_physical_rtd.nc"
  ice_rtdfileo="sc_${runid}_${year_rtdiag_start}${mon1}_${yearmo}${mon2}_nemo_ice_rtd.nc"
  carbon_rtdfileo="sc_${runid}_${year_rtdiag_start}${mon1}_${yearmo}${mon2}_nemo_carbon_rtd.nc"

# Access old RTD files from last year. Note, if these files exist, the RTD programs
# below will automatically append to them.

  if [ $yr1 -gt ${year_rtdiag_start} ] ; then
    access nemo_physical_rtd_old.nc $physical_rtdfile1 nocp=off ; cp nemo_physical_rtd_old.nc nemo_physical_rtd.nc ; chmod +w nemo_physical_rtd.nc
    access nemo_ice_rtd_old.nc $ice_rtdfile1 nocp=off ; cp nemo_ice_rtd_old.nc nemo_ice_rtd.nc ; chmod +w nemo_ice_rtd.nc
    access nemo_carbon_rtd_old.nc $carbon_rtdfile1 nocp=off na ; [ -s nemo_carbon_rtd_old.nc ] && cp nemo_carbon_rtd_old.nc nemo_carbon_rtd.nc ; chmod +w nemo_carbon_rtd.nc
  fi

# Get the RTD executables

  [ -z "$nemo_physical_rtd_exe" ] && bail "nemo_physical_rtd_exe is not defined."
  cp ${EXEC_STORAGE_DIR}/${nemo_physical_rtd_exe} .

  [ -z "$nemo_ice_rtd_exe" ] && bail "ice_rtd_exe is not defined."
  cp ${EXEC_STORAGE_DIR}/${nemo_ice_rtd_exe} .

  [ -z "$nemo_carbon_rtd_exe" ] && bail "nemo_carbon_rtd_exe is not defined."
  cp ${EXEC_STORAGE_DIR}/${nemo_carbon_rtd_exe} .

  ######## Compute the rt diagnostics #######

    # Invoke the nemo run time diagnostics
  echo "Doing  NEMO RTD for files with months ${nemo_rtd_mons}"

  # Month loop
  for mon in $nemo_rtd_mons ; do
      echo $mon
      if [ $mon1 -eq 1 ] ; then
        yearm=`echo $year | awk '{printf "%04d",$1}'`
      else
        if [ $mon -ge $mon1 ] ; then
          yearm=`echo $year | awk '{printf "%04d",$1-1}'`
        else
          yearm=`echo $year | awk '{printf "%04d",$1}'`
        fi
      fi
                   # Physical run time diagnostics

      # Attempt to access the "physical" history files
      rtd_hist1="mc_${runid}_${yearm}_m${mon}_1m_grid_t.nc"
      rtd_hist2="mc_${runid}_${yearm}_m${mon}_1m_grid_u.nc"
      rtd_hist3="mc_${runid}_${yearm}_m${mon}_1m_grid_v.nc"
      rtd_hist4="mc_${runid}_${yearm}_m${mon}_1m_grid_w.nc"
      rtd_hist5="mc_${runid}_${yearm}_m${mon}_mesh_mask.nc"

      access grid_t $rtd_hist1 na #|| bail "NEMO rdt cannot access $rtd_hist1"
      access grid_u $rtd_hist2 na #|| bail "NEMO rdt cannot access $rtd_hist2"
      access grid_v $rtd_hist3 na #|| bail "NEMO rdt cannot access $rtd_hist3"
      access grid_w $rtd_hist4 na #|| bail "NEMO rdt cannot access $rtd_hist4"
      access orca_mesh_mask $rtd_hist5 na

      # Create run time diagnostics for physical ocean variables
      $nemo_physical_rtd_exe ${yearm} ${mon}

                   # Sea-ice run time diagnostics

      # Access additional annual history files containing PISCES related variables
      rtd_hist7="mc_${runid}_${yearm}_m${mon}_1m_icemod.nc"
      access icemod $rtd_hist7 na

      # Create run time diagnostics for ice variables
      [ -s icemod ] && $nemo_ice_rtd_exe ${yearm} ${mon}

                   # Carbon run time diagnostics

      # Access additional annual history files containing PISCES related variables
      rtd_hist8="mc_${runid}_${yearm}_m${mon}_1m_ptrc_t.nc"
      rtd_hist9="mc_${runid}_${yearm}_m${mon}_1m_diad_t.nc"
      access ptrc_t $rtd_hist8 na
      access diad_t $rtd_hist9 na

      [ -s ptrc_t ] && [ -s diad_t ] && $nemo_carbon_rtd_exe ${yearm} ${mon}

      # Clean up
      release grid_t
      release grid_u
      release grid_v
      release grid_w
      release icemod
      release ptrc_t
      release diad_t
      release orca_mesh_mask
   done # month loop end

   # Cleanup
   release orca_mesh_mask

#   Move RTD files for staging back to Victoria
    rtd_staging nemo_physical_rtd.nc $physical_rtdfile
    rtd_staging nemo_ice_rtd.nc      $ice_rtdfile
    rtd_staging nemo_carbon_rtd.nc   $carbon_rtdfile

  # Save new RTD files.
    save nemo_physical_rtd.nc $physical_rtdfile
    save nemo_ice_rtd.nc      $ice_rtdfile
    save nemo_carbon_rtd.nc   $carbon_rtdfile

# ************************************ delete older rtd files
    if [ "$keep_old_rtdiag" != "on" ] ; then
      release old oldi oldc
      access  old  $physical_rtdfileo na
      access  oldi $ice_rtdfileo      na
      access  oldc $carbon_rtdfileo   na

      delete old  na
      delete oldi na
      delete oldc na
    fi   