#!/bin/bash
#============================================================================
#
#  NEMO PRELUDE
#
#  Gets the files required to run the model and sets up namelists etc.
#  Implemented here for running as a maestro task within the CanNEMO suite.
#  
#  Derived from XNEMO/nemo_prelude.
#
#  NCS, May 2017.  
#
#  Enable cmip6_diag option when with_nemo_diag=on (D. YANG, OCT 2017)
#
#  This script is peeled off from nemo_pre.tsk and adapted to facilitate migrations
#  (D. Yang, Aug 2020).
#============================================================================

  # nemo_rs is the restart file name that will be used below
  nemo_rs=${nemo_restart}

  # Total number of time steps in the current job
  steps_in_job=${nemo_steps_in_job}
  
  # Currently the only valid values for nemo_forcing are "bulk" or "flux"
  case $nemo_forcing in
    bulk|bulk_iaf|flux|flux_iaf) ;; # These are valid values
    *) bail "Invalid value for nemo_forcing = $nemo_forcing" ;;
  esac

  # The following parameters are used by nemo to define restart file names

# Configure namelist parameters etc appropriate for continuing from a restart or not
    # Only do this on the first loop iteration
    if [ $nemo_from_rest -eq 1 ] && [ "${nemo_loop_index}" -eq "${loop_start}" ] ; then
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
      #=#=#=#    Start from rest    #=#=#=#
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
      ln_rstart=.false.
      ln_rsttr=.false.

      [ -n "$nemo_restart" ] &&
	bail "nemo_from_rest = on and nemo_restart = ${nemo_restart}. These are mutually exclusive."
      echo "###############################"
      echo "###   Starting from rest.   ###"
      echo "###############################"

      # first time step
      # This will correspond with the beginning of ${start_year}:${start_mon}
      #[ -z "$days_since_run_start" ] && bail "Cannot determine days_since_run_start"
      #nn_it000=`echo $steps_per_day $days_since_run_start|awk '{printf "%d",1+$1*$2}' -`
      nn_it000=${chunk_nn_it000}

      # Use input namelist files from a user supplied location

      # The input ocean namelists are hard coded as "namelist_cfg" &
      # "namelist_ref" in nemo source
      [ -z "$nemo_namelist_cfg" ] && bail "nemo_namelist_cfg is not defined"
      acc_cp namelist_cfg $nemo_namelist_cfg
      [ -z "$nemo_namelist_ref" ] && bail "nemo_namelist_ref is not defined"
      acc_cp namelist_ref $nemo_namelist_ref

      # The input ice namelists are hard coded as "namelist_ice_cfg" & 
      # "namelist_ice_ref" in nemo source
      if [ $pisces_offline -ne  1 ]; then
        [ -z "$nemo_namelist_ice_cfg" ] && bail "nemo_namelist_ice_cfg is not defined"
        acc_cp namelist_ice_cfg $nemo_namelist_ice_cfg
        [ -z "$nemo_namelist_ice_ref" ] && bail "nemo_namelist_ice_ref is not defined"
        acc_cp namelist_ice_ref $nemo_namelist_ice_ref
      fi

      # The input top namelist is hard coded as "namelist_top_cfg" & 
      # "namelist_top_ref" in nemo source
      if [ $nemo_trc -eq 1 ] || [ $nemo_carbon -eq 1 ]; then
        [ -z "$nemo_namelist_top_cfg" ] && bail "nemo_namelist_top_cfg is not defined"
        acc_cp namelist_top_cfg $nemo_namelist_top_cfg
        [ -z "$nemo_namelist_top_ref" ] && bail "nemo_namelist_top_ref is not defined"
        acc_cp namelist_top_ref $nemo_namelist_top_ref
      fi

      if [ $nemo_carbon -eq 1 ]; then
	# The input pisces namelist is hard coded as "namelist_pisces_cfg" & 
	# "namelist_pisces_ref" in nemo source
	[ -z "$nemo_namelist_pisces_cfg" ] && bail "nemo_namelist_pisces_cfg is not defined"
	acc_cp namelist_pisces_cfg $nemo_namelist_pisces_cfg
        [ -z "$nemo_namelist_pisces_ref" ] && bail "nemo_namelist_pisces_ref is not defined"
        acc_cp namelist_pisces_ref $nemo_namelist_pisces_ref
        
        if [ $nemo_cmoc -eq 1 ]; then
          # The input cmoc namelist is hard coded as "namelist_cmoc_cfg" &
	  # "namelist_cmoc_ref" in nemo source
          [ -z "$nemo_namelist_cmoc_cfg" ] && bail "nemo_namelist_cmoc_cfg is not defined"
          acc_cp namelist_cmoc_cfg $nemo_namelist_cmoc_cfg 
	  [ -z "$nemo_namelist_cmoc_ref" ] && bail "nemo_namelist_cmoc_ref is not defined"
          acc_cp namelist_cmoc_ref $nemo_namelist_cmoc_ref
        fi
      fi

      # Get a local copy of the XIOS executable
      [ -z "$xios_exec" ] && bail "xios_exec must be defined."
      cp ${storage_dir}/executables/${xios_exec} xios_server.exe

      # Get a local copy of the NEMO executable
      [ -z "$nemo_exec" ] && bail "nemo_exec must be defined."
      cp ${storage_dir}/executables/${nemo_exec} nemo.exe

      # If CFCs are enabled we need to override namelist_top_cfg & namelist_top_ref,
      # and namelist_cfc_cfg & namelist_cfc_ref 
      if [ $nemo_cfc -eq 1 ]; then
        [ -z "$config_dir/EXP00/namelist_cfc_cfg" ] && bail "namelist_cfc_cfg not found in configuration directory"
	[ -z "$config_dir/EXP00/namelist_cfc_ref" ] && bail "namelist_cfc_ref not found in configuration directory"
        [ -z "$config_dir/EXP00/namelist_top_cfc_cfg" ] && bail "namelist_top_cfc_cfg not found in configuration directory"
	[ -z "$config_dir/EXP00/namelist_top_cfc_ref" ] && bail "namelist_top_cfc_ref not found in configuration directory"
        acc_cp namelist_top_cfg $config_dir/EXP00/namelist_top_cfc_cfg
	acc_cp namelist_top_ref $config_dir/EXP00/namelist_top_cfc_ref
        acc_cp namelist_cfc_cfg $config_dir/EXP00/namelist_cfc_cfg
	acc_cp namelist_cfc_ref $config_dir/EXP00/namelist_cfc_ref
      fi

    elif [ "${nemo_loop_index}" -eq "${loop_start}" ]; then
    # Only do a restart from an arhive if we are not on the first
    # loop of a chunk. If we are beyond the first loop iteration,
    # then we suck in the restarts from the previous tmp dir (below).

      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
      #=#=#=#   Start from a restart  #=#=#=#
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
      ln_rstart=.true.
      ln_rsttr=.true.

      # Extract restart and namelist files from the tar archive
      [ -z "$nemo_rs" ] && bail "No restart file name provided. Define nemo_restart."
      echo "##########################"
      echo "###   Restarting from file $nemo_rs"
      echo "##########################"
      rm -f RESTART_ARC
      acc_cp RESTART_ARC $nemo_rs
      tar xf RESTART_ARC || bail "Cannot extract files from $nemo_rs"
      #release RESTART_ARC
      rm RESTART_ARC  
      # Rename each task specific output restart file from the previous run
      # as the equivalent task specific input restart for the current run

      # When starting a new run from an existing restart the ocean and ice restart
      # file names found in the restart tar file will be prefixed with the runid used
      # in the previous run.

      # Ocean restart files
      if [ $pisces_offline -ne 1 ]; then
        found_rs=`(ls -1 *_${cn_ocerst_out}.nc || : ) 2>/dev/null`
        [ -z "$found_rs" ] && bail "No ocean restart files were found in $nemo_rs "
        rm -f input_ocn_restart_file_names
        touch input_ocn_restart_file_names
        for rsfile in $found_rs; do
	  mv $rsfile ${cn_ocerst_in}.nc
	  echo "${cn_ocerst_in}$sfx" >> input_ocn_restart_file_names
        done

        # Ice restart files
        rm -f input_ice_restart_file_names
        touch input_ice_restart_file_names
        found_rs=`(ls -1 *_${cn_icerst_out}.nc || : ) 2>/dev/null`
        if [ -z "$found_rs" ]; then
	  [ $nn_ice -ge 2 ] && bail "No ice restart files were found in $nemo_rs "
	  with_ice=off
        else
	  for rsfile in $found_rs; do
	    mv $rsfile ${cn_icerst_in}.nc
	    echo "${cn_icerst_in}$sfx" >> input_ice_restart_file_names
	  done
        fi
      fi

      if [ $nemo_trc -eq 1 ] || [ $nemo_carbon -eq 1 ]; then
        # Pisces restart files
        found_rs=`(ls -1 *_${cn_trcrst_out}.nc || : ) 2>/dev/null`
        [ -z "$found_rs" ] && bail "No TRC restart files were found in $nemo_rs"
        rm -f input_trc_restart_file_names
        touch input_trc_restart_file_names
        for rsfile in $found_rs; do
          mv $rsfile ${cn_trcrst_in}.nc
          echo "${cn_trcrst_in}$sfx" >> input_trc_restart_file_names
        done
      fi

      # Get a local copy of the NEMO executable
      rm -f nemo.exe
      if [ -z "$nemo_exec" ]; then
        # If nemo_exec is not defined then use the executable found in the restart archive
        [ ! -s rs_nemo_exec ] && bail "rs_nemo_exec is missing from $nemo_rs"
        mv rs_nemo_exec nemo.exe
      else
        # Otherwise use the executable pointed to by nemo_exec
        cp ${storage_dir}/executables/${nemo_exec} nemo.exe
      fi
      [ ! -s nemo.exe ] && bail "Unable to find nemo executable."

      if [ $nemo_force_namelists -eq 1 ]; then
	# Override the namelist files found in the restart archive with
	# files supplied by the user

	# The input ocean namelist is hard coded as "namelist_cfg" &
        # "namelist_ref" in nemo source
	[ -z "$nemo_namelist_cfg" ] && bail q"nemo_namelist_cfg is not defined"
	acc_cp namelist_cfg $nemo_namelist_cfg nocp=no
        [ -z "$nemo_namelist_ref" ] && bail q"nemo_namelist_ref is not defined"
        acc_cp namelist_ref $nemo_namelist_ref nocp=no

	if [ $pisces_offline -ne 1 ]; then
	  # The input ice namelist is hard coded as "namelist_ice_cfg" &
	  # "namelist_ice_ref" in nemo source
	  [ -z "$nemo_namelist_ice_cfg" ] && bail "nemo_namelist_ice_cfg is not defined"
	  acc_cp namelist_ice_cfg $nemo_namelist_ice_cfg nocp=no
	  [ -z "$nemo_namelist_ice_ref" ] && bail "nemo_namelist_ice_ref is not defined"
          acc_cp namelist_ice_ref $nemo_namelist_ice_ref nocp=no
        fi

	if [ $nemo_trc -eq 1 ] || [ $nemo_carbon -eq 1 ]; then
          # The input top namelist is hard coded as "namelist_top_cfg" &
	  # "namelist_top_ref" in nemo source
          [ -z "$nemo_namelist_top_cfg" ] && bail "nemo_namelist_top_cfg is not defined"
          acc_cp namelist_top_cfg $nemo_namelist_top_cfg nocp=no
	  [ -z "$nemo_namelist_top_ref" ] && bail "nemo_namelist_top_ref is not defined"
          acc_cp namelist_top_ref $nemo_namelist_top_ref nocp=no
        fi

	if [ $nemo_carbon -eq 1 ]; then
	  # The input pisces namelist is hard coded as "namelist_pisces_cfg" & 
	  # "namelist_pisces_ref" in nemo source
	  [ -z "$nemo_namelist_pisces_cfg" ] && bail "nemo_namelist_pisces_cfg is not defined"
	  acc_cp namelist_pisces_cfg $nemo_namelist_pisces_cfg nocp=no
	  [ -z "$nemo_namelist_pisces_ref" ] && bail "nemo_namelist_pisces_ref is not defined"
          acc_cp namelist_pisces_ref $nemo_namelist_pisces_ref nocp=no
        
          if [ $nemo_cmoc -eq 1 ]; then
	    # The input cmoc namelist is hard coded as "namelist_cmoc_cfg" &
	    # "namelist_cmoc_ref" in nemo source
	    [ -z "$nemo_namelist_cmoc_cfg" ] && bail "nemo_namelist_cmoc_cfg is not defined"
	    acc_cp namelist_cmoc_cfg $nemo_namelist_cmoc_cfg nocp=no
	    [ -z "$nemo_namelist_cmoc_ref" ] && bail "nemo_namelist_cmoc_ref is not defined"
            acc_cp namelist_cmoc_ref $nemo_namelist_cmoc_ref nocp=no
          fi
          # If CFCs are enabled we need to override namelist_top_cfg & namelist_top_ref,
	  #  and namelist_cfc_cfg & namelist_cfc_ref
          if [ $nemo_cfc -eq 1 ]; then
            [ -z "$config_dir/EXP00/namelist_cfc_cfg" ] && bail "namelist_cfc_cfg not found in configuration directory"
	    [ -z "$config_dir/EXP00/namelist_cfc_ref" ] && bail "namelist_cfc_ref not found in configuration directory"
            [ -z "$config_dir/EXP00/namelist_top_cfc_cfg" ] && bail "namelist_top_cfc_cfg not found in configuration directory"
	    [ -z "$config_dir/EXP00/namelist_top_cfc_ref" ] && bail "namelist_top_cfc_ref not found in configuration directory"
            acc_cp namelist_top_cfg $config_dir/EXP00/namelist_top_cfc_cfg
            acc_cp namelist_top_ref $config_dir/EXP00/namelist_top_cfc_ref
            acc_cp namelist_cfc_cfg $config_dir/EXP00/namelist_cfc_cfg
            acc_cp namelist_cfc_ref $config_dir/EXP00/namelist_cfc_ref
          fi
	fi
      else
	# Use the namelist files found in the restart archive
	[ ! -s rs_namelist_cfg ]     && bail "Missing rs_namelist_cfg from $nemo_rs"
	[ ! -s rs_namelist_ref ]     && bail "Missing rs_namelist_ref from $nemo_rs"

	if [ $pisces_offline -ne 1 ]; then
	  [ ! -s rs_namelist_ice_cfg ] && bail "Missing rs_namelist_ice_cfg from $nemo_rs"
	  [ ! -s rs_namelist_ice_ref ] && bail "Missing rs_namelist_ice_ref from $nemo_rs"
        fi

	if [ $nemo_trc -eq 1 ] || [ $nemo_carbon -eq 1 ]; then
	  [ ! -s rs_namelist_top_cfg ] && bail "Missing rs_namelist_top_cfg from $nemo_rs"
	  [ ! -s rs_namelist_top_ref ] && bail "Missing rs_namelist_top_ref from $nemo_rs"
        fi 

	if [ $nemo_carbon -eq 1 ]; then
	  [ ! -s rs_namelist_pisces_cfg ] && bail "Missing rs_namelist_pisces_cfg from $nemo_rs"
	  [ ! -s rs_namelist_pisces_ref ] && bail "Missing rs_namelist_pisces_ref from $nemo_rs"

          if [ $nemo_cmoc -eq 1 ]; then
	    [ ! -s rs_namelist_cmoc_cfg ] && bail "Missing rs_namelist_cmoc_cfg from $nemo_rs"
	    [ ! -s rs_namelist_cmoc_ref ] && bail "Missing rs_namelist_cmoc_ref from $nemo_rs"
          fi
	fi

	# The input ocean namelist_cfg is hard coded as "namelist_cfg" &
        # "namelist_ref" in nemo source
	cp -f rs_namelist_cfg namelist_cfg
	cp -f rs_namelist_ref namelist_ref

	# The input ice namelist is hard coded as "namelist_ice_cfg" 
	# "namelist_ice_ref" in nemo source
	if [ $pisces_offline -ne 1 ]; then
	  cp -f rs_namelist_ice_cfg namelist_ice_cfg
	  cp -f rs_namelist_ice_ref namelist_ice_ref
        fi
 
	# The input top namelist is hard coded as "namelist_top_cfg" &
	# "namelist_top_ref" in nemo source
	if [ $nemo_trc -eq 1 ] || [ $nemo_carbon -eq 1 ]; then
	  cp -f rs_namelist_top_cfg namelist_top_cfg
	  cp -f rs_namelist_top_ref namelist_top_ref
        fi

	if [ $nemo_carbon -eq 1 ]; then

	  # The input pisces namelist is hard coded as "namelist_pisces_cfg" &
	  # "namelist_pisces_ref" in nemo source
	  cp -f rs_namelist_pisces_cfg namelist_pisces_cfg
	  cp -f rs_namelist_pisces_ref namelist_pisces_ref

          if [ $nemo_cmoc -eq 1 ]; then
    	    # The input cmoc namelist is hard coded as "namelist_cmoc_cfg" &
	    # "namelist_cmoc_ref" in nemo source
	    cp -f rs_namelist_cmoc_cfg namelist_cmoc_cfg
	    cp -f rs_namelist_cmoc_ref namelist_cmoc_ref
          fi
	fi
      fi 

      # Current beginning time step
      # get this from the time.step file that was saved in the restart
      #prev_itend=`cat rs_time.step`
      #nn_it000=`expr $prev_itend + 1`   
      #echo "nn_it000 " : ${nn_it000}
      nn_it000=${chunk_nn_it000}

    elif [ "${nemo_loop_index}" -gt "${loop_start}" ]; then
      # In this case we are beyond the start of the loop. We want to use
      # restarts, but not those from "disk" / storage. Instead we use restarts
      # from the previous iteration of the loop, which should be in the previous
      # tmp directory and can be linked in.
            ln_rstart=.true.
            ln_rsttr=.true.
            nn_no=${nemo_loop_index}
            # Link in restarts from previous loop. Needs to be more robust.
            lm1=$((${nemo_loop_index} - 1))
	    # Define output directory for nemo_run+${lm1}
            outdir_lm1=${nemo_seq_workbase}/${nemo_seq_container}/nemo_run+${lm1}/output

            # Ocean restart files
            found_rs=`(ls -1 ${outdir_lm1}/*_${cn_ocerst_out}_[0-9][0-9][0-9][0-9].nc || : ) 2>/dev/null`
            for rsfile in $found_rs; do
                sfx=`echo $rsfile|sed 's/^.*\(_[0-9][0-9][0-9][0-9].nc\).*$/\1/'`
                #[ "x$sfx" = "x$rsfile" ] && sfx=`echo $rsfile|sed 's/^.*\(\.nc\).*$/\1/'`
                ln -s $rsfile ${cn_ocerst_in}$sfx
            done

            # Ice restart files
            found_rs=`(ls -1 ${outdir_lm1}/*_${cn_icerst_out}_[0-9][0-9][0-9][0-9].nc || : ) 2>/dev/null`
            for rsfile in $found_rs; do
                sfx=`echo $rsfile|sed 's/^.*\(_[0-9][0-9][0-9][0-9].nc\).*$/\1/'`
                ln -s $rsfile ${cn_icerst_in}$sfx
            done

            # trc restart files
            found_rs=`(ls -1 ${outdir_lm1}/*_${cn_trcrst_out}_[0-9][0-9][0-9][0-9].nc || : ) 2>/dev/null`
            for rsfile in $found_rs; do
                sfx=`echo $rsfile|sed 's/^.*\(_[0-9][0-9][0-9][0-9].nc\).*$/\1/'`
                ln -s $rsfile ${cn_trcrst_in}$sfx
            done

            cp ${nemo_seq_workbase}/${nemo_seq_container}/nemo_run+${lm1}/work/namelist* .
            cp ${nemo_seq_workbase}/${nemo_seq_container}/nemo_run+${lm1}/work/nemo.exe .

            nn_it000=${chunk_nn_it000}
    else
        bail "Don't know how to start."
    fi

    # The user may override the initial time step value by setting nemo_nn_it000
    [ -n "$nemo_nn_it000" ] && nn_it000=$nemo_nn_it000

    # Time-step information
    # last  time step
    # This will correspond with the end of ${stop_year}:${stop_mon}
    nn_itend=`echo $nn_it000 $steps_in_job|awk '{printf "%d",$1+$2-1}' -`

    # initial calendar date yymmdd
    nn_date0=${chunk_start_date}
    #`echo ${NEMO_CHUNK_START_YEAR} ${NEMO_CHUNK_START_MONTH} ${NEMO_CHUNK_START_DAY} |awk '{printf "%2.2d%2.2d%2.2d",$1,$2,$3}' -`

    # frequency of creation of a restart file (modulo referenced to 1)
    nn_stock=$nn_itend

    # frequency of write in the output file (modulo referenced to nn_it000)
    nn_write=$steps_in_job

    # Modify specified parameter values in the input ocean namelist
    mod_nl namelist_cfg nn_it000 nn_itend jpni jpnj \
      rn_rdt cn_exp nn_date0 ln_rstart nn_rstctl cn_ocerst_in cn_ocerst_out \
      ln_ctl nn_print nn_ice sn_rnf sn_cnf

    # Revise the ocean namelist flags to output on-the-fly diagnostics when with_nemo_diag = on.
    if [ ${with_nemo_diag} -eq 1 ] ; then
      ln_diaptr=.true.
      # ln_diaznl=.false.
      # nn_fwri=`echo $rn_rdt | awk '{printf "%d",86400/$1}'` # Frequency of ptr outputs (time steps in a day)
      mod_nl namelist_cfg ln_diaptr
    fi

    # Modify specified parameter values in the input ice namelist
    if [ $pisces_offline -ne 1 ]; then
      mod_nl namelist_ice_cfg cn_icerst_in cn_icerst_out
    fi

    # Modify the input top namelist
    if [ $nemo_trc -eq 1 ] || [ $nemo_carbon -eq 1 ]; then
	mod_nl namelist_top_cfg ln_rsttr cn_trcrst_in cn_trcrst_out
    fi

    # Modify cfc namelist and copy over input file
    if [ $nemo_cfc -eq 1 ]; then
        mod_nl namelist_cfc offset_cfc_year
        acc_cp nemo_cfc_atmhist.nc $nemo_cfc_atmhist 
    fi
    # Access or copy files required to run NEMO
    acc_cp iodef.xml              $nemo_iodef
    acc_cp context_nemo.xml       $nemo_context_nemo
    acc_cp field_def_nemo-oce.xml $nemo_field_def_nemo_oce
    acc_cp field_def_nemo-ice.xml $nemo_field_def_nemo_ice
    acc_cp field_def_nemo-pisces.xml $nemo_field_def_nemo_pisces
    acc_cp file_def_nemo-oce.xml  $nemo_file_def_nemo_oce
    acc_cp file_def_nemo-ice.xml  $nemo_file_def_nemo_ice
    acc_cp file_def_nemo-pisces.xml  $nemo_file_def_nemo_pisces
    acc_cp axis_def_nemo.xml      $nemo_axis_def_nemo
    acc_cp domain_def_nemo.xml    $nemo_domain_def_nemo
    acc_cp grid_def_nemo.xml      $nemo_grid_def_nemo
    #acc_cp xmlio_server.def       $nemo_xmlio_server_def
    acc_cp mixing_power_bot.nc    $nemo_mixing_power_bot
    acc_cp mixing_power_pyc.nc    $nemo_mixing_power_pyc
    acc_cp mixing_power_cri.nc    $nemo_mixing_power_cri
    acc_cp decay_scale_bot.nc     $nemo_decay_scale_bot
    acc_cp decay_scale_cri.nc     $nemo_decay_scale_cri
    acc_cp eddy_viscosity_3D.nc   $nemo_eddy_viscosity_3D
    acc_cp resto.nc               $nemo_resto
    acc_cp EMPave_old.dat         $nemo_empave_old
    acc_cp ORCA_R2_zps_domcfg.nc  $nemo_coordinates
    #acc_cp ahmcoef                $nemo_ahmcoef
    #acc_cp coordinates.nc         $nemo_coordinates
    #acc_cp bathy_meter.nc         $nemo_bathy_meter
    #acc_cp mask_itf.nc            $nemo_mask_itf
    #acc_cp M2rowdrg.nc            $nemo_M2rowdrg
    #acc_cp K1rowdrg.nc            $nemo_K1rowdrg
    #acc_cp Eddyengf.nc            $nemo_Eddyengf
    acc_cp geothermal_heating.nc  $nemo_geothermal_heating
    acc_cp runoff_core_monthly.nc $nemo_runoff_core_monthly
    acc_cp sss_data.nc            $nemo_sss_data
    acc_cp sst_data.nc            $nemo_sst_data
    acc_cp chlorophyll.nc         $nemo_chlorophyll
    acc_cp subbasins.nc           $nemo_subbasins
    acc_cp data_1m_potential_temperature_nomask.nc $nemo_data_1m_potential_temperature_nomask
    acc_cp data_1m_salinity_nomask.nc              $nemo_data_1m_salinity_nomask
    acc_cp weights_bic2.nc        $nemo_weights_bicubic2   # weights for bulk mode
    acc_cp weights_bil2.nc        $nemo_weights_bilinear2  # weights for bulk mode
    #acc_cp weights_bic3.nc        $nemo_weights_bicubic3   # weights for flux mode
    #acc_cp weights_bil3.nc        $nemo_weights_bilinear3  # weights for flux mode

    # Update NEMO priority level of output variables: output_level 
    sed -i "s/output_level=\"\([0-9]\+\)\"/output_level=\"${output_level}\"/" iodef.xml

    # Copy in surface forcing files depending on the option chosen.
    if [ x"$nemo_forcing" = "xbulk" ]; then
      acc_cp tair10m.nc          $nemo_t
      acc_cp uwnd10m.nc          $nemo_u
      acc_cp vwnd10m.nc          $nemo_v
      acc_cp humi10m.nc          $nemo_q
      acc_cp qsr.nc              $nemo_qsr
      acc_cp qlw.nc              $nemo_qlw
      acc_cp precip.nc           $nemo_precip
      acc_cp snow.nc             $nemo_snow
      acc_cp slp.nc              $nemo_slp

    elif [ x"$nemo_forcing" = "xbulk_iaf" ]; then
      # Get forcing data, whcih may vary by year
      for current_year in ${chunk_years}; do

        yearc=`expr ${current_year} | awk  '{printf "%04d",$1}'`
        fy=`expr ${yearc} + ${iaf_year_offset} | awk  '{printf "%04d",$1}'`
        yd=`expr ${iaf_loop_year} - ${iaf_year_offset}`
        fyl=`expr ${fy} - '(' '(' ${yearc} - 1 ')' / ${yd} ')' \* ${yd}`           
        forcing_year=`expr ${fyl} | awk  '{printf "%04d",$1}'`
      
        acc_cp tair10m_y${yearc}.nc          ${nemo_t}_y${forcing_year}.nc
        acc_cp uwnd10m_y${yearc}.nc          ${nemo_u}_y${forcing_year}.nc
        acc_cp vwnd10m_y${yearc}.nc          ${nemo_v}_y${forcing_year}.nc
        acc_cp humi10m_y${yearc}.nc          ${nemo_q}_y${forcing_year}.nc
        acc_cp qsw_y${yearc}.nc              ${nemo_qsw}_y${forcing_year}.nc
        acc_cp qlw_y${yearc}.nc              ${nemo_qlw}_y${forcing_year}.nc
        acc_cp precip_y${yearc}.nc           ${nemo_precip}_y${forcing_year}.nc
        acc_cp snow_y${yearc}.nc             ${nemo_snow}_y${forcing_year}.nc
      done
    elif [ x"$nemo_forcing" = "xflux" ]; then
      acc_cp utau.nc             $nemo_utau
      acc_cp vtau.nc             $nemo_vtau
      acc_cp qtot.nc             $nemo_qtot
      acc_cp qsr.nc              $nemo_qsr
      acc_cp emp.nc              $nemo_emp

    elif [ x"$nemo_forcing" = "xflux_iaf" ]; then
      # Get forcing data, whcih may vary by year
      for current_year in ${chunk_years}; do

        yearc=`expr ${current_year} | awk  '{printf "%04d",$1}'`
        fy=`expr ${yearc} + ${iaf_year_offset} | awk  '{printf "%04d",$1}'`
        yd=`expr ${iaf_loop_year} - ${iaf_year_offset}`
        fyl=`expr ${fy} - '(' '(' ${yearc} - 1 ')' / ${yd} ')' \* ${yd}`           
        forcing_year=`expr ${fyl} | awk  '{printf "%04d",$1}'`

        acc_cp utau_y${yearc}.nc             ${nemo_utau}_y${forcing_year}.nc
        acc_cp vtau_y${yearc}.nc             ${nemo_vtau}_y${forcing_year}.nc
        acc_cp qtot_y${yearc}.nc             ${nemo_qtot}_y${forcing_year}.nc
        acc_cp qsr_y${yearc}.nc              ${nemo_qsr}_y${forcing_year}.nc
        acc_cp emp_y${yearc}.nc              ${nemo_emp}_y${forcing_year}.nc
      done
    fi

    # Modify the namelists to suite the forcing, as specified in the makefile
    # Bulk forcing
    if [ x"$nemo_forcing" = "xbulk" ] || [ x"$nemo_forcing" = "xbulk_iaf" ] ; then
      #ln_ana=".false.   !  analytical formulation                    (T => fill namsbc_ana )"
      #ln_flx=".false.   !  flux formulation                          (T => fill namsbc_flx )"
      #ln_blk_clio=".false.   !  CLIO bulk formulation                     (T => fill namsbc_clio)"
      #ln_blk_core=".true.    !  CORE bulk formulation                     (T => fill namsbc_core)"
      ln_blk=".true.    !  CORE bulk formulation                     (T => fill namsbc_core)"
      #ln_blk_mfs=".false.   !  MFS bulk formulation                      (T => fill namsbc_mfs )"
      #ln_cpl=".false.   !  Coupled formulation                       (T => fill namsbc_cpl )"    
      mod_nl namelist_cfg ln_blk sn_wndi sn_wndj sn_qsr sn_qlw sn_tair sn_humi \
                      sn_prec sn_snow sn_slp sn_tdif
   
    elif [ x"$nemo_forcing" = "xflux" ] || [ x"$nemo_forcing" = "xflux_iaf" ] ; then
      #ln_ana=".false.   !  analytical formulation                    (T => fill namsbc_ana )"
      ln_flx=".true.    !  flux formulation                          (T => fill namsbc_flx )"
      #ln_blk_clio=".false.   !  CLIO bulk formulation                     (T => fill namsbc_clio)" 
      #ln_blk_core=".false.    !  CORE bulk formulation                     (T => fill namsbc_core)" 
      #ln_blk_mfs=".false.   !  MFS bulk formulation                      (T => fill namsbc_mfs )"
      #ln_cpl=".false.   !  Coupled formulation                       (T => fill namsbc_cpl )"    
      mod_nl namelist_cfg ln_flx sn_utau sn_vtau sn_qtot sn_qsr sn_emp
    fi

    if [ $pisces_offline -eq 1 ]; then
      acc_cp dyna_grid_t.nc $dyna_grid_t
      acc_cp dyna_grid_u.nc $dyna_grid_u
      acc_cp dyna_grid_v.nc $dyna_grid_v
      acc_cp dyna_grid_w.nc $dyna_grid_w
      acc_cp mesh_mask.nc   $offline_mash_mask
    fi 

    if [ $nemo_carbon -eq 1 ]; then
      # Access PISCES files
      acc_cp presatm.nc           $nemo_presatm
      acc_cp dust.orca.nc         $nemo_dust_orca
      acc_cp river.orca.nc        $nemo_river_orca
      acc_cp ndeposition.orca.nc  $nemo_ndeposition_orca
      acc_cp bathy.orca.nc        $nemo_bathy_orca
      acc_cp par.orca.nc          $nemo_par_orca
      acc_cp solubility.orca.nc   $nemo_solubility_orca
      acc_cp hydrofe.orca.nc      $nemo_hydrofe_orca

      if [ $nemo_cmoc -eq 1 ]; then
        acc_cp fermask.orca.nc $nemo_fermask_orca
      fi

      acc_cp data_DIC_nomask.nc        $nemo_data_dic_nomask
      acc_cp data_Alkalini_nomask.nc   $nemo_data_alkalini_nomask
      acc_cp data_O2_nomask.nc         $nemo_data_o2_nomask
      acc_cp data_PO4_nomask.nc        $nemo_data_po4_nomask
      acc_cp data_Si_nomask.nc         $nemo_data_si_nomask
      acc_cp data_DOC_nomask.nc        $nemo_data_doc_nomask
      acc_cp data_Fer_nomask.nc        $nemo_data_fer_nomask
      acc_cp data_NO3_nomask.nc        $nemo_data_no3_nomask
      acc_cp atcco2.nc                 $nemo_data_atcco2

      # Abiotic biogeochemistry
      acc_cp data_at_abio_nomask.nc    $nemo_data_at_abio_nomask
      acc_cp data_c14t_abio_nomask.nc  $nemo_data_c14t_abio_nomask
      acc_cp data_ct_abio_nomask.nc    $nemo_data_ct_abio_nomask
      acc_cp data_o2abio_nomask.nc     $nemo_data_o2_abio_nomask
    
      # Modify the pisces namelist to read CO2 from file, if specified in makefile
      ln_co2int=$ln_co2int # read atm pco2 from a file (T) or constant (F) 
      mod_nl namelist_pisces_cfg ln_co2int
    fi 
