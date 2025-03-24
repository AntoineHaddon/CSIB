#!/bin/bash
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Local function/variable defs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  stamp=$(date "+%j%H%M%S"$$)

  # this_host will simply be the output from uname
  this_host=$(uname -n|awk -F. '{print $1}' -)

  # this_mach will be a known alias (or possibly the actual machine name)
  this_mach=$this_host
  ROUTE_SSH=''
  case $this_mach in
             xc3*) this_mach=banting; ROUTE_SSH='ssh banting';;
             xc4*) this_mach=daley;   ROUTE_SSH='ssh daley'  ;;
      hpcr3*|cs3*|ppp3*|eccc*-ppp3*) this_mach=ppp3;   ROUTE_SSH='ssh ppp3'  ;;
      hpcr4*|cs4*|ppp4*|eccc*-ppp4*) this_mach=ppp4;   ROUTE_SSH='ssh ppp4'  ;;
  esac

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

  # bail is a simple error exit routine
  bail_prefix="NEMO RTD"
  bail(){
    echo `date`" $runid --- ${bail_prefix}: $*"
    exit 1
  }

  # function to extract files types from plot yaml files
  get_plot_filetypes() {
    plt_fts=`cat $1 | grep 'ifile' | awk '{print $2}'`
    fts=""
    for ft in $plt_fts; do
      if [[ ${fts} != *${ft}* ]]; then
        fts=${fts}" "${ft}
      fi
    done
    echo $fts
  }

#~~~~~~~~~~~~~~~~~~~
#   Main Script
#~~~~~~~~~~~~~~~~~~~

  if [[ $runmode == *"CanTODS"* ]]; then
  # JGI, Oct 2024
  # CanTODS RTD using a simple Python script
  
    # activate Python environment
    orig_env=$CONDA_PREFIX # current environment just in case different
    source activate /home/scrd102/cccma_conda/envs/py3_analysis_v2

    # year referenced to run_start_date
    rtdyear=$(expr $year + 0)           # strip leading zeros
    rtdrst=$(expr $run_start_time + 0)
    rtdYR=$(( rtdyear - rtdrst ))

    # run script
    python3 ${CCRNSRC}/CanESM/CanNEMO/rtd/cantods_rtd.py -r $runid -p  ${RUNPATH%%${runid_env}*} -o ${RUNPATH}/RTD -x [${rtdYR},${rtdYR}] -y ${run_start_time} -P True -s ${CCRNSRC}/CanESM/CanNEMO/rtd

    # revert to original environment
    source deactivate
    if [ ! -z "${orig_env}" ]; then
      source activate $orig_env
    fi

  else
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
    year0=$(pad_integer $start_rtdiag 4)
    physical_rtdfile="sc_${runid}_${year0}${mon1}_${year}${mon2}_nemo_physical_rtd.nc"
    ice_rtdfile="sc_${runid}_${year0}${mon1}_${year}${mon2}_nemo_ice_rtd.nc"
    carbon_rtdfile="sc_${runid}_${year0}${mon1}_${year}${mon2}_nemo_carbon_rtd.nc"

    # previous year rtd file names
    physical_rtdfile1="sc_${runid}_${year0}${mon1}_${yearm1}${mon2}_nemo_physical_rtd.nc"
    ice_rtdfile1="sc_${runid}_${year0}${mon1}_${yearm1}${mon2}_nemo_ice_rtd.nc"
    carbon_rtdfile1="sc_${runid}_${year0}${mon1}_${yearm1}${mon2}_nemo_carbon_rtd.nc"

    # older rtd file to be deleted (depends on keep_old_rtdiag_number)
    physical_rtdfileo="sc_${runid}_${year0}${mon1}_${yearmo}${mon2}_nemo_physical_rtd.nc"
    ice_rtdfileo="sc_${runid}_${year0}${mon1}_${yearmo}${mon2}_nemo_ice_rtd.nc"
    carbon_rtdfileo="sc_${runid}_${year0}${mon1}_${yearmo}${mon2}_nemo_carbon_rtd.nc"

  # Access old RTD files from last year. Note, if these files exist, the RTD programs
  # below will automatically append to them.
    if [ $yr1 -gt ${start_rtdiag} ] ; then
      access nemo_physical_rtd_old.nc $physical_rtdfile1 nocp=off ; cp nemo_physical_rtd_old.nc nemo_physical_rtd.nc ; chmod +w nemo_physical_rtd.nc
      access nemo_ice_rtd_old.nc $ice_rtdfile1 nocp=off ; cp nemo_ice_rtd_old.nc nemo_ice_rtd.nc ; chmod +w nemo_ice_rtd.nc
      if [ "$nemo_carbon" = "on" ] ; then
         access nemo_carbon_rtd_old.nc $carbon_rtdfile1 nocp=off na 
         [ -s nemo_carbon_rtd_old.nc ] && cp nemo_carbon_rtd_old.nc nemo_carbon_rtd.nc 
         chmod +w nemo_carbon_rtd.nc
      fi
    fi

  # Get the RTD executables
    nemo_physical_rtd_exe=nemo_physical_rtd.exe   
    nemo_ice_rtd_exe=nemo_ice_rtd.exe            
    nemo_carbon_rtd_exe=nemo_carbon-cmoc_rtd.exe 
    [[ $nemo_config == *"CMOC" ]] && nemo_carbon_rtd_exe=nemo_carbon-cmoc_rtd.exe
    [[ $nemo_config == *"CANOE" ]] && nemo_carbon_rtd_exe=nemo_carbon-canoe_rtd.exe

    cp ${EXEC_STORAGE_DIR}/${nemo_physical_rtd_exe} .

    cp ${EXEC_STORAGE_DIR}/${nemo_ice_rtd_exe} .

    [ "$nemo_carbon" = "on" ] &&cp ${EXEC_STORAGE_DIR}/${nemo_carbon_rtd_exe} .

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
        rtd_hist6="mc_${runid}_${yearm}_m${mon}_1m_scalar.nc"

        access grid_t $rtd_hist1 na #|| bail "NEMO rdt cannot access $rtd_hist1"
        access grid_u $rtd_hist2 na #|| bail "NEMO rdt cannot access $rtd_hist2"
        access grid_v $rtd_hist3 na #|| bail "NEMO rdt cannot access $rtd_hist3"
        access grid_w $rtd_hist4 na #|| bail "NEMO rdt cannot access $rtd_hist4"
        access orca_mesh_mask $rtd_hist5 na
        access scalar $rtd_hist6 na

        # Access additional annual history files containing ICE related variables (needed for qsr_ice anbd qns_ice)
        rtd_hist7="mc_${runid}_${yearm}_m${mon}_1m_icemod.nc"
        access icemod $rtd_hist7 na

        # Create run time diagnostics for physical ocean variables
        $nemo_physical_rtd_exe ${yearm} ${mon} 

                     # Sea-ice run time diagnostics

        # Create run time diagnostics for ice variables
        $nemo_ice_rtd_exe ${yearm} ${mon}

        if [ "$nemo_carbon" = "on" ]; then
                     # Carbon run time diagnostics
          # Access additional annual history files containing PISCES related variables
          rtd_hist8="mc_${runid}_${yearm}_m${mon}_1m_btrc_t.nc"
          rtd_hist9="mc_${runid}_${yearm}_m${mon}_1m_diad_t.nc"
          access ptrc_t $rtd_hist8 na
          access diad_t $rtd_hist9 na

          $nemo_carbon_rtd_exe ${yearm} ${mon}
        fi

        # Clean up
        release grid_t
        release grid_u
        release grid_v
        release grid_w
        release icemod
        release ptrc_t
        release diad_t
        release orca_mesh_mask
        release scalar
     done # month loop end

     # Cleanup
     release orca_mesh_mask

  #   Modify attributes
      [ -s nemo_physical_rtd.nc ] && ncatted -a RTD_version,global,o,c,"${physical_rtd_exe}\n" nemo_physical_rtd.nc
      [ -s nemo_ice_rtd.nc ] && ncatted  -a RTD_version,global,o,c,"${ice_rtd_exe}\n" nemo_ice_rtd.nc
      [ "$nemo_carbon" = "on" ] && [ -s nemo_carbon_rtd.nc ] && ncatted -a RTD_version,global,o,c,"${carbon_rtd_exe}\n" nemo_carbon_rtd.nc

  #   Move RTD files for staging back to Victoria
      rtd_staging nemo_physical_rtd.nc $physical_rtdfile
      rtd_staging nemo_ice_rtd.nc      $ice_rtdfile
      [ "$nemo_carbon" = "on" ] && rtd_staging nemo_carbon_rtd.nc   $carbon_rtdfile

    # Save new RTD files.
      save nemo_physical_rtd.nc $physical_rtdfile
      save nemo_ice_rtd.nc      $ice_rtdfile
      [ "$nemo_carbon" = "on" ] && save nemo_carbon_rtd.nc   $carbon_rtdfile

    # same operation for nemo_carbon
      if [ "$nemo_carbon" = "on" ]; then
        echo ---------------------
      fi

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

  ############# NemoVal diagnostics #############
      if [ $with_nemoval_diag -eq 1 ] ; then
          
      # determine if plotting/merging should be done
      nemoval_mod=$((year % nemoval_repeat_period)) 
      noplt_period=$((nemoval_repeat_period - nemoval_avg_period)) 
      if [ $nemoval_mod -eq 0 ] || [ $nemoval_mod -gt $noplt_period ]; then
            export SKIP_SAME_TIME=1

        # copy the NView settings files and source
        cp ${nemoval_root}/data/NViewConfig.sh .
        . ./NViewConfig.sh

        # determine plot years
        if [ $nemoval_mod -eq 0 ]; then
          pltey=$year
        else
          nval_tmp=$((year + (nemoval_repeat_period - nemoval_mod)))
          pltey=$(printf "%04d" $nval_tmp)
        fi
        # note: the + 1 is required since we want plots to INCLUDE the final year
        pltsy=$(printf "%04d" $((pltey - nemoval_avg_period + 1))) 

        # make sure data should actually exists by confirming the run start year 
        #    is less than the plot start year
        rsy=$start_rtdiag

              # output useful messages that will help for debugging
              echo "Trying to performing nemoval file handling for a plot for years ${pltey} to ${pltsy}"
              echo "Simulation at year ${year} and started at ${rsy}"
        if [ $pltsy -ge $rsy ]; then
                  echo "In a valid nemoval plotting period"
             
          # get previous end year 
            prey=$(printf "%04d" $((year - 1)))

            # set tmp merged file suffix
            new_int_mrg_file_suf="model_rtime_nval_tmp_${runid}_${pltsy}_${year}_"
            prev_int_mrg_file_suf="model_rtime_nval_tmp_${runid}_${pltsy}_${prey}_"

                  # get filetypes
                  filetypes=`get_plot_filetypes $nemoval_pltyml`

                  # check if this is a plot year
            if [ $year -eq $pltey ]; then
                      echo "Performing final file handling and launching plot job to front end machine"

                  # get template for nemoval job
                  cp ${nemoval_root}/data/model_rtime_files/model_rtime_template .

                  # perform final .nc handling, time average, and store in nemoval_root
                  for ft in $filetypes; do
                          # make filetype lowercase to match cpld naming convention
                          ft=$(echo $ft | tr '[:upper:]' '[:lower:]')

                          # determine year's input files and access them
                          yr_fls=""
                          for mon in $nemo_rtd_mons; do
                              fl="mc_${runid}_${year}_m${mon}_1m_${ft}.nc"
                              access $fl $fl
                              yr_fls=${yr_fls}${fl}" "
                          done

                          if [ $year -eq $pltsy ]; then
                              # then no merged file has been created yet and plots are just for this year
                              
                              # merge all year's files into one
                              outpt_fl=model_rtime_nval_tmp_${runid}_${pltsy}_${pltey}_${ft}.nc
                              ${CDO_PTH}/cdo mergetime $yr_fls $outpt_fl

                              # clean up
                              release $yr_fls
                          else
                              # perform final merger

                              # access previously merged file
                              inpt_fl="${prev_int_mrg_file_suf}${ft}.nc"
                              access $inpt_fl $inpt_fl

                              # merge
                              outpt_fl="${new_int_mrg_file_suf}${ft}.nc"
                              ${CDO_PTH}/cdo mergetime $inpt_fl $yr_fls $outpt_fl

                              # clean up
                              delete $inpt_fl    # delete previous merge file
                              release $yr_fls     # release files from this year
                          fi

                          # time average
                          inpt_fl="$outpt_fl"
                          outpt_fl="${MRG_DATA_DIR}/mc_${runid}_1m_${pltsy}0101_${pltey}1231_${ft}.nc"
                          ${CDO_PTH}/cdo timmean $inpt_fl $outpt_fl
                          rm $inpt_fl
                      done

                      #---- setup/launch plotting job 
                      # replace name of python environment with its absolute path
                      PY_ENV="/home/$SCI_ID/.conda/envs/$PY_ENV"
                      sed -i "s#PY_ENV=.*#PY_ENV=\"$PY_ENV\"#g" NViewConfig.sh

                      # replace activate exectubale with its full/absolute path
                      ACTV_EXE="$PY_ENV/bin/activate"
                      sed -i "s#ACTV_EXE=.*#ACTV_EXE=\"$ACTV_EXE\"#g" NViewConfig.sh
                          
                      # replace SCI_ID with user running the model
                      SCI_ID=$(whoami)
                      sed -i "s#SCI_ID=.*#SCI_ID=\"$SCI_ID\"#g" NViewConfig.sh

                      # set plotting machine
                      RUN_MACH=`echo $HDNODE | sed 's/eccc-//'`

                      # set mask to be used by the plotting scripts
                      if [ $nemoval_res == "2d" ]; then
                          MASK="orca2_tmask_miss.nc"
                      elif [ $nemoval_res == "1d" ]; then
                          MASK="orca1_tmask_miss.nc"
                      elif [ $nemoval_res == "025d" ]; then
                          MASK="orca025_tmask_miss.nc"
                      else
                          bail "you selected a resolution not supported by NemoVal plots..either decrease resolution or turn off NemoVal plots"
                      fi
                  
                      # get unique PID/DATE string to ID the plotting job 
                      plt_ID=$$_`date +%Y-%m-%d%H%M%S`

                      # produce job script
                      cat model_rtime_template | sed "s#RUNID#$runid#g;s#START#$pltsy#g;\
                      s#END#$pltey#g;s#SCI_ID#$SCI_ID#g;s#NV_USER#canesm_nemoval#g;s#NVAL_ROOT#$nemoval_root#g;\
                      s#PID#$plt_ID#g;s#RUN_MACH#$RUN_MACH#g;s#VERSION#$DEF_NVAL_VER#g;\
                      s#OBS_DIR#${nemoval_root}/data/${DEF_OBS_DIR}#g;s#MASKFILE#$MASK#g;s#PLTYML_PATH/PLOTS#$nemoval_pltyml#g;\
                      s#PLOTS#`basename $nemoval_pltyml`#g;s#OBSERVED#$DEF_OBS_LST#g;s#PY_ENV#$PY_ENV#g;s#ACTV_EXE#$ACTV_EXE#g;\
                      s#FILETYPES#$filetypes#g"\
                      >> canesm_nvaljob_${runid}_${pltsy}_${pltey}

                      # remove useless template and place NViewConfig.sh in NemoVal's staging directory
                      if [ $? -eq 0 ]; then
                          mv NViewConfig.sh ${nemoval_root}/stg/NViewConfig_${plt_ID}_.sh
                          rm model_rtime_template
                      fi
                      
                      # submit job
                      rsub $RUN_MACH canesm_nvaljob_${runid}_${pltsy}_${pltey} > /dev/null 2>&1
                      if [ $? -eq 0 ]; then
                          echo "NemoVal plotting job submitted"
                      fi
                  elif [ $year -eq $pltsy ]; then
                      echo "Creating new merged file"
                      
                      for ft in $filetypes; do
                          # translate ft to lower case
                          ft=$(echo $ft | tr '[:upper:]' '[:lower:]')

                          # determine year's input files and access them
                          yr_fls=""
                          for mon in $nemo_rtd_mons; do
                              fl="mc_${runid}_${year}_m${mon}_1m_${ft}.nc"
                              access $fl $fl
                              yr_fls=${yr_fls}${fl}" "
                          done

                          # merge year's files into one
                          outpt_fl="${new_int_mrg_file_suf}${ft}.nc"
                          ${CDO_PTH}/cdo mergetime $yr_fls $outpt_fl

                          # clean up
                          release $yr_fls
                          save $outpt_fl $outpt_fl
                      done
                  else
                      echo "Appending to existing merged file"
                      for ft in $filetypes; do
                          # translate ft to lower case
                          ft=$(echo $ft | tr '[:upper:]' '[:lower:]')

                          # determine year's input files and access them
                          yr_fls=""
                          for mon in $nemo_rtd_mons; do
                              fl="mc_${runid}_${year}_m${mon}_1m_${ft}.nc"
                              access $fl $fl
                              yr_fls=${yr_fls}${fl}" "
                          done

                          # set files
                          inpt_fl="${prev_int_mrg_file_suf}${ft}.nc" 
                          access $inpt_fl $inpt_fl
                          outpt_fl="${new_int_mrg_file_suf}${ft}.nc"

                          ${CDO_PTH}/cdo mergetime $inpt_fl $yr_fls $outpt_fl

                          # clean up
                          release $yr_fls
                          delete $inpt_fl
                          save $outpt_fl $outpt_fl; rm $outpt_fl
                      done
                  fi
              fi
          fi
       fi
    fi
