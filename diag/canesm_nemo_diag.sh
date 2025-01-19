#########################################################
# nemo diagnostics & time mean (1d -> 1m) & time series
# D. Yang, Nov 2018, A. Shao, S.Kharin
#
# This script is sourced in CanESM/CCCma_tools/cccjob_dir/
# lib/jobdefs/canesm_nemo_diag_jobdef
#########################################################

set -e

# Note that nemo_rtd_mons used below is first month of the time chunk. 
# nemo_rtd_mons=1 for a run starting from January in a single 12-month chunk;
# nemo_rtd_mons=6 for a run starting from June in a single 12-month chunk;
# nemo_rtd_mons='1 7' for a run starting from January in two 6-month chunks;
# nemo_rtd_mons='6 12' for a run starting from June in two 6-month chunks.

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

# Previous year
  if [[ "$year" == "$run_start_year" ]] && [[ $nemo_from_rest == 'on' ]]; then
    # use the current year because output.init.nc is used in that case (below)
    yearm1=`echo $year | awk '{printf "%04d", $1}'`
    file_state="initial"
    t_state="votemper"
    s_state="vosaline"
  else
    yearm1=`echo $year | awk '{printf "%04d", $1 - 1}'`
    file_state="restart"
    t_state="tn"
    s_state="sn"
  fi

# copy in the nemo diag executable
  diag_exe=nemo_diag.exe
  cp ${EXEC_STORAGE_DIR}/${diag_exe} .

# Access file containing grid information
  mask_mon=$(echo $nemo_rtd_mons | awk '{printf "%02d", $1}')  # get last element of nemo_rtd_mons, printed as 2 digit number
  orca_grid_info=mc_${runid}_${fyear}_m${mask_mon}_mesh_mask.nc
  [ -s orca_mesh_mask ] || access orca_mesh_mask $orca_grid_info || bail "Failed to access $orca_grid_info"

# Access file containing mfo line mask
  [ -s mfo_line_mask ] || access mfo_line_mask mfo_line_mask || bail "Failed to access $mfo_line_mask"

# suffix list for sub-yearly nemo historical files.
  nemo_diag_file_suffix_list=${nemo_diag_file_suffix_list}

# suffix list for yearly nemo historical files.
  nemo_diag_file_1y_suffix_list=${nemo_diag_file_1y_suffix_list}

# Access the history files
  for sfx in $nemo_diag_file_suffix_list $nemo_diag_file_1y_suffix_list ; do
    yr=$fyear
    mp=0
    for mm in $nemo_rtd_mons ; do
      if [ $mm -lt $mp ] ; then
        # increment year by 1 if the current month is smaller than the previous month
        yr=`echo $yr | awk '{printf "%04d", $1 + 1}'`;
      fi
      diag_hist="mc_${runid}_${yr}_m${mm}_${sfx}.nc"
      access ${sfx}_${mm} $diag_hist 
      mp=$mm
    done
# Merge sub-yearly files
    if [ $nmon -gt 1 ] ; then
      cdo mergetime ${sfx}_?? ${sfx}_m$fmon
      rm -f ${sfx}_??
      mv ${sfx}_m$fmon ${sfx}_$fmon
    fi
  done

##########################################################################
# 1. Access input files/variables                                        #
# 2. Run the Fortran executable to compute the CMIP6 offline diagnostics #
# 3. Process *diaptr* files                                              #
##########################################################################
#
# Note that the following two "for and case" structures translates 
# "output_leve" to "level" that represent the variable priority levels. 
# For example, output_leve=3 means variable priority levels 1, 2 and 3.
#
# Execute the following lines when output_level -ge 1;
  for level in $(seq 1 $output_level); do	  
    case $level in
      # output_level=1 and only if starting from January
      1)
        # access input variables for computing vars with priority level 1
        
        if [[ $fmon -eq 1 && $nemo_calc_diag == 1 ]] ; then
          ln -sf 1m_grid_t_${fmon} grid_t  || bail "Link to grid_t failed"
          ln -sf 1m_grid_u_${fmon} grid_u  || bail "Link to grid_u failed"
          ln -sf 1m_grid_v_${fmon} grid_v  || bail "Link to grid_v failed"

          ################################################################
          # Run the CMIP6 nemo offline diagnostics executable: $diag_exe #
          ################################################################
          # make sure inputs exist and run!
          if [[ -L grid_t && $output_level -eq 1 ]]; then
            $diag_exe
          elif [ ! -L grid_t ]; then
            bail "Inputs for $diag_exe (grid_t) don't exist!"
          fi
        fi
        ######################################
        # Time mean (1d_diaptr -> 1m_diaptr) #
        ######################################
        [[ -L 1d_diaptr_${fmon} || -s 1d_diaptr_${fmon} ]] && cdo -b F64 monmean 1d_diaptr_${fmon} 1m_diaptr_${fmon}
        # Replace 1d_diaptr with 1m_diaptr after doing time mean
        nemo_diag_file_suffix_list=`echo $nemo_diag_file_suffix_list | sed -e "s/1d_diaptr/1m_diaptr/"`
        ;;
      # output_level=3 and only if starting from January
      3)
        if [ $nemo_calc_diag -eq 1 ] ; then
          if [ $fmon -eq 1 ] ; then
            # Run offline computation of tendency terms only if starting from January and yearly chunk
            # Access the nemo restart files
            diag_rs1="mc_${runid}_${yearm1}_m${lmon}_nemors" # previous year
            diag_rs2="mc_${runid}_${year}_m${lmon}_nemors"   # current year
            access rsp $diag_rs1 || bail "Failed to access $diag_rs1"
            access rsc $diag_rs2 || bail "Failed to access $diag_rs2"

            # Get tn and sn from the last step of previous year
            if [ -L rsp ] ; then
              work_dir=$(pwd)
              cd rsp
              cdo select,name=$t_state,timestep=-1 *_$file_state.nc ${work_dir}/tnp.nc
              cdo select,name=$s_state,timestep=-1 *_$file_state.nc ${work_dir}/snp.nc
              cd $work_dir
              release rsp
            fi

            # Get tn and sn from the last step of current year
            if [ -L rsc ] ; then
              work_dir=$(pwd)
              cd rsc
              cdo select,name=tn,timestep=-1 ${runid}_*_restart.nc ${work_dir}/tnc.nc
              cdo select,name=sn,timestep=-1 ${runid}_*_restart.nc ${work_dir}/snc.nc
              cd $work_dir
              release rsc
            fi
            ################################################################
            # Run the CMIP6 nemo offline diagnostics executable: $diag_exe #
            ################################################################
            # make sure inputs exist and run!
            if [[ -L grid_t && -s tnp.nc ]]; then
              $diag_exe
            else
              bail "Inputs for $diag_exe (grid_t and tnp.nc) don't exist!"
            fi
          fi
        fi
        ;;
    esac
  done

  if [ $nemo_calc_diag -eq 1 ] ; then
    for level in $(seq 1 $output_level); do
      case $level in
        1)
          if [ $fmon -eq 1 ] ; then
          # Append mfo.nc to 1m_scalar_ar6_${fmon} if existing
            if [ -s mfo.nc ]; then
              chmod u+w mfo.nc
              cp 1m_scalar_ar6_${fmon} 1m_scalar_ar6.nc && chmod u+w 1m_scalar_ar6.nc || bail "1m_scalar_ar6_${fmon} does not exist"
              ncks -A mfo.nc 1m_scalar_ar6.nc
              rm -f 1m_scalar_ar6_${fmon}
              mv 1m_scalar_ar6.nc 1m_scalar_ar6_${fmon}
            else
              bail "mfo.nc does not exist"
            fi
            # Append msftbarot.nc to 1m_grid_u_ar6_${fmon} if existing
            if [ -s msftbarot.nc ]; then
              chmod u+w msftbarot.nc
              cp 1m_grid_u_ar6_${fmon} 1m_grid_u_ar6.nc && chmod u+w 1m_grid_u_ar6.nc || bail "1m_grid_u_ar6_${fmon} does not exist"
              ncks -A msftbarot.nc 1m_grid_u_ar6.nc
              rm -f 1m_grid_u_ar6_${fmon}
              mv 1m_grid_u_ar6.nc 1m_grid_u_ar6_${fmon}
            else
              bail "msftbarot.nc does not exist"
            fi
          fi 
          ;;
        3)
          if [[ ${nmon} -eq 1 && $fmon -eq 1 ]] ; then
            # Append tstend.nc to 1y_grid_t_ar6_${fmon} if existing
            if [ -s tstend.nc ]; then
              chmod u+w tstend.nc
              #cp 1y_grid_t_ar6_${fmon} 1y_grid_t_ar6.nc && chmod u+w 1y_grid_t_ar6.nc || bail "1y_grid_t_ar6_${fmon} does not exist"
              #ncks -A tstend.nc 1y_grid_t_ar6.nc
              #rm -f 1y_grid_t_ar6_${fmon}
              #mv 1y_grid_t_ar6.nc 1y_grid_t_ar6_${fmon}
              # Append yearly diagnostics suffix list
              chmod u+w 1y_grid_t_ar6_${fmon}
              ncks -A tstend.nc 1y_grid_t_ar6_${fmon}
            else
              bail "tstend.nc does not exist"
            fi
          fi
          ;;
      esac 
    done
  fi

##################################################
# Split historical files to time series and save #
##################################################
  # split to time series
  for sfx in $nemo_diag_file_suffix_list $nemo_diag_file_1y_suffix_list ; do
    [ ! -e ${sfx}_${fmon} ] && continue
    # Get the list of variable (use "coordinates" as the key word because it won't apply to the coordinates variables, e.g. nav_lat)
    vars_list=$( ncdump -h ${sfx}_${fmon} | grep coordinates | awk -F: '{print $1}' )
    for var in $vars_list; do
        # produce one file per variables
        ncks -v $var ${sfx}_${fmon} xxx-${sfx}_${var}.nc || ncks -v $var ${sfx}_${fmon} xxx-${sfx}_${var}.nc
        # Save time series
        save xxx-${sfx}_${var}.nc sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_${sfx/_ar6/}_${var}.nc
        release xxx-${sfx}_${var}.nc
    done
    rm  ${sfx}_${fmon}
  done

  # Save orca grid mask with consistent name as TS files
  save orca_mesh_mask sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_mesh_mask.nc

exit 0
