#########################################################
# nemo diagnostics & time mean (1d -> 1m) & time series
# D. Yang, Nov 2018, A. Shao, S.Kharin
#
# This script is sourced in CanESM/CCCma_tools/cccjob_dir/
# lib/jobdefs/canesm_nemo_diag_jobdef
#########################################################

set -x

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
  #last month of the first chunk (for istate)
  if [ $nmon -eq 1 ];then
    lmon0=$lmon
  else
    lmon0=`echo $fmon | awk '{printf "%02d", $2 - 1}'`
  fi

# Previous year
  yearm1=`echo $year | awk '{printf "%04d", $1 - 1}'`

# copy in the nemo diag executable
  diag_exe=nemo_diag.exe
  cp ${EXEC_STORAGE_DIR}/${diag_exe} .

# Access file containing grid information
  mask_mon=$(echo $nemo_rtd_mons | awk '{printf "%02d", $NF}')  # get last element of nemo_rtd_mons, printed as 2 digit number
  orca_grid_info=mc_${runid}_${fyear}_m${mask_mon}_mesh_mask.nc
  [ -s orca_mesh_mask ] || access orca_mesh_mask $orca_grid_info nocp=no

# Access file containing mfo line mask
  [ -s mfo_line_mask ] || access mfo_line_mask mfo_line_mask

# suffix list for sub-yearly nemo historical files.
  nemo_diag_file_suffix_list=${nemo_diag_file_suffix_list}

# suffix list for yearly nemo historical files.
  nemo_diag_file_1y_suffix_list=${nemo_diag_file_1y_suffix_list}

# Append yearly diagnostics
  if [ $nmon -eq 1 -a $fmon -eq 1 ] ; then
    nemo_diag_file_suffix_list="$nemo_diag_file_suffix_list $nemo_diag_file_1y_suffix_list"
  fi

# access to time series
  for sfx in $nemo_diag_file_suffix_list ; do
    diag_hist="mc_${runid}_${fyear}_m${fmon}_${sfx}.nc"
    access ${sfx}_${fmon} $diag_hist na
  done

# Execute the following lines when output_level -ge 1
  if [ $output_level -ge 1 ] ; then

      # Run offline computation only if starting from January 
      if [ $fmon -eq 1 ] ; then 
# Access the nemo restart files
        diag_rs1="mc_${runid}_${fyear}_m${lmon0}_nemors" # First restart of that year (for init.nc)
        diag_rs2="mc_${runid}_${lyear}_m${lmon}_nemors" # last restart of that year (for restart.nc)
        access rsp $diag_rs1 || ( echo "$diag_rs1 does not exist" ; exit 1 )
        access rsc $diag_rs2 || ( echo "$diag_rs2 does not exist" ; exit 1 )

# Get tn and sn from the last step of previous year
        if [ -L rsp ] ; then
          cd rsp
          cdo select,name=votemper,timestep=-1 *_istate.nc tnp.nc
          cdo select,name=vosaline,timestep=-1 *_istate.nc snp.nc
          ncrename -O -v votemper,tn tnp.nc tnp.nc
          ncrename -O -v vosaline,sn snp.nc snp.nc
          cd ..
          mv rsp/tnp.nc rsp/snp.nc ./
          release rsp
          rm -f -r dir_rsp
        fi

# Get tn and sn from the last step of current year
        if [ -L rsc ] ; then
          cd rsc
          cdo select,name=tn,timestep=-1 ${runid}_*_restart.nc tnc.nc
          cdo select,name=sn,timestep=-1 ${runid}_*_restart.nc snc.nc
          cd ..
          mv rsc/tnc.nc rsc/snc.nc ./
          release rsc
          rm -f -r dir_rsc
        fi

##############################
# Run CMIP6 nemo diagnostics #
##############################
        ln -s 1m_grid_t_${fmon} grid_t  || bail "Link to grid_t failed"
        ln -s 1m_grid_u_${fmon} grid_u  || bail "Link to grid_u failed"
        ln -s 1m_grid_v_${fmon} grid_v  || bail "Link to grid_v failed"

        # make sure inputs exist, and run!
        if [[ -L grid_t ]] && [[ -s tnp.nc ]]; then
          $diag_exe
        else
          bail "Inputs for $diag_exe (grid_t and tnp.nc) don't exist!"
        fi

        [ -s mfo.nc ] && chmod u+w mfo.nc || bail "mfo.nc does not exist"
        [ -s msftbarot.nc ] && chmod u+w msftbarot.nc || bail "msftbarot.nc does not exist"
        [ -s tstend.nc ] && chmod u+w tstend.nc || bail "tstend.nc does not exist"

# Append mfo.nc to 1m_scalar_ar6_${fmon}
        if [ -e 1m_scalar_ar6_${fmon} ]; then
          cp 1m_scalar_ar6_${fmon} 1m_scalar_ar6.nc && chmod u+w 1m_scalar_ar6.nc || bail "1m_scalar_ar6_${fmon} does not exist"
          ncks -A mfo.nc 1m_scalar_ar6.nc
          rm -f 1m_scalar_ar6_${fmon}
          mv 1m_scalar_ar6.nc 1m_scalar_ar6_${fmon}
        else
          echo  "1m_scalar_ar6_${fmon} does not exist"
        fi

# Append msftbarot.nc to 1m_grid_u_ar6_${fmon}
        if [ -e 1m_grid_u_ar6_${fmon} ]; then
          cp 1m_grid_u_ar6_${fmon} 1m_grid_u_ar6.nc && chmod u+w 1m_grid_u_ar6.nc || bail "1m_grid_u_ar6_${fmon} does not exist"
          ncks -A msftbarot.nc 1m_grid_u_ar6.nc
          rm -f 1m_grid_u_ar6_${fmon}
          mv 1m_grid_u_ar6.nc 1m_grid_u_ar6_${fmon}
        else
          echo  "1m_grid_u_ar6_${fmon} does not exist, created with msftbarot.nc"
          mv msftbarot.nc 1m_grid_u_ar6_${fmon}
        fi

# Append tstend.nc to 1y_grid_t_ar6_${fmon}
        if [ -e 1y_grid_t_ar6_${fmon} ]; then
          if [ ${nmon} -eq 1 -a $fmon -eq 1 ] ; then
            cp 1y_grid_t_ar6_${fmon} 1y_grid_t_ar6.nc && chmod u+w 1y_grid_t_ar6.nc || bail "1y_grid_t_ar6_${fmon} does not exist"
            ncks -A tstend.nc 1y_grid_t_ar6.nc
            rm -f 1y_grid_t_ar6_${fmon}
            mv 1y_grid_t_ar6.nc 1y_grid_t_ar6_${fmon}
          fi
        else
          echo "1y_grid_t_ar6_${fmon} does not exist"
        fi
      fi

  fi # end of "output_level -ge 1"

#########################################
# Split historical files to time series #
#########################################


# Split to time series
  for sfx in $nemo_diag_file_suffix_list ; do
    [ ! -e ${sfx}_${fmon} ] && continue
    ncks -O -C -x -v time_centered_bounds,time_centered ${sfx}_${fmon} ${sfx}_${fmon} 
    cdo splitname ${sfx}_${fmon} xxx-${sfx}_ || true
    # UGLY PATCH : Spetial treatments for diaptr (5D-variables not suported) || true to not cause error if no variable with that name (nil001, july 2023)
    if [ "${sfx}_${fmon}" == "1m_diaptr_01" ]; then
      ncks -v znltem  ${sfx}_${fmon} xxx-${sfx}_znltem.nc && ncrename -O -v time_counter_bounds,time_counter_bnds xxx-${sfx}_znltem.nc xxx-${sfx}_znltem.nc || true
      ncks -v znlsal  ${sfx}_${fmon} xxx-${sfx}_znlsal.nc && ncrename -O -v time_counter_bounds,time_counter_bnds xxx-${sfx}_znlsal.nc xxx-${sfx}_znlsal.nc || true
      ncks -v znlsrf  ${sfx}_${fmon} xxx-${sfx}_znlsrf.nc && ncrename -O -v time_counter_bounds,time_counter_bnds xxx-${sfx}_znlsrf.nc xxx-${sfx}_znlsrf.nc || true
      ncks -v msftyz  ${sfx}_${fmon} xxx-${sfx}_msftyz.nc && ncrename -O -v time_counter_bounds,time_counter_bnds xxx-${sfx}_msftyz.nc xxx-${sfx}_msftyz.nc || true

    fi
    rm  ${sfx}_${fmon}
  done

# Save time series
  tslist=`ls -1 xxx-*`
  for ts in $tslist ; do
    tssfx=`echo $ts |cut -f 2 -d '-' |sed 's/_ar6//'`
    save $ts sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_${tssfx}
    release $ts
  done

# Save orca grid mask with consistent name as TS files
  save orca_mesh_mask sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_mesh_mask.nc
