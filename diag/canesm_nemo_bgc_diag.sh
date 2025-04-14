#########################################################
# nemo diagnostics & time mean (1d -> 1m) & time series
# D. Yang, Nov 2018, A. Shao, S.Kharin
#
# Change annual mean computation from cdo yearmean to
# cdo yearmonmean (weighted) - D. Yang, 23/NOV/2020
#
# This script is sourced in CanESM/CCCma_tools/cccjob_dir/
# lib/jobdefs/canesm_nemo_bgc_diag_jobdef
#########################################################

set -e

# NEMO priority level
  output_level=${output_level}

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
  if [[ $lmon -eq 12 ]] && [[ "$year" == "$run_start_year" ]] && [[ $nemo_from_rest == 'on' ]]; then
    # use the current year because output.init.nc is used in that case (below)
    yearm1=`echo $year | awk '{printf "%04d", $1}'`
    file_state="initial_trc"
  else
    yearm1=`echo $year | awk '{printf "%04d", $1 - 1}'`
    file_state="restart_trc"
  fi

# Access file containing grid information
  mask_mon=$(echo $nemo_rtd_mons | awk '{printf "%02d", $1}')  # get first element of nemo_rtd_mons, printed as 2 digit number
  orca_grid_info=mc_${runid}_${fyear}_m${mask_mon}_mesh_mask.nc
  [ -s orca_mesh_mask ] || access orca_mesh_mask $orca_grid_info

# Access fiels for BGCM diag 
  sfx='1m_grid_t'
  diag_hist="mc_${runid}_${fyear}_m${mask_mon}_${sfx}.nc"
  access ${sfx}_${mask_mon} $diag_hist na
  sfx='1m_diad_t'
  diag_hist="mc_${runid}_${fyear}_m${mask_mon}_${sfx}.nc"
  access ${sfx}_${mask_mon} $diag_hist na
  sfx='1m_btrc_t'
  diag_hist="mc_${runid}_${fyear}_m${mask_mon}_${sfx}.nc"
  access ${sfx}_${mask_mon} $diag_hist na


##########################
# CMIP6 nemo diagnostics #
##########################
  ln -s 1m_grid_t_${fmon} grid_t  || ( echo "Link to grid_t failed" ; exit 1 )
  ln -s 1m_diad_t_${fmon} diad_t  || ( echo "Link to diad_t failed" ; exit 1 )
  ln -s 1m_btrc_t_${fmon} ptrc_t  || ( echo "Link to ptrc_t failed" ; exit 1 ) # OR Jan 09 '23 changed from ptrc_t 

  if [[ $CanNEMO_CONFIG == *'CMOC'* && ${output_level} -gt 0 ]]; then

    process_abio=0       # this flag needs to be set both here and inside nemo_diag_cmoc.F90 (process_abio = .false./.true.)

    # Expected outputs from CMOC or CanOE offline diagnostics
    if [ $process_abio -gt 0 ]; then
      cmoc_outvars_l1="Zsat_A Zsat_C o2min zo2min o2sol pH3D pHabio pHnat"
      cmoc_outvars_l2="CO3 CO3abio CO3nat CO3sata CO3satc"
      cmoc_outvars_l8="Omega_A Omega_A_abio Omega_A_nat Omega_C Omega_C_abio Omega_C_nat"
    else
      cmoc_outvars_l1="Zsat_A Zsat_C o2min zo2min o2sol pH3D pHnat"
      cmoc_outvars_l2="CO3 CO3nat CO3sata CO3satc"
      cmoc_outvars_l8="Omega_A Omega_A_nat Omega_C Omega_C_nat"
    fi
    case ${output_level} in
         1) cmoc_outvars="${cmoc_outvars_l1}"                          ;;
         2) cmoc_outvars="${cmoc_outvars_l1} ${cmoc_outvars_l2}"       ;;
     [3-7]) cmoc_outvars="${cmoc_outvars_l1} ${cmoc_outvars_l2}"       ;;
         8) cmoc_outvars="${cmoc_outvars_l1} ${cmoc_outvars_l2} ${cmoc_outvars_l8}" ;;
         *) cmoc_outvars="${cmoc_outvars_l1} ${cmoc_outvars_l2} ${cmoc_outvars_l8}" ;;
    esac

    # copy in diagnostics exec
    cp ${EXEC_STORAGE_DIR}/nemo_diag_cmoc.exe .

    # Get all auxiliary files needed before running the offline diagnostics
    access si.nc $nemo_data_si_nomask 
    # Get globally averaged surface salinity from previous year
    if [ $process_abio -gt 0 ]; then
      # Access the nemo restart files
      if [ $lmon -eq 12 ]; then
        # previous year
        diag_rs1="mc_${runid}_${yearm1}_m${lmon}_nemors" # previous year
      else
        diag_rs1="mc_${runid}_${year}_m${lmon}_nemors" # previous year
      fi	      
      access rsp $diag_rs1 || ( echo "$diag_rs1 does not exist" ; exit 1 )
    fi
    if [ -L rsp ] ; then
       work_dir=$(pwd)
       cd rsp
       ncks -v sss_glob_avg *_$file_state.nc ${work_dir}/sss_glob_avg.nc
       cd $work_dir
       release rsp
    fi
    # Run the offline diagnostics
    ./nemo_diag_cmoc.exe

    cmoc_destfile="1m_diad_t"

# Extract some grid variables from 1m_diad_t
    ncks -v nav_lat,nav_lon 1m_diad_t_${fmon} 1m_diad_t_header

# Save time series
    for f in $cmoc_outvars; do
      ncks -A 1m_diad_t_header $f.nc
      save $f.nc sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_${cmoc_destfile}_${f}.nc
    done

# Calculate annual averages for select variables
    if [ $fmon -eq 1 ] ; then
      cmoc_annual_ptrc=""
      if [[ ${output_level} -gt 3 ]]; then
        cmoc_annual_ptrc="cfc11 cfc12 sf6"
      fi
      cmoc_annual_ptrc+=" nchl di14c dic dicabio dicnat no3 o2 alkalini poc zoo phy"
      cmoc_annual_diad="cflx_14c cflx cflx_abio cflx_nat ph3d phabio phnat co3 co3sata co3satc ppphy co3abio co3nat o2sol"

      cmoc_src_file=sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_1m_btrc_t  # OR Jan 09 '23 changed from ptrc_t
      cmoc_dest_file=sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_1y_btrc_t # OR Jan 09 '23 changed from ptrc_t
      for f in ${cmoc_annual_ptrc}; do
        release tmp.nc
        access  tmp.nc ${cmoc_src_file}_$f.nc na
        if [ -s tmp.nc ] ; then
          cdo yearmonmean tmp.nc ann_$f.nc
          save ann_$f.nc ${cmoc_dest_file}_$f.nc
          release tmp.nc
        fi
      done

      cmoc_src_file=sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_1m_diad_t
      cmoc_dest_file=sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_1y_diad_t
      for f in ${cmoc_annual_diad}; do
        release tmp.nc
        access  tmp.nc ${cmoc_src_file}_$f.nc na
        if [ -s tmp.nc ] ; then
          cdo yearmonmean tmp.nc ann_$f.nc
          save ann_$f.nc ${cmoc_dest_file}_$f.nc
          release tmp.nc
        fi
      done
    fi

  # Similar but for CANOE configurations
  elif [[ $CanNEMO_CONFIG == *'CANOE'* && ${output_level} -gt 0 ]]; then
    # Expected outputs from CMOC or CanOE offline diagnostics
    canoe_outvars_l1="Zsat_A Zsat_C o2min zo2min o2sol pH3D"
    canoe_outvars_l2="CO3 CO3sata CO3satc"
    canoe_outvars_l8="Omega_C Omega_A"
    case ${output_level} in
         1) canoe_outvars="${canoe_outvars_l1}"                          ;;
         [2-7]) canoe_outvars="${canoe_outvars_l1} ${canoe_outvars_l2}"       ;;
         *) canoe_outvars="${canoe_outvars_l1} ${canoe_outvars_l2} ${canoe_outvars_l8}" ;;
    esac
    canoe_destfile="1m_diad_t"

    # copy in executable
    cp ${EXEC_STORAGE_DIR}/nemo_diag_canoe.exe .

    # Get all auxiliary files needed before running the offline diagnostics
    access si.nc $nemo_data_si_nomask
    # Run the offline diagnostics
    ./nemo_diag_canoe.exe

# Extract some grid variables from 1m_diad_t
    ncks -v nav_lat,nav_lon 1m_diad_t_${fmon} 1m_diad_t_header

# Save time series
    for f in $canoe_outvars; do
      ncks -A 1m_diad_t_header $f.nc
      save $f.nc sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_${canoe_destfile}_${f}.nc
    done

    if [ $fmon -eq 1 ] ; then
      canoe_annual_ptrc="nchl dchl dfe dic no3 o2 talk caco3 nh4 zoo zoo2 phy2c phyc phyn phy2n phyfe phy2fe poc goc"
      canoe_annual_diad="cflx ph3d co3 co3sata co3satc dcal graz1 graz2 pfen pfed pcal ppphy ppphy2 o2sol irondep"

      canoe_src_file=sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_1m_btrc_t  # OR Jan 09 '23 changed from ptrc_t
      canoe_dest_file=sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_1y_btrc_t # OR Jan 09 '23 changed from ptrc_t
      for f in ${canoe_annual_ptrc}; do
        release tmp.nc
        access  tmp.nc ${canoe_src_file}_$f.nc na
        if [ -s tmp.nc ] ; then
          cdo yearmonmean tmp.nc ann_$f.nc
          save ann_$f.nc ${canoe_dest_file}_$f.nc
          release tmp.nc
        fi
      done

      canoe_src_file=sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_1m_diad_t
      canoe_dest_file=sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_1y_diad_t
      for f in ${canoe_annual_diad}; do
        release tmp.nc
        access  tmp.nc ${canoe_src_file}_$f.nc na
        if [ -s tmp.nc ] ; then
          cdo yearmonmean tmp.nc ann_$f.nc
          save ann_$f.nc ${canoe_dest_file}_$f.nc
          release tmp.nc
        fi
      done
    fi

  fi
