#########################################################
# nemo diagnostics & time mean (1d -> 1m) & time series
# D. Yang, Nov 2018, A. Shao, S.Kharin
#
# This script is sourced in CanESM/CCCma_tools/cccjob_dir/
# lib/jobdefs/canesm_nemo_bgc_diag_jobdef
#########################################################

set -x
set -e

# NEMO priority level
  output_level=${output_level}

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
  yearm1=`echo $year | awk '{printf "%04d", $1 - 1}'`

# Access file containing grid information
  mask_mon=$(echo $nemo_rtd_mons | awk '{printf "%02d", $NF}')  # get last element of nemo_rtd_mons, printed as 2 digit number
  orca_grid_info=mc_${runid}_${year}_m${mask_mon}_mesh_mask.nc
  [ -s orca_mesh_mask ] || access orca_mesh_mask $orca_grid_info

# sfxlst is a suffix list for some nemo historical files.
  sfxlst="1m_grid_t 1m_diad_t 1m_ptrc_t"

# Access the history files
  for sfx in $sfxlst ; do
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
    if [ $nmon -gt 1 ] ; then
      cdo mergetime ${sfx}_?? ${sfx}_m$fmon
      # mergetime changes tbnds to bnds, which causes nemo_diag.exe to crash.
      ncrename -d bnds,tbnds ${sfx}_m$fmon
      rm -f ${sfx}_??
      mv ${sfx}_m$fmon ${sfx}_$fmon
    fi
  done

##########################
# CMIP6 nemo diagnostics #
##########################
  ln -s 1m_grid_t_${fmon} grid_t  || ( echo "Link to grid_t failed" ; exit 1 )
  ln -s 1m_diad_t_${fmon} diad_t  || ( echo "Link to diad_t failed" ; exit 1 )
  ln -s 1m_ptrc_t_${fmon} ptrc_t  || ( echo "Link to ptrc_t failed" ; exit 1 )

  if [[ $nemo_config == *'CMOC'* && ${output_level} -gt 0 ]]; then
    # Expected outputs from CMOC or CanOE offline diagnostics
    cmoc_outvars_l1="Zsat_A Zsat_C o2min zo2min o2sol pH3D pHabio pHnat"
    cmoc_outvars_l2="CO3 CO3abio CO3nat CO3sata CO3satc"
    cmoc_outvars_l8="Omega_A Omega_A_abio Omega_A_nat Omega_C Omega_C_abio Omega_C_nat"
    case ${output_level} in
         1) cmoc_outvars="${cmoc_outvars_l1}"                          ;;
         2) cmoc_outvars="${cmoc_outvars_l1} ${cmoc_outvars_l2}"       ;;
     [3-7]) cmoc_outvars="${cmoc_outvars_l1} ${cmoc_outvars_l2}"       ;;
         8) cmoc_outvars="${cmoc_outvars_l1} ${cmoc_outvars_l2} ${cmoc_outvars_l8}" ;;
         *) cmoc_outvars="${cmoc_outvars_l1} ${cmoc_outvars_l2} ${cmoc_outvars_l8}" ;;
    esac

    # Compile the diagnostic program
    WRKDIR=$PWD
    ( cd $CCRNSRC/CanESM/CanNEMO/diag ;
      ifort -o $WRKDIR/nemo_diag_cmoc.exe nemo_diag_glovars_cmoc.F90 nemo_diag_cal_cmoc.F90 nemo_diag_cmoc.F90 uvic_netcdf.f \
               `nc-config --fflags` `nc-config --flibs`
    )

    # Get all auxiliary files needed before running the offline diagnostics
    access si.nc uncs_orca1_data_si_nomask.nc
    # Get globally averaged surface salinity from previous year
    if [ -L rsp ] ; then
       mkdir dir_rsp; cd dir_rsp
       tar -xvf ../rsp
       ncks -v sss_glob_avg *_restart_trc.nc ../sss_glob_avg.nc
       cd ..
       release rsp
       rm -f -r dir_rsp
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
    cmoc_annual_ptrc=""
    if [[ ${output_level} -gt 3 ]]; then
        cmoc_annual_ptrc="cfc11 cfc12 sf6"
    fi
    cmoc_annual_ptrc+=" nchl di14c dic dicabio dicnat no3 o2 alkalini poc zoo phy"
    cmoc_annual_diad="cflx_14c cflx cflx_abio cflx_nat ph3d phabio phnat co3 co3sata co3satc ppphy co3abio co3nat o2sol"

    cmoc_src_file=sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_1m_ptrc_t
    cmoc_dest_file=sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_1y_ptrc_t
    for f in ${cmoc_annual_ptrc}; do
      release tmp.nc
      access  tmp.nc ${cmoc_src_file}_$f.nc na
      if [ -s tmp.nc ] ; then
        cdo yearmean tmp.nc ann_$f.nc
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
        cdo yearmean tmp.nc ann_$f.nc
        save ann_$f.nc ${cmoc_dest_file}_$f.nc
        release tmp.nc
      fi
    done

  # Similar but for CANOE configurations
  elif [[ $nemo_config == *'CANOE'* && ${output_level} -gt 0 ]]; then
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

    # Compile the diagnostic program
    WRKDIR=$PWD
    ( cd $CCRNSRC/CanESM/CanNEMO/diag ;
      ifort -o $WRKDIR/nemo_diag_canoe.exe nemo_diag_glovars_canoe.F90 nemo_diag_cal_canoe.F90 nemo_diag_canoe.F90 uvic_netcdf.f \
                `nc-config --fflags` `nc-config --flibs`
    )

    # Get all auxiliary files needed before running the offline diagnostics
    access si.nc uncs_orca1_data_si_nomask.nc
    # Run the offline diagnostics
    ./nemo_diag_canoe.exe

# Extract some grid variables from 1m_diad_t
    ncks -v nav_lat,nav_lon 1m_diad_t_${fmon} 1m_diad_t_header

# Save time series
    for f in $canoe_outvars; do
      ncks -A 1m_diad_t_header $f.nc
      save $f.nc sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_${canoe_destfile}_${f}.nc
    done

    canoe_annual_ptrc="nchl dchl dfe dic no3 o2 talk caco3 nh4 zoo zoo2 phy2c phyc phyn phy2n phyfe phy2fe poc goc"
    canoe_annual_diad="cflx ph3d co3 co3sata co3satc dcal graz1 graz2 pfen pfed pcal ppphy ppphy2 o2sol irondep"

    canoe_src_file=sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_1m_ptrc_t
    canoe_dest_file=sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_1y_ptrc_t
    for f in ${canoe_annual_ptrc}; do
      release tmp.nc
      access  tmp.nc ${canoe_src_file}_$f.nc na
      if [ -s tmp.nc ] ; then
        cdo yearmean tmp.nc ann_$f.nc
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
        cdo yearmean tmp.nc ann_$f.nc
        save ann_$f.nc ${canoe_dest_file}_$f.nc
        release tmp.nc
      fi
    done

  fi
