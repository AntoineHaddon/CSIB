#########################################################
# nemo diagnostics & time mean (1d -> 1m) & time series
# D. Yang, Nov 2018, A. Shao, S.Kharin
#
# This script is sourced in CanESM/CCCma_tools/cccjob_dir/
# lib/jobdefs/canesm_nemo_bgc_diag_jobdef
#########################################################

set -x
set -e
# Get CDO / TEMPORARY!
export PATH=/fs/ssm/hpco/exp/mib002/anaconda/anaconda-4.4.0/anaconda_4.4.0_ubuntu-14.04-amd64-64/envs/cdo-1.9.0/bin:$PATH
export LD_LIBRARY_PATH=/fs/ssm/hpco/tmp/eccc/201402/04/intel-2016.1.150/ubuntu-14.04-amd64-64/lib/:$LD_LIBRARY_PATH

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
  orca_grid_info=nemo_mesh_mask_rc3.nc
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
      . ssmuse-sh -d /fs/ssm/hpco/tmp/eccc/201402/03/base  -d main/opt/intelcomp/intelcomp-2016.1.156 ;
      ifort -o $WRKDIR/nemo_diag_cmoc.exe nemo_diag_glovars_cmoc.F90 nemo_diag_cal_cmoc.F90 nemo_diag_cmoc.F90 uvic_netcdf.f \
            -I/fs/ssm/hpco/tmp/eccc/201402/04/intel-2016.1.150/ubuntu-14.04-amd64-64/include/                                \
            -L/fs/ssm/hpco/tmp/eccc/201402/04/intel-2016.1.150/ubuntu-14.04-amd64-64/lib/ -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl;
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

# # SK: comment out until someone tells it's needed
#     # Merge all CMOC variables into a single file and save it
#     for f in $cmoc_outvars; do
#       ncks -A $f.nc $cmoc_destfile
#     done
#     save    $cmoc_destfile mc_${runid}_${year}_m${fmon}_${cmoc_destfile}.nc
#     release $cmoc_destfile

# Extract some grid variables from 1m_diad_t
    ncks -v nav_lat,nav_lon 1m_diad_t_${fmon} 1m_diad_t_header

# Save time series
    for f in $cmoc_outvars; do
      ncks -A 1m_diad_t_header $f.nc
      save $f.nc sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_${cmoc_destfile}_${f}.nc
    done

  # Similar but for CANOE configurations
  elif [[ $nemo_config == *'CANOE'* && ${output_level} -gt 0 ]]; then
    # Expected outputs from CMOC or CanOE offline diagnostics
    canoe_outvars_l1="o2sol pH3D"
    canoe_outvars_l2="CO3 CO3sata CO3satc"
    case ${output_level} in
         1) canoe_outvars="${canoe_outvars_l1}"                          ;;
         2) canoe_outvars="${canoe_outvars_l1} ${canoe_outvars_l2}"       ;;
         *) canoe_outvars="${canoe_outvars_l1} ${canoe_outvars_l2} ${canoe_outvars_l8}" ;;
    esac

    # Compile the diagnostic program
    WRKDIR=$PWD
    ( cd $CCRNSRC/CanESM/CanNEMO/diag ;
      . ssmuse-sh -d /fs/ssm/hpco/tmp/eccc/201402/03/base  -d main/opt/intelcomp/intelcomp-2016.1.156 ;
      ifort -o $WRKDIR/nemo_diag_cmoc.exe nemo_diag_glovars_cmoc.F90 nemo_diag_cal_cmoc.F90 nemo_diag_canoe.F90 uvic_netcdf.f \
            -I/fs/ssm/hpco/tmp/eccc/201402/04/intel-2016.1.150/ubuntu-14.04-amd64-64/include/                                \
            -L/fs/ssm/hpco/tmp/eccc/201402/04/intel-2016.1.150/ubuntu-14.04-amd64-64/lib/ -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl;
    )

    # Get all auxiliary files needed before running the offline diagnostics
    access si.nc uncs_orca1_data_si_nomask.nc
    # Run the offline diagnostics
    ./nemo_diag_canoe.exe

# Extract some grid variables from 1m_diad_t
    ncks -v nav_lat,nav_lon 1m_diad_t_${fmon} 1m_diad_t_header

# Save time series
    for f in $cmoc_outvars; do
      ncks -A 1m_diad_t_header $f.nc
      save $f.nc sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_${cmoc_destfile}_${f}.nc
    done


  fi
