#!/bin/bash

set -e

# get copy of the domain.cfg file
is_defined $nemo_coordinates || bail "The variable nemo_coordinates must be defined in the configuration file."
access domain_cfg.nc $nemo_coordinates

# activate correct Python environment
source activate /home/scrd102/cccma_conda/envs/py3_analysis_v2

# determine loop parameters/iteration info
python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -y $(($NEMO_CHUNK_START_YEAR - 1)) -y $(($NEMO_CHUNK_END_YEAR + 1)) -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -t -1 -R ${run_start_year} -r ${run_start_month} -f ${months}
# Source the file with this info
. nemo_counter_info.cfg

# initial conditions
if [[ $ctds_dnscl != 0 ]] && [[ $NEMO_CHUNK_START_DATE == $run_start_date ]] ; then
  if [[ $runmode == *"CanTODS"* ]]; then
    if [[ $nemo_from_rest == "on" ]] ; then
      # Start from rest
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -y ${dada_ic_year} -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -t 0 -i $(( $dada_ic_month - 1 )) > ic_status
    else
      # hot start; interpolate restart files from parent run
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -y ${dada_ic_year} -P ${dada_parent_path%%nc_output*} -p ${dada_parent_name/CanESM5-NEMO4-/} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -t 0 -i $(( $dada_ic_month - 1 )) -H 1 -V 1 > ic_status
    fi
  else
    # eORCA grids slightly different, so interpolate to an existing IC file
    access tmp_ic.nc $nemo_data_1m_temperature_rest
    if [[ $nemo_from_rest == "on" ]] ; then
      # Start from rest
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -y ${dada_ic_year} -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m tmp_ic.nc -t 0 -i $(( $dada_ic_month - 1 )) > ic_status
    else
      # hot start; interpolate restart files
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -y ${dada_ic_year} -P ${dada_parent_path%%nc_output*} -p ${dada_parent_name/CanESM5-NEMO4-/} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m tmp_ic.nc -t 0 -i $(( $dada_ic_month - 1 )) -H 1 > ic_status
    fi
  fi
  if [[ -z "$(ls ./data_1m_*_nomask.nc)" && -z "$(ls ./restart*.nc)" ]] ; then
    echo "ERROR: No IC files generated!"
    exit 29
  fi
fi

# boundary conditions
if [[ $ctds_dnscl != 0 ]] && [[ $runmode == *"CanTODS"* ]] ; then
  if [[ -z "${dada_year_offset}" ]] && [[ -z "${dada_loop_year}" ]] ; then
      # use current year for forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o obc_BDY_${dada_outfield}_yYYYY.nc -y $(($NEMO_CHUNK_START_YEAR - 1)) -y $(($NEMO_CHUNK_END_YEAR + 1)) -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -t 1 > bc_status
  else
      # use cyclical forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o obc_BDY_${dada_outfield}_yYYYY.nc -y $(($NEMO_CHUNK_START_YEAR - 1)) -y $(($NEMO_CHUNK_END_YEAR + 1)) -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -t 1 -A ${dada_year_offset} -a ${dada_loop_year} > bc_status
  fi
  if [ -z "$(ls ./obc_*_${dada_outfield}_y*.nc)" ] ; then
    echo "ERROR: No OBC files generated!"
    exit 10
  fi
fi

# atmospheric forcing if downscaling
# if not downscaling (i.e., forced CanTODS run, or using OMIP forcing),
# then forcing is obtained in nemo_prelude and weighting is used
if [[ $ctds_dnscl != 0 ]] || [[ $runmode != *"CanTODS"* ]]; then
  # CanESM atmospheric forcing
  if [[ -z "${iaf_year_offset}" ]] && [[ -z "${iaf_loop_year}" ]] ; then
      # use current year for forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o VAR_${dada_outfield}_yYYYY -y $(($NEMO_CHUNK_START_YEAR - 1)) -y $(($NEMO_CHUNK_END_YEAR + 1)) -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -M ${dada_forcing_freq} -t 2 > frc_status
  else
      # use cyclical forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o VAR_${dada_outfield}_yYYYY -y $(($NEMO_CHUNK_START_YEAR - 1)) -y $(($NEMO_CHUNK_END_YEAR + 1)) -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -M ${dada_forcing_freq} -t 2 -A ${iaf_year_offset} -a ${iaf_loop_year} > frc_status
  fi
  if [ -z "$(ls ./tair_${dada_outfield}_y*)" ] ; then
    echo "ERROR: Not all forcing files generated!"
    exit 92
  fi
fi

# sea-surface salinity
if [[ $nemo_nn_sssr != 0 ]] ; then
  if [[ -z "${dada_year_offset}" ]] && [[ -z "${dada_loop_year}" ]] ; then
      # use current year for forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o sss_${dada_outfield}_yYYYY -y $(($NEMO_CHUNK_START_YEAR - 1)) -y $(($NEMO_CHUNK_END_YEAR + 1)) -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -M ${dada_forcing_freq} -t 4 > sos_status
  else
      # use cyclical forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o sss_${dada_outfield}_yYYYY -y $(($NEMO_CHUNK_START_YEAR - 1)) -y $(($NEMO_CHUNK_END_YEAR + 1)) -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -M ${dada_forcing_freq} -t 4 -A ${dada_year_offset} -a ${dada_loop_year} > sos_status
  fi
  if [ -z "$(ls ./sss_${dada_outfield}_y*)" ] ; then
    echo "ERROR: SSS file not generated!"
    exit 2007
  fi
fi

# river forcing
# TODO: rivers will need to be remapped properly
if [[ $nemo_ln_rnf == "on" ]] && [[ $runmode != *"CanTODS"* ]]; then
  if [[ -z "${nemo_river_remap}" ]] ; then
   # don't remap/rescale rivers
    if [[ -z "${iaf_year_offset}" ]] && [[ -z "${iaf_loop_year}" ]] ; then
      # use current year for forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o rvr_${dada_outfield} -y $NEMO_CHUNK_START_YEAR -y $NEMO_CHUNK_END_YEAR -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble} -m domain_cfg.nc -t 3 > rvr_status
    else
      # use cyclical forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o rvr_${dada_outfield} -y $NEMO_CHUNK_START_YEAR -y $NEMO_CHUNK_END_YEAR -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble} -m domain_cfg.nc -t 3 -A ${iaf_year_offset} -a ${iaf_loop_year} > rvr_status
    fi
  else
    access nemo_river_remap_${NEMO_CHUNK_START_YEAR}-${NEMO_CHUNK_END_YEAR}.nc $nemo_river_remap
    # remap/rescale rivers
    access dpg.nc $dada_parent_grid || ln -s $dada_parent_grid dpg.nc # get parent grid for remapping
    if [[ -z "${iaf_year_offset}" ]] && [[ -z "${iaf_loop_year}" ]] ; then
      # use current year for forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o rvr_${dada_outfield} -y $NEMO_CHUNK_START_YEAR -y $NEMO_CHUNK_END_YEAR -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble} -m domain_cfg.nc -t 3 -v nemo_river_remap_${NEMO_CHUNK_START_YEAR}-${NEMO_CHUNK_END_YEAR}.nc -g dpg.nc > rvr_status
    else
      # use cyclical forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o rvr_${dada_outfield} -y $NEMO_CHUNK_START_YEAR -y $NEMO_CHUNK_END_YEAR -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble} -m domain_cfg.nc -t 3 -A ${iaf_year_offset} -a ${iaf_loop_year} -v nemo_river_remap_${NEMO_CHUNK_START_YEAR}-${NEMO_CHUNK_END_YEAR}.nc -g dpg.nc > rvr_status
    fi
    rm -f nemo_river_remap_${NEMO_CHUNK_START_YEAR}-${NEMO_CHUNK_END_YEAR}.nc
  fi
  if [ -z "$(ls ./rvr_${dada_outfield}).nc" ] ; then
    echo "ERROR: River files not generated!"
    exit 28
  fi
fi

# delete domain file (will be re-copied in model portion of loop)
rm -f domain_cfg.nc
rm -f tmp_ic.nc
