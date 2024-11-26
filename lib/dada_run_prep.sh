#!/bin/bash

# get copy of the domain.cfg file
is_defined $nemo_coordinates || bail "The variable nemo_coordinates must be defined in the configuration file."
acc_cp domain_cfg.nc $nemo_coordinates
if [[ $ctds_dnscl != 0 ]] && [[ $NEMO_CHUNK_START_DATE == $run_start_date ]] ; then
  acc_cp tmp_ic.nc $nemo_data_1m_temperature_rest
fi

# activate correct Python environment
source activate /home/scrd102/cccma_conda/envs/py3_analysis_v2

# determine loop parameters/iteration info
python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py.py -y $(($NEMO_CHUNK_START_YEAR - 1)) -y $(($NEMO_CHUNK_END_YEAR + 1)) -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -t -1 -R ${run_start_year} -r ${run_start_month} -f ${months}
# Source the file with this info
. nemo_counter_info.cfg

# initial conditions
if [[ $ctds_dnscl != 0 ]] && [[ $NEMO_CHUNK_START_DATE == $run_start_date ]] ; then
  #python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py.py -y ${dada_ic_year} -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -t 0 -i $(( $dada_ic_month - 1 )) > ic_status
  python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py.py -y ${dada_ic_year} -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m tmp_ic.nc -t 0 -i $(( $dada_ic_month - 1 )) > ic_status
  if [ -z "$(ls ./data_1m_*_nomask.nc)" ] ; then
    echo "ERROR: No IC files generated!"
    exit 29
  fi
fi

# boundary conditions
if [[ $ctds_dnscl != 0 ]] && [[ $runmode == *"dada"* ]] ; then
  python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py.py -o obc_BDY_${dada_outfield}_yYYYY -y $(($NEMO_CHUNK_START_YEAR - 1)) -y $(($NEMO_CHUNK_END_YEAR + 1)) -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -t 1 > bc_status
  if [ -z "$(ls ./obc_*_${dada_outfield}_y*.nc)" ] ; then
    echo "ERROR: No OBC files generated!"
    exit 10
  fi
fi

# atmospheric forcing
if [[ $ctds_dnscl == 1 ]] && [[ $runmode == *"CanTODS"* ]] ; then
  # OMIP atmospheric forcing
  python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py.py -o VAR_${dada_outfield}_yYYYY -y $(($NEMO_CHUNK_START_YEAR - 1)) -y $(($NEMO_CHUNK_END_YEAR + 1)) -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -t 2 -f OMIP -A ${iaf_year_offset} -a ${iaf_loop_year} > frc_status
elif [[ $ctds_dnscl != 0 ]] ; then
  # CanESM atmospheric forcing
  if [[ -z "${iaf_year_offset}" ]] && [[ -z "${iaf_loop_year}" ]] ; then
      # use current year for forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py.py -o VAR_${dada_outfield}_yYYYY -y $(($NEMO_CHUNK_START_YEAR - 1)) -y $(($NEMO_CHUNK_END_YEAR + 1)) -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -M ${dada_forcing_freq} -t 2 > frc_status
  else
      # use cyclical forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py.py -o VAR_${dada_outfield}_yYYYY -y $(($NEMO_CHUNK_START_YEAR - 1)) -y $(($NEMO_CHUNK_END_YEAR + 1)) -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble}  -m domain_cfg.nc -M ${dada_forcing_freq} -t 2 -A ${iaf_year_offset} -a ${iaf_loop_year} > frc_status
  fi
  if [ -z "$(ls ./tair_${dada_outfield}_y*)" ] ; then
    echo "ERROR: Not all forcing files generated!"
    exit 92
  fi
fi

# river forcing
if [[ $runmode != *"dada"* ]] && [[ $nemo_ln_rnf == "on" ]] ; then
  if [[ -z "${nemo_river_remap}" ]] ; then
   # don't remap/rescale rivers
    if [[ -z "${iaf_year_offset}" ]] && [[ -z "${iaf_loop_year}" ]] ; then
      # use current year for forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py.py -o rvr_${dada_outfield} -y $NEMO_CHUNK_START_YEAR -y $NEMO_CHUNK_END_YEAR -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble} -m domain_cfg.nc -t 3 > rvr_status
    else
      # use cyclical forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py.py -o rvr_${dada_outfield} -y $NEMO_CHUNK_START_YEAR -y $NEMO_CHUNK_END_YEAR -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble} -m domain_cfg.nc -t 3 -A ${iaf_year_offset} -a ${iaf_loop_year} > rvr_status
    fi
  else
    # remap/rescale rivers
    if [[ -z "${iaf_year_offset}" ]] && [[ -z "${iaf_loop_year}" ]] ; then
      # use current year for forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py.py -o rvr_${dada_outfield} -y $NEMO_CHUNK_START_YEAR -y $NEMO_CHUNK_END_YEAR -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble} -m domain_cfg.nc -t 3 -v ${nemo_river_remap} > rvr_status
    else
      # use cyclical forcing
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py.py -o rvr_${dada_outfield} -y $NEMO_CHUNK_START_YEAR -y $NEMO_CHUNK_END_YEAR -P ${dada_parent_path} -p ${dada_parent_name} -x ${dada_parent_experiment} -e ${dada_parent_ensemble} -m domain_cfg.nc -t 3 -A ${iaf_year_offset} -a ${iaf_loop_year} -v ${nemo_river_remap} > rvr_status
    fi
  fi
  if [ -z "$(ls ./rvr_${dada_outfield}).nc" ] ; then
    echo "ERROR: River files not generated!"
    exit 28
  fi
fi

# delete domain file (will be re-copied in model portion of loop)
rm -f domain_cfg.nc
rm -f tmp_ic.nc
