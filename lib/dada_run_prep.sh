#!/bin/bash

set -e

# get copy of the domain.cfg file
is_defined $nemo_coordinates || bail "The variable nemo_coordinates must be defined in the configuration file."
access domain_cfg.nc $nemo_coordinates

# copy remapping JSON file
cp ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.json remap_canesm.json

# overwrite default JSON file values
sed -i "/\"parent_path\":/c\ \"parent_path\": \"${dada_parent_path}\"," remap_canesm.json
sed -i "/\"parent_name\":/c\ \"parent_name\": \"${dada_parent_name}\"," remap_canesm.json
sed -i "/\"parent_ensemble\":/c\ \"parent_ensemble\": \"${dada_parent_ensemble}\"," remap_canesm.json
sed -i "/\"parent_experiment\":/c\ \"parent_experiment\": \"${dada_parent_experiment}\"," remap_canesm.json
sed -i "/\"meshfile\":/c\ \"meshfile\": \"domain_cfg.nc\"," remap_canesm.json
sed -i "/\"ic_ind\":/c\ \"ic_ind\": $(( $dada_ic_month - 1 ))," remap_canesm.json
sed -i "/\"run_start_year\":/c\ \"run_start_year\": ${run_start_year}," remap_canesm.json
sed -i "/\"run_start_month\":/c\ \"run_start_month\": ${run_start_month}," remap_canesm.json
sed -i "/\"nemo_freq_months\":/c\ \"nemo_freq_months\": ${months}," remap_canesm.json
if [[ ! -z "${dada_year_offset}" ]] && [[ ! -z "${dada_loop_year}" ]] ; then
  sed -i "/\"iaf_year_offset\":/c\ \"iaf_year_offset\": \"${dada_year_offset}\"," remap_canesm.json
  sed -i "/\"iaf_loop_year\":/c\ \"iaf_loop_year\": \"${dada_loop_year}\"," remap_canesm.json
fi
sed -i "/\"mor\":/c\ \"mor\": \"${dada_forcing_freq}\"," remap_canesm.json

# activate correct Python environment
source activate /home/scrd102/cccma_conda/envs/py3_analysis_v2

# determine loop parameters/iteration info
python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -t -1 
# Source the file with this info
. nemo_counter_info.cfg

# overwrite model years
sed -i "/\"years\":/c\ \"years\": [$(($NEMO_CHUNK_START_YEAR - 1)), $(($NEMO_CHUNK_START_YEAR +1))]," remap_canesm.json

# initial conditions
if [[ $ctds_dnscl != 0 ]] && [[ $NEMO_CHUNK_START_DATE == $run_start_date ]] ; then
  if [[ $runmode == *"CanTODS"* ]]; then
    if [[ $nemo_from_rest == "on" ]] ; then
      # Start from rest
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -y ${dada_ic_year} -t 0 > ic_status
    else
      # hot start; interpolate restart files from parent run
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -y ${dada_ic_year} -P ${dada_parent_path%%nc_output*} -p ${dada_parent_name/CanESM5-NEMO4-/} -t 0 -H 1 -V 1 > ic_status
    fi
  else
    # eORCA grids slightly different, so interpolate to an existing IC file
    access tmp_ic.nc $nemo_data_1m_temperature_rest
    if [[ $nemo_from_rest == "on" ]] ; then
      # Start from rest
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -y ${dada_ic_year} -m tmp_ic.nc -t 0 > ic_status
    else
      # hot start; interpolate restart files
      python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -y ${dada_ic_year} -P ${dada_parent_path%%nc_output*} -m tmp_ic.nc -t 0 -H 1 > ic_status
    fi
  fi
  if [[ -z "$(ls ./data_1m_*_nomask.nc)" && -z "$(ls ./restart*.nc)" ]] ; then
    echo "ERROR: No IC files generated!"
    exit 29
  fi
fi

# boundary conditions
if [[ $ctds_dnscl != 0 ]] && [[ $runmode == *"CanTODS"* ]] ; then
  python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o obc_BDY_${dada_outfield}_yYYYY.nc -t 1 > bc_status
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
  python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o VAR_${dada_outfield}_yYYYY -t 2 > frc_status
  if [ -z "$(ls ./tair_${dada_outfield}_y*)" ] ; then
    echo "ERROR: Not all forcing files generated!"
    exit 92
  fi
fi

# sea-surface salinity
if [[ $nemo_nn_sssr != 0 ]] ; then
  python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o sss_${dada_outfield}_yYYYY -t 4 > sos_status
  if [ -z "$(ls ./sss_${dada_outfield}_y*)" ] ; then
    echo "ERROR: SSS file not generated!"
    exit 2007
  fi
fi

# river forcing
# TODO: rivers will need to be remapped properly
if [[ $nemo_ln_rnf == "on" ]] && [[ $runmode != *"CanTODS"* ]]; then
  if [[ -z "${nemo_river_remap}" ]] ; then
    # don't remap/rescale rivers; simply copy
    python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o rvr_${dada_outfield} -y $NEMO_CHUNK_START_YEAR -y $NEMO_CHUNK_END_YEAR -t 3 > rvr_status
  else
    # remap/rescale rivers; need to get files for that
    access nemo_river_remap_${NEMO_CHUNK_START_YEAR}-${NEMO_CHUNK_END_YEAR}.nc $nemo_river_remap
    access dpg.nc $dada_parent_grid || ln -s $dada_parent_grid dpg.nc # get parent grid for remapping
    # rescale
    python3 ${CANESM_SRC_ROOT}/CanNEMO/lib/remap_canesm.py -o rvr_${dada_outfield} -y $NEMO_CHUNK_START_YEAR -y $NEMO_CHUNK_END_YEAR -t 3 -v nemo_river_remap_${NEMO_CHUNK_START_YEAR}-${NEMO_CHUNK_END_YEAR}.nc -g dpg.nc > rvr_status
    # remove temporary files
    rm -f nemo_river_remap_${NEMO_CHUNK_START_YEAR}-${NEMO_CHUNK_END_YEAR}.nc
    rm -f dpg.nc
  fi
  if [ -z "$(ls ./rvr_${dada_outfield}).nc" ] ; then
    echo "ERROR: River files not generated!"
    exit 28
  fi
fi

# delete domain file (will be re-copied in model portion of loop)
rm -f domain_cfg.nc
rm -f tmp_ic.nc
