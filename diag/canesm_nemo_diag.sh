#########################################################
# nemo diagnostics & time mean (1d -> 1m)
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
  if [[ $year == $run_start_time ]] && [[ $nemo_from_rest == 'on' ]]; then
    # use the current year because output.init.nc is used in that case (below)
    yearm1=`echo $year | awk '{printf "%04d", $1}'`
    file_state="initial"
    mon1=`echo $fmon $months | awk '{printf "%02d", $1 - 1 + $2 }'` 
    t_state="votemper"
    s_state="vosaline"
  else
    yearm1=`echo $year | awk '{printf "%04d", $1 - 1}'`
    mon1=$lmon
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
  [ -s orca_mesh_mask ] || access orca_mesh_mask $orca_grid_info

# Access file containing mfo line mask
  [ -s mfo_line_mask ] || access mfo_line_mask mfo_line_mask || bail "Failed to access $mfo_line_mask"

##########################################################################
# 1. Access input files/variables                                        #
# 2. Run the Fortran executable to compute the CMIP6 offline diagnostics #
##########################################################################
#
# Execute the following lines when output_level -ge 1;
if [[ $output_level -ge 1 ]]; then
   if [[ $fmon -eq 1 ]] ; then
   # access input variables for computing vars with priority level 1
     access grid_t mc_${runid}_${year}_m${fmon}_1m_grid_t.nc || bail "Failed to access $mc_${runid}_${year}_m${fmon}_1m_grid_t.nc"
     access grid_u mc_${runid}_${year}_m${fmon}_1m_grid_u.nc || bail "Failed to access $mc_${runid}_${year}_m${fmon}_1m_grid_t.nc"
     access grid_v mc_${runid}_${year}_m${fmon}_1m_grid_v.nc || bail "Failed to access $mc_${runid}_${year}_m${fmon}_1m_grid_t.nc"

     ################################################################
     # Run the CMIP6 nemo offline diagnostics executable: $diag_exe #
     ################################################################
     # make sure inputs exist and run!
     $diag_exe

     release grid_t
     release grid_u
     release grid_v
     # Append mfo.nc to 1m_scalar_ar6_${fmon} if existing
     if [ -s mfo.nc ]; then
       chmod u+w mfo.nc
       access 1m_scalar_ar6.nc mc_${runid}_${year}_m${fmon}_1m_scalar_ar6.nc || bail "Failed to access $mc_${runid}_${year}_m${fmon}_1m_scalar_ar6.nc"
       chmod u+w $(readlink -f 1m_scalar_ar6.nc)
       ncks -A mfo.nc $(readlink -f 1m_scalar_ar6.nc)
       chmod u-w $(readlink -f 1m_scalar_ar6.nc)
       release 1m_scalar_ar6.nc
     else
       bail "mfo.nc does not exist"
     fi
     # Append msftbarot.nc to 1m_grid_u_ar6_${fmon} if existing
     if [ -s msftbarot.nc ]; then
       chmod u+w msftbarot.nc
       access 1m_grid_u_ar6.nc mc_${runid}_${year}_m${fmon}_1m_grid_u_ar6.nc || bail "Failed to access $mc_${runid}_${year}_m${fmon}_1m_grid_u_ar6.nc"
       chmod u+w $(readlink -f 1m_grid_u_ar6.nc)
       ncks -A msftbarot.nc $(readlink -f 1m_grid_u_ar6.nc)
       chmod u-w $(readlink -f 1m_grid_u_ar6.nc)
       release 1m_grid_u_ar6.nc
     else
       bail "msftbarot.nc does not exist"
     fi
   fi
fi
# Execute the following lines when output_level -ge 3;
if [[ $output_level -ge 3 ]]; then
   if [ $fmon -eq 1 ] ; then
     access grid_t mc_${runid}_${year}_m${fmon}_1m_grid_t.nc || bail "Failed to access $mc_${runid}_${year}_m${fmon}_1m_grid_t.nc"
     access grid_u mc_${runid}_${year}_m${fmon}_1m_grid_u.nc || bail "Failed to access $mc_${runid}_${year}_m${fmon}_1m_grid_t.nc"
     access grid_v mc_${runid}_${year}_m${fmon}_1m_grid_v.nc || bail "Failed to access $mc_${runid}_${year}_m${fmon}_1m_grid_t.nc"
     # Run offline computation of tendency terms only if starting from January and yearly chunk
     # Access the nemo restart files
     diag_rs1="mc_${runid}_${yearm1}_m${mon1}_nemors" # previous year
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
     $diag_exe

     release grid_t
     release grid_u
     release grid_v
     release rsp
     release rsc
     # Append tstend.nc to 1y_grid_t_ar6_${fmon} if existing
     if [ -s tstend.nc ]; then
       access 1y_grid_t_ar6.nc mc_${runid}_${year}_m${fmon}_1y_grid_t_ar6.nc na
       if [ -s 1y_grid_t_ar6.nc ]; then
         chmod u+w $(readlink -f 1y_grid_t_ar6.nc)
         ncks -A tstend.nc $(readlink -f 1y_grid_t_ar6.nc)
         chmod u-w $(readlink -f 1y_grid_t_ar6.nc)
         release 1y_grid_t_ar6.nc
       else # 1y_grid_t_ar6 does not exist if jobs shorter than 12mo
         save tstend.nc mc_${runid}_${year}_m${fmon}_1y_grid_t_ar6.nc 
       fi
     else
       bail "tstend.nc does not exist"
     fi
   fi
fi
release orca_mesh_mask 
release mfo_line_mask 
exit 0
