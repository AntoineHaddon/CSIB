#!/bin/bash
#=========================================================
# This script contains three functions:
# * CMIP6 nemo diagnostics (offline)
# * Time mean (1d_diaptr -> 1m_diaptr)
# * Generation of time series
# (requires with_nemo_diag = on)
# 
# This script is peeled off from nemo_diag.tsk and adapted 
# to facilitate migrations (D. Yang, Sep 2020).
#=========================================================

####################
## Pre-processing ##
####################

if [ $with_nemo_diag -eq 1 ]; then

   # Access the nemo diag executable
   diag_exe=${nemo_diag_exe:=nemo_diag.exe}
   [ -z "$diag_exe" ] && ( echo "diag_exe is not defined." ; exit 1 )
   diag_exe_path=${storage_dir}/executables/${diag_exe}
   if [[ -e "$diag_exe_path" ]]; then 
       cp $diag_exe_path .
   else
       echo "$diag_exe does not exist in $(dirname $diag_exe_path)! Has it been compiled?"
       exit 1
   fi

   # Access file containing grid information
   [ -s orca_mesh_mask ] || access orca_mesh_mask $orca_grid_info

   # Access file containing mfo line mask
   [ -s mfo_line_mask ] || access mfo_line_mask mfo_line_mask

   # Number of months in a chunk
   nmon=`echo $chunk_start_month $chunk_end_month | awk '{printf "%2.2d", $2-$1+1}'`

   # suffix list for sub-yearly nemo historical files.
   sfxlst=("grid_T" "grid_U" "grid_V" "grid_W" "icemod")
   frqlst=("1m" "1m" "1m" "1m" "1m")
   if [ ${nemo_carbon} = 1 -a ${pisces_offline} = 0 ] ; then
     sfxlst=("${sfxlst[@]}" "btrc_T" "diad_T")
     frqlst=("${frqlst[@]}" "1m" "1m")
   fi
   if [ $output_level -ge 1 ] ; then
     sfxlst=("${sfxlst[@]}" "grid_T_ar6" "grid_U_ar6" "grid_V_ar6" "grid_W_ar6"  \
             "scalar_ar6"          \
             "grid_T_ar6" "grid_U_ar6" "grid_V_ar6" "icemod"  \
             "diaptr")
     frqlst=("${frqlst[@]}" "1m" "1m" "1m" "1m" "1m" "1d" "1d" "1d" "1d" "1d")
     if [ ${nmon} = 12 ] ; then
       sfxlst=("${sfxlst[@]}" "grid_T_ar6")
       frqlst=("${frqlst[@]}" "1y")
     fi
   fi

   # Check that the lists are the same length. Use this number to loop below..
   n_sfx=${#sfxlst[@]}
   n_frq=${#frqlst[@]}
   if [ "$n_sfx" -ne "$n_frq" ]; then
      bail "ERROR: sfxlst and frqlst MUST have the same number of elements."
   fi

   # Make links to the history files
   for i in $(seq 0 $(($n_sfx-1))); do 
     sfx=${sfxlst[$i]}
     frq=${frqlst[$i]}
     fhis=$OUTPUT_PATH/${RUNID}_${frq}_${chunk_start_date}_${chunk_end_date}_${sfx}.nc
     # Merge sub-yearly files 
     if [ $nmon -lt 12 ] ; then
       fyhis=$OUTPUT_PATH/${RUNID}_${frq}_${chunk_start_year}????_*_${sfx}.nc 
       [ -s $fyhis ] && cdo mergetime $fyhis ${frq}_${sfx}_${chunk_start_year}
     elif [ $nmon -eq 12 ] ; then
       [ -s $fhis ] && ln -s $fhis ${frq}_${sfx}_${chunk_start_year}
     else
       echo "Number of months in a chunk must be less or equal to 12"
     fi
   done

   # Make links to the nemo restart files
   yearm1=`echo $chunk_start_year | awk '{printf "%4.4d",$1-1}' -`
   if [ $yearm1 -ge 1 -a ${output_level} -ge 1 ] ; then
     # Make link to the previous year restart file
     frsp=$OUTPUT_PATH/${RUNID}_${yearm1}1231_restart.tar
     [ -s $frsp ] && ln -s $frsp rsp || ( echo "$frsp does not exist" ; exit 1 )
     # Make link to the current year restart file
     frs=$OUTPUT_PATH/${RUNID}_${chunk_start_year}1231_restart.tar
     [ -s $frs ] && ln -s $frs rsc || ( echo "$frs does not exist" ; exit 1 )

     # Get tn and sn from the last step of previous year
     if [ -L rsp ] ; then
       mkdir dir_rsp; cd dir_rsp
       tar -xvf ../rsp
       cdo select,name=tn,timestep=-1 *_restart.nc ../tnp.nc
       cdo select,name=sn,timestep=-1 *_restart.nc ../snp.nc
       cd ..
       rm -f -r rsp dir_rsp
     fi

     # Get tn and sn from the last step of current year
     if [ -L rsc ] ; then
       mkdir dir_rsc; cd dir_rsc
       tar -xvf ../rsc
       cdo select,name=tn,timestep=-1 ${runid}_*_restart.nc ../tnc.nc
       cdo select,name=sn,timestep=-1 ${runid}_*_restart.nc ../snc.nc
       cd ..
       rm -f -r rsc dir_rsc
     fi

     ##########################
     # CMIP6 nemo diagnostics #
     ##########################
     ln -s 1m_grid_T_${chunk_start_year} grid_t || bail "Link to grid_t failed"
     ln -s 1m_grid_U_${chunk_start_year} grid_u || bail "Link to grid_u failed"
     ln -s 1m_grid_V_${chunk_start_year} grid_v || bail "Link to grid_v failed"

     [ -L grid_t -a -s tnp.nc ] && $diag_exe || bail "grid_t or tnp.nc does not exist"
   
     # Change permisions
     [ -s mfo.nc ] && chmod u+w mfo.nc || bail "mfo.nc does not exist"
     [ -s msftbarot.nc ] && chmod u+w msftbarot.nc || bail "msftbarot.nc does not exist"
     [ -s tstend.nc ] && chmod u+w tstend.nc || bail "tstend.nc does not exist"

     # Append mfo.nc to 1m_scalar_ar6_${chunk_start_year}
     cp 1m_scalar_ar6_${chunk_start_year} 1m_scalar_ar6.nc && chmod u+w 1m_scalar_ar6.nc || bail "1m_scalar_ar6_${chunk_start_year} does not exist"
     ncks -A mfo.nc 1m_scalar_ar6.nc
     rm -f 1m_scalar_ar6_${chunk_start_year}
     mv 1m_scalar_ar6.nc 1m_scalar_ar6_${chunk_start_year}

     # Append msftbarot.nc to 1m_grid_u_ar6_${chunk_start_year}
     cp 1m_grid_U_ar6_${chunk_start_year} 1m_grid_u_ar6.nc && chmod u+w 1m_grid_u_ar6.nc || bail "1m_grid_U_ar6_${chunk_start_year} does not exist"
     ncks -A msftbarot.nc 1m_grid_u_ar6.nc
     rm -f 1m_grid_U_ar6_${chunk_start_year}
     mv 1m_grid_u_ar6.nc 1m_grid_U_ar6_${chunk_start_year}

     # Append tstend.nc to 1y_grid_t_ar6_${chunk_start_year}
     if [ $nmon -eq 12 ] ; then
       cp 1y_grid_T_ar6_${chunk_start_year} 1y_grid_t_ar6.nc && chmod u+w 1y_grid_t_ar6.nc || bail "1y_grid_T_ar6_${chunk_start_year} does not exist"
       ncks -A tstend.nc 1y_grid_t_ar6.nc
       rm -f 1y_grid_T_ar6_${chunk_start_year}
       mv 1y_grid_t_ar6.nc 1y_grid_T_ar6_${chunk_start_year}
     fi
   fi # end of $yearm1 -ge 1 -a $output_level -ge 1

   #####################################
   # Time mean (1d_diaptr -> 1m_diaptr)#
   #####################################
   #
   [ -L 1d_diaptr_${chunk_start_year} -o -s 1d_diaptr_${chunk_start_year} ] && cdo -b F64 monmean 1d_diaptr_${chunk_start_year} 1m_diaptr_${chunk_start_year}
 
   #########################################
   # Split historical files to time series #
   #########################################

   # Split to time series
   for i in $(seq 0 $(($n_sfx-1))); do
     sfx=${sfxlst[$i]}
     frq=${frqlst[$i]}
     # Replace 1d_diaptr with 1m_diaptr after doing time mean
     if [ $sfx == "diaptr" ] ; then
       frq=1m
     fi
     fhis=${frq}_${sfx}_${chunk_start_year}
     cdo splitname $fhis xxx-${frq}_${sfx}_
   done

   # Save time series
   tslist=`ls -1 xxx-*`
   RUNPATH=$OUTPUT_PATH
   for ts in $tslist ; do
     tssfx=`echo $ts |cut -f 2 -d '-' |sed 's/_ar6//' | tr '[:upper:]' '[:lower:]'`
     save $ts sc_${runid}_${chunk_start_year}_${chunk_end_year}_${tssfx}
     release $ts
   done
fi # with_nemo_diag
