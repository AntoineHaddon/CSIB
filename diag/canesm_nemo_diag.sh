#########################################################
# nemo diagnostics & time mean (1d -> 1m) & time series
# D. Yang, Nov 2018, A. Shao, S.Kharin
#
# This script is sourced in CanESM/CCCma_tools/cccjob_dir/
# lib/jobdefs/canesm_nemo_diag_jobdef
#########################################################

set -x

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

# Access the nemo diag executable
  diag_exe=${nemo_diag_exe:=nemo_diag.exe}
  [ -z "$diag_exe" ] && ( echo "diag_exe is not defined." ; exit 1 )
  [ -s "$diag_exe" ] || access $diag_exe $nemo_diag_exe

# Access file containing grid information
  orca_grid_info=nemo_mesh_mask_rc3.nc
  [ -s orca_mesh_mask ] || access orca_mesh_mask $orca_grid_info

# Access file containing mfo line mask
  [ -s mfo_line_mask ] || access mfo_line_mask mfo_line_mask

# suffix list for sub-yearly nemo historical files.
  sfxlst="1m_grid_t 1m_grid_u 1m_grid_v 1m_grid_w 1m_icemod 1m_ptrc_t 1m_diad_t"
  if [ $output_level -ge 1 ] ; then
      sfxlst="$sfxlst 1m_grid_t_ar6 1m_grid_u_ar6 1m_grid_v_ar6 1m_grid_w_ar6     \
              1m_scalar_ar6          \
              1d_grid_t_ar6 1d_grid_u_ar6 1d_grid_v_ar6 1d_icemod  \
              3h_grid_t_ar6 1d_diaptr"
  fi

# suffix list for yearly nemo historical files.
  if [ $output_level -ge 1 ] ; then
      sfxlst_1y="1y_grid_t_ar6"
  fi

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
      rm -f ${sfx}_??
      mv ${sfx}_m$fmon ${sfx}_$fmon
    fi
  done

# Execute the following lines when output_level -ge 1
  if [ $output_level -ge 1 ] ; then
      if [ $nmon -eq 1 ] ; then
        for sfx in $sfxlst_1y ; do
          diag_hist="mc_${runid}_${fyear}_m${fmon}_${sfx}.nc"
          access ${sfx}_${fmon} $diag_hist na
        done
      fi

# Access the nemo restart files
      diag_rs1="mc_${runid}_${yearm1}_m${lmon}_nemors.tar" # previous year
      diag_rs2="mc_${runid}_${year}_m${lmon}_nemors.tar"   # current year
      access rsp $diag_rs1 || ( echo "$diag_rs1 does not exist" ; exit 1 )
      access rsc $diag_rs2 || ( echo "$diag_rs2 does not exist" ; exit 1 )

# Get tn and sn from the last step of previous year
      if [ -L rsp ] ; then
        mkdir dir_rsp; cd dir_rsp
        tar -xvf ../rsp
        cdo select,name=tn,timestep=-1 *_restart.nc ../tnp.nc
        cdo select,name=sn,timestep=-1 *_restart.nc ../snp.nc
        cd ..
        release rsp
        rm -f -r dir_rsp
      fi

# Get tn and sn from the last step of current year
      if [ -L rsc ] ; then
        mkdir dir_rsc; cd dir_rsc
        tar -xvf ../rsc
        cdo select,name=tn,timestep=-1 ${runid}_*_restart.nc ../tnc.nc
        cdo select,name=sn,timestep=-1 ${runid}_*_restart.nc ../snc.nc
        cd ..
        release rsc
        rm -f -r dir_rsc
      fi

##########################
# CMIP6 nemo diagnostics #
##########################
      ln -s 1m_grid_t_${fmon} grid_t  || bail "Link to grid_t failed"
      ln -s 1m_grid_u_${fmon} grid_u  || bail "Link to grid_u failed"
      ln -s 1m_grid_v_${fmon} grid_v  || bail "Link to grid_v failed"

      # Compile the diagnostic program
      WRKDIR=$PWD
      ( cd $CCRNSRC/CanESM/CanNEMO/diag ;
      ifort -o $WRKDIR/nemo_diag.exe nemo_diag_glovars.F90 nemo_diag_cal.F90 nemo_diag.F90 uvic_netcdf.f \
               `nc-config --fflags` `nc-config --flibs`
    )

      [ -L grid_t -a -s tnp.nc ] && $diag_exe || bail "grid_t or tnp.nc does not exist"

######################################
# Time mean (1d_diaptr -> 1m_diaptr) #
######################################
      [ -L 1d_diaptr_${fmon} -o -s 1d_diaptr_${fmon} ] && cdo -b F64 monmean 1d_diaptr_${fmon} 1m_diaptr_${fmon}

# Remove attributes FillValue & missing_value
      [ -s mfo.nc ] && chmod u+w mfo.nc || bail "mfo.nc does not exist"
      ncatted -h -a FillValue,time_counter,d,,, mfo.nc
      ncatted -h -a missing_value,time_counter,d,,, mfo.nc
      ncatted -h -a FillValue,time_counter_bnds,d,,, mfo.nc
      ncatted -h -a missing_value,time_counter_bnds,d,,, mfo.nc
      ncatted -h -a missing_value,mfo,d,,, mfo.nc

      [ -s msftbarot.nc ] && chmod u+w msftbarot.nc || bail "msftbarot.nc does not exist"
      ncatted -h -a FillValue,time_counter,d,,, msftbarot.nc
      ncatted -h -a missing_value,time_counter,d,,, msftbarot.nc
      ncatted -h -a FillValue,time_counter_bnds,d,,, msftbarot.nc
      ncatted -h -a missing_value,time_counter_bnds,d,,, msftbarot.nc
      ncatted -h -a missing_value,msftbarot,d,,, msftbarot.nc

      [ -s tstend.nc ] && chmod u+w tstend.nc || bail "tstend.nc does not exist"
      ncatted -h -a FillValue,time_counter,d,,, tstend.nc
      ncatted -h -a missing_value,time_counter,d,,, tstend.nc
      ncatted -h -a FillValue,time_counter_bnds,d,,, tstend.nc
      ncatted -h -a missing_value,time_counter_bnds,d,,, tstend.nc
      
      ncatted -h -a missing_value,opottemptend,d,,, tstend.nc
      ncatted -h -a missing_value,osalttend,d,,, tstend.nc

# Append mfo.nc to 1m_scalar_ar6_${fmon}
      cp 1m_scalar_ar6_${fmon} 1m_scalar_ar6.nc && chmod u+w 1m_scalar_ar6.nc || bail "1m_scalar_ar6_${fmon} does not exist"
      ncks -A mfo.nc 1m_scalar_ar6.nc
      rm -f 1m_scalar_ar6_${fmon}
      mv 1m_scalar_ar6.nc 1m_scalar_ar6_${fmon}

# Append msftbarot.nc to 1m_grid_u_ar6_${fmon}
      cp 1m_grid_u_ar6_${fmon} 1m_grid_u_ar6.nc && chmod u+w 1m_grid_u_ar6.nc || bail "1m_grid_u_ar6_${fmon} does not exist"
      ncks -A msftbarot.nc 1m_grid_u_ar6.nc
      rm -f 1m_grid_u_ar6_${fmon}
      mv 1m_grid_u_ar6.nc 1m_grid_u_ar6_${fmon}

# Append tstend.nc to 1y_grid_t_ar6_${fmon}
      if [ ${nmon} -eq 1 ] ; then
        cp 1y_grid_t_ar6_${fmon} 1y_grid_t_ar6.nc && chmod u+w 1y_grid_t_ar6.nc || bail "1y_grid_t_ar6_${fmon} does not exist"
        ncks -A tstend.nc 1y_grid_t_ar6.nc
        rm -f 1y_grid_t_ar6_${fmon}
        mv 1y_grid_t_ar6.nc 1y_grid_t_ar6_${fmon}
      fi

#########################################
# Split historical files to time series #
#########################################

# Replace 1d_diaptr with 1m_diaptr after doing time mean
      sfxlst=`echo $sfxlst | sed -e "s/1d_diaptr/1m_diaptr/"`

# Append yearly diagnostics
      if [ $nmon -eq 1 ] ; then
        sfxlst="$sfxlst $sfxlst_1y"
      fi
  fi # end of "output_level -ge 1"

# Split to time series
  for sfx in $sfxlst ; do
    cdo splitname ${sfx}_${fmon} xxx-${sfx}_
  done

# Save time series
  tslist=`ls -1 xxx-*`
  for ts in $tslist ; do
    tssfx=`echo $ts |cut -f 2 -d '-' |sed 's/_ar6//'`
    save $ts sc_${runid}_${fyear}${fmon}_${lyear}${lmon}_${tssfx}
    release $ts
  done
