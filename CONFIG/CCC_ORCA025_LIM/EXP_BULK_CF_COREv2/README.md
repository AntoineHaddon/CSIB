nga is the run forced by CanESM2 historical forcing and reduced CanESM2 rout around sub-polar Atlantic (as runoff), having tidal mixing, double diffusion and solar penetration.

Restarted from mc_nfo_00101231_restart.tar.001

nemo_rs nfo nga 00101231 y

=========================================
code: nemo3.4.1

executable: 
  nemo_3.4.1_orca025_lim_zdftmx_exe

cpp_ORCA025_LIM_zdftmx.fcm:
bld::tool::fppkeys key_trabbl  key_lim2 key_dynspg_flt key_diaeiv key_ldfslp key_traldf_c2d key_ traldf_eiv key_dynldf_c3d key_zdftke key_zdfddm key_zdftmx key_iomput key_mpp_mpi key_orca_r025 key_mpp_rep

Initial conditions: 
  interpolated from orca1 IC:
  nemo_3.4_orca025_potemp_1m_z46_nomask_orca1jan.nc
  nemo_3.4_orca025_salin_1m_z46_nomask_orca1jan.nc

forcing 
  nemo_3.4_orca1_tas_day_r1i1p1_1979-2005-mean.nc
  nemo_3.4_orca1_uas_day_r1i1p1_1979-2005-mean.nc
  nemo_3.4_orca1_vas_day_r1i1p1_1979-2005-mean.nc
  nemo_3.4_orca1_huss_day_r1i1p1_1979-2005-mean.nc
  nemo_3.4_orca1_rsds_day_r1i1p1_1979-2005-mean.nc
  nemo_3.4_orca1_rlds_day_r1i1p1_1979-2005-mean.nc
  nemo_3.4_orca1_pr_day_r1i1p1_1979-2005-mean.nc
  nemo_3.4_orca1_prsn_day_r1i1p1_1979-2005-mean.nc
  nemo_3.4_weights_canam4grid_bicubic_orca1.nc
  nemo_3.4_weights_canam4grid_bilinear_orca1.nc
  nemo_3.4.1_orca025_rout_socoefr_month_r1i1p1_1979-2005-mean2.nc

output file:
  nemo_3.4.1_orca025_lim_iodef_5d_xml

30/APR/2013. D.YANG
