#!/bin/bash
# O Riche Aug 22nd 2022
# variant of the previous code
# Add a mask for eorca1 grid
# it's a square of 0s everywhere
# else is filled with 1s.
# lines# 
# # Caspian Sea  (spread over the globe)
# ncsnr[1]   = 1    ; ncstt[1]   = 0           
# ncsi1[1]   = 332  ; ncsj1[1]   = 243
# ncsi2[1]   = 344  ; ncsj2[1]   = 275
# These are python numbering so add 1
# to all integers for Fortran numbering
#
# O Riche Aug 15th 2022
# Additionally missed this reference in closea.F90/dom_clo s/routine:
#       !!                definitions is available at TOOLS/DOMAINcfg/make_closea_masks.py
# which refers to TOOLS folder in CanESM5 code/repo
# O Riche Aug 12th 2022
# 1. copy tmask from source file
#    to a temporary file
#    
# 2. set closea_mask to 0 in the 
#    Caspian Sea using CanESM5 code
#    in ./nemo/NEMO/OPA_SRC/DOM/closea.F90
#     IF( cp_cfg == "orca" ) THEN
#        !
#        SELECT CASE ( jp_cfg )
#        !                                           ! =======================
#        CASE ( 1 )                                  ! ORCA_R1 configuration
#           !                                        ! =======================
#           ncsnr(1)   = 1    ; ncstt(1)   = 0           ! Caspian Sea
#           ncsi1(1)   = 332  ; ncsj1(1)   = 203
#           ncsi2(1)   = 344  ; ncsj2(1)   = 235
#
# =============================================================== 
#    Below is wrong but possibly relevant to the 
#    Caspian Sea although not to its global grid
#    definition:
#    in ./nemo/NEMO/OPA_SRC/DOM/dommsk.F90
#    defined in ORCA1.
#    in the square defined by:
#    ij0 =  87; ij1 =  88
#    ii0 = 160; ii1 = 161
#    tmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) = 0._wp
#    mi0/mi1/mj0/mj1 are used to convert global indices
#    to local indices
#    Note from NEMO_manual.pdf (2019 Ed.):
#    The 1-d arrays mig(1 : jpi) and mjg(1 : jpj), defined in
#    dom_glo routine ( domain.F90 module), should be
#    used to get global domain indices from local domain indices.
#    The 1-d arrays, mi0(1 : jpiglo), mi1(1 : jpiglo) and
#    mj0(1 : jpjglo), mj1(1 : jpjglo) have the reverse purpose and
#    should be used to define loop indices expressed in
#    global domain indices (see examples in dtastd.F90 module).
#    NEMO Reference Manual Page 96 of 273
#    Chap. 7 Lateral Boundary Condition (LBC)
# =============================================================== 
# 3. regrid over eORCA1
#
# 4. add it to targe file as closea_mask
#
#  Alternatively find the positions and the lat/lon
#  with ncview or panoply or python and directly set to 0s
#  tmask/closea_mask at these locations

set -e

# target file is on eorca1 grid
# targetfile original db file:
# /space/hall6/sitestore/eccc/crd/ccrn/forcing/data/nemo/nemo_inputs/initial_conditions/
# nemo4_orca1_ice_jstart/urdy_eorca1_domcfg_eorca1v2.2x_orca_jstart.nc.001
dblocation=/space/hall6/sitestore/eccc/crd/ccrn/forcing/data/nemo/nemo_inputs/initial_conditions/nemo4_orca1_ice_jstart
targetfile=urdy_eorca1_domcfg_eorca1v2.2x_orca_jstart.nc.001
echo rm ${targetfile}                             \;
echo rm ${targetfile:0:-4}                        \;  #### read by ncview only if there is no .001 at the end
echo cp ${dblocation}/${targetfile} ${targetfile} \;
echo chmod u+wx ${targetfile}                     \;
echo ln -s ${targetfile} ${targetfile:0:-4}       \;
echo targetfile=${targetfile}                     \;

# echo "ncap2 -F -O -s 'closea_bgc_mask=array(1.,0,gdept_0) ' " ${targetfile} " closea_mask_tmp.nc;"           # 0s in the Caspian Sea box 1s everywhere else.
# echo "ncap2 -F -O -s 'closea_bgc_mask(1,1:75,244:276,332:344)= 0.' closea_mask_tmp.nc closea_mask_tmp00.nc;" # -F to be able to use indices starting at 1 not 0
echo "ncap2 -F -O -s 'closea_bgc_mask=array(1.,0,bathy_metry) ' " ${targetfile} " closea_mask_tmp.nc;"           # 0s in the Caspian Sea box 1s everywhere else.
echo "ncap2 -O -s 'closea_bgc_mask(0,243:275,332:344)= 0.' closea_mask_tmp.nc closea_mask_tmp00.nc;" # -F to be able to use indices starting at 1 not 0
echo "ncap2 -F -O -s 'closea_bgc_mask=float(closea_bgc_mask)' closea_mask_tmp00.nc closea_mask_tmp01.nc;"

echo mv closea_mask_tmp01.nc ${targetfile} \;

echo rm closea_mask_tmp00.nc \;
# echo rm closea_mask_tmp01.nc \;
echo rm closea_mask_tmp.nc   \;

echo cp ${targetfile} urdy_domcfg_urormod_closeamask.nc  \;
echo save urdy_domcfg_urormod_closeamask.nc urdy_domcfg_urormod_closeamask.nc    \;

set +e

# echo release urdy_domcfg_urormod_closeamask.nc urdy_domcfg_urormod_closeamask.nc \;
# echo fdb delete urdy_domcfg_urormod_closeamask.nc         \;
# echo rm -f ../../ urdy_domcfg_urormod_closeamask.nc       \;