#!/bin/bash
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

# sourcefile original db file: 
#/space/hall6/sitestore/eccc/crd/ccrn/forcing/data/nemo/nemo_inputs/initial_conditions/nemo4_orca1_ice_jstart/
# nemo4-eorca1-blk4-mesh_mask.nc.001 
set -e
gridtarget=urdy_eorca1_en4_v1.1.1995_2014.monthlymean_eorca1t_nemo_l75_jstart.nc.001
echo gridtarget=$gridtarget

sourcefile=nemo4-eorca1-blk4-mesh_mask.nc.001
echo sourcefile=$sourcefile

echo ncks -A -v tmask           ${sourcefile} -o closea_mask_tmp.nc    \;
echo ncks -A -v nav_lon,nav_lat ${sourcefile} -o closea_mask_tmp.nc    \;

echo ncrename -v tmask,closea_mask closea_mask_tmp.nc \;
echo ncrename -v nav_lev,depth     closea_mask_tmp.nc \;
echo ncrename -d nav_lev,depth     closea_mask_tmp.nc \;

echo ncatted -O -a standard_name,nav_lon,o,c,"longitude"           closea_mask_tmp.nc\;   #### need this attribute to work properly
echo ncatted -O -a long_name,nav_lon,o,c,"Longitude"               closea_mask_tmp.nc\;   #### need this attribute to work properly
echo ncatted -O -a _CoordinateAxisType,nav_lon,o,c,"Lon"           closea_mask_tmp.nc\;   #### need this attribute to work properly
echo ncatted -O -a standard_name,nav_lat,o,c,"latitude"            closea_mask_tmp.nc\;   #### need this attribute to work properly
echo ncatted -O -a long_name,nav_lat,o,c,"Latitude"                closea_mask_tmp.nc\;   #### need this attribute to work properly
echo ncatted -O -a _CoordinateAxisType,nav_lat,o,c,"lat"           closea_mask_tmp.nc\;   #### need this attribute to work properly
echo ncatted -O -a coordinates,closea_mask,o,c,\"nav_lon nav_lat\" closea_mask_tmp.nc\;   #### need this attribute to work properly

echo "ncap2 -F -O -s 'closea_mask(1,1:75,203:235,332:344)= 0.' closea_mask_tmp.nc closea_mask_tmp00.nc;" # -F to be able to use indices starting at 1 not 0
echo "ncap2 -F -O -s 'closea_mask=float(closea_mask)' closea_mask_tmp00.nc closea_mask_tmp01.nc;"

echo cdo -remapnn,${gridtarget} closea_mask_tmp01.nc  closea_mask_tmp.nc \;

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

echo ncrename -d time_counter,t closea_mask_tmp.nc                  \;
echo ncks -A -v closea_mask     closea_mask_tmp.nc -o ${targetfile} \;

echo rm closea_mask_tmp00.nc \;
echo rm closea_mask_tmp01.nc \;
echo rm closea_mask_tmp.nc   \;

echo cp ${targetfile} urdy_domcfg_urormod_closeamask.nc  \;
echo save urdy_domcfg_urormod_closeamask.nc urdy_domcfg_urormod_closeamask.nc    \;

set +e

# echo release urdy_domcfg_urormod_closeamask.nc urdy_domcfg_urormod_closeamask.nc \;
# echo fdb delete urdy_domcfg_urormod_closeamask.nc         \;
# echo rm -f ../../ urdy_domcfg_urormod_closeamask.nc       \;
