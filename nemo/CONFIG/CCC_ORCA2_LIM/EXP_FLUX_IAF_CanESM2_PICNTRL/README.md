CCC_ORCA2_LIM_CMOC/EXP_FLUX_IAF_CanESM2_PICNTRL
===============================================

`CCC_ORCA2_LIM_CMOC` is the CCCma configuration of ORCA2-CMOC.

The `EXP_FLUX_IAF_CanESM2_PICNTRL` experiment is a flux-forced run
with interannual forcing derived from the CanESM2 pre-industrial
control experiment (run IGA).

NOTE: 

There is restoring the SSS and SST. At the time of writing, the
SSS was from Levitus. In the future SST and SSS should come from 
the CanESM2 runs. 

Sometimes runs crash due to negative salinities (normally associated with
Arctic rivers). One quick fix is to modify `stpctl.F90` to include:

             WHERE ( tsn(:,:,1,jp_sal) < 0.0 ) tsn(:,:,1,jp_sal) = 0.1
             kindic = 0

near line 140.


CO2 is not dealt with specifically. It is set by default to 280 ppm in
namelist_pisces.

