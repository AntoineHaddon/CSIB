CCC_ORCA2_LIM_PISCES/EXP_FLUX_IAF_CanESM2_HIST
==============================================

`CCC_ORCA2_LIM_PISCES` is the CCCma configuration of ORCA2-PISCES

The `EXP_FLUX_IAF_CanESM2_HIST` experiment is a flux-forced run
with interannual forcing derived from the CanESM2 historical experiment 
(run IGM).

NOTE: 

There is restoring the SSS and SST. At the time of writing, the
SSS was from Levitus. In the future SST and SSS should come 
from the CanESM2 runs. 

CO2 IS NOT DEALT WITH EXPLICITLY at this time, but is set to the default 280 ppm
specified in namelist_pisces. A txt file of valid co2 concentrations should be 
provided for PISCES in the future.

Sometimes runs crash due to negative salinities (normally associated with
Arctic rivers). One quick fix is to modify `stpctl.F90` to include:

             WHERE ( tsn(:,:,1,jp_sal) < 0.0 ) tsn(:,:,1,jp_sal) = 0.1
             kindic = 0

near line 140.

