# DMS model

Notes on implementation in NEMO4.2 (CanNEMo version) of DMS model (in ocean and in sea ice)

Ocean DMS is implemented as a TOP module, with main code in folder TOP/DMS.

Sea ice DMS is added as a sea ice tracer like the other variables of CSIB. See CSIB.md for calling sequence and files modified.


## Calling sequence of ocean DMS

Init in `trc_ini`
- `trc_ini_sms` calls `trc_ini_dms` which allocates arrays, calls `trc_nam_dms` to read ocean DMS namelist in CanOe namelist files
- `trc_ini_state` calls `trc_rst_read` reading of restart files

TOP main routine `trcstp.F90/trc_stp` calls:
- `trc_wri` output, calls `trc_wri_dms`: ocean DMS outputs
- `trc_sms` calls `trc_sms_dms`: compute DMS variables tendencies 

Transport done by TOP (like CanOE variables) `trc_trp`
- Correction of negative concentrations `trc_rad`




## Files modidfied for ocean DMS

(In addition to files in folder `TOP/DMS`)

TOP files modified 
- `par_trc.F90` declares DMS flag and number of tracers 
- `trcini.F90` stores DMS tracers info into arrays and calls `trc_ini_dms`
- `trcnam.F90` reads DMS tracers info and flag from TOP namelist, sets total number of tracers to account for dms tracers
- `trcrst.F90` add a check for DMS variables in restart file, in case of restart from a simulation without DMS
- `trcsbc` check for each variable that it is in restart (was only done for first)
- `trcsms.F90` calls `trc_sms_dms`
- `trcwri.F90` calls `trc_wri_dms`
- `TRP/trcrad.F90` correction of negative concentrations 

CanOE files modified
- `/CANBGC/CANOE/sms_canoe.F90` declaring canoe variables (that were local) to make them available for DMS computations
- `/CANBGC/CANOE/canoeprod.F90` change local variable to be global: `zprocn` and `zprocd` 
- `/CANBGC/CANOE/canoenzd.F90` change local varaible to be global: `grazing1`, `grazing2`, `grazing3`, `zmortpn`, `zmortpd`




## Changes from Tessa's version
- par_trc: moved DMS flag and number of tracers here
- trcrst: csib restart moved to trcini_csib
    + automatic spin up of ice BGC if no ice BGC variable in ice restart file
- trc_ini : csib ini moved to SI3
- trcice: not needed
- DMS/trcnam_dms: DMS parameters in CanOE namelist
- ice DMS 
    - parameters in CSIB namelist in sea ice model namelist
    - didn't include DMSP(D) molecular diffusion as it wasn't in previous model version
    - moved checks for consistency in flags to trcsms_csib because CSIB init is in sea ice model init which is before CanOE init
    - jp_csib (mumber of sea ice tracers) set in csib init
    - corrected output of release rates (descriptions in xml file where also incorrect) and remove output of 'npp' (already there with icenpp)
    - exudation f_ei not f_e1 (but they have the same value)
    - for exudation it is indeed NPP and not specific NPP. actually you need (specific NPP * ice algal C biomass * dmsp-to-C) which is simply (NPP * dmsp-to-C)
- added ice dms processes output (only for nonlinear processes)





## Parameters

In `namelist_canoe_ref`:

    !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    &namdmsoce     !   Ocean DMS module
    !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
    q_p1    = 12.           ! intracellular DMSPp-to-Carbon ratio for p1 (umol S:mmol C)
    q_p2    = 4.            ! intracellular DMSPp-to-Carbon ratio for p2 (umol S:mmol C)
    f_z1    = 0.3           ! sloppy feeding fraction for z1 (-)
    f_z2    = 0.3           ! sloppy feeding fraction for z2 (-)
    f_e1    = 0.05          ! exudation fraction for p1 (-)
    f_e2    = 0.05          ! exudation fraction for p2 (-)
    f_yield = 0.2           ! DMS yield (-)
    k_dmspd = 5.            ! bacterial DMSPd consumption rate constant (d-1)
    k_dms   = 0.5           ! bacterial DMS consumption rate constant (d-1)
    k_free  = 0.02          ! free DMSPd-lyase rate constant (d-1)
    k_photo = 0.1           ! photolysis rate constant (d-1)
    /

In `namelist_ice_ref`:

    !-----------------------------------------------------------------------    
    &namicetra    !   Sea ice tracers (CSIB2)
    !-----------------------------------------------------------------------    
    ! Ice DMS
        ln_dmsice       = .false.        ! Run with ice dms
        q_pi            = 4.             ! intracellular DMSPp-to-Carbon ratio (umol S:mmol C)
        f_zi            = 0.3            ! sloppy feeding fraction (-)
        f_ei            = 0.05           ! exudation fraction (-)
        f_yieldi        = 0.2            ! bacterial conversion fraction (-)
        k_dmspdi        = 1.             ! bacterial dmspd consumption rate constant (d-1)
        k_dmsi          = 0.2            ! bacterial dms consumption rate constat (d-1)
        k_freei         = 0.02           ! free lyase rate constant (d-1)
        k_photoi        = 0.1            ! photolysis rate constant (d-1)
        h_ni            = 1.0            ! half-saturation constant used in the Monod equation
    /





## Output

For file `field_def_nemo-canoe.xml`

        <!-- Ocean DMS variables and diagnostics -->
       <field id="dmspd"        long_name="Dissolved dimethylsulfoniopropionate"  unit="umolS/m3"             grid_ref="grid_T_3D" />     
       <field id="dms"          long_name="Dimethylsulfide "                      unit="umolS/m3"             grid_ref="grid_T_3D" /> 
       <field id="sfcdmspd"    long_name="Suface ocean DMSPD"                     unit="umolS/m3"             grid_ref="grid_T_2D" />
       <field id="sfcdms"      long_name="Suface ocean DMS"                       unit="umolS/m3"             grid_ref="grid_T_2D" />
       
       <field id="dmsflux"     long_name="Sea-to-air DMS flux"                    unit="umolS/m2/s"           grid_ref="grid_T_2D" />
       <field id="smsdmspd"    long_name="Uppermost DMSPd SMS"                    unit="umolS/m2/s"           grid_ref="grid_T_2D" />
       <field id="smsdms"      long_name="Uppermost DMS SMS"                      unit="umolS/m2/s"           grid_ref="grid_T_2D" />

        <!-- Ice DMS variables and diagnostics -->
       <field id="icedmspd"     long_name="Bottom ice DMSPd"                      unit="umolS/m3"             grid_ref="grid_T_2D" />
       <field id="icedms"       long_name="Bottom ice DMS"                        unit="umolS/m3"             grid_ref="grid_T_2D" />
       <field id="icedmspd_cat" long_name="Bottom ice DMSPd per ice category"     unit="umolS/m3"             grid_ref="grid_T_ncatice" />
       <field id="icedms_cat"   long_name="Bottom ice DMS per ice category"       unit="umolS/m3"             grid_ref="grid_T_ncatice" />
       <field id="icedmspd_gca" long_name="Bottom ice DMSPd grid cell avergae"    unit="umolS/m3"             grid_ref="grid_T_ncatice" />
       <field id="icedms_gca"   long_name="Bottom ice DMS grid cell avergae"      unit="umolS/m3"             grid_ref="grid_T_ncatice" />
       
       <field id="icedmspdrls"  long_name="DMSPd release rate"                    unit="umolS/m3/s"           grid_ref="grid_T_ncatice" />
       <field id="icedmsrls"    long_name="DMS release rate"                      unit="umolS/m3/s"           grid_ref="grid_T_ncatice" />
       <field id="dmsp_exud"    long_name="DMSPd production from exudation"       unit="umolS/m3/s"           grid_ref="grid_T_ncatice" />
       <field id="dmsp_lysis"   long_name="DMSPd production from lysis"           unit="umolS/m3/s"           grid_ref="grid_T_ncatice" />
       <field id="dms_phot"     long_name="DMS photolysis"                        unit="umolS/m3/s"           grid_ref="grid_T_ncatice" />
