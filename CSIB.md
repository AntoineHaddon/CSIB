# CSIB 

Notes on implementation in NEMO4.2 (CanNEMo version) of  Canadian Sea Ice Biogeochemistry model (CSIB) 



## Calling sequence of CSIB

Init of CSIB in ice model init: `icestep.F90/ice_init` 
- calls `trcini_csib.F90/trc_ini_csib`
    - calls `trc_sms_csib.f90/trc_sms_csib_alloc` allocate ice tracer variables
    - calls `trcnam_csib.F90/trc_nam_csib` reads `icetra` namelist from ice namelist files `namelist_ice_ref/cfg`
    - set ice tracer variables from ice restart file or set equal to ocean variables  
- calls `ice_dyn_init` calls `ice_dyn_adv_init` calls `icedyn_adv_pra.F90/adv_pra_init` : read moments from ice tracer restart file or set to 0 in case of spinup

Time stepping of CSIB
- Horizontal transport: by sea ice module
    - surface boundary conditions `sbc` calls sea ice model SI3, 
        - `icedyn_adv` sea ice variables transport, does also CSIB variables horizontal transport (only `icedyn_adv_pra`, which is the default advection scheme)
        - redistribution between sea ice categories of sea ice tracers from thermodynamics and ridging 
            - ridging/rafting: `icedyn_rdgrft`
            - thermodynamics: `itd_shiftice`
        - `ice_rst_write` writes CSIB variables in ice restart file at end of run
- Biogeochemistry and sea ice ocean exchanges: with passive tracers module TOP 
    - TOP main routine `trcstp.F90/trc_stp` calls:
        - `trc_wri` output, calls `trc_wri_csib` CSIB outputs
        - `trc_sms` calls `trc_sms_csib` 
            - first time step:
                - checks that CanOE (ocean BGC model) is active, if not turns off CSIB
                - init of CSIB variables from CanOE variables in case of spin-up or if no CSIB data in ice restart file
            - compute CSIB variables tendencies and time steping
            - compute sea ice-ocean exchanges, modifies tendencies of CanOE varariables, ocean surface only (i.e tr(ji,jj,1,...,Krhs))




## Changes to sea ice model SI3 

- `icestp.F90` 
    - `ice_init` ice init does CSIB init (see above) 
    - `ice_stp` calls `ice_rst_write` which writes CSIB variables in ice restart file
    - `diag_set0` reset diagnostics to 0, add reseting of diagnostics added for CSIB 

- ice dynamics
    - `icedynadv_pra.F90` advection of sea ice trancer with ice motion
        - `ice_dyn_adv_pra` computes advection
        - `adv_pra_init` allocates and initiates momments for transport of sea ice variables, with reading of ice restart files
        - called by `icedyn_adv` which passes `icetra_gca` array of CSIB variables to `ice_dyn_adv_pra`
    - `icedyn_rdgrft` ridging/rafting , leads to redistribution of ice area between sea ice thickness categories -> transport of ice BGC between  sea ice thickness categories

- ice thermodynamics
    - `iceitd.f90/itd_shiftice` does the redistribution after thermodynamics. add redistribution of ice tracers with area of ice redistributed
        - called by `ice_thd_rem`  for remapping thickness distribution
        - also called by `ice_cor` which checks ice variables are in bounds and leads to corrections to ice thicknesses, area, etc
    
- diagnostics added, mainly for sea ice - ocean exchanges of ice tracers. All 
    - ice growth/melt: `icethd_dh`: add diagnostics for recording of rate of ice thickness change from bottom and surface melt, bottom growth; snow thickness change for bottom and surface melt
    - ice area melt: `icethd_da`:  add diagnostics for recording of rate of ice area change from lateral melt
    - melt ponds `ice_thd_pnd`:  add diagnostics for recording of rate of melt pond drainage volume
    -  ice growth in leads `icethd_do`: add diagnostics for recording of rate of ice area change from lateral growth
    - diagnostics declared and allocated in `ice.F90` and `ice1d.F90`


## Changes to tracer module TOP

- `trc_wri` calls `trc_wri_csib` CSIB outputs
- `trc_sms` calls `trc_sms_csib` main CSIB routine
- `trcopt_canbgc` computes photosynthetic active radiation (PAR) in ocean
    - changed PAR computation to have under-ice PAR
    - add attenuation of PAR by ice algae


## Parameters

In `namelist_ice_ref`:

    &namicetra    !   Sea ice tracers (CSIB2)
    !-----------------------------------------------------------------------
        ln_csib         = .true.         ! Run CSIB 
    ! Ice algae 
        z_ia            = 0.03           ! height of skeletal layer (m)
        qnidiamin       = 0.04           ! Mininum ice algae N/C (gN gC-1)
        cn_fct          = 2.61           ! Factor for target ice algae C:N as function of C:Chl (-)
        cn_pow          = 0.33           ! Power exponent for target ice algae C:N as function of C:Chl (-)    
        qchidiaref      = 0.03           ! Reference ice algae Chl:C for uptake from ice growth (gChl gC-1)
        alrefidia       = 202.44         ! Reference initial slope of photosynthesis-irradiance curve (gC gChl-1 (W m-2)-1 d-1)
        pcrefidia       = 1.964          ! Reference photosynthesis rate (gC gChl-1 d-1)    
        betaidia        = 1.             ! photo-inhibition for ice algae (gC gChl-1 (W m-2)-1 d-1) 
        cigr            = 0.015          ! Critical ice growth rate for ice algal limitation  (m d-1)  
        vnref           = 0.6            ! Reference N uptake rate (gN gC d-1)
        knh4            = 0.05           ! NH4 limitation half saturation constant (mmolN m-3)
        kno3            = 1.0            ! NO3 limitation half saturation constant (mmolN m-3)
        etares          = 2.0            ! Respiratory cost of biosynthesis (gC gN-1)
        ch2nmax         = 0.35           ! Maximum Chl synthesis rate to N uptake rate (gChl gN-1)
        min_icedia      = 500.0          ! Mortality threshold for ice algae (mgC m-3)
        t_ia            = 0.0633         ! Temperature sensitivity coefficient for the ice algae (C)-1
        r_m1            = 0.008          ! Linear Mortality rate for ice algae (d-1)  
        r_m2            = 1.E-6          ! Quadratic Mortality rate for ice algae ((mgC m-3)-1 d-1)  
        f_p2            = 0.1            ! Seeding fraction (-)
        f_flsh          = 0.8            ! Flushing fraction (-)
        f_slgh          = 0.8            ! Sloughing fraction (-)
        dt_mo           = 0.15           ! Heating export sea ice warming threshold (deg C d-1)
        t_mo            = -5.0           ! Heating export sea ice temp trheshold (deg C)
        d_mo            = 0.0006         ! Heating export coeffecient ((deg C mg m-3)-1)
    ! Ice nitrogen
        f_rm            = 0.3            ! Remineralization fraction (-)
        r_ni            = 0.01           ! Nitrification rate (d-1 W m-2)  
        c_di            = 1E-9           ! Molecular diffusion coefficient for dissolved nutrients at the ice-water 
        c_nu            = 1.85E-6        ! Kinematic viscosity of seawater (m2/s)
    ! Ice DIC
        sicpump         = .false.        ! Flag for activation of sea ice C pump 
        icedicref       = 0.0004         ! Sea ice reference DIC (mol L-1)
        icetalref       = 0.0005         ! Sea ice reference TA  (mol L-1)
        f_dicsw         = 0.99           ! fraction of DIC rejected into seawater during growth (-)
        f_dicsw_melt    = 0.975          ! fraction of DIC rejected into seawater during melt (-)
    ! Ice DMS
        ln_dmsice       = .false.        ! Run with ice dms - requires ocean DMS (ln_dmsoce in TOP namelist + parameters in CanOE namelists)
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

      <!-- CSIB variables -->
       <field id="icediac"      long_name="Ice algae C biomass"                         unit="mg C m-3"             grid_ref="grid_T_2D" />
       <field id="icedian"      long_name="Ice algae N biomass"                         unit="mg N m-3"             grid_ref="grid_T_2D" />
       <field id="icediach"     long_name="Ice algae Chl biomass"                       unit="mg Chl m-3"           grid_ref="grid_T_2D" />
       <field id="iceno3"       long_name="Ice nitrate concentration"                   unit="mmol N m-3"           grid_ref="grid_T_2D" />
       <field id="icenh4"       long_name="Ice ammonium concentration"                  unit="mmol N m-3"           grid_ref="grid_T_2D" />

       <field id="icediac_cat"  long_name="Ice algae C biomass per ice category"        unit="mg C m-3"             grid_ref="grid_T_ncatice" />
       <field id="icedian_cat"  long_name="Ice algae N biomass per ice category"        unit="mg N m-3"             grid_ref="grid_T_ncatice" />
       <field id="icediach_cat" long_name="Ice algae Chl biomass per ice category"      unit="mg Chl m-3"           grid_ref="grid_T_ncatice" />
       <field id="iceno3_cat"   long_name="Ice nitrate concentration per ice category"  unit="mmol N m-3"           grid_ref="grid_T_ncatice" />
       <field id="icenh4_cat"   long_name="Ice ammonium concentration per ice category" unit="mmol N m-3"           grid_ref="grid_T_ncatice" />

       <field id="icediac_gca"  long_name="Ice algae C biomass grid cell average"       unit="mg C m-3"             grid_ref="grid_T_ncatice" />
       <field id="icedian_gca"  long_name="Ice algae N biomass grid cell average"       unit="mg N m-3"             grid_ref="grid_T_ncatice" />
       <field id="icediach_gca" long_name="Ice algae Chl biomass grid cell average"     unit="mg Chl m-3"           grid_ref="grid_T_ncatice" />
       <field id="iceno3_gca"   long_name="Ice nitrate grid cell average"               unit="mmol N m-3"           grid_ref="grid_T_ncatice" />
       <field id="icenh4_gca"   long_name="Ice ammonium grid cell average"              unit="mmol N m-3"           grid_ref="grid_T_ncatice" />

      <!-- CSIB processes -->
       <field id="flushrate"    long_name="Flushrate per category"                    unit="m s-1"                 grid_ref="grid_T_ncatice" />
       <field id="bogup"        long_name="Water uptake from ice bottom growth"       unit="m s-1"                 grid_ref="grid_T_ncatice" />
       
       <field id="flush_dia"    long_name="Ice diatoms flushrate"                       unit="mg m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="slough_dia"   long_name="Ice diatoms sloughing rate"                  unit="mg m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="lamloss_dia"  long_name="Ice diatoms loss from lateral melt"          unit="mg m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="heatexp_dia"  long_name="Ice diatoms loss from heating export"        unit="mg m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="dt_i"         long_name="Sea ice temperature tendency"                unit="deg C s-1"            grid_ref="grid_T_ncatice" />
       <field id="bogup_dia"    long_name="Ice diatoms uptake from ice bottom growth"   unit="mg m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="lagup_dia"    long_name="Ice diatoms uptake from ice latteral growth" unit="mg m-3 s-1"           grid_ref="grid_T_ncatice" />
       
       <field id="nxsicedia"    long_name="N excess from ice algal export directly remineralized" unit="mg N m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="cxsicedia"    long_name="C excess from ice algal export directly remineralized" unit="mg C m-3 s-1"           grid_ref="grid_T_ncatice" />

       <field id="flush_no3"    long_name="Ice NO3 flushrate"                                 unit="mmol N m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="slough_no3"   long_name="Ice NO3 sloughrate"                                unit="mmol N m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="lamloss_no3"  long_name="Ice NO3 loss from lateral melt"                    unit="mmol N m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="lagup_no3"    long_name="Ice NO3 uptake from ice latteral growth"           unit="mmol N m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="bogup_no3"    long_name="Ice NO3 uptake from ice bottom growth"             unit="mmol N m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="moldif_no3"   long_name="Ocean to ice NO3 flux from molecular diffusion"    unit="mmol N m-3 s-1"           grid_ref="grid_T_ncatice" />

       <field id="flush_nh4"    long_name="Ice NH4 flushrate"                                 unit="mmol N m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="slough_nh4"   long_name="Ice NH4 sloughrate"                                unit="mmol N m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="lamloss_nh4"  long_name="Ice NH4 loss from lateral melt"                    unit="mmol N m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="lagup_nh4"    long_name="Ice NH4 uptake from ice latteral growth"           unit="mmol N m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="bogup_nh4"    long_name="Ice NH4 uptake from ice bottom growth"             unit="mmol N m-3 s-1"           grid_ref="grid_T_ncatice" />
       <field id="moldif_nh4"   long_name="Ocean to ice NH4 flux from molecular diffusion"    unit="mmol N m-3 s-1"           grid_ref="grid_T_ncatice" />

       <field id="phot_dia"     long_name="Ice diatoms photosynthesis rate"            unit="mg C m-3 s-1"              grid_ref="grid_T_ncatice" />
       <field id="lim_PAR"      long_name="Ice diatoms light limitation factor"        unit="-"                         grid_ref="grid_T_ncatice" />
       <field id="lim_nut"      long_name="Ice diatoms nutrient limitation factor"     unit="-"                         grid_ref="grid_T_ncatice" />
       <field id="lim_ice"      long_name="Ice growth limitation factor"               unit="-"                         grid_ref="grid_T_ncatice" />
       <field id="diaupn"       long_name="N uptake by ice diatoms in ice"             unit="mg N m-3 s-1"              grid_ref="grid_T_ncatice" />
       <field id="chlsyn"       long_name="Chl synsthesis by ice diatoms in ice"       unit="mg Chl m-3 s-1"            grid_ref="grid_T_ncatice" />
       <field id="mortlin_dia"  long_name="Ice diatoms linear mortality rate"          unit="mg C m-3 s-1"              grid_ref="grid_T_ncatice" />
       <field id="mortquad_dia" long_name="Ice diatoms quadratic mortality rate"       unit="mg C m-3 s-1"              grid_ref="grid_T_ncatice" />
       <field id="remin_dia"    long_name="Remineralization rate"                      unit="mmol N m-3 s-1"            grid_ref="grid_T_ncatice" />
       <field id="nitri"        long_name="Nitrification rate"                         unit="mmol N m-3 s-1"            grid_ref="grid_T_ncatice" />
      
      <!--CSIB diagnostics -->
       <field id="icenpp"       long_name="Sea ice net primary production"             unit="mg C m-3 s-1"              grid_ref="grid_T_2D" />
       <field id="icenpp_cat"   long_name="Sea ice net primary production per ice category"   unit="mg C m-3 s-1"       grid_ref="grid_T_ncatice" />
       <field id="qnidia"       long_name="Ice diatoms N:C"                            unit="g N g C-1"                 grid_ref="grid_T_ncatice" />
       <field id="qnidiamax"    long_name="Ice diatoms target N:C"                     unit="g N g C-1"                 grid_ref="grid_T_ncatice" />
       <field id="qchidia"      long_name="Ice diatoms Chl:C"                          unit="g Chl g C-1"               grid_ref="grid_T_ncatice" />

      <!-- SI3 variables added for CSIB -->
       <field id="dh_bom_cat"     long_name="ice thickness change from ice bottom melt per category"                  unit="m/s"   grid_ref="grid_T_ncatice" />
       <field id="dh_sum_cat"     long_name="ice thickness change from ice surface melt per category"                 unit="m/s"   grid_ref="grid_T_ncatice" />
       <field id="da_lam_cat"     long_name="fraction of ice concentration loss from ice lateral melt per category"   unit="1/s"   grid_ref="grid_T_ncatice" />
       <field id="dh_snw_sum_cat" long_name="snow thickness change from surface melt per ice category"                unit="m/s"   grid_ref="grid_T_ncatice" />
       <field id="dh_mpdrn_cat"   long_name="flowrate of melt pond drainge volume per sea ice area per ice category"  unit="m/s"   grid_ref="grid_T_ncatice" />
       <field id="dh_bog_cat"     long_name="ice thickness change from ice bottom growth per category"                unit="m/s"   grid_ref="grid_T_ncatice" />
       <field id="da_lag_cat"     long_name="fraction of ice concentration gain from ice lateral melt per category"   unit="1/s"   grid_ref="grid_T_ncatice" />

      <field id="qtr_ice_bot_cat" long_name="solar heat flux transmitted through the ice per category"   unit="W/m2"    grid_ref="grid_T_ncatice" />
      <field id="fric_vel"        long_name="ice-ocean friction velocity"                                unit="m/s"     grid_ref="grid_T_2D" />

      <!-- CanOE variables added for CSIB -->
      <field id="PHY2c_os"        long_name="Ocean surface diatoms C concentration "       unit="mmol/m3"       grid_ref="grid_T_2D" />
      <field id="PHY2N_os"        long_name="Ocean surface diatoms N concentration "       unit="mmol/m3"       grid_ref="grid_T_2D" />
      <field id="PHY2CHL_os"      long_name="Ocean surface diatoms Chlorophyll concentration "       unit="mmol/m3"       grid_ref="grid_T_2D" />
      <field id="NO3_os"          long_name="Ocean surface NO3 concentration "             unit="mmol/m3"       grid_ref="grid_T_2D" />
      <field id="NH4_os"          long_name="Ocean surface NH4 concentration "             unit="mmol/m3"       grid_ref="grid_T_2D" />
      <field id="GOC_os"          long_name="Ocean surface large POC concentration "             unit="mmol/m3"       grid_ref="grid_T_2D" />
      <field id="DIC_os"          long_name="Ocean surface DIC concentration "             unit="mmol/m3"       grid_ref="grid_T_2D" />

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
