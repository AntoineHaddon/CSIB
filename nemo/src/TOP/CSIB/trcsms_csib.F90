MODULE trcsms_csib
   !!======================================================================
   !!                         ***  MODULE trcsms_csib  ***
   !! TOP :   Main module of the CSIB tracers
   !!======================================================================
   !! History :      !  2007  (C. Ethe, G. Madec)  Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_sms_csib       : CSIB model main routine
   !! trc_sms_csib_alloc : allocate arrays specific to CSIB sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables   
   USE dom_oce         ! ocean domain
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc

   USE ice              ! ice variables
   USE phycst           ! physical constants: rhoi, rhos
   USE sbc_oce , ONLY : ssu_m, ssv_m         ! sea surface velocity, for computation of friction velocity
   USE sbc_oce , ONLY : tprecip, sprecip     ! total and solid precipitation
   USE sbc_oce , ONLY : sst_m                ! sea surface temperature (Celsius)
   USE zdfmxl  , ONLY : nmln                 ! level of mixed layer depth for dic/tak fluxes

   USE par_canoe        ! indices of CanOE model variables, e.g. jrdia: diatoms

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_csib       ! called by trcsms.F90 module
   PUBLIC   trc_sms_csib_alloc ! called by trcini_csib.F90 module

!! * Substitutions
#  include "do_loop_substitute.h90"

   ! Defined HERE the arrays specific to CSIB sms and ALLOCATE them in trc_sms_csib_alloc
   REAL(wp), PUBLIC, PARAMETER ::   epsi10 = 1.e-10_wp  !: small number


   !! Sea ice tracers are like ice model variables and have 2 equivalent variables: 
   !!    - one extenisve ("global") for dynamics, ice transport
   !!    - one intensive ("equivalent") for biogeochemistry
   !! 
   !! For extensive variables we use grid cell averages:
   !! = in situ concentration * sea ice concentration
   !! = mass/sea ice area     * sea ice area / cell area  / height
   !! = mass per unit area / height
   !!
   !! ** variables 
   !! name        | indice                                                 |
   !!-------------|-------------|---------------------------------|------------|
   !! icedia      | jridiac     | Ice diatoms C  per ice area     | mg C /m3   |
   !! icedia      | jridian     | Ice diatoms N per ice area      | mmol N /m3 |
   !! icedia      | jridiach    | Ice diatoms Chl per ice area    | mg Chl /m3 |
   !! iceno3      | jrino3      | Ice NO3 per ice area            | mmol N /m3 |
   !! icenh4      | jrinh4      | Ice NH4 per ice area            | mmol N /m3 |
   !!-------------|-------------|---------------------------------|------------|
   !!

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)     :: icetra            !  Ice tracer per ice area (4d: 2d horizontal * ice cat * ice tracers)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)     :: icetra_gca        !  Ice tracer grid cell average
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)       :: icetragca_2d      !  Ice tracer grid cell average, 2d version for ice model

   !! Biomass ratios (used to convert C fluxes to N or Chl fluxes)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)       :: qnidia             !  Ice diatoms N/C 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)       :: qchidia            !  Ice diatoms Chl/C


   !!
   !! Sources and sinks
   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flushrate        !  Flushrate per ice category (m/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup            !  Flowrate of water uptake from bottom ice growth per ice category (m/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup            !  Flowrate of water uptake per ice area from lateral ice growth (m/s)
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flush_dia        !  Loss rate of ice diatoms from flushing per ice category (mg C/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lamloss_dia      !  Loss rate from lateral melt per ice category (mg C/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup_dia        !  Diatoms uptake rate from bottom ice growth per ice category (mg C/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_dia        !  Diatoms uptake rate from lateral ice growth  per ice category (mg C/m3/s)
   
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flush_no3        !  Loss rate of ice no3 from flushing per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flush_nh4        !  Loss rate of ice nh4 from flushing per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lamloss_no3      !  Loss rate from lateral melt per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lamloss_nh4      !  Loss rate from lateral melt per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_no3        !  NO3 uptake rate from lateral ice growth  per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_nh4        !  NH4 uptake rate from lateral ice growth  per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup_no3        !  No3 uptake rate from bottom ice growth per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup_nh4        !  Nh4 uptake rate from bottom ice growth per ice category (mmol/m3/s)
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fric_vel         !  Ice frictional velocity (m/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: moldif_no3       !  Molecular diffusion ratea at ice ocean interface for NO3 per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: moldif_nh4       !  Molecular diffusion ratea at ice ocean interface for NH4  per ice category (mmol/m3/s)
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: phot_dia         !  Diatoms growth rate per ice category (mg C/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lim_PAR          !  Diatoms light limitation factor per ice category (-)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lim_nut          !  Diatoms nutrients limitation factor per ice category (-)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: diaupn           !  N uptake by ice diatoms per ice category (mg N/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: chlsyn           !  Ice diatoms Chl synthesis rate per ice category (mg Chl/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: mortlin_dia      !  Linear mortality rate of ice diatoms per ice category (mg C/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: mortquad_dia     !  Quadratic mortality rate of ice diatoms per ice category (mg C/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: remin_dia        !  N remineralization rate in ice per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nitri            !  Nitrification rate in ice per ice category (mmol/m3/s)


   ! model parameters
   REAL(wp), PUBLIC, SAVE ::   z_ia                 ! height of skeletal layer
   ! ice diatoms
   REAL(wp), PUBLIC, SAVE ::   qnidiamin            ! Mininum ice diatom N/C (mmolN mmolC-1)
   REAL(wp), PUBLIC, SAVE ::   qnidiamax            ! Maximum ice diatom N/C (mmolN mmolC-1)
   REAL(wp), PUBLIC, SAVE ::   ear=4498._wp         ! activation energy / R gas constant (deg K)
   REAL(wp), PUBLIC, SAVE ::   tempref=298.15_wp    ! Reference temperature for photosynthesis (deg K)
   REAL(wp), PUBLIC, SAVE ::   alphaidia            ! Initial slope of photosynthesis-irradiance curve (mgC mgChl-1 (W m-2)-1 d-1)) / molC/gC
   REAL(wp), PUBLIC, SAVE ::   pcrefidia            ! Reference photosynthesis rate (d-1) / sec per day
   REAL(wp), PUBLIC, SAVE ::   vnref                ! Reference N uptake rate (gN gC d-1) / sec per day
   REAL(wp), PUBLIC, SAVE ::   knh4                 ! NH4 limitation half saturation constant (mmol N m-3)
   REAL(wp), PUBLIC, SAVE ::   kno3                 ! NO3 limitation half saturation constant (mmol N m-3)
   REAL(wp), PUBLIC, SAVE ::   etares               ! Respiratory cost of biosynthesis (gC gN-1)
   REAL(wp), PUBLIC, SAVE ::   ch2nmax              ! Maximum Chl synthesis rate to N uptake rate (gChl gN-1)
   REAL(wp), PUBLIC, SAVE ::   min_icedia           ! Mortality threshold for ice diatoms (mmol C m-3)
   REAL(wp), PUBLIC, SAVE ::   t_ia                 ! Temperature sensitivity coefficient for ice diatoms (C)-1
   REAL(wp), PUBLIC, SAVE ::   r_m1                 ! Linear Mortality rate for ice diatoms (d-1)  / sec per day
   REAL(wp), PUBLIC, SAVE ::   r_m2                 ! Quadratic Mortality rate for ice diatoms (mmol C m-3 d-1)  / sec per day

   REAL(wp), PUBLIC, SAVE ::   f_p2                 ! Seeding fraction (-)
   ! ice N
   REAL(wp), PUBLIC, SAVE ::   f_rm                 ! Remineralization fraction (-)
   REAL(wp), PUBLIC, SAVE ::   r_ni                 ! Nitrification rate (d-1 W m-2)  / sec per day
   REAL(wp), PUBLIC, SAVE ::   c_di                 ! Molecular diffusion coefficient for dissolved nutrients at the ice-water interface (m/s2)
   REAL(wp), PUBLIC, SAVE ::   c_nu                 ! Kinematic viscosity of seawater (m2/s)
  ! sea ice C pump
   LOGICAL , PUBLIC, SAVE ::   sicpump              ! Flag for activation of sea ice C pump
   REAL(wp), PUBLIC, SAVE ::   icedicref            ! Sea ice reference DIC (mmol C m-3?)
   REAL(wp), PUBLIC, SAVE ::   icetalref            ! Sea ice reference TA (mmol C m-3?)
   REAL(wp), PUBLIC, SAVE ::   f_dicsw              ! fraction of DIC rejected into seawater during growth (-)
   REAL(wp), PUBLIC, SAVE ::   f_dicsw_melt         ! fraction of DIC rejected into seawater during melt (-)

  
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms_my_trc.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_csib( kt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_csib  ***
      !!
      !! ** Purpose :   main routine of CSIB model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt               ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs   ! time level indices
      
      INTEGER  :: ji,jj,jl,jn          ! dummy loop index
      REAL(wp) :: mmc=12._wp           ! Molar mass of Carbon (g mol-1)
      REAL(wp) :: mmn=14._wp           ! Molar mass of Nitrogen (g mol-1)
      REAL(wp) :: zscale               ! scale factor between sea ice skeletal layer and ocean surface layer
      REAL(wp) :: if_below_min=0._wp   ! mortality switch
      REAL(wp) :: zdtDif               ! temp variable for molecular diff
      REAL(wp) :: zmaxia               ! for diagnostics/debug
      REAL(wp) :: zln2 = 0.693147_wp   ! natural log ln(2)
      REAL(wp) :: zsigup_tot           ! total water uptake from sea ice growth per ice area
      REAL(wp) :: ztemp                ! Temperature factor
      REAL(wp) :: znut                 ! N switch
      REAL(wp) :: zlim_nh4, zlim_no3   ! NH4 and NO3 limitation factors
      REAL(wp) :: zrhoch               ! for Chl synthesis

      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_csib')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_csib:  CSIB model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)
      ! IF(lwp) WRITE(numout,*) ' time step kt ', kt
      ! IF(lwp) WRITE(numout,*) ' rDt_trc ', rDt_trc
      ! IF(lwp) WRITE(numout,*) ' rn_Dt ', rn_Dt
      ! IF(lwp) WRITE(numout,*) ' nsec_day ', nsec_day
      ! IF(lwp) WRITE(numout,*)
      ! IF(lwp) WRITE(numout,*)

      ! Initiation from ocean surface concentrations (need to do it here and not in trcini_csib because CanOE initiation occurs after?)
      IF ( (kt == 1) .AND. (.NOT. ln_rsttr) ) THEN
         IF(lwp) WRITE(numout,*) 'Init from ocean surface'
         DO jl = 1, jpl ! loop ice categories
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( a_i(ji,jj,jl) > 1.e-4_wp ) THEN ! if ice
                     icetra(ji,jj,jl,jridiac) = tr(ji,jj,1,jrdia,Kmm) *mmc ! convert to mass units
                     ! icetra(ji,jj,jl,jridian) = tr(ji,jj,1,jrdn,Kmm) *mmn ! convert to mass units
                     icetra(ji,jj,jl,jridian) = icetra(ji,jj,jl,jridiac) /8._wp ! init at C:N=8
                     ! icetra(ji,jj,jl,jridiach) = tr(ji,jj,1,jrdch,Kmm)
                     icetra(ji,jj,jl,jridiach) = icetra(ji,jj,jl,jridiac) /10._wp ! init at C:Chl=10
                     icetra(ji,jj,jl,jrino3) = tr(ji,jj,1,jqno3,Kmm)
                     icetra(ji,jj,jl,jrinh4) = tr(ji,jj,1,jrnh4,Kmm)
                  END IF
               ENDDO
            ENDDO
         ENDDO
         DO jn = 1,jp_csib
            icetra_gca(:,:,:,jn) = icetra(:,:,:,jn) * a_i(:,:,:)
         ENDDO
      END IF

      ! Conversion from global (extensive used by ice transport model) to equivalent (intensive, used by bgc model) variables
      ! and set to 0 if low ice concentration 
      DO jl = 1, jpl ! loop ice categories
         DO jj = 1, jpj
            DO ji = 1, jpi
               
               IF( a_i(ji,jj,jl) > 1.e-4_wp ) THEN ! if ice
                  DO jn = 1,jp_csib
                     icetra(ji,jj,jl,jn) = icetra_gca(ji,jj,jl,jn) / a_i(ji,jj,jl)
                  ENDDO
                  
               ELSE ! low ice
                  ! set tracers to 0 
                  icetra(ji,jj,jl,:)=0._wp
                  ! send what was there to the ocean
                  zscale = z_ia / e3t_0(ji,jj,1) ! dilution ratio = skeletal layer height/ surface ocean layer height
                  tr(ji,jj,1,jrdia,Kmm) = tr(ji,jj,1,jrdia,Kmm) + icetra_gca(ji,jj,jl,jridiac) * zscale /mmc ! convert to mmol
                  tr(ji,jj,1,jrdn ,Kmm) = tr(ji,jj,1,jrdn ,Kmm) + icetra_gca(ji,jj,jl,jridian) * zscale /mmn ! convert to mmol
                  tr(ji,jj,1,jrdch,Kmm) = tr(ji,jj,1,jrdch,Kmm) + icetra_gca(ji,jj,jl,jridiach) * zscale
                  tr(ji,jj,1,jqno3,Kmm) = tr(ji,jj,1,jqno3,Kmm) + icetra_gca(ji,jj,jl,jrino3) * zscale
                  tr(ji,jj,1,jrnh4,Kmm) = tr(ji,jj,1,jrnh4,Kmm) + icetra_gca(ji,jj,jl,jrinh4) * zscale
               ENDIF ! if ice
               
            ENDDO ! loop jpi
         ENDDO ! loop jpj
      ENDDO ! loop jpl ice categories
      
      
      ! reset rates: so that they are 0 where there is no ice
      flushrate(:,:,:) = 0._wp
      bogup(:,:,:) = 0._wp
      lagup(:,:,:) = 0._wp
      
      flush_dia(:,:,:) = 0._wp
      lamloss_dia(:,:,:) = 0._wp
      bogup_dia(:,:,:) = 0._wp
      lagup_dia(:,:,:) = 0._wp
      
      flush_no3(:,:,:) = 0._wp
      lamloss_no3(:,:,:) = 0._wp
      moldif_no3(:,:,:) = 0._wp
      lagup_no3(:,:,:) = 0._wp
      bogup_no3(:,:,:) = 0._wp
      
      flush_nh4(:,:,:) = 0._wp
      lamloss_nh4(:,:,:) = 0._wp
      moldif_nh4(:,:,:) = 0._wp
      lagup_nh4(:,:,:) = 0._wp
      bogup_nh4(:,:,:) = 0._wp
      
      phot_dia(:,:,:) = 0._wp
      lim_PAR(:,:,:) = 0._wp
      lim_nut(:,:,:) = 0._wp
      diaupn(:,:,:) = 0._wp
      chlsyn(:,:,:) = 0._wp
      mortlin_dia(:,:,:) = 0._wp
      mortquad_dia(:,:,:) = 0._wp
      remin_dia(:,:,:) = 0._wp
      nitri(:,:,:) = 0._wp
      
      qnidia(:,:,:) = 0._wp
      qchidia(:,:,:) = 0._wp

      ! Compute friction velocity, for molecular diffusion at sea ice ocean interface
      CALL ice_friction_velocity
      
      
      
      DO jj = 1, jpj
         DO ji = 1, jpi
            
            ! Bottom ice variables
            DO jl = 1, jpl ! loop ice cat
               
               IF( a_i(ji,jj,jl) > epsi10 ) THEN ! precence of ice

                  qnidia(ji,jj,jl) = icetra(ji,jj,jl,jridian) /  (icetra(ji,jj,jl,jridiac)+rtrn )
                  qchidia(ji,jj,jl) = icetra(ji,jj,jl,jridiach) /  (icetra(ji,jj,jl,jridiac)+rtrn )
                  
               ! Loss of ice tracers from ice-ocean exchanges

                  ! flushrate: water flowrate per ice area (m/s) 
                  flushrate(ji,jj,jl) = &
                           ! change of ice thickness from bottom+surface melt (m/s) * ice density (kg/m3) / freshwater density (kg/m3)
                     &     ( dh_bom_cat(ji,jj,jl) + dh_sum_cat(ji,jj,jl) ) * rhoi / rhow    &
                           ! change of snow thickness from surface melt (m/s) * snow density (kg/m3) / freshwater density (kg/m3)
                     &     + dh_snw_sum_cat(ji,jj,jl) * rhos / rhow     &
                           ! flowrate of melt pond drainage volume per ice area (m/s)
                     &     + dh_mpdrn_cat(ji,jj,jl)      &
                           ! SIC * (total precipation (Kg/m2/s) - solid precipitation (Kg/m2/s) ) / freshwater density (kg/m3)
                   ! &     + a_i(ji,jj,jl) * MAX( 0._wp, tprecip(ji,jj) - sprecip(ji,jj) ) / rhow
                     &     + a_i(ji,jj,jl) * tprecip(ji,jj) / rhow ! dev run : tprecip seems to be only rain

                  ! flushing of ice tracers: flushrate * ice tracers concentration/ height of skeletal layer  = mass flux (mg/m3/s)
                  flush_dia(ji,jj,jl) = flushrate(ji,jj,jl)/z_ia  * icetra(ji,jj,jl,jridiac)  ! ice diatoms   
                  flush_no3(ji,jj,jl) = flushrate(ji,jj,jl)/z_ia  * icetra(ji,jj,jl,jrino3)  ! ice no3   
                  flush_nh4(ji,jj,jl) = flushrate(ji,jj,jl)/z_ia  * icetra(ji,jj,jl,jrinh4)  ! ice nh4   

                  ! loss of ice tracers from lateral melt : fraction of ice concentration lost (1/s) * ice tracers concentration (mg/m3)
                  lamloss_dia(ji,jj,jl) = da_lam_cat(ji,jj,jl) * icetra(ji,jj,jl,jridiac)   ! ice diatoms
                  lamloss_no3(ji,jj,jl) = da_lam_cat(ji,jj,jl) * icetra(ji,jj,jl,jrino3)   ! ice no3
                  lamloss_nh4(ji,jj,jl) = da_lam_cat(ji,jj,jl) * icetra(ji,jj,jl,jrinh4)   ! ice nh4


               ! Uptake of ice tracers from ice growth

                  ! bogup: water uptake flowrate per ice area from bottom ice growth (m/s) 
                  ! = rate of ice thickness change (m/s) * ice density (kg/m3) / freshwater density (kg/m3)
                  bogup(ji,jj,jl) = dh_bog_cat(ji,jj,jl) * rhoi / rhow

                  ! Uptake of ice tracers from bottom ice growth 
                  ! = flowrate per ice area (m/s) * ocean surface concentration (mmol/m3) /skeletal layer (m) (*molar mass)
                  bogup_dia(ji,jj,jl) = bogup(ji,jj,jl) * tr(ji,jj,1,jrdia,Kbb) /z_ia *mmc
                  bogup_no3(ji,jj,jl) = bogup(ji,jj,jl) * tr(ji,jj,1,jqno3,Kbb) /z_ia
                  bogup_nh4(ji,jj,jl) = bogup(ji,jj,jl) * tr(ji,jj,1,jrnh4,Kbb) /z_ia

                  ! lagup: water uptake flowrate per ice area from lateral ice growth (m/s) 
                  ! = (sic increase rate / sic) (1/s) * (height new ice) (s) * ice density (kg/m3) / freshwater density (kg/m3)
                  lagup(ji,jj,jl) = da_lag_cat(ji,jj,jl) * ht_i_new(ji,jl) * rhoi / rhow

                  ! Uptake of ice diatoms from lateral ice growth  
                  lagup_dia(ji,jj,jl) = &
                        ! water uptake flowrate per ice area * skeletal layer * ocean surface diatoms concentration *molar mass    
                        & lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrdia,Kbb) *mmc     &
                        ! - ice tracer conc * (sic increase rate / sic)
                        & - icetra(ji,jj,jl,jridiac) * da_lag_cat(ji,jj,jl)
                  ! Uptake of ice N from lateral ice growth  
                  lagup_no3(ji,jj,jl) = lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jqno3,Kbb) - icetra(ji,jj,jl,jrino3) * da_lag_cat(ji,jj,jl)
                  lagup_nh4(ji,jj,jl) = lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrnh4,Kbb) - icetra(ji,jj,jl,jrinh4) * da_lag_cat(ji,jj,jl)


                  
               ! Biogeochemical processes
                  
                  ! N limitation factor (for photosynthesis) (-)
                  lim_nut(ji,jj,jl) = MIN( MAX(0._wp, (qnidia(ji,jj,jl) - qnidiamin) / (qnidiamax - qnidiamin) ) , 1._wp)

                  ! Temperature factor (-)
                  ztemp = EXP(-ear*( 1._wp/(sst_m(ji,jj)+273.15_wp) - 1._wp/tempref) ) 

                  ! Light limitation factor (-), from qtr_ice_bot: shortwave radiation transmitted through ice (W/m2) (convert to PAR?)
                  lim_PAR(ji,jj,jl) = 1._wp - exp(- alphaidia/(pcrefidia * ztemp * lim_nut(ji,jj,jl) +rtrn) * qchidia(ji,jj,jl) * qtr_ice_bot(ji,jj,jl) ) 
                  
                  ! Photosynthesis rate (mg C s-1)
                  phot_dia(ji,jj,jl) = pcrefidia * ztemp * lim_nut(ji,jj,jl) * lim_PAR(ji,jj,jl) * icetra(ji,jj,jl,jridiac)
                  
                  ! N uptake switch (-)
                  znut = ( MIN( MAX(0._wp, (qnidiamax - qnidia(ji,jj,jl)) / (qnidiamax - qnidiamin) ), 1._wp) )**0.05

                  ! nh4 and no3 limitation factors (-)
                  zlim_nh4 = icetra(ji,jj,jl,jrinh4) / (icetra(ji,jj,jl,jrinh4) + knh4)
                  zlim_no3 = icetra(ji,jj,jl,jrino3) / (icetra(ji,jj,jl,jrino3) + kno3)

                  ! N uptake by ice diatoms (mg N s-1)
                  diaupn(ji,jj,jl) = vnref * ztemp * znut * (zlim_nh4 + (1._wp - zlim_nh4) * zlim_no3) * icetra(ji,jj,jl,jridiac)

                  ! Chl synthesis rate (mg Chl m-3 s-1)
                  zrhoch = pcrefidia * ztemp * lim_nut(ji,jj,jl) * lim_PAR(ji,jj,jl) / ( alphaidia * qchidia(ji,jj,jl) * qtr_ice_bot(ji,jj,jl) + rtrn )
                  chlsyn(ji,jj,jl) = zrhoch * ch2nmax * diaupn(ji,jj,jl)

                  ! Ice diatom mortality, off if below threshold min_icedia
                  if_below_min = MAX( 0._wp , SIGN(1._wp, icetra(ji,jj,jl,jridiac) - min_icedia) )
                  ! Linear mortality (mg C m-3 s-1)
                  mortlin_dia(ji,jj,jl)  = if_below_min * r_m1 * exp(t_ia* sst_m(ji,jj)) * icetra(ji,jj,jl,jridiac)
                  ! Quadratic mortality (mg C m-3 s-1)
                  mortquad_dia(ji,jj,jl) = if_below_min * r_m2 * icetra(ji,jj,jl,jridiac) * icetra(ji,jj,jl,jridiac)

                  ! Remineralization (mmol N m-3 s-1)
                  remin_dia(ji,jj,jl) = f_rm * mortlin_dia(ji,jj,jl) * qnidia(ji,jj,jl) /mmn

                  ! Nitrification, reduced by light
                  nitri(ji,jj,jl) = r_ni / (1.0_wp+qtr_ice_bot(ji,jj,jl)) * icetra(ji,jj,jl,jrinh4)




                  ! Ice diatoms C biomass dynamics
                  icetra(ji,jj,jl,jridiac) = icetra(ji,jj,jl,jridiac) + rn_Dt * (      &
                              &        - flush_dia(ji,jj,jl)                           & ! loss from flushing 
                              &        - lamloss_dia(ji,jj,jl)                         & ! loss from lateral melting of ice
                              &        + bogup_dia(ji,jj,jl)                           & ! Uptake from bottom ice growth
                              &        + lagup_dia(ji,jj,jl)                           & ! Uptake from lateral ice growth
                              &        + phot_dia(ji,jj,jl)                            & ! Photosynthesis
                              &        - etares * diaupn(ji,jj,jl)                     & ! Respiration
                              &        - mortlin_dia(ji,jj,jl)                         & ! linear mortality
                              &        - mortquad_dia(ji,jj,jl)                        & ! quadratic mortality
                              )
                  ! guarantee positive concentration
                  icetra(ji,jj,jl,jridiac) = MAX(0._wp, icetra(ji,jj,jl,jridiac) )
                  ! icetra(ji,jj,jl,jridiac) = MIN(1000._wp, icetra(ji,jj,jl,jridiac) ) 


                  ! Ice diatoms N biomass dynamics
                  icetra(ji,jj,jl,jridian) = icetra(ji,jj,jl,jridian) + rn_Dt * (            &
                              &    - flush_dia(ji,jj,jl) * qnidia(ji,jj,jl)                  & ! loss from flushing 
                              &    - lamloss_dia(ji,jj,jl) * qnidia(ji,jj,jl)                & ! loss from lateral melting of ice 
                              &    + bogup(ji,jj,jl) * tr(ji,jj,1,jrdn,Kbb) /z_ia *mmn       & ! Uptake from bottom ice growth
                              &    + lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrdn,Kbb) *mmn - icetra(ji,jj,jl,jridian) * da_lag_cat(ji,jj,jl)   & ! Uptake from lateral ice growth
                              &    + diaupn(ji,jj,jl)                                        & ! NO3+NH4 uptake
                              &    - mortlin_dia(ji,jj,jl) * qnidia(ji,jj,jl)                & ! linear mortality
                              &    - mortquad_dia(ji,jj,jl) * qnidia(ji,jj,jl)               & ! quadratic mortality
                              )
                  ! guarantee positive concentration
                  icetra(ji,jj,jl,jridian) = MAX(0._wp, icetra(ji,jj,jl,jridian) )


                  ! Ice diatoms Chl biomass dynamics
                  icetra(ji,jj,jl,jridiach) = icetra(ji,jj,jl,jridiach) + rn_Dt * (          &
                              &    - flush_dia(ji,jj,jl) * qchidia(ji,jj,jl)                 & ! loss from flushing
                              &    - lamloss_dia(ji,jj,jl) * qchidia(ji,jj,jl)               & ! loss from lateral melting of ice 
                              &    + bogup(ji,jj,jl) * tr(ji,jj,1,jrdch,Kbb) /z_ia           & ! Uptake from bottom ice growth 
                              &    + lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrdch,Kbb) - icetra(ji,jj,jl,jridiach) * da_lag_cat(ji,jj,jl)   & ! Uptake from lateral ice growth
                              &    + chlsyn(ji,jj,jl)                                        & ! Chl synthesis
                              &    - mortlin_dia(ji,jj,jl) * qchidia(ji,jj,jl)               & ! linear mortality
                              &    - mortquad_dia(ji,jj,jl) * qchidia(ji,jj,jl)              & ! quadratic mortality
                              )
                  ! guarantee positive concentration
                  icetra(ji,jj,jl,jridiach) = MAX(0._wp, icetra(ji,jj,jl,jridiach) )


                  ! Ice NO3 dynamics
                  icetra(ji,jj,jl,jrino3) = icetra(ji,jj,jl,jrino3) + rn_Dt * (        &
                              &        - flush_no3(ji,jj,jl)                           & ! loss from flushing 
                              &        - lamloss_no3(ji,jj,jl)                         & ! loss from lateral melting of ice
                              &        + bogup_no3(ji,jj,jl)                           & ! Uptake from bottom ice growth
                              &        + lagup_no3(ji,jj,jl)                           & ! Uptake from lateral ice growth
                              &        - diaupn(ji,jj,jl)/mmn * (1._wp - zlim_nh4) * zlim_no3 / (zlim_nh4 + (1._wp - zlim_nh4) * zlim_no3 +rtrn)   & ! uptake by ice diatoms
                              &        + nitri(ji,jj,jl)                               & ! nitrification     
                              )
                  ! guarantee positive concentration
                  icetra(ji,jj,jl,jrino3) = MAX(0._wp, icetra(ji,jj,jl,jrino3) )


                  ! Ice NH4 dynamics
                  icetra(ji,jj,jl,jrinh4) = icetra(ji,jj,jl,jrinh4) + rn_Dt * (        &
                              &        - flush_nh4(ji,jj,jl)                           & ! loss from flushing 
                              &        - lamloss_nh4(ji,jj,jl)                         & ! loss from lateral melting of ice
                              &        + bogup_nh4(ji,jj,jl)                           & ! Uptake from bottom ice growth
                              &        + lagup_nh4(ji,jj,jl)                           & ! Uptake from lateral ice growth
                              &        - diaupn(ji,jj,jl)/mmn * zlim_nh4 / (zlim_nh4 + (1._wp - zlim_nh4) * zlim_no3 +rtrn)    & ! uptake by ice diatoms
                              &        - nitri(ji,jj,jl)                               & ! nitrification                
                              &        + remin_dia(ji,jj,jl)                           & ! remineraliztion                
                              )
                  ! guarantee positive concentration
                  icetra(ji,jj,jl,jrinh4) = MAX(0._wp, icetra(ji,jj,jl,jrinh4) )


               ! Diffusion of N at ice ocean interface: 
                  ! can cause numerical problems if large time step and explicit euler scheme. 
                  ! instead implicit euler scheme after all other calculations
                  zdtDif = rn_Dt * c_di / c_nu * abs(fric_vel(ji,jj)) / z_ia
                  zscale = a_i(ji,jj,jl) * z_ia / e3t_0(ji,jj,1)
                  ! NO3
                  icetra(ji,jj,jl,jrino3) = ( (1._wp+zdtDif*zscale)*icetra(ji,jj,jl,jrino3) + zdtDif*tr(ji,jj,1,jqno3,Kmm) ) / ( 1._wp + zdtDif*(1._wp+zscale) )
                  tr(ji,jj,1,jqno3,Kmm) = ( zdtDif*zscale*icetra(ji,jj,jl,jrino3) + (1._wp+zdtDif)*tr(ji,jj,1,jqno3,Kmm) ) / ( 1._wp + zdtDif*(1._wp+zscale) )
                  ! NH4
                  icetra(ji,jj,jl,jrinh4) = ( (1._wp+zdtDif*zscale)*icetra(ji,jj,jl,jrinh4) + zdtDif*tr(ji,jj,1,jrnh4,Kmm) ) / ( 1._wp + zdtDif*(1._wp+zscale) )
                  tr(ji,jj,1,jrnh4,Kmm) = ( zdtDif*zscale*icetra(ji,jj,jl,jrinh4) + (1._wp+zdtDif)*tr(ji,jj,1,jrnh4,Kmm) ) / ( 1._wp + zdtDif*(1._wp+zscale) )

                  ! Molecular diffusion flux, for ouput 
                  ! D / (nu / |friction velocty|) * ( N_ocean - N_ice ) / skeletal layer
                  moldif_no3(ji,jj,jl) = c_di / c_nu * abs(fric_vel(ji,jj)) * ( tr(ji,jj,1,jqno3,Kmm) - icetra(ji,jj,jl,jrino3) ) /z_ia
                  moldif_nh4(ji,jj,jl) = c_di / c_nu * abs(fric_vel(ji,jj)) * ( tr(ji,jj,1,jrnh4,Kmm) - icetra(ji,jj,jl,jrinh4) ) /z_ia



               ! Sea ice C pump
                  IF( sicpump ) THEN
                     ! total water uptake from ice growth: sea ice conc * ( seawater volume uptake from bottom + lateral ice growth per ice area) (m/s)
                     zsigup_tot = a_i(ji,jj,jl) * ( bogup(ji,jj,jl) + lagup(ji,jj,jl) ) 
                     
                     ! DIC flux from ice growth
                     tr(ji,jj,nmln(ji,jj)-1,jqdic, Krhs) = tr(ji,jj,nmln(ji,jj)-1,jqdic, Krhs) + zsigup_tot * ( tr(ji,jj,nmln(ji,jj)-1,jqdic, Kbb) - icedicref ) * f_dicsw / e3t_0(ji,jj,nmln(ji,jj)-1)
                     ! DIC flux from ice melt
                     tr(ji,jj,1,jqdic, Krhs) = tr(ji,jj,1,jqdic, Krhs) + a_i(ji,jj,jl) * flushrate(ji,jj,jl) * (tr(ji,jj,1,jqdic, Kbb) - icedicref ) * f_dicsw_melt / e3t_0(ji,jj,1)
                     
                     ! TA flux from ice growth
                     tr(ji,jj,nmln(ji,jj)-1,jqtal, Krhs) = tr(ji,jj,nmln(ji,jj)-1,jqtal, Krhs) + zsigup_tot * ( tr(ji,jj,nmln(ji,jj)-1,jqtal, Kbb) - icetalref ) / e3t_0(ji,jj,nmln(ji,jj)-1)
                     ! TA flux from ice melt
                     tr(ji,jj,1,jqtal, Krhs) = tr(ji,jj,1,jqtal, Krhs) + a_i(ji,jj,jl) * flushrate(ji,jj,jl) * (tr(ji,jj,1,jqtal, Kbb) - icetalref ) / e3t_0(ji,jj,1)
                  ENDIF ! if sea ice C pump

               ENDIF ! if ice
               
            ENDDO ! loop jpl ice categories


         ! ocean surface variables affected by exchanges between sea ice and ocean
            ! done after fluxes have been calculated for all ice categories 
            DO jl = 1, jpl ! loop ice cat
         
               IF( a_i(ji,jj,jl) > epsi10 ) THEN ! precence of ice

                  ! scaling factor: conversion of flowrate per sea ice area to flowrate per unit volume
                  ! = sea ice concentration / ocean surface layer height
                  ! = sea ice area / (cell area * ocean surface layer height)
                  zscale = a_i(ji,jj,jl) / e3t_0(ji,jj,1)

               ! Ocean surface large phytoplankton C biomass
                  tr(ji,jj,1,jrdia,Krhs) = tr(ji,jj,1,jrdia,Krhs) + zscale * (    &
                  &        + f_p2 * flush_dia(ji,jj,jl) * z_ia                    & ! flushing of ice diatoms
                  &        + f_p2 * lamloss_dia(ji,jj,jl) * z_ia                  & ! ice diatoms from lateral melting of ice
                  &        - bogup_dia(ji,jj,jl) * z_ia                           & ! Uptake from bottom ice growth
                  &        - lagup(ji,jj,jl) * tr(ji,jj,1,jrdia,Kbb)*mmc          & ! Uptake from lateral ice growth
                  ) /mmc ! convert to mmol

               ! Ocean surface large phytoplankton N biomass
                  tr(ji,jj,1,jrdn,Krhs) = tr(ji,jj,1,jrdn,Krhs) + zscale * (              &
                  &        + f_p2 * flush_dia(ji,jj,jl) * qnidia(ji,jj,jl) * z_ia /mmc    & ! flushing of ice diatoms 
                  &        + f_p2 * lamloss_dia(ji,jj,jl) * qnidia(ji,jj,jl) * z_ia /mmc  & ! ice diatoms from lateral melting of ice
                  &        - bogup(ji,jj,jl) * tr(ji,jj,1,jrdn,Kbb)                       & ! Uptake from bottom ice growth
                  &        - lagup(ji,jj,jl) * tr(ji,jj,1,jrdn,Kbb)                       & ! Uptake from lateral ice growth
                  )

               ! Ocean surface large phytoplankton Chl biomass
                  tr(ji,jj,1,jrdch,Krhs) = tr(ji,jj,1,jrdch,Krhs) + zscale * (         &
                  &        + f_p2 * flush_dia(ji,jj,jl) * qchidia(ji,jj,jl) * z_ia     & ! flushing of ice diatoms 
                  &        + f_p2 * lamloss_dia(ji,jj,jl) * qchidia(ji,jj,jl) * z_ia   & ! ice diatoms from lateral melting of ice
                  &        - bogup(ji,jj,jl) * tr(ji,jj,1,jrdch,Kbb)                   & ! Uptake from bottom ice growth
                  &        - lagup(ji,jj,jl) * tr(ji,jj,1,jrdch,Kbb)                   & ! Uptake from lateral ice growth
                  )

               ! Ocean surface large POC
                  tr(ji,jj,1,jrgoc,Krhs) = tr(ji,jj,1,jrgoc,Krhs) + zscale * (            &
                  &        + (1.0_wp - f_p2) * flush_dia(ji,jj,jl) * z_ia                 & ! flushing of ice diatoms
                  &        + (1.0_wp - f_p2) * lamloss_dia(ji,jj,jl) * z_ia               & ! ice diatoms from lateral melting of ice
                  ! &        + (1.0_wp - f_rm) * mortlin_dia(ji,jj,jl) * z_ia               & ! linear mortality 
                  ! &        + mortquad_dia(ji,jj,jl) * z_ia                                & ! quadratic mortality 
                  ) /mmc ! convert to mmol

               ! Ocean surface NO3 
                  tr(ji,jj,1,jqno3,Krhs) = tr(ji,jj,1,jqno3,Krhs) + zscale * (   &
                  &        + flush_no3(ji,jj,jl) * z_ia                          & ! flushing of ice NO3
                  &        + lamloss_no3(ji,jj,jl) * z_ia                        & ! ice NO3 from lateral melting of ice
                  &        - bogup_no3(ji,jj,jl) * z_ia                          & ! Uptake from bottom ice growth
                  &        - lagup(ji,jj,jl) * tr(ji,jj,1,jqno3,Kbb)             & ! Uptake from lateral ice growth
                  )

               ! Ocean surface NH4 
                  tr(ji,jj,1,jrnh4,Krhs) = tr(ji,jj,1,jrnh4,Krhs) + zscale * (   &
                  &        + flush_nh4(ji,jj,jl) * z_ia                          & ! flushing of ice NO3
                  &        + lamloss_nh4(ji,jj,jl) * z_ia                        & ! ice NO3 from lateral melting of ice
                  &        - bogup_nh4(ji,jj,jl) * z_ia                          & ! Uptake from bottom ice growth
                  &        - lagup(ji,jj,jl) * tr(ji,jj,1,jrnh4,Kbb)             & ! Uptake from lateral ice growth
                  )

                  
               ENDIF ! if ice
               
            ENDDO ! loop jpl ice categories
               
         ENDDO ! loop jpi
      ENDDO ! loop jpj





      ! ! For debug/diagnostics: print max ice diatoms
      ! IF(lwp) WRITE(numout,*) 
      ! IF(lwp) WRITE(numout,*) 'max ice diatoms N hemisphere : '
      ! DO jl = 1, jpl
      !    zmaxia = MAXVAL( icedia(:,:,jl), MASK= gphit(:,:) > 0._wp )
      !    CALL mpp_max( "trc_sms_csib", zmaxia )
      !    IF(lwp) WRITE(numout,*) 'ice category ', jl , ' : ' , zmaxia
      ! ENDDO
      ! IF(lwp) WRITE(numout,*) 
      ! IF(lwp) WRITE(numout,*) 'max N hemisphere : '
      ! zmaxia = MAXVAL( iceno3(:,:,1), MASK= gphit(:,:) > 0._wp )
      ! CALL mpp_max( "trc_sms_csib", zmaxia )
      ! IF(lwp) WRITE(numout,*) 'max ice no3 cat 1' , zmaxia
      ! IF(lwp) WRITE(numout,*) 


      ! Conversion from intensive/equivalent to extensive/global variables
      DO jn = 1,jp_csib
         icetra_gca(:,:,:,jn) = icetra(:,:,:,jn) * a_i(:,:,:)
      ENDDO

      IF( ln_timing )   CALL timing_stop('trc_sms_csib')
      !
   END SUBROUTINE trc_sms_csib



   SUBROUTINE ice_friction_velocity
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_friction_velocity ***
      !!
      !! ** Purpose :   compute friction velocity, for molecular diffusion at sea ice- ocean interface
      !!                code adapted from ICE/icesbc.F90/ice_flx_other
      !!
      !! ** Inputs  :   u_ice, v_ice, ssu_m, ssv_m
      !! ** Outputs :   fric_vel
      !!-----------------------------------------------------------------------
      INTEGER  ::   ji, jj             ! dummy loop indices
      REAL(wp) ::   zu_io, zv_io, zu_iom1, zv_iom1
      !!-----------------------------------------------------------------------
      DO_2D( 0, 0, 0, 0 )
         zu_io   = u_ice(ji  ,jj  ) - ssu_m(ji  ,jj  )
         zu_iom1 = u_ice(ji-1,jj  ) - ssu_m(ji-1,jj  )
         zv_io   = v_ice(ji  ,jj  ) - ssv_m(ji  ,jj  )
         zv_iom1 = v_ice(ji  ,jj-1) - ssv_m(ji  ,jj-1)
         !
         fric_vel(ji,jj) = SQRT( rn_cio * 0.5_wp * ( zu_io*zu_io + zu_iom1*zu_iom1 + zv_io*zv_io + zv_iom1*zv_iom1 ) ) * tmask(ji,jj,1)
      END_2D
   END SUBROUTINE ice_friction_velocity


   INTEGER FUNCTION trc_sms_csib_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_csib_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to CSIB
      ! ALLOCATE( tab(...) , STAT=trc_sms_csib_alloc )
      trc_sms_csib_alloc = 0      ! set to zero if no array to be allocated
      
      ALLOCATE( &
      ! variables
         &     icetra(jpi,jpj,jpl,jp_csib) , icetra_gca (jpi,jpj,jpl,jp_csib) , icetragca_2d(jpij,jpl,jp_csib) , &
         &     qnidia(jpi,jpj,jpl), qchidia(jpi,jpj,jpl),   &
      ! sources and sinks
         &     flushrate (jpi,jpj,jpl) , flush_dia  (jpi,jpj,jpl) , lamloss_dia (jpi,jpj,jpl) , &
         &     bogup     (jpi,jpj,jpl) , bogup_dia  (jpi,jpj,jpl) , lagup       (jpi,jpj,jpl) , &
         &     lagup_dia (jpi,jpj,jpl) , lagup_no3  (jpi,jpj,jpl) , lagup_nh4   (jpi,jpj,jpl) , &
         &     phot_dia  (jpi,jpj,jpl) , mortlin_dia(jpi,jpj,jpl) , mortquad_dia(jpi,jpj,jpl) , &
         &     lim_PAR   (jpi,jpj,jpl) , lim_nut    (jpi,jpj,jpl) ,                             &
         &     diaupn    (jpi,jpj,jpl) , chlsyn     (jpi,jpj,jpl) ,                             &
         &     remin_dia (jpi,jpj,jpl) , nitri      (jpi,jpj,jpl) ,                             &
         &     flush_no3 (jpi,jpj,jpl) , lamloss_no3(jpi,jpj,jpl) , moldif_no3  (jpi,jpj,jpl) , &
         &     flush_nh4 (jpi,jpj,jpl) , lamloss_nh4(jpi,jpj,jpl) , moldif_nh4  (jpi,jpj,jpl) , &
         &     bogup_no3 (jpi,jpj,jpl) , bogup_nh4  (jpi,jpj,jpl) ,                             &
         &     fric_vel  (jpi,jpj)     ,                                                        &
         &     STAT=trc_sms_csib_alloc)

      IF( trc_sms_csib_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_csib_alloc : failed to allocate arrays' )
      !
   END FUNCTION trc_sms_csib_alloc

   !!======================================================================
END MODULE trcsms_csib
