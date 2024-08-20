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
   USE sbc_oce , ONLY :  tprecip, sprecip    ! total and solid precipitation
   USE par_canoe        ! indices of CanOE model variables, e.g. jrdia: diatoms

   USE lbclnk         ! lateral boundary conditions (or mpp links)

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
   !! *************************************************************************|
   !! ***         Category dependent state variables (prognostic)           ***|
   !! *************************************************************************|
   !!                                                                          |
   !! ** Global variables                                                      |
   !!-------------|-------------|---------------------------------|------------|
   !! icedia_gca  |      -      | Ice diatoms grid cell average   | mmol C /m3 |
   !! iceno3_gca  |      -      | Ice NO3 grid cell average       | mmol N /m3 |
   !! icenh4_gca  |      -      | Ice NH4 grid cell average       | mmol N /m3 |
   !!                                                                          |
   !!-------------|-------------|---------------------------------|------------|
   !!                                                                          | 
   !! ** Equivalent variables                                                  |
   !!-------------|-------------|---------------------------------|------------|
   !! icedia      | -           | Ice diatoms per ice area        | mmol C /m3 |
   !! iceno3      | -           | Ice NO3 per ice area            | mmol N /m3 |
   !! icenh4      | -           | Ice NH4 per ice area            | mmol N /m3 |
   !!                                                                          |
   !!-------------|-------------|---------------------------------|------------|
   !!

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: icedia            !  Ice diatoms per ice area
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: icedia_gca        !  Ice diatoms grid cell average
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: icediagca_2d      !  Ice diatoms grid cell average, 2d version for ice model

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: iceno3            !  Ice NO3 per ice area
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: iceno3_gca        !  Ice NO3 grid cell average
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: iceno3gca_2d      !  Ice NO3 grid cell average, 2d version for ice model

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: icenh4            !  Ice NH4 per ice area
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: icenh4_gca        !  Ice NH4 grid cell average
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: icenh4gca_2d      !  Ice NH4 grid cell average, 2d version for ice model

   !!
   !! Process
   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flushrate        !  Flushrate per ice category (m/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup            !  Flowrate of water uptake from bottom ice growth per ice category (m/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup            !  Flowrate of water uptake per ice area from lateral ice growth (m/s)
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flush_dia        !  Loss rate of ice diatoms from flushing per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lamloss_dia      !  Loss rate from lateral melt per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup_dia        !  Diatoms uptake rate from bottom ice growth per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_dia        !  Diatoms uptake rate from lateral ice growth  per ice category (mmol/m3/s)
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: growth_dia       !  Diatoms growth rate per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lim_lig          !  Diatoms light limitation factor per ice category (-)

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flush_no3        !  Loss rate of ice no3 from flushing per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lamloss_no3      !  Loss rate from lateral melt per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flush_nh4        !  Loss rate of ice nh4 from flushing per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lamloss_nh4      !  Loss rate from lateral melt per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_no3        !  NO3 uptake rate from lateral ice growth  per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_nh4        !  NH4 uptake rate from lateral ice growth  per ice category (mmol/m3/s)

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fric_vel         ! Ice frictional velocity (m/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: moldif_no3       ! Molecular diffusion ratea at ice ocean interface for NO3 per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: moldif_nh4       ! Molecular diffusion ratea at ice ocean interface for NH4  per ice category (mmol/m3/s)


   ! model parameters
   REAL(wp), SAVE ::   z_ia = 0.03_wp                    ! height of skeletal layer
   REAL(wp), SAVE ::   r_pp = 2.0_wp                     ! ratio of photosynthetic parameters (W m-2)-1
   REAL(wp), SAVE ::   mu_max = 0.85_wp / 86400._wp      ! Maximum specific growth rate  (d)-1
   ! REAL(wp), SAVE ::   c_di = 4.7e-8_wp                  ! Molecular diffusion coefficient for dissolved nutrients at the ice-water interface (m/s2)
   REAL(wp), SAVE ::   c_di = 4.7e-11_wp                  ! Molecular diffusion coefficient for dissolved nutrients at the ice-water interface (m/s2) [reduced for large time steps, otherwise numerical problems]
   REAL(wp), SAVE ::   c_nu = 1.85e-6_wp                 ! Kinematic viscosity of seawater (m2/s)

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
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      
      INTEGER ::   ji,jj,jl   ! dummy loop index
      REAL(wp) :: zscale ! scale factor between sea ice skeletal layer and ocean surface layer
      REAL(wp) :: zmax !for diagnostics/debug
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_csib')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_csib:  CSIB model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      ! IF( .NOT. ln_rsttr .AND. kt==nit000 ) THEN
      !    ! initiate here nitrogen, because during spin up  big changes can occur during first timestep in ocean N, then molecular diffusion becones very big and causes numerical instabilities
      !    iceno3(:,:,:)=0._wp
      !    icenh4(:,:,:)=0._wp
         
      !    ! init from ocean surface concentration
      !    DO jl = 1, jpl
      !       WHERE( a_i(:,:,jl) > epsi10 )
      !          iceno3(:,:,jl) = tr(:,:,1,jqno3,Kmm)
      !          icenh4(:,:,jl) = tr(:,:,1,jrnh4,Kmm)
      !       END WHERE
      !    ENDDO

      !    iceno3_gca(:,:,:) = iceno3(:,:,:) * a_i(:,:,:)
      !    icenh4_gca(:,:,:) = icenh4(:,:,:) * a_i(:,:,:)
         
         
      !    ! For debug/diagnostics: print max 
      !    IF(lwp) WRITE(numout,*) 
      !    IF(lwp) WRITE(numout,*) 'init, max N hemisphere : '
         
      !    zmax = MAXVAL( iceno3(:,:,1), MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,1) > epsi10 )
      !    CALL mpp_max( "trc_sms_csib", zmax )
      !    IF(lwp) WRITE(numout,*) 'ice category ', 1 , ' no3_i : ' , zmax
         
      !    zmax = MAXVAL( tr(:,:,1,jqno3,Kmm) - iceno3(:,:,1), MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,1) > epsi10 )
      !    CALL mpp_max( "trc_sms_csib", zmax )
      !    IF(lwp) WRITE(numout,*) 'no3_o(Kmm) - no3_i : ' , zmax

      !    zmax = MAXVAL( tr(:,:,1,jqno3,Kmm) , MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,1) > epsi10 )
      !    CALL mpp_max( "trc_sms_csib", zmax )
      !    IF(lwp) WRITE(numout,*) 'no3_o(Kmm): ' , zmax

      ! ENDIF

      ! Conversion from global to equivalent variables
      ! and set to 0 if low ice concentration 
      DO jl = 1, jpl ! loop ice categories
         DO jj = 1, jpj
            DO ji = 1, jpi
               
               IF( a_i(ji,jj,jl) > 1e-4 ) THEN ! if ice
                  icedia(ji,jj,jl) = icedia_gca(ji,jj,jl) / a_i(ji,jj,jl)
                  iceno3(ji,jj,jl) = iceno3_gca(ji,jj,jl) / a_i(ji,jj,jl)
                  icenh4(ji,jj,jl) = icenh4_gca(ji,jj,jl) / a_i(ji,jj,jl)

               ELSE ! very low ice
                  ! set tracers to 0 
                  icedia(ji,jj,jl)=0._wp
                  iceno3(ji,jj,jl)=0._wp
                  icenh4(ji,jj,jl)=0._wp
                  
                  ! send what was in the ice to the ocean
                  zscale = z_ia / e3t_0(ji,jj,1)
                  tr(ji,jj,1,jrdia,Kmm) = tr(ji,jj,1,jrdia,Kmm) + icedia_gca(ji,jj,jl) * zscale
                  ! tr(ji,jj,1,jqno3,Kmm) = tr(ji,jj,1,jqno3,Kmm) + iceno3_gca(ji,jj,jl) * zscale
                  ! tr(ji,jj,1,jrnh4,Kmm) = tr(ji,jj,1,jrnh4,Kmm) + icenh4_gca(ji,jj,jl) * zscale

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
      
      growth_dia(:,:,:) = 0._wp
      lim_lig(:,:,:) = 0._wp

      flush_no3(:,:,:) = 0._wp
      lamloss_no3(:,:,:) = 0._wp
      moldif_no3(:,:,:) = 0._wp
      lagup_no3(:,:,:) = 0._wp

      flush_nh4(:,:,:) = 0._wp
      lamloss_nh4(:,:,:) = 0._wp
      moldif_nh4(:,:,:) = 0._wp
      lagup_nh4(:,:,:) = 0._wp

      ! Compute friction velocity, for molecular diffusion at sea ice ocean interface
      CALL ice_friction_velocity


      ! For debug/diagnostics: print max 
      ! IF(lwp) WRITE(numout,*) 
      ! IF(lwp) WRITE(numout,*) 'before computations, max N hemisphere, time step number : ' , kt
      ! IF(lwp) WRITE(numout,*) 'tracer time step : ' , rDt_trc
      ! IF(lwp) WRITE(numout,*) 'ice category ', 1

      ! zmax = MAXVAL( iceno3(:,:,1), MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,1) > 1e-4 )
      ! CALL mpp_max( "trc_sms_csib", zmax )
      ! IF(lwp) WRITE(numout,*) '   no3_i : ' , zmax

      ! zmax = MAXVAL( tr(:,:,1,jqno3,Kmm) - iceno3(:,:,1), MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,1) > 1e-4 )
      ! CALL mpp_max( "trc_sms_csib", zmax )
      ! IF(lwp) WRITE(numout,*) '   no3_o(Kmm) - no3_i : ' , zmax
      
      ! zmax = MAXVAL( tr(:,:,1,jqno3,Kmm) , MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,1) > 1e-4 )
      ! CALL mpp_max( "trc_sms_csib", zmax )
      ! IF(lwp) WRITE(numout,*) 'no3_o(Kmm): ' , zmax


      
      DO jl = 1, jpl
         DO jj = 1, jpj
            DO ji = 1, jpi
         
               IF( a_i(ji,jj,jl) > epsi10 ) THEN ! precence of ice
         
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

                  ! flushing of ice tracers: flushrate * ice tracers concentration/ height of skeletal layer  = mass flux (mmol/m3/s)
                  flush_dia(ji,jj,jl) = flushrate(ji,jj,jl)/z_ia  * icedia(ji,jj,jl)  ! ice diatoms   
                  flush_no3(ji,jj,jl) = flushrate(ji,jj,jl)/z_ia  * iceno3(ji,jj,jl)  ! ice no3   
                  flush_nh4(ji,jj,jl) = flushrate(ji,jj,jl)/z_ia  * icenh4(ji,jj,jl)  ! ice nh4   

                  ! loss of ice tracers from lateral melt : fraction of ice concentration lost (1/s) * ice tracers concentration (mmol/m3)
                  lamloss_dia(ji,jj,jl) = da_lam_cat(ji,jj,jl) * icedia(ji,jj,jl)   ! ice diatoms
                  lamloss_no3(ji,jj,jl) = da_lam_cat(ji,jj,jl) * iceno3(ji,jj,jl)   ! ice no3
                  lamloss_nh4(ji,jj,jl) = da_lam_cat(ji,jj,jl) * icenh4(ji,jj,jl)   ! ice nh4


            ! Uptake of ice tracers from ice-ocean exchanges

                  ! bogup: water uptake flowrate per ice area from bottom ice growth (m/s) 
                  ! = rate of ice thickness change (m/s) * ice density (kg/m3) / freshwater density (kg/m3)
                  bogup(ji,jj,jl) = dh_bog_cat(ji,jj,jl) * rhoi / rhow

                  ! Uptake of ice diatoms from bottom ice growth : flowrate per ice area (m/s) * ocean surface concentration (mmol/m3) /skeletal layer (m)
                  bogup_dia(ji,jj,jl) = bogup(ji,jj,jl) * tr(ji,jj,1,jrdia,Kmm) /z_ia

                  ! lagup: water uptake flowrate per ice area from lateral ice growth (m/s) 
                  ! = (sic increase rate / sic) (1/s) * (height new ice) (s) * ice density (kg/m3) / freshwater density (kg/m3)
                  lagup(ji,jj,jl) = da_lag_cat(ji,jj,jl) * ht_i_new(ji,jl) * rhoi / rhow

                  ! Uptake of ice diatoms from lateral ice growth  
                  lagup_dia(ji,jj,jl) = &
                              ! water uptake flowrate per ice area * skeletal layer * ocean surface diatoms concentration     
                        &     lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrdia,Kmm)     &
                              ! - ice tracer conc * (sic increase rate / sic)
                        &     - icedia(ji,jj,jl) * da_lag_cat(ji,jj,jl)
                  
                  ! Uptake of ice N from lateral ice growth  
                  lagup_no3(ji,jj,jl) = lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jqno3,Kmm) - iceno3(ji,jj,jl) * da_lag_cat(ji,jj,jl)
                  lagup_nh4(ji,jj,jl) = lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrnh4,Kmm) - icenh4(ji,jj,jl) * da_lag_cat(ji,jj,jl)
                  
                  
                  ! Diffusion of N at ice ocean interface
                  ! = D / (nu / |friction velocity|) * ( N_ocean - N_ice ) / skeletal layer
                  IF( a_i(ji,jj,jl) > 1e-4 ) THEN
                     moldif_no3(ji,jj,jl) = 0._wp !c_di / c_nu * abs(fric_vel(ji,jj)) * ( tr(ji,jj,1,jqno3,Kmm) - iceno3(ji,jj,jl) ) /z_ia
                     moldif_nh4(ji,jj,jl) = 0._wp !c_di / c_nu * abs(fric_vel(ji,jj)) * ( tr(ji,jj,1,jrnh4,Kmm) - icenh4(ji,jj,jl) ) /z_ia
                     ! limit molecular diffusion to +/- 10%  per timestep
                     ! moldif_no3(ji,jj,jl) = MAX( MIN( moldif_no3(ji,jj,jl), iceno3(ji,jj,jl)*0.1_wp/rDt_trc ), -iceno3(ji,jj,jl)*0.1_wp/rDt_trc)
                     ! moldif_nh4(ji,jj,jl) = MAX( MIN( moldif_nh4(ji,jj,jl), icenh4(ji,jj,jl)*0.1_wp/rDt_trc ), -icenh4(ji,jj,jl)*0.1_wp/rDt_trc)
                  ENDIF 

                  ! scaling factor for conversion of ice-ocean exchanges rates for ocean side 
                  ! conversion of flowrate per sea ice area to flowrate per unit volume
                  ! = sea ice area / (cell area * ocean surface layer height)
                  ! = sea ice concentration / ocean surface layer height
                  zscale = a_i(ji,jj,jl) / e3t_0(ji,jj,1)
                  ! zscale = a_i(ji,jj,jl)  / e3t(ji,jj,1,Kmm) ! (time dependent scale factor) compilation fails: e3t  only defined if not using key_qco


            ! Ice diatom growth

                  ! Light limitation factor, from qtr_ice_bot: shortwave radiation transmitted through ice (W/m2)
                  lim_lig(ji,jj,jl) = tanh( r_pp * qtr_ice_bot(ji,jj,jl) ) 

                  ! Ice diatom growth rate 
                  growth_dia(ji,jj,jl) = mu_max * lim_lig(ji,jj,jl) * icedia(ji,jj,jl)


            ! Ice diatoms dynamics
                  icedia(ji,jj,jl) = icedia(ji,jj,jl) + rDt_trc * (         &
                                       ! sink: flushing from bottom and surface ice melt, snow melt, rain on ice, melt pond drainage
                              &        - flush_dia(ji,jj,jl)                  & 
                                       ! sink: loss from lateral melting of ice
                              &        - lamloss_dia(ji,jj,jl)                & 
                                       ! source: Uptake from bottom ice growth
                              &        + bogup_dia(ji,jj,jl)                  &
                                       ! source: Uptake from lateral ice growth
                              &        + lagup_dia(ji,jj,jl)                  &
                                       ! source: growth
                              &        + growth_dia(ji,jj,jl)                  &
                              )
                  ! guarantee positive concentration
                  icedia(ji,jj,jl) = MAX(0._wp, icedia(ji,jj,jl) )
                  icedia(ji,jj,jl) = MIN(1000._wp, icedia(ji,jj,jl) ) 
                  

            ! Ocean surface phytoplankton dynamics
                  tr(ji,jj,1,jrdia,Kmm) = tr(ji,jj,1,jrdia,Kmm) + rDt_trc * (             &
                           ! source: flushing of ice diatoms from bottom and surface ice melt, snow melt, rain on ice, melt pond drainage
                  &        + flushrate(ji,jj,jl) * zscale * icedia(ji,jj,jl)              & 
                           ! source: ice diatoms from lateral melting of ice
                  &        + da_lam_cat(ji,jj,jl) * z_ia * zscale * icedia(ji,jj,jl)      & 
                           ! sink: Uptake from bottom ice growth
                  &        - bogup(ji,jj,jl) * zscale * tr(ji,jj,1,jrdia,Kmm)             &
                           ! sink: Uptake from lateral ice growth
                  &        - lagup(ji,jj,jl) * zscale * tr(ji,jj,1,jrdia,Kmm)             &
                  )
                  ! guarantee positive concentration
                  tr(ji,jj,1,jrdia,Kmm) = MAX(0._wp, tr(ji,jj,1,jrdia,Kmm) )


            ! Ice NO3 dynamics
                  ! iceno3(ji,jj,jl) = iceno3(ji,jj,jl) + rDt_trc * (           &
                                       ! sink: flushing from bottom and surface ice melt, snow melt, rain on ice, melt pond drainage
                              ! &        - flush_no3(ji,jj,jl)                  & 
                                       ! sink: loss from lateral melting of ice
                              ! &        - lamloss_no3(ji,jj,jl)                &           
                                       ! source: Uptake from lateral ice growth
                              ! &        + lagup_no3(ji,jj,jl)                  &
                                       ! source/sink: diffusion at ocean interface
                              ! &        + moldif_no3(ji,jj,jl)                 &
                              ! )
                  ! guarantee positive concentration
                  iceno3(ji,jj,jl) = MAX(0._wp, iceno3(ji,jj,jl) )
                  iceno3(ji,jj,jl) = MIN(1000._wp, iceno3(ji,jj,jl) )

 
            ! Ocean surface NO3 dynamics
                  ! tr(ji,jj,1,jqno3,Kmm) = tr(ji,jj,1,jqno3,Kmm) + rDt_trc * (             &
                           ! source: flushing of ice NO3 from bottom and surface ice melt, snow melt, rain on ice, melt pond drainage
                  ! &        + flushrate(ji,jj,jl) * zscale * iceno3(ji,jj,jl)              & 
                           ! source: ice NO3 from lateral melting of ice
                  ! &        + da_lam_cat(ji,jj,jl) * z_ia * zscale * iceno3(ji,jj,jl)      &
                           ! sink: Uptake from lateral ice growth
                  ! &        - lagup(ji,jj,jl) * zscale * tr(ji,jj,1,jqno3,Kmm)             &
                           ! sink/source: diffusion at ocean interface
                  ! &        - moldif_no3(ji,jj,jl) * z_ia * zscale                         &
                  ! )
                  ! guarantee positive concentration
                  ! tr(ji,jj,1,jqno3,Kmm) = MAX(0._wp, tr(ji,jj,1,jqno3,Kmm) )
                  ! tr(ji,jj,1,jqno3,Kmm) = MIN(1000._wp, tr(ji,jj,1,jqno3,Kmm) )


            ! Ice NH4 dynamics
                  icenh4(ji,jj,jl) = icenh4(ji,jj,jl) + rDt_trc * (           &
                                       ! sink: flushing from bottom and surface ice melt, snow melt, rain on ice, melt pond drainage
                              &        - flush_nh4(ji,jj,jl)                  & 
                                       ! sink: loss from lateral melting of ice
                              &        - lamloss_nh4(ji,jj,jl)                &
                                       ! source: Uptake from lateral ice growth
                              &        + lagup_nh4(ji,jj,jl)                  & 
                                       ! source/sink: diffusion at ocean interface
                              &        + moldif_nh4(ji,jj,jl)                 &
                              )
                  ! guarantee positive concentration
                  icenh4(ji,jj,jl) = MAX(0._wp, icenh4(ji,jj,jl) )
                  icenh4(ji,jj,jl) = MIN(1000._wp, icenh4(ji,jj,jl) )
               

            ! Ocean surface NH4 dynamics
                  tr(ji,jj,1,jrnh4,Kmm) = tr(ji,jj,1,jrnh4,Kmm) + rDt_trc * (             &
                           ! source: flushing of ice NO3 from bottom and surface ice melt, snow melt, rain on ice, melt pond drainage
                  &        + flushrate(ji,jj,jl) * zscale * icenh4(ji,jj,jl)              & 
                           ! source: ice NO3 from lateral melting of ice
                  &        + da_lam_cat(ji,jj,jl) * z_ia * zscale * icenh4(ji,jj,jl)      & 
                           ! sink: Uptake from lateral ice growth
                  &        - lagup(ji,jj,jl) * zscale * tr(ji,jj,1,jrnh4,Kmm)             &
                           ! sink/source: diffusion at ocean interface
                  &        - moldif_nh4(ji,jj,jl) * z_ia * zscale                         &
                  )
                  ! guarantee positive concentration
                  tr(ji,jj,1,jrnh4,Kmm) = MAX(0._wp, tr(ji,jj,1,jrnh4,Kmm) )
                  tr(ji,jj,1,jrnh4,Kmm) = MAX(1000._wp, tr(ji,jj,1,jrnh4,Kmm) )


               ENDIF ! if ice

            ENDDO ! loop jpi
         ENDDO ! loop jpj
      ENDDO ! loop jpl ice categories

      ! lateral boundary conditions (or mpp links) needed??
      ! CALL lbc_lnk( 'trc_sms_csib', icedia , 'T',  1._wp)
      ! CALL lbc_lnk( 'trc_sms_csib', iceno3 , 'T',  1._wp)
      ! CALL lbc_lnk( 'trc_sms_csib', icenh4 , 'T',  1._wp)

      ! For debug/diagnostics: print max 
      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) 'after computations, max N hemisphere : '
      ! DO jl = 1, jpl
         ! IF(lwp) WRITE(numout,*) 'ice category ', jl 
         ! zmax = MAXVAL( icedia(:,:,jl), MASK= gphit(:,:) > 0._wp )
         ! CALL mpp_max( "trc_sms_csib", zmax )
         ! IF(lwp) WRITE(numout,*) ' icedia : ' , zmax
      ! ENDDO

      ! zmax = MAXVAL( iceno3(:,:,1), MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,1) > 1e-4 )
      ! CALL mpp_max( "trc_sms_csib", zmax )
      ! IF(lwp) WRITE(numout,*) '   no3 : ' , zmax

      ! zmax = MAXVAL( moldif_no3(:,:,1), MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,1) > 1e-4 )
      ! CALL mpp_max( "trc_sms_csib", zmax )
      ! IF(lwp) WRITE(numout,*) '   moldif_no3 : ' , zmax
      
      ! zmax = MAXVAL( tr(:,:,1,jqno3,Kmm) - iceno3(:,:,1), MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,1) > 1e-4 )
      ! CALL mpp_max( "trc_sms_csib", zmax )
      ! IF(lwp) WRITE(numout,*) '   no3_o(Kmm) - no3_i : ' , zmax
      
      ! zmax = MAXVAL( tr(:,:,1,jqno3,Kmm) , MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,1) > 1e-4 )
      ! CALL mpp_max( "trc_sms_csib", zmax )
      ! IF(lwp) WRITE(numout,*) 'no3_o(Kmm): ' , zmax

      ! zmax = MAXVAL( fric_vel(:,:), MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,1) > 1e-4 )
      ! CALL mpp_max( "trc_sms_csib", zmax )
      ! IF(lwp) WRITE(numout,*) 'fric_vel : ' , zmax


      ! Conversion from intensive/equivalent to extensive/global variables
      icedia_gca(:,:,:) = icedia(:,:,:) * a_i(:,:,:)
      iceno3_gca(:,:,:) = iceno3(:,:,:) * a_i(:,:,:)
      icenh4_gca(:,:,:) = icenh4(:,:,:) * a_i(:,:,:)


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
         IF( at_i(ji,jj) > epsi10 ) THEN ! precence of ice
            zu_io   = u_ice(ji  ,jj  ) - ssu_m(ji  ,jj  )
            zv_io   = v_ice(ji  ,jj  ) - ssv_m(ji  ,jj  )
            ! zu_iom1 = u_ice(ji-1,jj  ) - ssu_m(ji-1,jj  )
            ! zv_iom1 = v_ice(ji  ,jj-1) - ssv_m(ji  ,jj-1)
            !
            ! fric_vel(ji,jj) = SQRT( rn_cio * 0.5_wp * ( zu_io*zu_io + zu_iom1*zu_iom1 + zv_io*zv_io + zv_iom1*zv_iom1 ) ) 
            fric_vel(ji,jj) = SQRT( rn_cio * 0.5_wp * ( zu_io*zu_io + zv_io*zv_io ) ) 
         ELSE ! no ice
            fric_vel(ji,jj) = 0.0_wp
         ENDIF
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
      
      ALLOCATE(icedia   (jpi,jpj,jpl) , icedia_gca (jpi,jpj,jpl) , icediagca_2d(jpij,jpl)    , &
         &     iceno3   (jpi,jpj,jpl) , iceno3_gca (jpi,jpj,jpl) , iceno3gca_2d(jpij,jpl)    , &
         &     icenh4   (jpi,jpj,jpl) , icenh4_gca (jpi,jpj,jpl) , icenh4gca_2d(jpij,jpl)    , &
         &     flushrate(jpi,jpj,jpl) , flush_dia  (jpi,jpj,jpl) , lamloss_dia (jpi,jpj,jpl) , &
         &     bogup    (jpi,jpj,jpl) , bogup_dia  (jpi,jpj,jpl) , lagup       (jpi,jpj,jpl) , &
         &     lagup_dia(jpi,jpj,jpl) , lagup_no3  (jpi,jpj,jpl) , lagup_nh4   (jpi,jpj,jpl) , &
         &     lim_lig  (jpi,jpj,jpl) , growth_dia (jpi,jpj,jpl) ,                             &
         &     flush_no3(jpi,jpj,jpl) , lamloss_no3(jpi,jpj,jpl) , moldif_no3(jpi,jpj,jpl)   , &
         &     flush_nh4(jpi,jpj,jpl) , lamloss_nh4(jpi,jpj,jpl) , moldif_nh4(jpi,jpj,jpl)   , &
         &     fric_vel (jpi,jpj)     ,                                                        &
         &     STAT=trc_sms_csib_alloc)

      IF( trc_sms_csib_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_csib_alloc : failed to allocate arrays' )
      !
   END FUNCTION trc_sms_csib_alloc

   !!======================================================================
END MODULE trcsms_csib
