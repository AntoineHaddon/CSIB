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
   USE sbc_oce , ONLY :  tprecip, sprecip    ! total and solid precipitation
   USE par_canoe        ! indices of CanOE model variables, e.g. jrdia: diatoms

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_csib       ! called by trcsms.F90 module
   PUBLIC   trc_sms_csib_alloc ! called by trcini_csib.F90 module

   ! Defined HERE the arrays specific to CSIB sms and ALLOCATE them in trc_sms_csib_alloc
   REAL(wp), PUBLIC, PARAMETER ::   epsi10 = 1.e-10_wp  !: small number

   !! Sea ice tracers are like ice model variables and have 2 equivalent variables: 
   !!    - one extenisve for dynamics
   !!    - one intensive for thermodynanics and biogeochemistry
   !!
   !! **********************************************************************|
   !! ***         Category dependent state variables (prognostic)        ***|
   !! **********************************************************************|
   !!                                                                       |
   !! ** Global variables                                                   |
   !!-------------|-------------|---------------------------------|---------|
   !! icedia_gca  |      -      |    Ice diatonms grid cell average  | mmol/m3 |
   !!                                                                       |
   !!-------------|-------------|---------------------------------|---------|
   !!                                                                       |
   !! ** Equivalent variables                                               |
   !!-------------|-------------|---------------------------------|---------|
   !! icedia      | -           |    Ice diatonms per ice area       | mmol/m3 |


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: icedia            !  Ice diatonms per ice area
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: icedia_gca        !  Ice diatonms grid cell average
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: icediagca_2d      !  Ice diatonms grid cell average, 2d version for ice model

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flushrate        !  Flushrate per ice category (m/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flush_dia        !  Loss rate of ice diatonms from flushing per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lamloss_dia      !  Loss rate from lateral melt per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup            !  Flowrate of water uptake from bottom ice growth per ice category (m/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup_dia        !  Diatoms uptake rate from bottom ice growth per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup            !  Flowrate uptake from lateral ice growth  per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_dia        !  Diatoms uptake rate from lateral ice growth  per ice category (mmol/m3/s)

   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: growth_dia       !  Diatoms growth rate per ice category (mmol/m3/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lim_lig          !  Diatoms light limitation factor per ice category (-)

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
      ! REAL(wp) :: zscale ! scale factor between sea ice skeletal layer and ocean surface layer
      ! REAL(wp) :: zmaxia !for diagnostics/debug
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_csib')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_csib:  CSIB model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'



      ! ! Conversion from global to equivalent variables
      ! ! and set to 0 if low ice concentration 
      ! DO jl = 1, jpl ! loop ice categories
      !    DO jj = 1, jpj
      !       DO ji = 1, jpi
               
      !          IF( a_i(ji,jj,jl) > 1.e-4_wp ) THEN ! if ice
      !             icedia(ji,jj,jl) = icedia_gca(ji,jj,jl) / a_i(ji,jj,jl)

      !          ELSE ! low ice
      !             ! set tracers to 0 
      !             icedia(ji,jj,jl)=0._wp
                  
      !             ! send what was there to the ocean
      !             zscale = z_ia / e3t_0(ji,jj,1)
      !             tr(ji,jj,1,jrdia,Kmm) = tr(ji,jj,1,jrdia,Kmm) + icedia_gca(ji,jj,jl) * zscale

      !          ENDIF ! if  ice
           
      !       ENDDO ! loop jpi
      !    ENDDO ! loop jpj
      ! ENDDO ! loop jpl ice categories


      ! reset fluxes: so that they are 0 where there is no ice
      ! flushrate(:,:,:) = 0._wp
      ! flush_dia(:,:,:)  = 0._wp
      ! lamloss_dia(:,:,:)   = 0._wp
      ! bogup(:,:,:) = 0._wp
      ! bogup_dia(:,:,:) = 0._wp
      ! lagup(:,:,:) = 0._wp
      ! lagup_dia(:,:,:) = 0._wp
      ! reset BGC processesd
      growth_dia(:,:,:)=0._wp
      lim_lig(:,:,:)=0._wp

      DO jl = 1, jpl
         DO jj = 1, jpj
            DO ji = 1, jpi
         
               IF( a_i(ji,jj,jl) > epsi10 ) THEN ! precence of ice
         
                  ! ! Flushing of ice tracers from ice-ocean exchanges

                  ! ! flushrate: water flowrate per ice area (m/s) 
                  ! flushrate(ji,jj,jl) = &
                  !          ! change of ice thickness from bottom+surface melt (m/s) * ice density (kg/m3) / freshwater density (kg/m3)
                  !    &     ( dh_bom_cat(ji,jj,jl) + dh_sum_cat(ji,jj,jl) ) * rhoi / rhow    &
                  !          ! change of snow thickness from surface melt (m/s) * snow density (kg/m3) / freshwater density (kg/m3)
                  !    &     + dh_snw_sum_cat(ji,jj,jl) * rhos / rhow     &
                  !          ! flowrate of melt pond drainage volume per ice area (m/s)
                  !    &     + dh_mpdrn_cat(ji,jj,jl)      &
                  !          ! SIC * (total precipation (Kg/m2/s) - solid precipitation (Kg/m2/s) ) / freshwater density (kg/m3)
                  !  ! &     + a_i(ji,jj,jl) * MAX( 0._wp, tprecip(ji,jj) - sprecip(ji,jj) ) / rhow
                  !    &     + a_i(ji,jj,jl) * tprecip(ji,jj) / rhow ! dev run : tprecip seems to be only rain

                  ! ! flushing of ice diatoms: flushrate * ice diatoms concentration/ height of skeletal layer  = biomass flux (mmol/m3/s)
                  ! flush_dia(ji,jj,jl) = flushrate(ji,jj,jl)/z_ia  * icedia(ji,jj,jl)     

                  ! ! loss of ice diatoms from lateral melt : fraction of ice concentration lost (1/s) * ice diatoms concentration (mmol/m3)
                  ! lamloss_dia(ji,jj,jl) = da_lam_cat(ji,jj,jl) * icedia(ji,jj,jl)


                  ! ! Uptake of ice tracers from ice growth

                  ! ! bogup: water uptake flowrate per ice area from bottom ice growth (m/s) 
                  ! ! = rate of ice thickness change (m/s) * ice density (kg/m3) / freshwater density (kg/m3)
                  ! bogup(ji,jj,jl) = dh_bog_cat(ji,jj,jl) * rhoi / rhow

                  ! ! Uptake of ice diatoms from bottom ice growth : flowrate per ice area (m/s) * ocean surface concentration (mmol/m3) /skeletal layer (m)
                  ! bogup_dia(ji,jj,jl) = bogup(ji,jj,jl) * tr(ji,jj,1,jrdia,Kmm) /z_ia

                  ! ! lagup: water uptake flowrate per ice area from lateral ice growth (m/s) 
                  ! ! = (sic increase rate / sic) (1/s) * (height new ice) (s) * ice density (kg/m3) / freshwater density (kg/m3)
                  ! lagup(ji,jj,jl) = da_lag_cat(ji,jj,jl) * ht_i_new(ji,jl) * rhoi / rhow

                  ! ! Uptake of ice diatoms from lateral ice growth  
                  ! lagup_dia(ji,jj,jl) = &
                  !       ! water uptake flowrate per ice area * skeletal layer * ocean surface diatoms concentration     
                  !       & lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrdia,Kmm)     &
                  !       ! - ice tracer conc * (sic increase rate / sic)
                  !       & - icedia(ji,jj,jl) * da_lag_cat(ji,jj,jl)

                  
                  ! Ice diatom growth

                  ! Light limitation factor, from qtr_ice_bot: shortwave radiation transmitted through ice (W/m2)
                  lim_lig(ji,jj,jl)  = tanh( r_pp * qtr_ice_bot(ji,jj,jl) ) 

                  ! Ice diatom growth rate 
                  growth_dia(ji,jj,jl) = lim_lig(ji,jj,jl) * icedia(ji,jj,jl)


                  ! Ice diatoms dynamics
                  icedia(ji,jj,jl) = icedia(ji,jj,jl) + rDt_trc * (         &
                              !          ! sink: flushing from bottom and surface ice melt, snow melt, rain on ice, melt pond drainage
                              ! &        - flush_dia(ji,jj,jl)                  & 
                              !          ! sink: loss from lateral melting of ice
                              ! &        - lamloss_dia(ji,jj,jl)                & 
                              !          ! source: Uptake from bottom ice growth
                              ! &        + bogup_dia(ji,jj,jl)                  &
                              !          ! source: Uptake from lateral ice growth
                              ! &        + lagup_dia(ji,jj,jl)                  &
                                       ! source: growth
                              &        growth_dia(ji,jj,jl)                  &
                              )
                  ! guarantee positive concentration
                  icedia(ji,jj,jl) = MAX(0._wp, icedia(ji,jj,jl) )


                  ! ! Ocean surface phytoplankton seeding and removal
                  
                  ! ! scaling factor: conversion of flowrate per sea ice area to flowrate per unit volume
                  ! ! = sea ice concentration / ocean surface layer height
                  ! ! = sea ice area / (cell area * ocean surface layer height)
                  ! zscale = a_i(ji,jj,jl) / e3t_0(ji,jj,1)
                  ! ! zscale = a_i(ji,jj,jl)  / e3t(ji,jj,1,Kmm) ! (time dependent scale factor) compilation fails: e3t  only defined if not using key_qco

                  ! ! ocean surface phytoplankton dynamics
                  ! tr(ji,jj,1,jrdia,Kmm) = tr(ji,jj,1,jrdia,Kmm) + rDt_trc * (             &
                  !          ! source: flushing of ice diatoms from bottom and surface ice melt, snow melt, rain on ice, melt pond drainage
                  ! &        + flushrate(ji,jj,jl) * zscale * icedia(ji,jj,jl)              & 
                  !          ! source: ice diatoms from lateral melting of ice
                  ! &        + da_lam_cat(ji,jj,jl) * z_ia * zscale * icedia(ji,jj,jl)      & 
                  !          ! sink: Uptake from bottom ice growth
                  ! &        - bogup(ji,jj,jl) * zscale * tr(ji,jj,1,jrdia,Kmm)             &
                  !          ! sink: Uptake from lateral ice growth
                  ! &        - lagup(ji,jj,jl) * zscale * tr(ji,jj,1,jrdia,Kmm)             &
                  ! )
                  ! ! guarantee positive concentration
                  ! tr(ji,jj,1,jrdia,Kmm) = MAX(0._wp, tr(ji,jj,1,jrdia,Kmm) )
               
               ENDIF ! if ice

            ENDDO ! loop jpi
         ENDDO ! loop jpj
      ENDDO ! loop jpl ice categories


      ! For debug/diagnostics: print max ice diatoms
      ! IF(lwp) WRITE(numout,*) 
      ! IF(lwp) WRITE(numout,*) 'max ice diatoms N hemisphere : '
      ! DO jl = 1, jpl
      !    zmaxia = MAXVAL( icedia(:,:,jl), MASK= gphit(:,:) > 0._wp )
      !    CALL mpp_max( "trc_sms_csib", zmaxia )
      !    IF(lwp) WRITE(numout,*) 'ice category ', jl , ' : ' , zmaxia
      ! ENDDO



      ! Conversion from intensive/equivalent to extensive/global variables
      ! icedia_gca(:,:,:) = icedia(:,:,:) * a_i(:,:,:)


      IF( ln_timing )   CALL timing_stop('trc_sms_csib')
      !
   END SUBROUTINE trc_sms_csib


   INTEGER FUNCTION trc_sms_csib_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_csib_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to CSIB
      ! ALLOCATE( tab(...) , STAT=trc_sms_csib_alloc )
      trc_sms_csib_alloc = 0      ! set to zero if no array to be allocated
      
      ALLOCATE(icedia   (jpi,jpj,jpl) , icedia_gca(jpi,jpj,jpl) , icediagca_2d(jpij,jpl)    , &
         &     flushrate(jpi,jpj,jpl) , flush_dia (jpi,jpj,jpl) , lamloss_dia (jpi,jpj,jpl) , &
         &     bogup    (jpi,jpj,jpl) , bogup_dia (jpi,jpj,jpl) , lagup       (jpi,jpj,jpl) , &
         &     lagup_dia(jpi,jpj,jpl) ,                                                       &
         &     lim_lig  (jpi,jpj,jpl) , growth_dia(jpi,jpj,jpl) ,                             &
         &     STAT=trc_sms_csib_alloc)

      IF( trc_sms_csib_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_csib_alloc : failed to allocate arrays' )
      !
   END FUNCTION trc_sms_csib_alloc

   !!======================================================================
END MODULE trcsms_csib
