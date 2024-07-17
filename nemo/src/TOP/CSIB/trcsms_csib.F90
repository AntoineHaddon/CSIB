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
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_dia        !  Diatoms uptake rate from lateral ice growth  per ice category (mmol/m3/s)


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

      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_csib')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_csib:  CSIB model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'



      ! Conversion from global to equivalent variables
      WHERE( a_i(:,:,:) > epsi10 )
         icedia(:,:,:) = icedia_gca(:,:,:) / a_i(:,:,:)
      ELSEWHERE
         icedia(:,:,:)=0._wp
      END WHERE


      

      DO jl = 1, jpl
         DO ji = 1, jpi
            DO jj = 1, jpj
         
               IF( a_i(ji,jj,jl) > epsi10 ) THEN ! precence of ice
         
                  ! Flushing of ice tracers from ice-ocean exchanges
                  ! flushrate: water flowrate per ice area (m/s) 
                  flushrate(ji,jj,jl) = &
                           ! change of ice thickness from bottom+surface melt (m/s) * ice density (kg/m3) / freshwater density (kg/m3)
                     &     ( dh_bom_cat(ji,jj,jl) + dh_sum_cat(ji,jj,jl) ) * rhoi / rhow    &
                           ! change of snow thickness from surface melt (m/s) * snow density (kg/m3) / freshwater density (kg/m3)
                     &     + dh_snw_sum_cat(ji,jj,jl) * rhos / rhow     &
                           ! SIC * (total precipation (Kg/m2/s) - solid precipitation (Kg/m2/s) ) / freshwater density (kg/m3)
                     &     + a_i(ji,jj,jl) * MAX( 0._wp, tprecip(ji,jj) - sprecip(ji,jj) ) / rhow     &
                           ! flowrate of melt pond drainage volume per ice area (m/s)
                     &     + dh_mpdrn_cat(ji,jj,jl)

                  ! flushing of ice algea: flushrate * ice algae concentration/ height of skeletal layer  = biomass flux (mmol/m3/s)
                  flush_dia(ji,jj,jl) = flushrate(ji,jj,jl)/z_ia  * icedia(ji,jj,jl)     

                  ! loss of ice algae from lateral melt : fraction of ice area lost (1/s) * ice algae concentration (mmol/m3)
                  lamloss_dia(ji,jj,jl) = da_lam_cat(ji,jj,jl) * icedia(ji,jj,jl)


                  ! Uptake of ice algae from ice growth
                  ! bogup: water uptake flowrate per ice area from bottom ice growth (m/s) 
                  ! = rate of ice thickness change (m/s) * ice density (kg/m3) / freshwater density (kg/m3)
                  bogup(ji,jj,jl) = dh_bog_cat(ji,jj,jl) * rhoi / rhow

                  ! Uptake of ice algae from bottom ice growth : flowrate per ice area (m/s) * ocean surface concentration (mmol/m3) /skeletal layer (m)
                  bogup_dia(ji,jj,jl) = bogup(ji,jj,jl) * tr(ji,jj,1,jrdia,Kmm) /z_ia


                  ! ice algae dynamics
                  icedia(ji,jj,jl) = icedia(ji,jj,jl) + rDt_trc * (           &
                                       ! sink: flushing from bottom and surface ice melt, snow melt, rain on ice, melt pond drainage
                              ! &        - flush_dia(ji,jj,jl)                  & 
                                       ! sink: loss from lateral melting of ice
                              ! &        - lamloss_dia(ji,jj,jl)                & 
                                       ! source: Uptake from bottom ice growth
                              &        + bogup_dia(ji,jj,jl)                  &
                              )

               ELSE ! no ice
                  flushrate(ji,jj,jl) = 0._wp
                  flush_dia(ji,jj,jl)  = 0._wp
                  lamloss_dia(ji,jj,jl)   = 0._wp
                  bogup(ji,jj,jl) = 0._wp
                  bogup_dia(ji,jj,jl) = 0._wp
                  lagup_dia(ji,jj,jl) = 0._wp
               ENDIF

            ENDDO ! loop jpj
         ENDDO ! loop jpi
      ENDDO ! loop jpl ice categories





      ! Conversion from equivalent to global variables
      icedia_gca(:,:,:) = icedia(:,:,:) * a_i(:,:,:)


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
         &     bogup    (jpi,jpj,jpl) , bogup_dia (jpi,jpj,jpl) , lagup_dia   (jpi,jpj,jpl) , &
         &     STAT=trc_sms_csib_alloc)

      IF( trc_sms_csib_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_csib_alloc : failed to allocate arrays' )
      !
   END FUNCTION trc_sms_csib_alloc

   !!======================================================================
END MODULE trcsms_csib
