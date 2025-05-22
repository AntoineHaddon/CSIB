MODULE zdftmx
   !!========================================================================
   !!                       ***  MODULE  zdftmx  ***
   !! Ocean physics: vertical tidal mixing coefficient
   !!========================================================================
   !! History :  1.0  !  2004-04  (L. Bessieres, G. Madec)  Original code
   !!             -   !  2006-08  (A. Koch-Larrouy) Indonesian strait
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            4.0.3!  2021-07  (D. Yang)   add tidal mixing scheme by Simmons et al 2004
   !!            4.0.3!  2021-07  (D. Yang)   Constrain tidal energy to be positive
   !!            4.0.3!  2021-08  (D. Yang)   Vertical diffusivity resulting from internal tide breaking is now capped by 20 cm2/s
   !!            4.0.3!  2022-07  (D. Yang)   Add two flags (ln_s2004 & ln_sm2005) to make computations of zav_tide from Simmons et al (2004) and Saenko and Merryfield (2005) optionally available
   !!            4.0.3!  2022-08  (D. Yang)   Add lee wave mixing scheme controlled by ln_leewmx.
   !!            4.2  !  2024-06  (G. Stanley) Fix multiplication by lee wave mixing factor (rn_lwm)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_zdftmx'                                  Tidal vertical mixing
   !!----------------------------------------------------------------------
   !!   zdf_tmx       : global     momentum & tracer Kz with tidal induced Kz
   !!   tmx_itf       : Indonesian momentum & tracer Kz with tidal induced Kz 
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE zdf_oce        ! ocean vertical physics variables
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE eosbn2         ! ocean equation of state
   USE phycst         ! physical constants
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE iom            ! I/O Manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_tmx         ! called in step module 
   PUBLIC   zdf_tmx_init    ! called in opa module 
   PUBLIC   zdf_tmx_alloc   ! called in nemogcm module

   LOGICAL, PUBLIC, PARAMETER ::   lk_zdftmx = .TRUE.    !: tidal mixing flag

   !                       !!* Namelist  namzdf_tmx : tidal mixing *
   REAL(wp) ::  rn_htmx     ! vertical decay scale for turbulence (meters)
   REAL(wp) ::  rn_n2min    ! threshold of the Brunt-Vaisala frequency (s-1)
   REAL(wp) ::  rn_tfe      ! tidal dissipation efficiency (St Laurent et al. 2002)
   REAL(wp) ::  rn_me       ! mixing efficiency (Osborn 1980)
   LOGICAL  ::  ln_tmx_itf  ! Indonesian Through Flow (ITF): Koch-Larrouy et al. (2007) parameterization
   REAL(wp) ::  rn_tfe_itf  ! ITF tidal dissipation efficiency (St Laurent et al. 2002)
   LOGICAL  ::  ln_s2004    ! Compute zav_tide following Simmons et al (2004)
   LOGICAL  ::  ln_sm2005   ! Compute zav_tide following Saenko and Merryfield (2005)
   LOGICAL  ::  ln_leewmx   ! Add (.TRUE.) lee wave mixing or not (.FALSE.).
   REAL(wp) ::  rn_lwm      ! lee wave mixing scaling

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   en_tmx     ! energy available for tidal mixing (W/m2)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   mask_itf   ! mask to use over Indonesian area
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   az_tmx     ! coefficient used to evaluate the tidal induced Kz

   !! * Substitutions
!#  include "vectopt_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdftmx.F90 8788 2017-11-22 18:01:02Z davestorkey $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_tmx_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION zdf_tmx_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE(en_tmx(jpi,jpj), mask_itf(jpi,jpj), az_tmx(jpi,jpj,jpk), STAT=zdf_tmx_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum ('zdftmx', zdf_tmx_alloc )
      IF( zdf_tmx_alloc /= 0 )   CALL ctl_warn('zdf_tmx_alloc: failed to allocate arrays')
   END FUNCTION zdf_tmx_alloc


   SUBROUTINE zdf_tmx( kt,Kmm, p_avm, p_avt, p_avs)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_tmx  ***
      !!                   
      !! ** Purpose :   add to the vertical mixing coefficients the effect of
      !!              tidal mixing (Simmons et al 2004).
      !!
      !! ** Method  : - tidal-induced vertical mixing is given by:
      !!                  Kz_tides = az_tmx / max( rn_n2min, N^2 )
      !!              where az_tmx is a coefficient that specified the 3D space 
      !!              distribution of the faction of tidal energy taht is used
      !!              for mixing. Its expression is set in zdf_tmx_init routine,
      !!              following Simmons et al. 2004.
      !!                NB: a specific bounding procedure is performed on av_tide
      !!              so that the input tidal energy is actually almost used. The
      !!              basic maximum value is 60 cm2/s, but values of 300 cm2/s 
      !!              can be reached in area where bottom stratification is too 
      !!              weak.
      !!
      !!              - update av_tide in the Indonesian Through Flow area
      !!              following Koch-Larrouy et al. (2007) parameterisation
      !!              (see tmx_itf routine).
      !!
      !!              - update the model vertical eddy viscosity and diffusivity: 
      !!                     avt  = avt  +    av_tides
      !!                     avm  = avm  +    av_tides
      !!
      !! ** Action  :   avt, avm   increased by tidal mixing
      !!
      !! References : Simmons et al. 2004, Ocean Modelling, 6, 3-4, 245-263.
      !!              Koch-Larrouy et al. 2007, GRL.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt,Kmm   ! ocean time-step 
      REAL(wp), DIMENSION(:,:,:) , INTENT(inout) ::   p_avm          ! momentum Kz (w-points)
      REAL(wp), DIMENSION(:,:,:) , INTENT(inout) ::   p_avt, p_avs   ! tracer   Kz (w-points)
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   ztpc         ! scalar workspace
      REAL(wp), DIMENSION(jpi,jpj) ::   zkz
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zav_tide
      !!----------------------------------------------------------------------
      !
      !                          ! ----------------------- !
      !                          !  Standard tidal mixing  !  (compute zav_tide)
      !                          ! ----------------------- !
      ! --------------------------------------------------------!
      ! original computation of zav_tide by Simmons et al. 2004 !
      ! --------------------------------------------------------!
      IF( ln_s2004 ) THEN           
      !                             !* First estimation (with n2 bound by rn_n2min) bounded by 60 cm2/s
         zav_tide(:,:,:) = MIN(  60.e-4, az_tmx(:,:,:) / MAX( rn_n2min, rn2(:,:,:) )  )

         zkz(:,:) = 0.e0               !* Associated potential energy consummed over the whole water column
         DO jk = 2, jpkm1
            zkz(:,:) = zkz(:,:) + e3w(:,:,jk,Kmm) * MAX( 0.e0, rn2(:,:,jk) ) * rho0 * zav_tide(:,:,jk) * wmask(:,:,jk)
         END DO

         DO jj = 1, jpj                !* Here zkz should be equal to en_tmx ==> multiply by en_tmx/zkz to recover en_tmx
            DO ji = 1, jpi
               IF( zkz(ji,jj) /= 0.e0 )   zkz(ji,jj) = en_tmx(ji,jj) / zkz(ji,jj)
            END DO
         END DO

         DO jk = 2, jpkm1     !* Mutiply by zkz to recover en_tmx, BUT bound by 30/6 ==> zav_tide bound by 300 cm2/s
            DO jj = 1, jpj                !* Here zkz should be equal to en_tmx ==> multiply by en_tmx/zkz to recover en_tmx
               DO ji = 1, jpi
                  zav_tide(ji,jj,jk) = zav_tide(ji,jj,jk) * MIN( zkz(ji,jj), 30./6. ) * wmask(ji,jj,jk)  !kz max = 300 cm2/s
               END DO
            END DO
         END DO
      ENDIF   

      ! ------------------------------------------------------- !
      ! Compute zav_tide following Saenko and Merryfield (2005) !
      ! ------------------------------------------------------- !

      IF( ln_sm2005 ) THEN
         zav_tide(:,:,:) = MIN(  20.e-4, az_tmx(:,:,:) / MAX( rn_n2min, rn2(:,:,:) )  )
      ENDIF

      IF( kt == nit000 ) THEN       !* check at first time-step: diagnose the energy consumed by zav_tide
         ztpc = 0.e0
         DO jk= 1, jpk
            DO jj= 1, jpj
               DO ji= 1, jpi
                  ztpc = ztpc + e3w(ji,jj,jk,Kmm) * e1t(ji,jj) * e2t(ji,jj)   &
                     &         * MAX( 0.e0, rn2(ji,jj,jk) ) * zav_tide(ji,jj,jk) * tmask(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
         ztpc= rho0 / ( rn_tfe * rn_me ) * ztpc
         IF(lwp) WRITE(numout,*) 
         IF(lwp) WRITE(numout,*) '          N Total power consumption by av_tide    : ztpc = ', ztpc * 1.e-12 ,'TW'
      ENDIF
       
      !                          ! ----------------------- !
      !                          !    ITF  tidal mixing    !  (update zav_tide)
      !                          ! ----------------------- !
      IF( ln_tmx_itf )   CALL tmx_itf( kt,Kmm, zav_tide )

      !                          ! ----------------------- !
      !                          !   Update  mixing coefs  !                          
      !                          ! ----------------------- !
      DO jk = 2, jpkm1              !* update momentum & tracer diffusivity with tidal mixing
         DO jj = 1, jpj                !* Here zkz should be equal to en_tmx ==> multiply by en_tmx/zkz to recover en_tmx
            DO ji = 1, jpi
               p_avt(ji,jj,jk) = p_avt(ji,jj,jk) + zav_tide(ji,jj,jk) * wmask(ji,jj,jk)
               p_avs(ji,jj,jk) = p_avs(ji,jj,jk) + zav_tide(ji,jj,jk) * wmask(ji,jj,jk)
               p_avm(ji,jj,jk) = p_avm(ji,jj,jk) + zav_tide(ji,jj,jk) * wmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !                             !* output tidal mixing coefficient
      CALL iom_put( "av_tmx", zav_tide )

      IF(sn_cfctl%l_prtctl)   CALL prt_ctl(tab3d_1=zav_tide , clinfo1=' tmx - av_tide: ', tab3d_2=p_avt, clinfo2=' p_avt: ', kdim=jpk)
      !
   END SUBROUTINE zdf_tmx


   SUBROUTINE tmx_itf( kt,Kmm, pav )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tmx_itf  ***
      !!                   
      !! ** Purpose :   modify the vertical eddy diffusivity coefficients 
      !!              (pav) in the Indonesian Through Flow area (ITF).
      !!
      !! ** Method  : - Following Koch-Larrouy et al. (2007), in the ITF defined
      !!                by msk_itf (read in a file, see tmx_init), the tidal
      !!                mixing coefficient is computed with :
      !!                  * q=1 (i.e. all the tidal energy remains trapped in
      !!                         the area and thus is used for mixing)
      !!                  * the vertical distribution of the tifal energy is a
      !!                    proportional to N above the thermocline (d(N^2)/dz > 0)
      !!                    and to N^2 below the thermocline (d(N^2)/dz < 0)
      !!
      !! ** Action  :   av_tide   updated in the ITF area (msk_itf)
      !!
      !! References :  Koch-Larrouy et al. 2007, GRL 
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   )                         ::   kt,Kmm   ! ocean time-step
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pav  ! Tidal mixing coef.
      !! 
      INTEGER  ::   ji, jj, jk    ! dummy loop indices
      REAL(wp) ::   zcoef, ztpc   ! temporary scalar
      REAL(wp), DIMENSION(jpi,jpj)   ::   zkz                        ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj)   ::   zsum1 , zsum2 , zsum       !  -      -
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zempba_3d_1, zempba_3d_2   ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zempba_3d  , zdn2dz        !  -      -
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zavt_itf                   !  -      -
      !!----------------------------------------------------------------------
      !
      !                             ! compute the form function using N2 at each time step
      zdn2dz     (:,:,jpk) = 0.e0
      zempba_3d_1(:,:,jpk) = 0.e0
      zempba_3d_2(:,:,jpk) = 0.e0
      DO jk = 1, jpkm1             
         zdn2dz     (:,:,jk) = rn2(:,:,jk) - rn2(:,:,jk+1)           ! Vertical profile of dN2/dz
!CDIR NOVERRCHK
         zempba_3d_1(:,:,jk) = SQRT(  MAX( 0.e0, rn2(:,:,jk) )  )    !    -        -    of N
         zempba_3d_2(:,:,jk) =        MAX( 0.e0, rn2(:,:,jk) )       !    -        -    of N^2
      END DO
      !
      zsum (:,:) = 0.e0
      zsum1(:,:) = 0.e0
      zsum2(:,:) = 0.e0
      DO jk= 2, jpk
         zsum1(:,:) = zsum1(:,:) + zempba_3d_1(:,:,jk) * e3w(:,:,jk,Kmm) * tmask(:,:,jk) * tmask(:,:,jk-1)
         zsum2(:,:) = zsum2(:,:) + zempba_3d_2(:,:,jk) * e3w(:,:,jk,Kmm) * tmask(:,:,jk) * tmask(:,:,jk-1)               
      END DO
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( zsum1(ji,jj) /= 0.e0 )   zsum1(ji,jj) = 1.e0 / zsum1(ji,jj)
            IF( zsum2(ji,jj) /= 0.e0 )   zsum2(ji,jj) = 1.e0 / zsum2(ji,jj)                
         END DO
      END DO

      DO jk= 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcoef = 0.5 - SIGN( 0.5, zdn2dz(ji,jj,jk) )       ! =0 if dN2/dz > 0, =1 otherwise 
               ztpc  = zempba_3d_1(ji,jj,jk) * zsum1(ji,jj) *        zcoef     &
                  &  + zempba_3d_2(ji,jj,jk) * zsum2(ji,jj) * ( 1. - zcoef )
               !
               zempba_3d(ji,jj,jk) =               ztpc 
               zsum     (ji,jj)    = zsum(ji,jj) + ztpc * e3w(ji,jj,jk,Kmm)
            END DO
         END DO
       END DO
       DO jj = 1, jpj
          DO ji = 1, jpi
             IF( zsum(ji,jj) > 0.e0 )   zsum(ji,jj) = 1.e0 / zsum(ji,jj)                
          END DO
       END DO

      !                             ! first estimation bounded by 10 cm2/s (with n2 bounded by rn_n2min) 
      zcoef = rn_tfe_itf / ( rn_tfe * rho0 )
      DO jk = 1, jpk
         zavt_itf(:,:,jk) = MIN(  10.e-4, zcoef * en_tmx(:,:) * zsum(:,:) * zempba_3d(:,:,jk)   &
            &                                      / MAX( rn_n2min, rn2(:,:,jk) ) * tmask(:,:,jk)  )
      END DO           

      zkz(:,:) = 0.e0               ! Associated potential energy consummed over the whole water column
      DO jk = 2, jpkm1
         zkz(:,:) = zkz(:,:) + e3w(:,:,jk,Kmm) * MAX( 0.e0, rn2(:,:,jk) ) * rho0 * zavt_itf(:,:,jk) * tmask(:,:,jk) * tmask(:,:,jk-1)
      END DO

      DO jj = 1, jpj                ! Here zkz should be equal to en_tmx ==> multiply by en_tmx/zkz to recover en_tmx
         DO ji = 1, jpi
            IF( zkz(ji,jj) /= 0.e0 )   zkz(ji,jj) = en_tmx(ji,jj) * rn_tfe_itf / rn_tfe / zkz(ji,jj)
         END DO
      END DO

      DO jk = 2, jpkm1              ! Mutiply by zkz to recover en_tmx, BUT bound by 30/6 ==> zavt_itf bound by 300 cm2/s
         zavt_itf(:,:,jk) = zavt_itf(:,:,jk) * MIN( zkz(:,:), 120./10. ) * tmask(:,:,jk) * tmask(:,:,jk-1)   ! kz max = 120 cm2/s
      END DO

      IF( kt == nit000 ) THEN       ! diagnose the nergy consumed by zavt_itf
         ztpc = 0.e0
         DO jk= 1, jpk
            DO jj= 1, jpj
               DO ji= 1, jpi
                  ztpc = ztpc + e1t(ji,jj) * e2t(ji,jj) * e3w(ji,jj,jk,Kmm) * MAX( 0.e0, rn2(ji,jj,jk) )   &
                     &                     * zavt_itf(ji,jj,jk) * tmask(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
         ztpc= rho0 * ztpc / ( rn_me * rn_tfe_itf )
         IF(lwp) WRITE(numout,*) '          N Total power consumption by zavt_itf: ztpc = ', ztpc * 1.e-12 ,'TW'
      ENDIF

      !                             ! Update pav with the ITF mixing coefficient
      DO jk = 2, jpkm1
         pav(:,:,jk) = pav     (:,:,jk) * ( 1.e0 - mask_itf(:,:) )   &
            &        + zavt_itf(:,:,jk) *          mask_itf(:,:) 
      END DO
      !
   END SUBROUTINE tmx_itf


   SUBROUTINE zdf_tmx_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_tmx_init  ***
      !!                     
      !! ** Purpose :   Initialization of the vertical tidal mixing, Reading
      !!              of M2 and K1 tidal energy in nc files
      !!
      !! ** Method  : - Read the namtmx namelist and check the parameters
      !!
      !!              - Read the input data in NetCDF files :
      !!              M2 and K1 tidal energy. The total tidal energy, en_tmx, 
      !!              is the sum of M2, K1 and S2 energy where S2 is assumed 
      !!              to be: S2=(1/2)^2 * M2
      !!              mask_itf, a mask array that determine where substituing 
      !!              the standard Simmons et al. (2005) formulation with the
      !!              one of Koch_Larrouy et al. (2007).
      !!
      !!              - Compute az_tmx, a 3D coefficient that allows to compute
      !!             the standard tidal-induced vertical mixing as follows:
      !!                  Kz_tides = az_tmx / max( rn_n2min, N^2 )
      !!             with az_tmx a bottom intensified coefficient is given by:
      !!                 az_tmx(z) = en_tmx / ( rho0 * rn_htmx ) * EXP( -(H-z)/rn_htmx )
      !!                                                  / ( 1. - EXP( - H   /rn_htmx ) ) 
      !!             where rn_htmx the characteristic length scale of the bottom 
      !!             intensification, en_tmx the tidal energy, and H the ocean depth
      !!
      !! ** input   :   - Namlist namtmx
      !!                - NetCDF file : M2_ORCA2.nc, K1_ORCA2.nc, and mask_itf.nc
      !!
      !! ** Action  : - Increase by 1 the nstop flag is setting problem encounter
      !!              - defined az_tmx used to compute tidal-induced mixing
      !!
      !! References : Simmons et al. 2004, Ocean Modelling, 6, 3-4, 245-263.
      !!              Koch-Larrouy et al. 2007, GRL.
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   inum         ! local integer
      INTEGER  ::   ios
      REAL(wp) ::   ztpc, ze_z   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj)   ::  zem2, zek1, zeef ! read M2 and K1 tidal energy and mesoscale eddy energy flux
      REAL(wp), DIMENSION(jpi,jpj)   ::  zkz          ! total M2, K1 and S2 tidal energy
      REAL(wp), DIMENSION(jpi,jpj)   ::  zfact        ! used for vertical structure function
      REAL(wp), DIMENSION(jpi,jpj)   ::  zhdep        ! Ocean depth 
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zpc        ! power consumption
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zav_tide   ! tidal mixing coefficient
      !!
      NAMELIST/namzdf_tmx/ rn_htmx, rn_n2min, rn_tfe, rn_me, ln_tmx_itf, rn_tfe_itf, ln_s2004, ln_sm2005, ln_leewmx, rn_lwm
      !!----------------------------------------------------------------------
      !
      
      ! Namelist namzdf_tmx in reference namelist : Tidal Mixing
      READ  ( numnam_ref, namzdf_tmx, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_tmx in reference namelist'  )

      ! Namelist namzdf_tmx in configuration namelist : Tidal Mixing
      READ  ( numnam_cfg, namzdf_tmx, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_tmx in configuration namelist' )
      IF(lwm) WRITE ( numond, namzdf_tmx )

      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_tmx_init : tidal mixing'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namzdf_tmx : set tidal mixing parameters'
         WRITE(numout,*) '      Vertical decay scale for turbulence   = ', rn_htmx 
         WRITE(numout,*) '      Brunt-Vaisala frequency threshold     = ', rn_n2min
         WRITE(numout,*) '      Tidal dissipation efficiency          = ', rn_tfe
         WRITE(numout,*) '      Mixing efficiency                     = ', rn_me
         WRITE(numout,*) '      ITF specific parameterisation         = ', ln_tmx_itf
         WRITE(numout,*) '      ITF tidal dissipation efficiency      = ', rn_tfe_itf
         WRITE(numout,*) '     zav_tide computed following Simmons et al =', ln_s2004
         WRITE(numout,*) '     zav_tide computed following Saenko and Merryfield =', ln_sm2005
         WRITE(numout,*) '     Lee wave mixing                       = ', ln_leewmx
         WRITE(numout,*) '     Lee wave mixing scaling               = ', rn_lwm
      ENDIF

      !                              ! allocate tmx arrays
      IF( zdf_tmx_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_tmx_init : unable to allocate tmx arrays' )

      IF( ln_tmx_itf ) THEN          ! read the Indonesian Through Flow mask
         CALL iom_open('mask_itf',inum)
         CALL iom_get (inum, jpdom_global, 'tmaskitf',mask_itf,1) ! 
         CALL iom_close(inum)
      ENDIF

      ! read M2 tidal energy flux : W/m2  ( zem2 < 0 )
      CALL iom_open('M2rowdrg',inum)
      CALL iom_get (inum, jpdom_global, 'field',zem2,1) ! 
      CALL iom_close(inum)

      ! read K1 tidal energy flux : W/m2  ( zek1 < 0 )
      CALL iom_open('K1rowdrg',inum)
      CALL iom_get (inum, jpdom_global, 'field',zek1,1) ! 
      CALL iom_close(inum)

      IF( ln_leewmx ) THEN ! read mesoscale eddy energy flux : W/m2  ( zeef < 0 )
         CALL iom_open('Eddyengf',inum)
         CALL iom_get (inum, jpdom_global, 'field',zeef,1) !
         CALL iom_close(inum)
      ENDIF
 
      ! Total tidal energy ( M2, S2 and K1  with S2=(1/2)^2 * M2 )
      ! only the energy available for mixing is taken into account,
      ! (mixing efficiency tidal dissipation efficiency)
      IF( ln_leewmx ) THEN ! include mesoscale eddy energy flux (zeef) in en_tmx
         en_tmx(:,:) = - rn_me * ( rn_tfe * (min(0.,zem2(:,:)) * 1.25 + min(0.,zek1(:,:))) + rn_lwm * zeef(:,:) ) * ssmask(:,:)
      ELSE   
         en_tmx(:,:) = - rn_tfe * rn_me * ( min(0.,zem2(:,:)) * 1.25 + min(0.,zek1(:,:)) ) * ssmask(:,:)
      ENDIF

!============
!TG: Bug for VVL? Should this section be moved out of _init and be updated at every timestep?
      ! Vertical structure (az_tmx)
      DO jj = 1, jpj                ! part independent of the level
         DO ji = 1, jpi
            zhdep(ji,jj) = gdepw_0(ji,jj,mbkt(ji,jj)+1)       ! depth of the ocean
            zfact(ji,jj) = rho0 * rn_htmx * ( 1. - EXP( -zhdep(ji,jj) / rn_htmx ) )
            IF( zfact(ji,jj) /= 0 )   zfact(ji,jj) = en_tmx(ji,jj) / zfact(ji,jj)
         END DO
      END DO
      DO jk= 1, jpk                 ! complete with the level-dependent part
         DO jj = 1, jpj
            DO ji = 1, jpi
               az_tmx(ji,jj,jk) = zfact(ji,jj) * EXP( -( zhdep(ji,jj)-gdepw_0(ji,jj,jk) ) / rn_htmx ) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
!===========

      IF( lwp ) THEN
         ! Control print
         ! Total power consumption due to vertical mixing
         ! zpc = rho0 * 1/rn_me * rn2 * zav_tide
         zav_tide(:,:,:) = 0.e0
         DO jk = 2, jpkm1
            zav_tide(:,:,jk) = az_tmx(:,:,jk) / MAX( rn_n2min, rn2(:,:,jk) )
         END DO

         ztpc = 0.e0
         zpc(:,:,:) = MAX(rn_n2min,rn2(:,:,:)) * zav_tide(:,:,:)
         DO jk= 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztpc = ztpc + e3w_0(ji,jj,jk) * e1t(ji,jj) * e2t(ji,jj) * zpc(ji,jj,jk) * wmask(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
         ztpc= rho0 * 1/(rn_tfe * rn_me) * ztpc

         WRITE(numout,*) 
         WRITE(numout,*) '          Total power consumption of the tidally driven part of Kz : ztpc = ', ztpc * 1.e-12 ,'TW'


         ! control print 2
         zav_tide(:,:,:) = MIN( zav_tide(:,:,:), 60.e-4 )   
         zkz(:,:) = 0.e0
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zkz(ji,jj) = zkz(ji,jj) + e3w_0(ji,jj,jk) * MAX(0.e0, rn2(ji,jj,jk)) * rho0 * zav_tide(ji,jj,jk) * wmask(ji,jj,jk)
               END DO
            END DO
         END DO
         ! Here zkz should be equal to en_tmx ==> multiply by en_tmx/zkz
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( zkz(ji,jj) /= 0.e0 )   THEN
                   zkz(ji,jj) = en_tmx(ji,jj) / zkz(ji,jj)
               ENDIF
            END DO
         END DO
         ztpc = 1.e50
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( zkz(ji,jj) /= 0.e0 )   THEN
                   ztpc = Min( zkz(ji,jj), ztpc)
               ENDIF
            END DO
         END DO
         WRITE(numout,*) '          Min de zkz ', ztpc, ' Max = ', maxval(zkz(:,:) )

         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zav_tide(ji,jj,jk) = zav_tide(ji,jj,jk) * MIN( zkz(ji,jj), 30./6. ) * wmask(ji,jj,jk)  !kz max = 300 cm2/s
               END DO
            END DO
         END DO
         ztpc = 0.e0
         zpc(:,:,:) = Max(0.e0,rn2(:,:,:)) * zav_tide(:,:,:)
         DO jk= 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztpc = ztpc + e3w_0(ji,jj,jk) * e1t(ji,jj) * e2t(ji,jj) * zpc(ji,jj,jk) * wmask(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
         ztpc= rho0 * 1/(rn_tfe * rn_me) * ztpc
         WRITE(numout,*) '          2 Total power consumption of the tidally driven part of Kz : ztpc = ', ztpc * 1.e-12 ,'TW'

         DO jk = 1, jpk
            ze_z =                  SUM( e1t(:,:) * e2t(:,:) * zav_tide(:,:,jk)     * tmask_i(:,:) )   &
               &     / MAX( 1.e-20, SUM( e1t(:,:) * e2t(:,:) * wmask (:,:,jk) * tmask_i(:,:) ) )
            ztpc = 1.E50
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( zav_tide(ji,jj,jk) /= 0.e0 )   ztpc =Min( ztpc, zav_tide(ji,jj,jk) )
               END DO
            END DO
            WRITE(numout,*) '            N2 min - jk= ', jk,'   ', ze_z * 1.e4,' cm2/s min= ',ztpc*1.e4,   &
               &       'max= ', MAXVAL(zav_tide(:,:,jk) )*1.e4, ' cm2/s'
         END DO

         WRITE(numout,*) '          e_tide : ', SUM( e1t*e2t*en_tmx ) / ( rn_tfe * rn_me ) * 1.e-12, 'TW'
         WRITE(numout,*) 
         WRITE(numout,*) '          Initial profile of tidal vertical mixing'
         DO jk = 1, jpk
            DO jj = 1,jpj
               DO ji = 1,jpi
                  zkz(ji,jj) = az_tmx(ji,jj,jk) /MAX( rn_n2min, rn2(ji,jj,jk) )
               END DO
            END DO
            ze_z =                  SUM( e1t(:,:) * e2t(:,:) * zkz(:,:)     * tmask_i(:,:) )   &
               &     / MAX( 1.e-20, SUM( e1t(:,:) * e2t(:,:) * wmask (:,:,jk) * tmask_i(:,:) ) )
            WRITE(numout,*) '                jk= ', jk,'   ', ze_z * 1.e4,' cm2/s'
         END DO
         DO jk = 1, jpk
            zkz(:,:) = az_tmx(:,:,jk) /rn_n2min
            ze_z =                  SUM( e1t(:,:) * e2t(:,:) * zkz(:,:)     * tmask_i(:,:) )   &
               &     / MAX( 1.e-20, SUM( e1t(:,:) * e2t(:,:) * wmask (:,:,jk) * tmask_i(:,:) ) )
            WRITE(numout,*) 
            WRITE(numout,*) '          N2 min - jk= ', jk,'   ', ze_z * 1.e4,' cm2/s min= ',MINVAL(zkz)*1.e4,   &
               &       'max= ', MAXVAL(zkz)*1.e4, ' cm2/s'
         END DO
         !
      ENDIF
      !
   END SUBROUTINE zdf_tmx_init

   !!======================================================================
END MODULE zdftmx
