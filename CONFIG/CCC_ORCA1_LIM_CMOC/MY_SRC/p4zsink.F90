MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  PISCES  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_init  :  Unitialisation of sinking speed parameters
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_sink_alloc

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio4   !: GOC sinking speed
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wscal    !: Calcite and BSi sinking speeds

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking ! <CMOC OR 05/06/2014> Removal of GOC tracer ! , sinking2  !: POC sinking fluxes 
   !                                                          !  (different meanings depending on the parameterization)
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkcal, sinksil   !: CaCO3 and BSi sinking fluxes
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer            !: Small BFe sinking fluxes
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   #if ! defined key_kriest
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer2           !: Big iron sinking fluxes
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   #endif

   INTEGER  :: iksed  = 10

#if  defined key_kriest
   REAL(wp) ::  xkr_sfact    = 250.     !: Sinking factor
   REAL(wp) ::  xkr_stick    = 0.2      !: Stickiness
   REAL(wp) ::  xkr_nnano    = 2.337    !: Nbr of cell in nano size class
   REAL(wp) ::  xkr_ndiat    = 3.718    !: Nbr of cell in diatoms size class
   REAL(wp) ::  xkr_nmeso    = 7.147    !: Nbr of cell in mesozoo  size class
   REAL(wp) ::  xkr_naggr    = 9.877    !: Nbr of cell in aggregates  size class

   REAL(wp) ::  xkr_frac 

   REAL(wp), PUBLIC ::  xkr_dnano       !: Size of particles in nano pool
   REAL(wp), PUBLIC ::  xkr_ddiat       !: Size of particles in diatoms pool
   REAL(wp), PUBLIC ::  xkr_dmeso       !: Size of particles in mesozoo pool
   REAL(wp), PUBLIC ::  xkr_daggr       !: Size of particles in aggregates pool
   REAL(wp), PUBLIC ::  xkr_wsbio_min   !: min vertical particle speed
   REAL(wp), PUBLIC ::  xkr_wsbio_max   !: max vertical particle speed

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   xnumm   !:  maximum number of particles in aggregates
#endif

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zsink.F90 3685 2012-11-27 15:39:02Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_kriest
   !!----------------------------------------------------------------------
   !!   'key_kriest'                                                    ???
   !!----------------------------------------------------------------------

   SUBROUTINE p4z_sink ( kt, jnt )
      !!---------------------------------------------------------------------
      !!                ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to
      !!              gravitational sinking - Kriest parameterization
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt, jnt
      !
      INTEGER  :: ji, jj, jk
      ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   REAL(wp) :: zagg1, zagg4, zagg5, zaggsi, zaggsh ! <CMOC OR 05/06/2014> Removal of GOC tracer !  zagg2, zagg3, zagg4, zagg5, zaggsi, zaggsh
      ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   REAL(wp) :: zagg ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !  , zaggdoc, znumdoc
      REAL(wp) :: znum , zeps, zfm, zgm, zsm
      REAL(wp) :: zdiv , zdiv1, zdiv2, zdiv3, zdiv4, zdiv5
      REAL(wp) :: zval1, zval2, zval3, zval4
      REAL(wp) :: zrfact2
      INTEGER  :: ik1
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: znum3d 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink')
      !
      CALL wrk_alloc( jpi, jpj, jpk, znum3d )
      !
      !     Initialisation of variables used to compute Sinking Speed
      !     ---------------------------------------------------------

      znum3d(:,:,:) = 0.e0
      zval1 = 1. + xkr_zeta
      zval2 = 1. + xkr_zeta + xkr_eta
      zval3 = 1. + xkr_eta

      !     Computation of the vertical sinking speed : Kriest et Evans, 2000
      !     -----------------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj,jk) /= 0.e0 ) THEN
                  znum = trn(ji,jj,jk,jppoc) / ( trn(ji,jj,jk,jpnum) + rtrn ) / xkr_massp
                  ! -------------- To avoid sinking speed over 50 m/day -------
                  znum  = MIN( xnumm(jk), znum )
                  znum  = MAX( 1.1      , znum )
                  znum3d(ji,jj,jk) = znum
                  !------------------------------------------------------------
                  zeps  = ( zval1 * znum - 1. )/ ( znum - 1. )
                  zfm   = xkr_frac**( 1. - zeps )
                  zgm   = xkr_frac**( zval1 - zeps )
                  zdiv  = MAX( 1.e-4, ABS( zeps - zval2 ) ) * SIGN( 1., ( zeps - zval2 ) )
                  zdiv1 = zeps - zval3
                  wsbio3(ji,jj,jk) = xkr_wsbio_min * ( zeps - zval1 ) / zdiv    &
                     &             - xkr_wsbio_max *   zgm * xkr_eta  / zdiv
                  wsbio4(ji,jj,jk) = xkr_wsbio_min *   ( zeps-1. )    / zdiv1   &
                     &             - xkr_wsbio_max *   zfm * xkr_eta  / zdiv1
                  IF( znum == 1.1)   wsbio3(ji,jj,jk) = wsbio4(ji,jj,jk)
               ENDIF
            END DO
         END DO
      END DO

      wscal(:,:,:) = MAX( wsbio3(:,:,:), 50._wp )

      !   INITIALIZE TO ZERO ALL THE SINKING ARRAYS
      !   -----------------------------------------

      sinking (:,:,:) = 0.e0
      ! <CMOC OR 05/06/2014> Removal of GOC tracer ! sinking2(:,:,:) = 0.e0
      ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   sinkcal (:,:,:) = 0.e0
      ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   sinkfer (:,:,:) = 0.e0
      ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   sinksil (:,:,:) = 0.e0

     !   Compute the sedimentation term using p4zsink2 for all the sinking particles
     !   -----------------------------------------------------

      CALL p4z_sink2( wsbio3, sinking , jppoc )
      ! <CMOC OR 05/06/2014> Removal of GOC tracer ! CALL p4z_sink2( wsbio4, sinking2, jpnum )
      ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   CALL p4z_sink2( wsbio3, sinkfer , jpsfe )
      ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   CALL p4z_sink2( wscal , sinksil , jpgsi )
      ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   CALL p4z_sink2( wscal , sinkcal , jpcal )

     !  Exchange between organic matter compartments due to coagulation/disaggregation
     !  ---------------------------------------------------

      zval1 = 1. + xkr_zeta
      zval2 = 1. + xkr_eta
      zval3 = 3. + xkr_eta
      zval4 = 4. + xkr_eta

      DO jk = 1,jpkm1
         DO jj = 1,jpj
            DO ji = 1,jpi
               IF( tmask(ji,jj,jk) /= 0.e0 ) THEN

                  znum = trn(ji,jj,jk,jppoc)/(trn(ji,jj,jk,jpnum)+rtrn) / xkr_massp
                  !-------------- To avoid sinking speed over 50 m/day -------
                  znum  = min(xnumm(jk),znum)
                  znum  = MAX( 1.1,znum)
                  !------------------------------------------------------------
                  zeps  = ( zval1 * znum - 1.) / ( znum - 1.)
                  zdiv  = MAX( 1.e-4, ABS( zeps - zval3) ) * SIGN( 1., zeps - zval3 )
                  zdiv1 = MAX( 1.e-4, ABS( zeps - 4.   ) ) * SIGN( 1., zeps - 4.    )
                  zdiv2 = zeps - 2.
                  zdiv3 = zeps - 3.
                  zdiv4 = zeps - zval2
                  zdiv5 = 2.* zeps - zval4
                  zfm   = xkr_frac**( 1.- zeps )
                  zsm   = xkr_frac**xkr_eta

                  !    Part I : Coagulation dependant on turbulence
                  !    ----------------------------------------------

                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zagg1 = ( 0.163 * trn(ji,jj,jk,jpnum)**2               &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      &            * 2.*( (zfm-1.)*(zfm*xkr_mass_max**3-xkr_mass_min**3)    &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      &            * (zeps-1)/zdiv1 + 3.*(zfm*xkr_mass_max-xkr_mass_min)    &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      &            * (zfm*xkr_mass_max**2-xkr_mass_min**2)                  &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      &            * (zeps-1.)**2/(zdiv2*zdiv3)) 
                  ! <CMOC OR 05/06/2014> Removal of GOC tracer ! zagg2 =  2*0.163*trn(ji,jj,jk,jpnum)**2*zfm*                       &
                  ! <CMOC OR 05/06/2014> Removal of GOC tracer !    &                   ((xkr_mass_max**3+3.*(xkr_mass_max**2          &
                  ! <CMOC OR 05/06/2014> Removal of GOC tracer !    &                    *xkr_mass_min*(zeps-1.)/zdiv2                 &
                  ! <CMOC OR 05/06/2014> Removal of GOC tracer !    &                    +xkr_mass_max*xkr_mass_min**2*(zeps-1.)/zdiv3)    &
                  ! <CMOC OR 05/06/2014> Removal of GOC tracer !    &                    +xkr_mass_min**3*(zeps-1)/zdiv1)                  &
                  ! <CMOC OR 05/06/2014> Removal of GOC tracer !    &                    -zfm*xkr_mass_max**3*(1.+3.*((zeps-1.)/           &
                  ! <CMOC OR 05/06/2014> Removal of GOC tracer !    &                    (zeps-2.)+(zeps-1.)/zdiv3)+(zeps-1.)/zdiv1))    

                  ! <CMOC OR 05/06/2014> Removal of GOC tracer ! zagg3 =  0.163*trn(ji,jj,jk,jpnum)**2*zfm**2*8. * xkr_mass_max**3  
                  
                 !    Aggregation of small into large particles
                 !    Part II : Differential settling
                 !    ----------------------------------------------

                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zagg4 =  2.*3.141*0.125*trn(ji,jj,jk,jpnum)**2*                       &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      &                 xkr_wsbio_min*(zeps-1.)**2                         &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      &                 *(xkr_mass_min**2*((1.-zsm*zfm)/(zdiv3*zdiv4)      &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      &                 -(1.-zfm)/(zdiv*(zeps-1.)))-                       &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      &                 ((zfm*zfm*xkr_mass_max**2*zsm-xkr_mass_min**2)     &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !     &                 *xkr_eta)/(zdiv*zdiv3*zdiv5) )   

                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zagg5 =   2.*3.141*0.125*trn(ji,jj,jk,jpnum)**2                         &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      &                 *(zeps-1.)*zfm*xkr_wsbio_min                        &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      &                 *(zsm*(xkr_mass_min**2-zfm*xkr_mass_max**2)         &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      &                 /zdiv3-(xkr_mass_min**2-zfm*zsm*xkr_mass_max**2)    &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !      &                 /zdiv)  
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zaggsi = ( zagg4 + zagg5 ) * xstep / 10.

                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zagg = 0.5 * xkr_stick * ( zaggsh + zaggsi )

                  !     Aggregation of DOC to small particles
                  !     --------------------------------------

                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !  zaggdoc = ( 0.4 * trn(ji,jj,jk,jpdoc)               &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !     &        + 1018.  * trn(ji,jj,jk,jppoc)  ) * xstep    &
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !     &        * xdiss(ji,jj,jk) * trn(ji,jj,jk,jpdoc)

! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   # if defined key_degrad
                   ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zagg1   = zagg1   * facvol(ji,jj,jk)                 
                   ! <CMOC OR 05/06/2014> Removal of GOC tracer ! zagg2   = zagg2   * facvol(ji,jj,jk)                 
                   ! <CMOC OR 05/06/2014> Removal of GOC tracer ! zagg3   = zagg3   * facvol(ji,jj,jk)                 
                   ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zagg4   = zagg4   * facvol(ji,jj,jk)                 
                   ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zagg5   = zagg5   * facvol(ji,jj,jk)                 
                   ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !  zaggdoc = zaggdoc * facvol(ji,jj,jk)                 
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   # endif
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zaggsh =   zagg1   * rfact2 * xdiss(ji,jj,jk) / 1000. ! <CMOC OR 05/06/2014> Removal of GOC tracer ! ( zagg1 + zagg2 + zagg3 ) * rfact2 * xdiss(ji,jj,jk) / 1000.
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zaggsi = ( zagg4 + zagg5 ) * xstep / 10.
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zagg = 0.5 * xkr_stick * ( zaggsh + zaggsi )
                  !
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !  znumdoc = trn(ji,jj,jk,jpnum) / ( trn(ji,jj,jk,jppoc) + rtrn )
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !  tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zaggdoc
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) - zagg ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !  + zaggdoc * znumdoc - zagg
                  ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zaggdoc

               ENDIF
            END DO
         END DO
      END DO

      IF( ln_diatrc ) THEN
         !
         ik1 = iksed + 1
         zrfact2 = 1.e3 * rfact2r
         IF( jnt == nrdttrc ) THEN
           CALL iom_put( "POCFlx"  , sinking (:,:,:)      * zrfact2 * tmask(:,:,:) )  ! POC export
           ! <CMOC OR 05/06/2014> Removal of GOC tracer ! CALL iom_put( "NumFlx"  , sinking2 (:,:,:)     * zrfact2 * tmask(:,:,:) )  ! Num export
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !              CALL iom_put( "SiFlx"   , sinksil (:,:,:)      * zrfact2 * tmask(:,:,:) )  ! Silica export
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !              CALL iom_put( "CaCO3Flx", sinkcal (:,:,:)      * zrfact2 * tmask(:,:,:) )  ! Calcite export
           CALL iom_put( "xnum"    , znum3d  (:,:,:)                * tmask(:,:,:) )  ! Number of particles in aggregats
           CALL iom_put( "W1"      , wsbio3  (:,:,:)                * tmask(:,:,:) )  ! sinking speed of POC
           CALL iom_put( "W2"      , wsbio4  (:,:,:)                * tmask(:,:,:) )  ! sinking speed of aggregats
           CALL iom_put( "PMO"     , sinking (:,:,ik1)    * zrfact2 * tmask(:,:,1) )  ! POC export at 100m
           ! <CMOC OR 05/06/2014> Removal of GOC tracer ! CALL iom_put( "PMO2"    , sinking2(:,:,ik1)    * zrfact2 * tmask(:,:,1) )  ! Num export at 100m
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !              CALL iom_put( "ExpFe1"  , sinkfer (:,:,ik1)    * zrfact2 * tmask(:,:,1) )  ! Export of iron at 100m
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !              CALL iom_put( "ExpSi"   , sinksil (:,:,ik1)    * zrfact2 * tmask(:,:,1) )  ! export of silica at 100m
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !              CALL iom_put( "ExpCaCO3", sinkcal (:,:,ik1)    * zrfact2 * tmask(:,:,1) )  ! export of calcite at 100m
         ENDIF
# if ! defined key_iomput
         trc2d(:,:  ,jp_pcs0_2d + 4)  = sinking (:,:,ik1)    * zrfact2 * tmask(:,:,1)
         ! <CMOC OR 05/06/2014> Removal of GOC tracer ! trc2d(:,:  ,jp_pcs0_2d + 5)  = sinking2(:,:,ik1)    * zrfact2 * tmask(:,:,1)
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !            trc2d(:,:  ,jp_pcs0_2d + 6)  = sinkfer (:,:,ik1)    * zrfact2 * tmask(:,:,1)
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !            trc2d(:,:  ,jp_pcs0_2d + 7)  = sinksil (:,:,ik1)    * zrfact2 * tmask(:,:,1)
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !            trc2d(:,:  ,jp_pcs0_2d + 8)  = sinkcal (:,:,ik1)    * zrfact2 * tmask(:,:,1)
         trc3d(:,:,:,jp_pcs0_3d + 11) = sinking (:,:,:)      * zrfact2 * tmask(:,:,:)
         ! <CMOC OR 05/06/2014> Removal of GOC tracer ! trc3d(:,:,:,jp_pcs0_3d + 12) = sinking2(:,:,:)      * zrfact2 * tmask(:,:,:)
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !            trc3d(:,:,:,jp_pcs0_3d + 13) = sinksil (:,:,:)      * zrfact2 * tmask(:,:,:)
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !            trc3d(:,:,:,jp_pcs0_3d + 14) = sinkcal (:,:,:)      * zrfact2 * tmask(:,:,:)
         trc3d(:,:,:,jp_pcs0_3d + 15) = znum3d  (:,:,:)                * tmask(:,:,:)
         trc3d(:,:,:,jp_pcs0_3d + 16) = wsbio3  (:,:,:)                * tmask(:,:,:)
         trc3d(:,:,:,jp_pcs0_3d + 17) = wsbio4  (:,:,:)                * tmask(:,:,:)
# endif
        !
      ENDIF
      
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, jpk, znum3d )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink


   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!
      !! ** Purpose :   Initialization of sinking parameters
      !!                Kriest parameterization only
      !!
      !! ** Method  :   Read the nampiskrs namelist and check the parameters
      !!      called at the first timestep 
      !!
      !! ** input   :   Namelist nampiskrs
      !!----------------------------------------------------------------------
      INTEGER  ::   jk, jn, kiter
      REAL(wp) ::   znum, zdiv
      REAL(wp) ::   zws, zwr, zwl,wmax, znummax
      REAL(wp) ::   zmin, zmax, zl, zr, xacc
      !
      NAMELIST/nampiskrs/ xkr_sfact, xkr_stick ,  &
         &                xkr_nnano, xkr_ndiat, xkr_nmeso, xkr_naggr
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink_init')
      !
      REWIND( numnatp )                     ! read nampiskrs
      READ  ( numnatp, nampiskrs )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist : nampiskrs'
         WRITE(numout,*) '    Sinking factor                           xkr_sfact    = ', xkr_sfact
         WRITE(numout,*) '    Stickiness                               xkr_stick    = ', xkr_stick
         WRITE(numout,*) '    Nbr of cell in nano size class           xkr_nnano    = ', xkr_nnano
         WRITE(numout,*) '    Nbr of cell in diatoms size class        xkr_ndiat    = ', xkr_ndiat
         WRITE(numout,*) '    Nbr of cell in mesozoo size class        xkr_nmeso    = ', xkr_nmeso
         WRITE(numout,*) '    Nbr of cell in aggregates size class     xkr_naggr    = ', xkr_naggr
      ENDIF


      ! max and min vertical particle speed
      xkr_wsbio_min = xkr_sfact * xkr_mass_min**xkr_eta
      xkr_wsbio_max = xkr_sfact * xkr_mass_max**xkr_eta
      WRITE(numout,*) ' max and min vertical particle speed ', xkr_wsbio_min, xkr_wsbio_max

      !
      !    effect of the sizes of the different living pools on particle numbers
      !    nano = 2um-20um -> mean size=6.32 um -> ws=2.596 -> xnum=xnnano=2.337
      !    diat and microzoo = 10um-200um -> 44.7 -> 8.732 -> xnum=xndiat=3.718
      !    mesozoo = 200um-2mm -> 632.45 -> 45.14 -> xnum=xnmeso=7.147
      !    aggregates = 200um-10mm -> 1414 -> 74.34 -> xnum=xnaggr=9.877
      !    doc aggregates = 1um
      ! ----------------------------------------------------------

      xkr_dnano = 1. / ( xkr_massp * xkr_nnano )
      xkr_ddiat = 1. / ( xkr_massp * xkr_ndiat )
      xkr_dmeso = 1. / ( xkr_massp * xkr_nmeso )
      xkr_daggr = 1. / ( xkr_massp * xkr_naggr )

      !!---------------------------------------------------------------------
      !!    'key_kriest'                                                  ???
      !!---------------------------------------------------------------------
      !  COMPUTATION OF THE VERTICAL PROFILE OF MAXIMUM SINKING SPEED
      !  Search of the maximum number of particles in aggregates for each k-level.
      !  Bissection Method
      !--------------------------------------------------------------------
      WRITE(numout,*)
      WRITE(numout,*)'    kriest : Compute maximum number of particles in aggregates'

      xacc     =  0.001_wp
      kiter    = 50
      zmin     =  1.10_wp
      zmax     = xkr_mass_max / xkr_mass_min
      xkr_frac = zmax

      DO jk = 1,jpk
         zl = zmin
         zr = zmax
         wmax = 0.5 * fse3t(1,1,jk) * rday / rfact2
         zdiv = xkr_zeta + xkr_eta - xkr_eta * zl
         znum = zl - 1.
         zwl =  xkr_wsbio_min * xkr_zeta / zdiv &
            & - ( xkr_wsbio_max * xkr_eta * znum * &
            &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
            & - wmax

         zdiv = xkr_zeta + xkr_eta - xkr_eta * zr
         znum = zr - 1.
         zwr =  xkr_wsbio_min * xkr_zeta / zdiv &
            & - ( xkr_wsbio_max * xkr_eta * znum * &
            &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
            & - wmax
iflag:   DO jn = 1, kiter
            IF    ( zwl == 0._wp ) THEN   ;   znummax = zl
            ELSEIF( zwr == 0._wp ) THEN   ;   znummax = zr
            ELSE
               znummax = ( zr + zl ) / 2.
               zdiv = xkr_zeta + xkr_eta - xkr_eta * znummax
               znum = znummax - 1.
               zws =  xkr_wsbio_min * xkr_zeta / zdiv &
                  & - ( xkr_wsbio_max * xkr_eta * znum * &
                  &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
                  & - wmax
               IF( zws * zwl < 0. ) THEN   ;   zr = znummax
               ELSE                        ;   zl = znummax
               ENDIF
               zdiv = xkr_zeta + xkr_eta - xkr_eta * zl
               znum = zl - 1.
               zwl =  xkr_wsbio_min * xkr_zeta / zdiv &
                  & - ( xkr_wsbio_max * xkr_eta * znum * &
                  &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
                  & - wmax

               zdiv = xkr_zeta + xkr_eta - xkr_eta * zr
               znum = zr - 1.
               zwr =  xkr_wsbio_min * xkr_zeta / zdiv &
                  & - ( xkr_wsbio_max * xkr_eta * znum * &
                  &     xkr_frac**( -xkr_zeta / znum ) / zdiv ) &
                  & - wmax
               !
               IF ( ABS ( zws )  <= xacc ) EXIT iflag
               !
            ENDIF
            !
         END DO iflag

         xnumm(jk) = znummax
         WRITE(numout,*) '       jk = ', jk, ' wmax = ', wmax,' xnum max = ', xnumm(jk)
         !
      END DO
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink_init')
      !
  END SUBROUTINE p4z_sink_init

#else

   SUBROUTINE p4z_sink ( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, jnt
      INTEGER  ::   ji, jj, jk
      ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   REAL(wp) ::   zagg1, zagg4 ! <CMOC OR 05/05/2014> Removal of GOC tracer ! zagg2, zagg3, zagg4
      ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   REAL(wp) ::   zagg , zaggfe ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !  , zaggdoc, zaggdoc3 ! <CMOC OR 05/05/2014> Removal of GOC tracer ! zaggdoc2, zaggdoc3
      REAL(wp) ::   zfact, zwsmax, zmax, zstep
      REAL(wp) ::   zrfact2
      INTEGER  ::   ik1
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink')
      !
      !    Sinking speeds of detritus is increased with depth as shown
      !    by data and from the coagulation theory
      !    -----------------------------------------------------------
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         DO jk = 1, jpkm1
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !            DO jj = 1, jpj
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !               DO ji = 1,jpi
      !         zmax  = MAX( heup(ji,jj), hmld(ji,jj) )
      !         zfact = MAX( 0., fsdepw(ji,jj,jk+1) - zmax ) / 5000._wp
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !                  zmax = hmld(ji,jj)
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !                  zfact = MAX( 0., fsdepw(ji,jj,jk+1) - zmax ) / 4000._wp
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !                  wsbio4(ji,jj,jk) = wsbio2 + ( 200.- wsbio2 ) * zfact
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !               END DO
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !            END DO
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         END DO

      ! limit the values of the sinking speeds to avoid numerical instabilities  
      wsbio3(:,:,:) = wsbio
      !
      ! OA Below, this is garbage. the ideal would be to find a time-splitting 
      ! OA algorithm that does not increase the computing cost by too much
      ! OA In ROMS, I have included a time-splitting procedure. But it is 
      ! OA too expensive as the loop is computed globally. Thus, a small e3t
      ! OA at one place determines the number of subtimesteps globally
      ! OA AWFULLY EXPENSIVE !! Not able to find a better approach. Damned !!

      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zwsmax = 0.8 * fse3t(ji,jj,jk) / xstep
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !                  wsbio4(ji,jj,jk) = MIN( wsbio4(ji,jj,jk), zwsmax )
               wsbio3(ji,jj,jk) = MIN( wsbio3(ji,jj,jk), zwsmax )
            END DO
         END DO
      END DO

! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         wscal(:,:,:) = wsbio4(:,:,:)

      !  Initializa to zero all the sinking arrays 
      !   -----------------------------------------

      sinking (:,:,:) = 0.e0
      ! <CMOC OR 05/06/2014> Removal of GOC tracer ! sinking2(:,:,:) = 0.e0
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         sinkcal (:,:,:) = 0.e0
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         sinkfer (:,:,:) = 0.e0
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         sinksil (:,:,:) = 0.e0
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         sinkfer2(:,:,:) = 0.e0

      !   Compute the sedimentation term using p4zsink2 for all the sinking particles
      !   -----------------------------------------------------

      CALL p4z_sink2( wsbio3, sinking , jppoc )
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         ! <CMOC OR 05/05/2014> Removal of GOC tracer ! CALL p4z_sink2( wsbio4, sinking2, jpgoc )
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         CALL p4z_sink2( wsbio4, sinkfer2, jpbfe )
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         CALL p4z_sink2( wsbio4, sinksil , jpgsi )
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         CALL p4z_sink2( wscal , sinkcal , jpcal )

      !  Exchange between organic matter compartments due to coagulation/disaggregation
      !  ---------------------------------------------------

! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         DO jk = 1, jpkm1
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !            DO jj = 1, jpj
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !               DO ji = 1, jpi
               !
               ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zstep = xstep 
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   # if defined key_degrad
               ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zstep = zstep * facvol(ji,jj,jk)
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   # endif
               ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zfact = zstep * xdiss(ji,jj,jk)
               !  Part I : Coagulation dependent on turbulence
               ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zagg1 = 354.  * zfact * trn(ji,jj,jk,jppoc) * trn(ji,jj,jk,jppoc)
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! zagg2 = 4452. * zfact * trn(ji,jj,jk,jppoc) * trn(ji,jj,jk,jpgoc)

               ! Part II : Differential settling

               !  Aggregation of small into large particles
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! zagg3 =  4.7 * zstep * trn(ji,jj,jk,jppoc) * trn(ji,jj,jk,jpgoc)
               ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zagg4 =  0.4 * zstep * trn(ji,jj,jk,jppoc) * trn(ji,jj,jk,jppoc)

               ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zagg   = zagg1 + zagg4 ! <CMOC OR 05/05/2014> Removal of GOC tracer ! + zagg2 + zagg3 + zagg4
               ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   zaggfe = zagg * trn(ji,jj,jk,jpsfe) / ( trn(ji,jj,jk,jppoc) + rtrn )

               ! Aggregation of DOC to small particles
               ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !  zaggdoc  = ( 0.83 * trn(ji,jj,jk,jpdoc) + 271. * trn(ji,jj,jk,jppoc) ) * zfact * trn(ji,jj,jk,jpdoc)
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! zaggdoc2 = 1.07e4 * zfact * trn(ji,jj,jk,jpgoc) * trn(ji,jj,jk,jpdoc)
               ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !  zaggdoc3 =   0.02 * ( 16706. * trn(ji,jj,jk,jppoc) + 231. * trn(ji,jj,jk,jpdoc) ) * zstep * trn(ji,jj,jk,jpdoc)

               !  Update the trends
               ! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !   tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) !- zagg + zaggdoc + zaggdoc3 <CMOC OR 11/06/2013> No aggregation allowed in the CMOC experiment
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zagg + zaggdoc2
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !                  tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zaggfe
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !                  tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zaggfe
               ! <CMOC OR 06/30/2014> Trimming code, tracers (jpdoc) !  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zaggdoc - zaggdoc3 ! <CMOC OR 05/05/2014> Removal of GOC tracer ! - zaggdoc2 - zaggdoc3
               !
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !               END DO
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !            END DO
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !         END DO

      IF( ln_diatrc ) THEN
         zrfact2 = 1.e3 * rfact2r
         ik1  = iksed + 1
         IF( lk_iomput ) THEN
           IF( jnt == nrdttrc ) THEN
              CALL iom_put( "EPC100"  ,   sinking(:,:,ik1)                       * zrfact2 * tmask(:,:,1) ) ! <CMOC OR 0115/2014> ( sinking(:,:,ik1) + sinking2(:,:,ik1) ) * zrfact2 * tmask(:,:,1) ) ! Export of carbon at 100m
! <CMOC OR 05/27/2014> trimming PISCES code  in p4zprod        CALL iom_put( "EPFE100" , ( sinkfer(:,:,ik1) + sinkfer2(:,:,ik1) ) * zrfact2 * tmask(:,:,1) ) ! Export of iron at 100m
              CALL iom_put( "EPCALC100",  wsbio3(:,:,11) / rday * trn(:,:,11,jppoc) * 1e3_wp * xrcico(:,:) * tmask(:,:,1) ) ! <CMOC OR 03/13/2014> reverse to 2D xrcico ! <CMOC OR 02/19/2014> replace xfracal(:,:,ik1) by xrcico(:,:) and POC flux by proxy used for alkalinity and dic ! <CMOC OR 01/27/2014> use CMOC definition of Calcite Export FPON * Rain Ratio  ! Export of calcite  at 100m 
! <CMOC OR 05/27/2014> trimming PISCES code  in p4zprod                CALL iom_put( "EPSI100" ,   sinksil(:,:,ik1)                       * zrfact2 * tmask(:,:,1) ) ! Export of biogenic silica at 100m
           ENDIF
         ELSE
           trc2d(:,:,jp_pcs0_2d + 4) = sinking (:,:,ik1) * zrfact2 * tmask(:,:,1)
           ! <CMOC OR 05/06/2014> Removal of GOC tracer ! trc2d(:,:,jp_pcs0_2d + 5) = sinking2(:,:,ik1) * zrfact2 * tmask(:,:,1)
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !              trc2d(:,:,jp_pcs0_2d + 6) = sinkfer (:,:,ik1) * zrfact2 * tmask(:,:,1)
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !              trc2d(:,:,jp_pcs0_2d + 7) = sinkfer2(:,:,ik1) * zrfact2 * tmask(:,:,1)
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !              trc2d(:,:,jp_pcs0_2d + 8) = sinksil (:,:,ik1) * zrfact2 * tmask(:,:,1)
! <CMOC OR 06/30/2014> Trimming code, tracers (jpsil .. jpgsi, jpdch, jpcal, jpfer .. jpdfe .. jpdfe, jpnum) !              trc2d(:,:,jp_pcs0_2d + 9) = sinkcal (:,:,ik1) * zrfact2 * tmask(:,:,1)
         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink

   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!----------------------------------------------------------------------

       ! <CMOC OR 01/28/2014> Initialize xfracal to 0 everywhere to avoid NaN values when calculating PIC export
! <CMOC OR 06/13/2014> Code trimming !        INTEGER :: ji, jj, jk
 
            ! <CMOC OR 03/08/2014> CMOC namelist
      NAMELIST/namcmocrr/  cnrr_cmoc, ncrr_cmoc
      ! <CMOC OR 03/08/2014> CMOC namelist end 
      !!----------------------------------------------------------------------

       ! <CMOC OR 01/28/2014> Initialize xfracal to 0 everywhere to avoid NaN values when calculating PIC export
! <CMOC OR 06/13/2014> Code trimming !      DO jk = 1, jpk
! <CMOC OR 06/13/2014> Code trimming !           DO ji = 1, jpi
! <CMOC OR 06/13/2014> Code trimming !             DO jj = 1, jpj
 
! <CMOC OR 06/13/2014> Code trimming !  ! <CMOC OR 06/13/2014> Code trimming !                xfracal(ji,jj,jk) = 0._wp
 
! <CMOC OR 06/13/2014> Code trimming !             ENDDO
! <CMOC OR 06/13/2014> Code trimming !           ENDDO
! <CMOC OR 06/13/2014> Code trimming !         ENDDO
       ! <CMOC OR 01/28/2014>

      REWIND( numcmoc )                     ! <CMOC OR 03/10/2014> ! read numcmoc, cmocrr
      READ  ( numcmoc, namcmocrr  )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters, Redfield ratios, namcmocrr'    
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    C:N ratio                                 cnrr_cmoc =', cnrr_cmoc
         WRITE(numout,*) '    N:C ratio                                 ncrr_cmoc =', ncrr_cmoc
      ENDIF

   END SUBROUTINE p4z_sink_init

#endif



   SUBROUTINE p4z_sink2( pwsink, psinkflx, jp_tra )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink2  ***
      !!
      !! ** Purpose :   Compute the sedimentation terms for the various sinking
      !!     particles. The scheme used to compute the trends is based
      !!     on MUSCL.
      !!
      !! ** Method  : - this ROUTINE compute not exactly the advection but the
      !!      transport term, i.e.  div(u*tra).
      !!---------------------------------------------------------------------
      !
      INTEGER , INTENT(in   )                         ::   jp_tra    ! tracer index index      
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj,jpk) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   psinkflx  ! sinking fluxe
      !!
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zigma,zew,zign, zflx, zstep
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztraz, zakz, zwsink2, ztrb 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink2')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, jpk, ztraz, zakz, zwsink2, ztrb )

      zstep = rfact2 / 2.

      ztraz(:,:,:) = 0.e0
      zakz (:,:,:) = 0.e0
      ztrb (:,:,:) = trn(:,:,:,jp_tra)

      DO jk = 1, jpkm1
         zwsink2(:,:,jk+1) = -pwsink(:,:,jk) / rday * tmask(:,:,jk+1) 
      END DO
      zwsink2(:,:,1) = 0.e0
      IF( lk_degrad ) THEN
         zwsink2(:,:,:) = zwsink2(:,:,:) * facvol(:,:,:)
      ENDIF


      ! Vertical advective flux
      DO jn = 1, 2
         !  first guess of the slopes interior values
         DO jk = 2, jpkm1
            ztraz(:,:,jk) = ( trn(:,:,jk-1,jp_tra) - trn(:,:,jk,jp_tra) ) * tmask(:,:,jk)
         END DO
         ztraz(:,:,1  ) = 0.0
         ztraz(:,:,jpk) = 0.0

         ! slopes
         DO jk = 2, jpkm1
            DO jj = 1,jpj
               DO ji = 1, jpi
                  zign = 0.25 + SIGN( 0.25, ztraz(ji,jj,jk) * ztraz(ji,jj,jk+1) )
                  zakz(ji,jj,jk) = ( ztraz(ji,jj,jk) + ztraz(ji,jj,jk+1) ) * zign
               END DO
            END DO
         END DO
         
         ! Slopes limitation
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zakz(ji,jj,jk) = SIGN( 1., zakz(ji,jj,jk) ) *        &
                     &             MIN( ABS( zakz(ji,jj,jk) ), 2. * ABS(ztraz(ji,jj,jk+1)), 2. * ABS(ztraz(ji,jj,jk) ) )
               END DO
            END DO
         END DO
         
         ! vertical advective flux
         DO jk = 1, jpkm1
            DO jj = 1, jpj      
               DO ji = 1, jpi    
                  zigma = zwsink2(ji,jj,jk+1) * zstep / fse3w(ji,jj,jk+1)
                  zew   = zwsink2(ji,jj,jk+1)
                  psinkflx(ji,jj,jk+1) = -zew * ( trn(ji,jj,jk,jp_tra) - 0.5 * ( 1 + zigma ) * zakz(ji,jj,jk) ) * zstep
               END DO
            END DO
         END DO
         !
         ! Boundary conditions
         psinkflx(:,:,1  ) = 0.e0
         psinkflx(:,:,jpk) = 0.e0
         
         DO jk=1,jpkm1
            DO jj = 1,jpj
               DO ji = 1, jpi
                  zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
                  trn(ji,jj,jk,jp_tra) = trn(ji,jj,jk,jp_tra) + zflx
               END DO
            END DO
         END DO

      ENDDO

      DO jk=1,jpkm1
         DO jj = 1,jpj
            DO ji = 1, jpi
               zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
               ztrb(ji,jj,jk) = ztrb(ji,jj,jk) + 2. * zflx
            END DO
         END DO
      END DO

      trn     (:,:,:,jp_tra) = ztrb(:,:,:)
      psinkflx(:,:,:)        = 2. * psinkflx(:,:,:) ! <CMOC OR 01/15/2014> aka "sinking" term
      !
      CALL wrk_dealloc( jpi, jpj, jpk, ztraz, zakz, zwsink2, ztrb )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink2')
      !
   END SUBROUTINE p4z_sink2


   INTEGER FUNCTION p4z_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( wsbio3 (jpi,jpj,jpk) , wsbio4  (jpi,jpj,jpk) , wscal(jpi,jpj,jpk) ,     &
         &      sinking(jpi,jpj,jpk) ,                                                  & ! <CMOC OR 05/06/2014> Removal of GOC tracer ! sinking2(jpi,jpj,jpk)                      ,     &                
 ! <CMOC OR 07/15/2014> Revert for consistency !          &      sinkcal(jpi,jpj,jpk) , sinksil (jpi,jpj,jpk)                      ,     &                
#if defined key_kriest
         &      xnumm(jpk)                                                        ,     &                
! <CMOC OR 07/15/2014> ! Removal of all the tracers !  #else
! <CMOC OR 07/15/2014> ! Removal of all the tracers !           &      sinkfer2(jpi,jpj,jpk)                                             ,     &                
#endif
         &                                                                          STAT=p4z_sink_alloc )    ! sinkfer(jpi,jpj,jpk)            ! <CMOC OR 07/15/2014> ! Removal of all the tracers !  
         !
      IF( p4z_sink_alloc /= 0 ) CALL ctl_warn('p4z_sink_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_sink_alloc
   
#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sink                    ! Empty routine
   END SUBROUTINE p4z_sink
#endif 

   !!======================================================================
END MODULE  p4zsink
