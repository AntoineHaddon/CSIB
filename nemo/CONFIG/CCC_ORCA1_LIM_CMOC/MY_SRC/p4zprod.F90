MODULE p4zprod
   !!======================================================================
   !!                         ***  MODULE p4zprod  ***
   !! TOP :  Growth Rate of the two phytoplanktons groups 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_prod       :   Compute the growth Rate of the two phytoplanktons groups
   !!   p4z_prod_init  :   Initialization of the parameters for growth
! <CMOC OR 05/26/2014> trimming PISCES code     !!   p4z_prod_alloc :   Allocate variables for growth
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zopt          !  optical model
! <CMOC OR 07/15/2014> ! Removal of all the tracers !     USE p4zlim          !  Co-limitations of differents nutrients
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_prod         ! called in p4zbio.F90
   PUBLIC   p4z_prod_init    ! called in trcsms_pisces.F90
! <CMOC OR 05/26/2014> trimming PISCES code    PUBLIC   p4z_prod_alloc

   !! * Shared module variables
! <CMOC OR 05/26/2014> trimming PISCES code    LOGICAL , PUBLIC ::  ln_newprod = .FALSE.
! <CMOC OR 05/26/2014> trimming PISCES code     REAL(wp), PUBLIC ::  pislope    = 3.0_wp            !:
!   REAL(wp), PUBLIC ::  pislope2   = 3.0_wp            !:
!   REAL(wp), PUBLIC ::  excret     = 10.e-5_wp         !:
!   REAL(wp), PUBLIC ::  excret2    = 0.05_wp           !:
!   REAL(wp), PUBLIC ::  bresp      = 0.00333_wp        !:
!   REAL(wp), PUBLIC ::  chlcnm     = 0.033_wp          !:
!   REAL(wp), PUBLIC ::  chlcdm     = 0.05_wp           !:
!   REAL(wp), PUBLIC ::  chlcmin    = 0.00333_wp        !:
!   REAL(wp), PUBLIC ::  fecnm      = 10.E-6_wp         !:
!   REAL(wp), PUBLIC ::  fecdm      = 15.E-6_wp         !:
! <CMOC OR 05/26/2014> trimming PISCES code     REAL(wp), PUBLIC ::  grosip     = 0.151_wp          !:

! <CMOC OR 05/26/2014> trimming PISCES code     REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prmax    !: optimal production = f(temperature)
!   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotan   !: proxy of N quota in Nanophyto
! <CMOC OR 05/26/2014> trimming PISCES code     REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotad   !: proxy of N quota in diatomee
   
   REAL(wp) :: r1_rday                !: 1 / rday
! <CMOC OR 05/26/2014> trimming PISCES code     REAL(wp) :: texcret                !: 1 - excret 
! <CMOC OR 05/26/2014> trimming PISCES code     REAL(wp) :: texcret2               !: 1 - excret2        
! <CMOC OR 05/26/2014> trimming PISCES code     REAL(wp) :: tpp                    !: Total primary production


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zprod.F90 3773 2013-02-07 11:06:58Z cbricaud $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_prod( kt , jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod  ***
      !!
      !! ** Purpose :   Compute the phytoplankton production depending on
      !!              light, temperature and nutrient availability
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt, jnt
      !
      INTEGER  ::   ji, jj, jk
! <CMOC OR 05/26/2014> trimming PISCES code       REAL(wp) ::   zsilfac, zfact, znanotot, zdiattot, zconctemp, zconctemp2
      REAL(wp) ::   zfact ! <CMOC OR 05/26/2014> trimming PISCES code 
! <CMOC OR 05/26/2014> trimming PISCES code       REAL(wp) ::   zratio, zmax, zsilim, ztn, zadap
      REAL(wp) ::   ztn, zadap
! <CMOC OR 05/26/2014> trimming PISCES code        REAL(wp) ::   zlim, zsilfac2, zsiborn, zprod, zproreg, zproreg2
      REAL(wp) ::   zprod
! <CMOC OR 05/26/2014> trimming PISCES code        REAL(wp) ::   zmxltst, zmxlday, zmaxday
! <CMOC OR 05/26/2014> trimming PISCES code        REAL(wp) ::   zpislopen  , zpislope2n
      REAL(wp) ::   zpislopen 
! <CMOC OR 05/26/2014> trimming PISCES code       REAL(wp) ::   zrum, zcodel, zargu, zval
      REAL(wp) ::   zrfact2
      CHARACTER (len=25) :: charout
! <CMOC OR 05/26/2014> trimming PISCES code        REAL(wp), POINTER, DIMENSION(:,:  ) :: zmixnano, zmixdiat, zstrn
! <CMOC OR 05/26/2014> trimming PISCES code        REAL(wp), POINTER, DIMENSION(:,:,:) :: zpislopead, zpislopead2, zprdia, zprbio, zprdch, zprnch, zysopt   
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zpislopead, zprbio, zprnch
! <CMOC OR 05/26/2014> trimming PISCES code        REAL(wp), POINTER, DIMENSION(:,:,:) :: zprorca, zprorcad, zprofed, zprofen, zprochln, zprochld, zpronew, zpronewd
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zprorca, zprochln
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zlimn,zliml ! <CMOC OR 02/21/2014> nitrogen and light limitation functions
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_prod')
      !
      !  Allocate temporary workspace
! <CMOC OR 05/26/2014> trimming PISCES code        CALL wrk_alloc( jpi, jpj,      zmixnano, zmixdiat, zstrn                                                  )
! <CMOC OR 05/26/2014> trimming PISCES code        CALL wrk_alloc( jpi, jpj, jpk, zpislopead, zpislopead2, zprdia, zprbio, zprdch, zprnch, zysopt            ) 
      CALL wrk_alloc( jpi, jpj, jpk, zpislopead, zprbio, zprnch            ) 
! <CMOC OR 05/26/2014> trimming PISCES code        CALL wrk_alloc( jpi, jpj, jpk, zprorca, zprorcad, zprofed, zprofen, zprochln, zprochld, zpronew, zpronewd )
      CALL wrk_alloc( jpi, jpj, jpk, zprorca, zprochln                     )
      CALL wrk_alloc( jpi, jpj, jpk, zlimn, zliml                                                               ) ! <CMOC OR 02/21/2014> 
      !
      zprorca (:,:,:) = 0._wp
! <CMOC OR 05/26/2014> trimming PISCES code        zprorcad(:,:,:) = 0._wp
! <CMOC OR 05/26/2014> trimming PISCES code        zprofed (:,:,:) = 0._wp
! <CMOC OR 05/26/2014> trimming PISCES code        zprofen (:,:,:) = 0._wp
      zprochln(:,:,:) = 0._wp
! <CMOC OR 05/26/2014> trimming PISCES code        zprochld(:,:,:) = 0._wp
! <CMOC OR 05/26/2014> trimming PISCES code        zpronew (:,:,:) = 0._wp
! <CMOC OR 05/26/2014> trimming PISCES code        zpronewd(:,:,:) = 0._wp
! <CMOC OR 05/26/2014> trimming PISCES code        zprdia  (:,:,:) = 0._wp
      zprbio  (:,:,:) = 0._wp
! <CMOC OR 05/26/2014> trimming PISCES code        zprdch  (:,:,:) = 0._wp
      zprnch  (:,:,:) = 0._wp
! <CMOC OR 05/26/2014> trimming PISCES code        zysopt  (:,:,:) = 0._wp
      zlimn   (:,:,:) = 0._wp
      zliml   (:,:,:) = 0._wp

! <CMOC OR 05/26/2014> trimming PISCES code      ! Computation of the optimal production
!      prmax(:,:,:) = 0.6_wp * r1_rday * tgfunc(:,:,:) 
!      IF( lk_degrad )  prmax(:,:,:) = prmax(:,:,:) * facvol(:,:,:) 
!
!      ! compute the day length depending on latitude and the day
!      zrum = REAL( nday_year - 80, wp ) / REAL( nyear_len(1), wp )
!      zcodel = ASIN(  SIN( zrum * rpi * 2._wp ) * SIN( rad * 23.5_wp )  )
!
!      ! day length in hours
!      zstrn(:,:) = 0.
!      DO jj = 1, jpj
!         DO ji = 1, jpi
!            zargu = TAN( zcodel ) * TAN( gphit(ji,jj) * rad )
!            zargu = MAX( -1., MIN(  1., zargu ) )
!            zstrn(ji,jj) = MAX( 0.0, 24. - 2. * ACOS( zargu ) / rad / 15. )
!         END DO
! <CMOC OR 05/26/2014> trimming PISCES code      END DO

! <CMOC OR 05/26/2014> trimming PISCES code
!      IF( ln_newprod ) THEN
!         ! Impact of the day duration on phytoplankton growth
!         DO jk = 1, jpkm1
!            DO jj = 1 ,jpj
!               DO ji = 1, jpi
!                  zval = MAX( 1., zstrn(ji,jj) )
!                  zval = 1.5 * zval / ( 12. + zval )
!                  zprbio(ji,jj,jk) = prmax(ji,jj,jk) * zval
!                  zprdia(ji,jj,jk) = zprbio(ji,jj,jk)
!               END DO
!            END DO
!         END DO
!      ENDIF
! <CMOC OR 05/26/2014> trimming PISCES code

! <CMOC OR 05/26/2014> trimming PISCES code
!      ! Maximum light intensity
!      WHERE( zstrn(:,:) < 1.e0 ) zstrn(:,:) = 24.
!      zstrn(:,:) = 24. / zstrn(:,:)
! <CMOC OR 05/26/2014> trimming PISCES code

! <CMOC OR 05/26/2014> trimming PISCES code
!      IF( ln_newprod ) THEN
!!CDIR NOVERRCHK
!         DO jk = 1, jpkm1
!!CDIR NOVERRCHK
!            DO jj = 1, jpj
!!CDIR NOVERRCHK
!               DO ji = 1, jpi
!
!                  ! Computation of the P-I slope for nanos and diatoms
!                  IF( etot(ji,jj,jk) > 1.E-3 ) THEN
!                      ztn    = MAX( 0., tsn(ji,jj,jk,jp_tem) - 15. )
!                      zadap  = ztn / ( 2.+ ztn )
!
!                      zconctemp   = MAX( 0.e0 , trn(ji,jj,jk,jpdia) - 5e-7 )
!                      zconctemp2  = trn(ji,jj,jk,jpdia) - zconctemp
!
!                      znanotot = enano(ji,jj,jk) * zstrn(ji,jj)
!                      zdiattot = ediat(ji,jj,jk) * zstrn(ji,jj)
!
!                      zfact  = EXP( -0.21 * znanotot )
!                      zpislopead (ji,jj,jk) = pislope  * ( 1.+ zadap  * zfact )  &
!                         &                   * trn(ji,jj,jk,jpnch) /( trn(ji,jj,jk,jpphy) * 12. + rtrn)
!
!                      zpislopead2(ji,jj,jk) = (pislope * zconctemp2 + pislope2 * zconctemp) / ( trn(ji,jj,jk,jpdia) + rtrn )   &
!                         &                   * trn(ji,jj,jk,jpdch) /( trn(ji,jj,jk,jpdia) * 12. + rtrn)
!
!                      ! Computation of production function for Carbon
!                      !  ---------------------------------------------
!                      zpislopen  = zpislopead (ji,jj,jk) / ( ( r1_rday + bresp * r1_rday / chlcnm ) * rday + rtrn)
!                      zpislope2n = zpislopead2(ji,jj,jk) / ( ( r1_rday + bresp * r1_rday / chlcdm ) * rday + rtrn)
!                      zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1.- EXP( -zpislopen  * znanotot )  )
!                      zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * ( 1.- EXP( -zpislope2n * zdiattot )  )
!
!                      !  Computation of production function for Chlorophyll
!                      !--------------------------------------------------
!                      zmaxday  = 1._wp / ( prmax(ji,jj,jk) * rday + rtrn )
!                      zprnch(ji,jj,jk) = prmax(ji,jj,jk) * ( 1.- EXP( -zpislopead (ji,jj,jk) * zmaxday * znanotot ) ) 
!                      zprdch(ji,jj,jk) = prmax(ji,jj,jk) * ( 1.- EXP( -zpislopead2(ji,jj,jk) * zmaxday * zdiattot ) )
!                  ENDIF
!               END DO
!            END DO
!         END DO
!      ELSE
! <CMOC OR 05/26/2014> trimming PISCES code

!CDIR NOVERRCHK
         DO jk = 1, jpkm1
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi

                  ! Computation of the P-I slope for nanos and diatoms
                  IF( etot(ji,jj,jk) > 1.E-3 ) THEN
                      ztn    = tsn(ji,jj,jk,jp_tem) + 273.15_wp                                  !MAX( 0., tsn(ji,jj,jk,jp_tem) - 15. ) <CMOC OR 09/19/13>
                      zadap  = ep_cmoc * 1e+3_wp / 8.31_wp * ( 1._wp / ( ztn + rtrn ) - 1._wp / ( tvm_cmoc + 273.15_wp) ) !ztn / ( 2.+ ztn )  <CMOC OR 03/09/2014> "softwired" parameters replace hardwired ones  <CMOC OR 09/19/13>
                      !<CMOC OR 09/25/13> zadap is the argument of the Arrhenius' function
                      
                      zfact  = EXP ( -zadap )                                                  !EXP( -0.21 * enano(ji,jj,jk) )         <CMOC OR 09/19/13>
                      zpislopead (ji,jj,jk) = vm_cmoc * r1_rday * zfact   ! <CMOC OR 03/09/2014> "softwired" parameters replace hardwired ones      !pislope  * ( 1.+ zadap  * zfact ) <CMOC OR 09/19/13>, this is CMOC's Vm, the photosynthetic rate at temperate ToC
! <CMOC OR 05/26/2014> trimming PISCES code                        zpislopead2(ji,jj,jk) = pislope2

                      zpislopen =  achl_cmoc * trn(ji,jj,jk,jpnch)                                & ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones  ! zpislopead(ji,jj,jk) * trn(ji,jj,jk,jpnch) <CMOC OR 09/19/13>
                        &          / ( trn(ji,jj,jk,jpphy) * 12._wp + rtrn)  & ! <CMOC OR 10/29/2013> Phytoplankton is converted from molC to gC
                        &          / ( zpislopead(ji,jj,jk) * rday  + rtrn )                            ! ( prmax(ji,jj,jk) * rday * xlimphy(ji,jj,jk) <CMOC OR 09/19/13> this is alpha x theta x PAR / Vm

! <CMOC OR 05/26/2014> trimming PISCES code                        zpislope2n = zpislopead2(ji,jj,jk) * trn(ji,jj,jk,jpdch)                &
!                        &          / ( trn(ji,jj,jk,jpdia) * 12._wp                  + rtrn )    &
! <CMOC OR 05/26/2014> trimming PISCES code                          &          / ( prmax(ji,jj,jk) * rday * xlimdia(ji,jj,jk) + rtrn )

                      ! Computation of production function for Carbon
                      !  ---------------------------------------------
                      zliml (ji,jj,jk) = 1.- EXP( -zpislopen  * etot(ji,jj,jk) )                                               ! <CMOC OR 02/21/2014> light limitation
                      zlimn (ji,jj,jk) = trn(ji,jj,jk,jpno3) / ( kn_cmoc * 1e-6_wp * cnrr_cmoc + trn(ji,jj,jk,jpno3)+ rtrn )              ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones ! <CMOC OR 02/21/2014> nutrient limitation
                      zprbio(ji,jj,jk) = zpislopead(ji,jj,jk) * min ( zliml(ji,jj,jk) , zlimn(ji,jj,jk) , xlimnfecmoc(ji,jj) ) ! <CMOC OR 02/21/2014> switch for zliml and zlimn ! <CMOC OR 01/21/2014> add the iron limitation !<CMOC OR 11/19/2013> enano (PAR for nanophyto) is replaced by etot (PAR for phyto) !<CMOC OR 10/29/13> this is Gamma
! <CMOC OR 05/26/2014> trimming PISCES code                        zprdia(ji,jj,jk) = prmax(ji,jj,jk)      *     ( 1.- EXP( -zpislope2n * ediat(ji,jj,jk) ) )

                      !  Computation of production function for Chlorophyll
                      !--------------------------------------------------
                      !zprnch(ji,jj,jk) = prmax(ji,jj,jk) * ( 1.- EXP( -zpislopen  * enano(ji,jj,jk) * zstrn(ji,jj) ) ) !<CMOC OR 09/20/13>
                      
                      ! <CMOC OR 09/29/13> Below is CMOC's Chl:C balanced ratio unit is gChl per molC 
                      zprnch(ji,jj,jk) = 12._wp *  thm_cmoc  * 2._wp  * zpislopead(ji,jj,jk)  /  ( 2._wp * zpislopead(ji,jj,jk)  +  achl_cmoc * thm_cmoc  * etot(ji,jj,jk) * r1_rday + rtrn )  ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones ! <CMOC OR 11/19/2013> replace enano with etot ! <CMOC OR 10/29/13> enano has to be a layer-averaged daily mean of the PAR otherwise the previous equation is invalid

                      ! <CMOC OR 05/26/2014> trimming PISCES code  zprdch(ji,jj,jk) = prmax(ji,jj,jk) * ( 1.- EXP( -zpislope2n * ediat(ji,jj,jk) * zstrn(ji,jj) ) )
                  ENDIF
               END DO
            END DO
         END DO

! <CMOC OR 05/26/2014> trimming PISCES code      ENDIF

! <CMOC OR 05/26/2014> trimming PISCES code       !  Computation of a proxy of the N/C ratio
!      !  ---------------------------------------
!!CDIR NOVERRCHK
!      DO jk = 1, jpkm1
!!CDIR NOVERRCHK
!         DO jj = 1, jpj
!!CDIR NOVERRCHK
!            DO ji = 1, jpi
!                zval = ( xnanonh4(ji,jj,jk) + xnanono3(ji,jj,jk) ) * prmax(ji,jj,jk) / ( zprbio(ji,jj,jk) + rtrn )
!                quotan(ji,jj,jk) = MIN( 1., 0.5 + 0.5 * zval )
!                zval = ( xdiatnh4(ji,jj,jk) + xdiatno3(ji,jj,jk) ) * prmax(ji,jj,jk) / ( zprdia(ji,jj,jk) + rtrn )
!                quotad(ji,jj,jk) = MIN( 1., 0.5 + 0.5 * zval )
!            END DO
!         END DO
!      END DO
!
!
!      DO jk = 1, jpkm1
!         DO jj = 1, jpj
!            DO ji = 1, jpi
!
!                IF( etot(ji,jj,jk) > 1.E-3 ) THEN
!                   !    Si/C of diatoms
!                   !    ------------------------
!                   !    Si/C increases with iron stress and silicate availability
!                   !    Si/C is arbitrariliy increased for very high Si concentrations
!                   !    to mimic the very high ratios observed in the Southern Ocean (silpot2)
!                  zlim  = trn(ji,jj,jk,jpsil) / ( trn(ji,jj,jk,jpsil) + xksi1 )
!                  zsilim = MIN( zprdia(ji,jj,jk) / ( prmax(ji,jj,jk) + rtrn ), xlimsi(ji,jj,jk) )
!                  zsilfac = 4.4 * EXP( -4.23 * zsilim ) * MAX( 0.e0, MIN( 1., 2.2 * ( zlim - 0.5 ) )  ) + 1.e0
!                  zsiborn = MAX( 0.e0, ( trn(ji,jj,jk,jpsil) - 15.e-6 ) )
!                  zsilfac2 = 1.+ 2.* zsiborn / ( zsiborn + xksi2 )
!                  zsilfac = MIN( 5.4, zsilfac * zsilfac2)
!                  zysopt(ji,jj,jk) = grosip * zlim * zsilfac
!              ENDIF
!            END DO
!         END DO
!      END DO
!
!      !  Computation of the limitation term due to a mixed layer deeper than the euphotic depth
!      DO jj = 1, jpj
!         DO ji = 1, jpi
!            zmxltst = MAX( 0.e0, hmld(ji,jj) - heup(ji,jj) )
!            zmxlday = zmxltst * zmxltst * r1_rday
!!            zmixnano(ji,jj) = 1. - zmxlday / ( 3. + zmxlday )             ! <CMOC OR 02/26/2014>
!            zmixdiat(ji,jj) = 1. - zmxlday / ( 4. + zmxlday )
!         END DO
!      END DO
! 
!      !  Mixed-layer effect on production                                                                               
!      DO jk = 1, jpkm1
!         DO jj = 1, jpj
!            DO ji = 1, jpi
!               IF( fsdepw(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
!!                  zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * zmixnano(ji,jj)   ! <CMOC OR 02/26/2014>
!                  zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * zmixdiat(ji,jj)
!               ENDIF
!            END DO
!         END DO
!      END DO
! <CMOC OR 05/26/2014> trimming PISCES code 

      ! Computation of the various production terms 
!CDIR NOVERRCHK
      DO jk = 1, jpkm1
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               IF( etot(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for nanophyto.
                  zprorca(ji,jj,jk) = zprbio(ji,jj,jk)  * trn(ji,jj,jk,jpphy) * rfact2 !* xlimphy(ji,jj,jk) * trn(ji,jj,jk,jpphy) * rfact2  <CMOC OR 09/19/13>
                  !zpronew(ji,jj,jk) = zprorca(ji,jj,jk) * xnanono3(ji,jj,jk) / ( xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) + rtrn ) <CMOC OR 09/25/13>
                  !
! <CMOC OR 05/26/2014> trimming PISCES code                    zratio = trn(ji,jj,jk,jpnfe) / ( trn(ji,jj,jk,jpphy) + rtrn )
! <CMOC OR 05/26/2014> trimming PISCES code                    zratio = zratio / fecnm 
! <CMOC OR 05/26/2014> trimming PISCES code                    zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) ) 
! <CMOC OR 05/26/2014> trimming PISCES code                    zprofen(ji,jj,jk) = fecnm * prmax(ji,jj,jk)  &
!                  &             * ( 4. - 4.5 * xlimnfe(ji,jj,jk) / ( xlimnfe(ji,jj,jk) + 0.5 ) )    &
!                  &             * trn(ji,jj,jk,jpfer) / ( trn(ji,jj,jk,jpfer) + concnfe(ji,jj,jk) )  &
! <CMOC OR 05/26/2014> trimming PISCES code                    &             * zmax * trn(ji,jj,jk,jpphy) * rfact2
                  !  production terms for diatomees
! <CMOC OR 05/26/2014> trimming PISCES code                    zprorcad(ji,jj,jk) = zprdia(ji,jj,jk) * xlimdia(ji,jj,jk) * trn(ji,jj,jk,jpdia) * rfact2
! <CMOC OR 05/26/2014> trimming PISCES code                    zpronewd(ji,jj,jk) = zprorcad(ji,jj,jk) * xdiatno3(ji,jj,jk) / ( xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) + rtrn )
                  !
! <CMOC OR 05/26/2014> trimming PISCES code                    zratio = trn(ji,jj,jk,jpdfe) / ( trn(ji,jj,jk,jpdia) + rtrn )
! <CMOC OR 05/26/2014> trimming PISCES code                    zratio = zratio / fecdm 
! <CMOC OR 05/26/2014> trimming PISCES code                    zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) ) 
! <CMOC OR 05/26/2014> trimming PISCES code                    zprofed(ji,jj,jk) = fecdm * prmax(ji,jj,jk)  &
!                  &             * ( 4. - 4.5 * xlimdfe(ji,jj,jk) / ( xlimdfe(ji,jj,jk) + 0.5 ) )    &
!                  &             * trn(ji,jj,jk,jpfer) / ( trn(ji,jj,jk,jpfer) + concdfe(ji,jj,jk) )  &
! <CMOC OR 05/26/2014> trimming PISCES code                    &             * zmax * trn(ji,jj,jk,jpdia) * rfact2
               ENDIF
            END DO
         END DO
      END DO

! <CMOC OR 05/26/2014> trimming PISCES code       IF( ln_newprod ) THEN
!!CDIR NOVERRCHK
!         DO jk = 1, jpkm1
!!CDIR NOVERRCHK
!            DO jj = 1, jpj
!!CDIR NOVERRCHK
!               DO ji = 1, jpi
!                  IF( fsdepw(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
!                     zprnch(ji,jj,jk) = zprnch(ji,jj,jk) * zmixnano(ji,jj)
!                     zprdch(ji,jj,jk) = zprdch(ji,jj,jk) * zmixdiat(ji,jj)
!                  ENDIF
!                  IF( etot(ji,jj,jk) > 1.E-3 ) THEN
!                     !  production terms for nanophyto. ( chlorophyll )
!                     znanotot = enano(ji,jj,jk) * zstrn(ji,jj)
!                     zprod    = rday * zprorca(ji,jj,jk) * zprnch(ji,jj,jk) * xlimphy(ji,jj,jk)
!                     zprochln(ji,jj,jk) = chlcmin * 12. * zprorca (ji,jj,jk)
!                     zprochln(ji,jj,jk) = zprochln(ji,jj,jk) + chlcnm * 12. * zprod /  &
!                                        & (  zpislopead(ji,jj,jk) * znanotot +rtrn)
!                     !  production terms for diatomees ( chlorophyll )
!                     zdiattot = ediat(ji,jj,jk) * zstrn(ji,jj)
!                     zprod = rday * zprorcad(ji,jj,jk) * zprdch(ji,jj,jk) * xlimdia(ji,jj,jk)
!                     zprochld(ji,jj,jk) = chlcmin * 12. * zprorcad(ji,jj,jk)
!                     zprochld(ji,jj,jk) = zprochld(ji,jj,jk) + chlcdm * 12. * zprod /  &
!                                        & ( zpislopead2(ji,jj,jk) * zdiattot +rtrn )
!                  ENDIF
!               END DO
!            END DO
!         END DO
! <CMOC OR 05/26/2014> trimming PISCES code       ELSE

!CDIR NOVERRCHK
         DO jk = 1, jpkm1
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi
                  IF( etot(ji,jj,jk) > 1.E-3 ) THEN
                     !  production terms for nanophyto. ( chlorophyll )
                     ! <CMOC OR 05/26/2014> trimming PISCES code  znanotot = enano(ji,jj,jk) * zstrn(ji,jj)
                     zprod = zprbio(ji,jj,jk)  * trn(ji,jj,jk,jpnch) * rfact2 !rday * zprorca(ji,jj,jk) * zprnch(ji,jj,jk) * trn(ji,jj,jk,jpphy) * xlimphy(ji,jj,jk) <CMOC OR 09/25/13> Chlorophyll production
                     !zprochln(ji,jj,jk) = chlcnm * 144. * zprod / &
                     !                   & (  zpislopead(ji,jj,jk) * trn(ji,jj,jk,jpnch) * znanotot +rtrn)
                     !<CMOC OR 09/20/13>

                     zprochln(ji,jj,jk) = zprod + ( zprnch (ji,jj,jk) * trn(ji,jj,jk,jpphy) - trn(ji,jj,jk,jpnch) ) * itau_cmoc * r1_rday * rfact2 ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones ! <CMOC OR 11/27/2013> convert 1/tau0 (0.5) from per day to per second  !<CMOC OR 09/20/13> chlorophyll production plus relaxation-to-balanced-ratio term
                     !  production terms for diatomees ( chlorophyll ) 
                     ! <CMOC OR 05/26/2014> trimming PISCES code  zdiattot = ediat(ji,jj,jk) * zstrn(ji,jj)
! <CMOC OR 05/26/2014> trimming PISCES code                       zprod = rday * zprorcad(ji,jj,jk) * zprdch(ji,jj,jk) * trn(ji,jj,jk,jpdia) * xlimdia(ji,jj,jk)
! <CMOC OR 05/26/2014> trimming PISCES code                       zprochld(ji,jj,jk) = chlcdm * 144._wp * zprod / &
! <CMOC OR 05/26/2014> trimming PISCES code                                          & ( zpislopead2(ji,jj,jk) * trn(ji,jj,jk,jpdch) * zdiattot +rtrn )
                  ENDIF
               END DO
            END DO
         END DO
! <CMOC OR 05/26/2014> trimming PISCES code       ENDIF

      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1 ,jpi
              !zproreg  = zprorca(ji,jj,jk) - zpronew(ji,jj,jk) <CMOC OR 09/25/13>
! <CMOC OR 05/26/2014> trimming PISCES code                zproreg2 = zprorcad(ji,jj,jk) - zpronewd(ji,jj,jk)
! <CMOC OR 05/26/2014> trimming PISCES code                tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - zprorca(ji,jj,jk) - zprorcad(ji,jj,jk)
              tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zprorca(ji,jj,jk) ! - zpronew(ji,jj,jk) - zpronewd(ji,jj,jk) <CMOC OR 09/25/13>
! <CMOC OR 05/26/2014> trimming PISCES code                tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zproreg2 !- zproreg - zproreg2 <CMOC OR 09/25/13>
              tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprorca(ji,jj,jk) ! * texcret <CMOC OR 09/19/13>, equivalent to set texcret to 1; NOTE: This is still in Carbon currency and it needs to be set in Nitrogen currency
              tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) + zprochln(ji,jj,jk)! * texcret <CMOC OR 09/19/13>, equivalent to set texcret to 1; NOTE: This is still in Carbon currency and it needs to be set in Nitrogen currency
! <CMOC OR 05/26/2014> trimming PISCES code                tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) + zprofen(ji,jj,jk) * texcret
! <CMOC OR 05/26/2014> trimming PISCES code                tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) + zprorcad(ji,jj,jk) * texcret2
! <CMOC OR 05/26/2014> trimming PISCES code                tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) + zprochld(ji,jj,jk) * texcret2
! <CMOC OR 05/26/2014> trimming PISCES code                tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) + zprofed(ji,jj,jk) * texcret2
! <CMOC OR 05/26/2014> trimming PISCES code                tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zprorcad(ji,jj,jk) * zysopt(ji,jj,jk) * texcret2
! <CMOC OR 05/26/2014> trimming PISCES code                tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + excret2 * zprorcad(ji,jj,jk) + excret * zprorca(ji,jj,jk)
              tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + zprorca(ji,jj,jk) ! * 1.3_wp!+ o2ut * ( zproreg + zproreg2 ) & <CMOC OR 11/14/13> dioxygen changes are mirroring DIC changes as far as PP is concerned
                ! &                + ( o2ut + o2nit ) * ( zpronew(ji,jj,jk) + zpronewd(ji,jj,jk) ) <CMOC OR 09/29/13>
! <CMOC OR 05/26/2014> trimming PISCES code                tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - texcret * zprofen(ji,jj,jk) - texcret2 * zprofed(ji,jj,jk)
! <CMOC OR 05/26/2014> trimming PISCES code                tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) - texcret2 * zprorcad(ji,jj,jk) * zysopt(ji,jj,jk)
              tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprorca(ji,jj,jk) !- zprorcad(ji,jj,jk) <CMOC OR 10/03/2013>
              tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + ncrr_cmoc * zprorca(ji,jj,jk) ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones !+ rno3 * ( zpronew(ji,jj,jk) + zpronewd(ji,jj,jk) ) & <CMOC OR 09/29/13>
                ! &                                      - rno3 * ( zproreg + zproreg2 ) <CMOC OR 09/29/13>
          END DO
        END DO
     END DO

! <CMOC OR 05/26/2014> trimming PISCES code       ! Total primary production per year
!     tpp = tpp + glob_sum( ( zprorca(:,:,:) + zprorcad(:,:,:) ) * cvol(:,:,:) )
!
!     IF( kt == nitend .AND. jnt == nrdttrc ) THEN
!        WRITE(numout,*) 'Total PP (Gtc) :'
!        WRITE(numout,*) '-------------------- : ',tpp * 12. / 1.E12
!        WRITE(numout,*) 
! <CMOC OR 05/26/2014> trimming PISCES code        ENDIF

! <CMOC OR 05/26/2014> trimming PISCES code       IF( ln_diatrc ) THEN
         !
         zrfact2 = 1.e3 * rfact2r  ! conversion from mol/L/timestep into mol/m3/s
         IF( lk_iomput ) THEN
           IF( jnt == nrdttrc ) THEN
              CALL iom_put( "PPPHY"   , zprorca (:,:,:) * zrfact2 * tmask(:,:,:) )  ! primary production by nanophyto
! <CMOC OR 05/26/2014> trimming PISCES code                CALL iom_put( "PPPHY2"  , zprorcad(:,:,:) * zrfact2 * tmask(:,:,:) )  ! primary production by diatom
              !CALL iom_put( "PPNEWN"  , zpronew (:,:,:) * zrfact2 * tmask(:,:,:) )  <CMOC OR 09/26/13> ! new primary production by nanophyto
! <CMOC OR 05/26/2014> trimming PISCES code                CALL iom_put( "PPNEWD"  , zpronewd(:,:,:) * zrfact2 * tmask(:,:,:) )  ! new primary production by diatom
! <CMOC OR 05/26/2014> trimming PISCES code                CALL iom_put( "PBSi"    , zprorcad(:,:,:) * zrfact2 * tmask(:,:,:) * zysopt(:,:,:) ) ! biogenic silica production
! <CMOC OR 05/26/2014> trimming PISCES code                CALL iom_put( "PFeD"    , zprofed (:,:,:) * zrfact2 * tmask(:,:,:) )  ! biogenic iron production by diatom
! <CMOC OR 05/26/2014> trimming PISCES code                CALL iom_put( "PFeN"    , zprofen (:,:,:) * zrfact2 * tmask(:,:,:) )  ! biogenic iron production by nanophyto
              CALL iom_put( "Mumax"   , zpislopead  (:,:,:) * rday * tmask(:,:,:) )  ! Maximum growth rate ! <CMOC OR 02/26/2014>
! <CMOC OR 05/26/2014> trimming PISCES code                CALL iom_put( "MuN"     , zprbio  (:,:,:) * xlimphy(:,:,:) * tmask(:,:,:) )  ! Realized growth rate for nanophyto
! <CMOC OR 05/26/2014> trimming PISCES code                CALL iom_put( "MuD"     , zprdia  (:,:,:) * xlimdia(:,:,:) * tmask(:,:,:) )  ! Realized growth rate for diatoms
! <CMOC OR 05/26/2014> trimming PISCES code                CALL iom_put( "MuNlight", zprbio  (:,:,:) * tmask(:,:,:) )  ! Light limited growth rate phytoplankton
! <CMOC OR 05/26/2014> trimming PISCES code                CALL iom_put( "MuDlight", zprdia  (:,:,:) * tmask(:,:,:) )  ! Light limited growth rate diatoms
              CALL iom_put( "LNnut"   , zlimn   (:,:,:) * tmask(:,:,:) )  ! <CMOC OR 02/21/2014> xlimphy (:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
! <CMOC OR 05/26/2014> trimming PISCES code                CALL iom_put( "LDnut"   , xlimdia (:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
              CALL iom_put( "LNFe"    , xlimnfecmoc (:,:) * tmask(:,:,1) )  ! <CMOC OR 02/21/2014> xlimnfecmoc instead of xlimnfe ! Iron limitation term
! <CMOC OR 05/26/2014> trimming PISCES code                CALL iom_put( "LDFe"    , xlimdfe (:,:,:) * tmask(:,:,:) )  ! Iron limitation term
              CALL iom_put( "LNlight" , zliml   (:,:,:) * tmask(:,:,:) )  ! <CMOC OR 02/21/2014> zprbio  (:,:,:) / (prmax(:,:,:) + rtrn) * tmask(:,:,:) )  ! light limitation term
! <CMOC OR 05/26/2014> trimming PISCES code                CALL iom_put( "LDlight" , zprdia  (:,:,:) / (prmax(:,:,:) + rtrn) * tmask(:,:,:) )  ! light limitation term
           ENDIF
! <CMOC OR 05/26/2014> trimming PISCES code           ELSE
!              trc3d(:,:,:,jp_pcs0_3d + 4)  = zprorca (:,:,:) * zrfact2 * tmask(:,:,:)
!              trc3d(:,:,:,jp_pcs0_3d + 5)  = zprorcad(:,:,:) * zrfact2 * tmask(:,:,:)
!              !trc3d(:,:,:,jp_pcs0_3d + 6)  = zpronew (:,:,:) * zrfact2 * tmask(:,:,:) <CMOC OR 09/26/13>
!              trc3d(:,:,:,jp_pcs0_3d + 7)  = zpronewd(:,:,:) * zrfact2 * tmask(:,:,:)
!              trc3d(:,:,:,jp_pcs0_3d + 8)  = zprorcad(:,:,:) * zrfact2 * tmask(:,:,:) * zysopt(:,:,:)
!              trc3d(:,:,:,jp_pcs0_3d + 9)  = zprofed (:,:,:) * zrfact2 * tmask(:,:,:)
!#  if ! defined key_kriest
!              trc3d(:,:,:,jp_pcs0_3d + 10) = zprofen (:,:,:) * zrfact2 * tmask(:,:,:)
!#  endif
          ENDIF
!         !
! <CMOC OR 05/26/2014> trimming PISCES code        ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
! <CMOC OR 05/26/2014> trimming PISCES code        CALL wrk_dealloc( jpi, jpj,      zmixnano, zmixdiat, zstrn                                                  )
! <CMOC OR 05/26/2014> trimming PISCES code        CALL wrk_dealloc( jpi, jpj, jpk, zpislopead, zpislopead2, zprdia, zprbio, zprdch, zprnch, zysopt            ) 
      CALL wrk_dealloc( jpi, jpj, jpk, zpislopead, zprbio, zprnch            ) 
! <CMOC OR 05/26/2014> trimming PISCES code        CALL wrk_dealloc( jpi, jpj, jpk, zprorca, zprorcad, zprofed, zprofen, zprochln, zprochld, zpronew, zpronewd )
      CALL wrk_dealloc( jpi, jpj, jpk, zprorca, zprochln                     )
      CALL wrk_dealloc( jpi, jpj, jpk, zlimn, zliml                                                               ) ! <CMOC OR 02/21/2014>
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_prod')
      !
   END SUBROUTINE p4z_prod


   SUBROUTINE p4z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      !
! <CMOC OR 05/26/2014> trimming PISCES code        NAMELIST/nampisprod/ pislope, pislope2, ln_newprod, bresp, excret, excret2,  &
! <CMOC OR 05/26/2014> trimming PISCES code           &                 chlcnm, chlcdm, chlcmin, fecnm, fecdm, grosip
      ! <CMOC OR 03/08/2014> CMOC namelist
      NAMELIST/namcmocphy/ achl_cmoc, thm_cmoc, tau_cmoc, itau_cmoc, ep_cmoc,     &
         &                 tvm_cmoc, vm_cmoc, kn_cmoc
      ! <CMOC OR 03/08/2014> CMOC namelist end 
         !!----------------------------------------------------------------------

      
      
! <CMOC OR 05/26/2014> trimming PISCES code        REWIND( numnatp )                     ! read numnatp
! <CMOC OR 05/26/2014> trimming PISCES code        READ  ( numnatp, nampisprod )

      REWIND( numcmoc )                     ! <CMOC OR 03/10/2014> ! read numcmoc, cmocphy
      READ  ( numcmoc, namcmocphy )


      
      IF(lwp) THEN                         ! control print
! <CMOC OR 05/26/2014> trimming PISCES code           WRITE(numout,*) ' '
!         WRITE(numout,*) ' Namelist parameters for phytoplankton growth, nampisprod'
!         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!         WRITE(numout,*) '    Enable new parame. of production (T/F)   ln_newprod   =', ln_newprod
!         WRITE(numout,*) '    mean Si/C ratio                           grosip       =', grosip
!         WRITE(numout,*) '    P-I slope                                 pislope      =', pislope
!         WRITE(numout,*) '    excretion ratio of nanophytoplankton      excret       =', excret
!         WRITE(numout,*) '    excretion ratio of diatoms                excret2      =', excret2
! <CMOC OR 05/26/2014> trimming PISCES code           WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton growth, namcmocphy'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Initial slope of the P-I curve            achl_cmoc    =', achl_cmoc
         WRITE(numout,*) '    Maximum chlorophyll relaxation time       thm_cmoc     =', thm_cmoc
         WRITE(numout,*) '    Chlorophyll relaxation time               tau_cmoc     =', tau_cmoc
         WRITE(numout,*) '    Inverse of chlorophyll relaxation time    itau_cmoc    =', itau_cmoc
         WRITE(numout,*) '    Activation energy for growth              ep_cmoc      =', ep_cmoc
         WRITE(numout,*) '    Reference ocean temperature               tvm_cmoc     =', tvm_cmoc
         WRITE(numout,*) '    Reference maximum photosynth. rate        vm_cmoc      =', vm_cmoc
         WRITE(numout,*) '    Half-saturation constant for N            kn_cmoc      =', kn_cmoc
         WRITE(numout,*) ' '
! <CMOC OR 05/26/2014> trimming PISCES code          IF( ln_newprod )  THEN
!            WRITE(numout,*) ' '
!            WRITE(numout,*) ' Namelist parameters for phytoplankton growth, nampisprod'
!            WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!            WRITE(numout,*) '    basal respiration in phytoplankton        bresp        =', bresp
!            WRITE(numout,*) '    Maximum Chl/C in phytoplankton            chlcmin      =', chlcmin
! <CMOC OR 05/26/2014> trimming PISCES code          ENDIF
! <CMOC OR 05/26/2014> trimming PISCES code           WRITE(numout,*) ' Namelist parameters for phytoplankton growth, nampisprod'
!         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!         WRITE(numout,*) '    P-I slope  for diatoms                    pislope2     =', pislope2
!         WRITE(numout,*) '    Minimum Chl/C in nanophytoplankton        chlcnm       =', chlcnm
!         WRITE(numout,*) '    Minimum Chl/C in diatoms                  chlcdm       =', chlcdm
!         WRITE(numout,*) '    Maximum Fe/C in nanophytoplankton         fecnm        =', fecnm
! <CMOC OR 05/26/2014> trimming PISCES code           WRITE(numout,*) '    Minimum Fe/C in diatoms                   fecdm        =', fecdm
      ENDIF
      !
      r1_rday   = 1._wp / rday 
! <CMOC OR 05/26/2014> trimming PISCES code        texcret   = 1._wp - excret
! <CMOC OR 05/26/2014> trimming PISCES code        texcret2  = 1._wp - excret2
! <CMOC OR 05/26/2014> trimming PISCES code        tpp       = 0._wp
      !
   END SUBROUTINE p4z_prod_init


! <CMOC OR 05/26/2014> trimming PISCES code     INTEGER FUNCTION p4z_prod_alloc()
!      !!----------------------------------------------------------------------
!      !!                     ***  ROUTINE p4z_prod_alloc  ***
!      !!----------------------------------------------------------------------
!      ALLOCATE( prmax(jpi,jpj,jpk), quotan(jpi,jpj,jpk), quotad(jpi,jpj,jpk), STAT = p4z_prod_alloc )
!      !
!      IF( p4z_prod_alloc /= 0 ) CALL ctl_warn('p4z_prod_alloc : failed to allocate arrays.')
!      !
! <CMOC OR 05/26/2014> trimming PISCES code     END FUNCTION p4z_prod_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_prod                    ! Empty routine
   END SUBROUTINE p4z_prod
#endif 

   !!======================================================================
END MODULE  p4zprod
