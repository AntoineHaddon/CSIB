MODULE canoetemp
   !!======================================================================
   !!                         ***  MODULE canoetemp  ***
   !! TOP :   PISCES interpolation and computation of various accessory fields
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
   !
   !!----------------------------------------------------------------------
   !!   'key_canoe'?                                       CanOE
   !!----------------------------------------------------------------------
   !!   canoe_temp        :  interpolation and computation of various accessory fields
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_canoe       !  CanOE basic arrays and parameters
   !
   USE sms_top_canbgc  ! shared arrays across BGCM code
   !
   USE in_out_manager  ! nn_timing integer
   USE timing          ! *_timing subroutines

   IMPLICIT NONE
   PRIVATE

   PUBLIC   canoe_temp  
   PUBLIC   canoe_temp_init    ! called in trcsms_canoe.F90 / canoeprod.F90

   REAL(wp), PUBLIC ::  T0C        = 273.15_wp         !:
   REAL(wp), PUBLIC ::  Tref       = 25._wp            !:
   REAL(wp), PUBLIC ::  AEP        = -4500._wp         !:
   REAL(wp), PUBLIC ::  AEZ        = -4500._wp         !:
   REAL(wp), PUBLIC ::  AEZ2       = -4500._wp         !:
   REAL(wp), PUBLIC ::  AER        = -4500._wp         !:

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: canoetemp.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE canoe_temp( Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE canoe_temp  ***
      !!
      !! ** Purpose :   interpolation and computation of various accessory fields
      !!
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   Kmm  ! time level indices
      IF( ln_timing )  CALL timing_start('canoe_temp')
      !
      ! Computation of phyto and zoo metabolic rate
      ! -------------------------------------------
      tgfuncp0(:,:,:)  = EXP(AEP*( 1./(ts(:,:,:,jp_tem,Kmm)+T0C)-1./(Tref+T0C)))
      tgfuncz0(:,:,:)  = EXP(AEZ*( 1./(ts(:,:,:,jp_tem,Kmm)+T0C)-1./(Tref+T0C)))
      tgfuncz20(:,:,:) = EXP(AEZ2*(1./(ts(:,:,:,jp_tem,Kmm)+T0C)-1./(Tref+T0C)))
      tgfuncr0(:,:,:)  = EXP(AER*( 1./(ts(:,:,:,jp_tem,Kmm)+T0C)-1./(Tref+T0C)))
      !
      IF( ln_timing )  CALL timing_stop('canoe_temp')
      !
   END SUBROUTINE canoe_temp

   SUBROUTINE canoe_temp_init
      !
      INTEGER ::   ios, ierr                     ! Local integers
      NAMELIST/namcanoetf/ AEP, AEZ, AEZ2, AER
      !!----------------------------------------------------------------------

      REWIND( numnatp_refb )              ! Namelist namcanoetf in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcanoetf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcanoetf in reference namelist_canoe' )
      REWIND( numnatp_cfgb )              ! Namelist namcanoetf in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcanoetf, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcanoetf in configuration namelist_canoe' )

      IF(lwm) WRITE( numonpb, namcanoetf ) 
      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' canoe_temp_init:'
         WRITE(numout,*) ' Namelist parameters for temperature, namcanoetf'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Activation energy for phytoplankton      AEP          =', AEP
         WRITE(numout,*) '    Activation energy for microzooplankton   AEZ          =', AEZ
         WRITE(numout,*) '    Activation energy for mesozooplankton    AEZ2         =', AEZ2
         WRITE(numout,*) '    Activation energy for remineralization   AER          =', AER
      ENDIF
      !
      ! Allocate mem to temperature dependency arrays
      ALLOCATE( tgfuncp0(jpi,jpj,jpk)  , tgfuncz0(jpi,jpj,jpk) ,     &
         &      tgfuncz20(jpi,jpj,jpk) , tgfuncr0(jpi,jpj,jpk) , STAT=ierr )
      !
      IF( ierr /= 0 )   CALL ctl_stop('STOP', 'canoe_temp_init : unable to allocate temp. dependency arrays')
      !
   END SUBROUTINE canoe_temp_init

   !!======================================================================
END MODULE  canoetemp

