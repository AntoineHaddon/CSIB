MODULE p4zint
   !!======================================================================
   !!                         ***  MODULE p4zint  ***
   !! TOP :   PISCES interpolation and computation of various accessory fields
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_int        :  interpolation and computation of various accessory fields
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_int  
   PUBLIC   p4z_int_init    ! called in trcsms_pisces.F90

   REAL(wp), PUBLIC ::  T0C        = 273.15_wp         !:
   REAL(wp), PUBLIC ::  Tref       = 25._wp            !:
   REAL(wp), PUBLIC ::  AEP        = -4500._wp         !:
   REAL(wp), PUBLIC ::  AEZ        = -4500._wp         !:
   REAL(wp), PUBLIC ::  AEZ2       = -4500._wp         !:
   REAL(wp), PUBLIC ::  AER        = -4500._wp         !:

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zint.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_int
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_int  ***
      !!
      !! ** Purpose :   interpolation and computation of various accessory fields
      !!
      !!---------------------------------------------------------------------
      INTEGER  ::   ji, jj                 ! dummy loop indices
      REAL(wp) ::   zvar                   ! local variable
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_int')
      !
      ! Computation of phyto and zoo metabolic rate
      ! -------------------------------------------
      tgfuncp(:,:,:) = EXP(AEP*(1./(tsn(:,:,:,jp_tem)+T0C)-1./(Tref+T0C)))
      tgfuncz(:,:,:) = EXP(AEZ*(1./(tsn(:,:,:,jp_tem)+T0C)-1./(Tref+T0C)))
      tgfuncz2(:,:,:) = EXP(AEZ2*(1./(tsn(:,:,:,jp_tem)+T0C)-1./(Tref+T0C)))
      tgfuncr(:,:,:) = EXP(AER*(1./(tsn(:,:,:,jp_tem)+T0C)-1./(Tref+T0C)))
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_int')
      !
   END SUBROUTINE p4z_int

   SUBROUTINE p4z_int_init
      !
      NAMELIST/nampistf/ AEP, AEZ, AEZ2, AER
      !!----------------------------------------------------------------------

      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampistf )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for temperature, nampistf'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Activation energy for phytoplankton      AEP          =', AEP
         WRITE(numout,*) '    Activation energy for microzooplankton   AEZ          =', AEZ
         WRITE(numout,*) '    Activation energy for mesozooplankton    AEZ2         =', AEZ2
         WRITE(numout,*) '    Activation energy for remineralization   AER          =', AER
      ENDIF
      !
   END SUBROUTINE p4z_int_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_int                   ! Empty routine
      WRITE(*,*) 'p4z_int: You should not have seen this print! error?'
   END SUBROUTINE p4z_int
#endif 

   !!======================================================================
END MODULE  p4zint
