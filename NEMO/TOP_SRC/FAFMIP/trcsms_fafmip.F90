MODULE trcsms_fafmip
   !!======================================================================
   !!                         ***  MODULE trcsms_fafmip  ***
   !! TOP :   Main module of the fafmip tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_fafmip
   !!----------------------------------------------------------------------
   !!   'key_fafmip'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_sms_fafmip       : fafmip model main routine
   !! trc_sms_fafmip_alloc : allocate arrays specific to fafmip sms
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trdmod_oce
   USE trdmod_trc
   USE phycst
   USE fldread         ! read input fields
   USE sbc_oce         ! surface boundary condition: ocean fields

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_fafmip       ! called by trcsms.F90 module
   PUBLIC   trc_sms_fafmip_alloc ! called by trcini_fafmip.F90 module
   

   !! * Substitution
#  include "domzgr_substitute.h90"

   ! Defined HERE the arrays specific to fafmip sms and ALLOCATE them in trc_sms_fafmip_alloc

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_fafmip.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_fafmip( kt )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_fafmip  ***
      !!
      !! ** Purpose :   main routine of fafmip model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER ::   ji, jn                       ! dummy loop index
      INTEGER  ::  ierror                       ! return error code
!!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_fafmip')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_fafmip:  fafmip model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      ! Apply fluxes to redistributed heat tracer
      tra(:,:,:,jpTr) = Tr_sbc(:,:,:)
      ! Reset the sbc flux array to 0.
      Tr_sbc(:,:,:) = 0.

      ! Apply fluxes to added heat tracer
      tra(:,:,1,jpTa) = sf_fafmip(jp_hflx)%fnow(:,:,1)

      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_fafmip')
      !
   END SUBROUTINE trc_sms_fafmip


   INTEGER FUNCTION trc_sms_fafmip_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_fafmip_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to fafmip
      ! ALLOCATE( tab(...) , STAT=trc_sms_fafmip_alloc )
      trc_sms_fafmip_alloc = 0      ! set to zero if no array to be allocated
      !
      IF( trc_sms_fafmip_alloc /= 0 ) CALL ctl_warn('trc_sms_fafmip_alloc : failed to allocate arrays')
      !
   END FUNCTION trc_sms_fafmip_alloc


#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No fafmip model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_fafmip( kt )             ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_fafmip: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_fafmip
#endif

   !!======================================================================
END MODULE trcsms_fafmip
