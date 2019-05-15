MODULE trcsms_fafmip
   !!======================================================================
   !!                         ***  MODULE trcsms_fafmip  ***
   !! TOP :   Main module of the fafmip tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_fafmip
   !!----------------------------------------------------------------------
   !!   'key_fafmip'                                         FAFMIP tracers
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
   USE sbcfaf, only : sf_fafmip, Tr_sbc, jp_hflx ! Surface boundary conditions for fafmip
   USE par_fafmip, only : jpTr, jpTa

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_fafmip       ! called by trcsms.F90 module

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
      INTEGER ::   ji, jj, jn                       ! dummy loop index
      INTEGER  ::  ierror                       ! return error code
!!----------------------------------------------------------------------
      REAL(wp) :: zfact, z1_e3t
!
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_fafmip')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_fafmip:  fafmip model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      !                                        ***********************************
      ! These next blocks are needed to deal with the fact that passive tracers
      ! and active tracers may be timestepped differently because PISCES assumes
      ! Euler timestepping (and thus need to multiple the tendency by 2)
# if defined key_pisces
     zfact = 2.
# else
     zfact = 1.
# endif

      ! Apply fluxes to redistributed heat tracer
      tra(:,:,:,jpTr) = zfact*Tr_sbc(:,:,:)  ! This is the heat flux from actual model components
      ! Reset the sbc flux array to 0.
      Tr_sbc(:,:,:) = 0.

      ! Apply fluxes to added heat tracer
      DO jj = 1, jpj
         DO ji = 1, jpi
            z1_e3t = zfact / fse3t(ji,jj,1)
            tra(ji,jj,1,jpTa) = (sf_fafmip(jp_hflx)%fnow(ji,jj,1) * ro0cpr) * z1_e3t
         ENDDO
      ENDDO

      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_fafmip')
      !
   END SUBROUTINE trc_sms_fafmip

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
