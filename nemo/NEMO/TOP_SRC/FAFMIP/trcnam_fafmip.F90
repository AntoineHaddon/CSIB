MODULE trcnam_fafmip
   !!======================================================================
   !!                      ***  MODULE trcnam_fafmip  ***
   !! TOP :   initialisation of some run parameters for LOBSTER bio-model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_fafmip
   !!----------------------------------------------------------------------
   !!   'key_fafmip'   :                                       fafmip model
   !!----------------------------------------------------------------------
   !! trc_nam_fafmip      : fafmip model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
!   USE sms_fafmip      ! my namelist variable

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_fafmip   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_fafmip.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_fafmip
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_fafmip  ***  
      !!
      !! ** Purpose :   read fafmip namelist
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER :: jl, jn
      !!
!      NAMELIST/nampisbio/ nrdttrc

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_fafmip : read fafmip namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~'

      !                               ! Open the namelist file
      !                               ! ----------------------
!      CALL ctl_opn( numnatp, 'namelist_fafmip', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )

!      REWIND( numnatp )
!      READ  ( numnatp, nampisbio )

!      IF(lwp) THEN                         ! control print
!         WRITE(numout,*) ' Namelist : nampisbio'
!         WRITE(numout,*) '    frequence pour la biologie                nrdttrc   =', nrdttrc
!      ENDIF


      !
   END SUBROUTINE trc_nam_fafmip
   
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                             No fafmip
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_fafmip                      ! Empty routine
   END  SUBROUTINE  trc_nam_fafmip
#endif  

   !!======================================================================
END MODULE trcnam_fafmip
