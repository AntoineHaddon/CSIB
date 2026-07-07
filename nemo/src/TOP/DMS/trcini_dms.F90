MODULE trcini_dms
   !!======================================================================
   !!                         ***  MODULE trcini_dms  ***
   !! TOP :   initialisation of the DMS  tracers
   !!======================================================================
   !! History :      ! 2026 (T. Sou, A. Haddon) Ocean DMS 
   !!----------------------------------------------------------------------
   !! trc_ini_dms   : DMS model initialisation
   !!----------------------------------------------------------------------
   USE par_kind   !: access wp kind
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE par_dms
   USE trcnam_dms     !  DMS namelist
   USE trcsms_dms

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_dms   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_dms( Kmm )
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_dms  ***  
      !!
      !! ** Purpose :   initialization for DMS model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kmm  ! time level indices
      INTEGER  ::   jp, jn          ! dummy loop indices
      CHARACTER (len=20)   :: cltra

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_dms:'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      !
      ! Allocate sms_DMS arrays
      IF( trc_sms_dms_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_dms: unable to allocate OCEAN DMS arrays' )

      CALL trc_nam_dms

      DO jp = 1,jp_dmsoce
         jn = jp + jp_bgc+jp_canoe
         write(numout,*) ctrcnm(jn)
         cltra = TRIM( ctrcnm(jn) )
         IF( cltra == 'dmspd'      ) jrdmspd = jn      !: Dissolved dimethylsulfoniopropionate
         IF( cltra == 'dms'      )   jrdms = jn      !: Dimethylsulfide
      ENDDO

      
      IF( .NOT. ln_rsttr ) THEN
         tr(:,:,:,jrdmspd,Kmm) = 0._wp
         tr(:,:,:,jrdms,Kmm) = 0._wp
      ENDIF

      zdmsflx(:,:) = 0._wp
      smsdmspd(:,:) = 0._wp
      smsdms(:,:) = 0._wp





      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_dms: passive tracer unit vector'
      !
   END SUBROUTINE trc_ini_dms

   !!======================================================================
END MODULE trcini_dms



