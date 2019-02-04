MODULE trcini_fafmip
   !!======================================================================
   !!                         ***  MODULE trcini_fafmip  ***
   !! TOP :   initialisation of the fafmip tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_fafmip
   !!----------------------------------------------------------------------
   !!   'key_fafmip'                                         FAFMIP tracers
   !!----------------------------------------------------------------------
   !! trc_ini_fafmip   : fafmip model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE par_fafmip, only : jpTr, jpTa

   IMPLICIT NONE
   PRIVATE
   PUBLIC   trc_ini_fafmip   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_fafmip.F90 2787 2011-06-27 09:54:00Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_fafmip( force_init )
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_fafmip  ***  
      !!
      !! ** Purpose :   initialization for fafmip model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
      LOGICAL, OPTIONAL, INTENT(IN) :: force_init ! Allow the tracer to initialize 
      ! Local variables
      LOGICAL :: initialize

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_fafmip: initialisation of fafmip model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
     
      initialize = .NOT. ln_rsttr
      IF ( PRESENT(force_init) ) initialize = force_init
      ! Initialize the redistributed heat tracer to the temperature variable (Section 2.4, paragraph 4)
      IF( initialize ) THEN
         trn(:,:,:,jpTr) = 0. 
         tra(:,:,:,jpTr) = 0. 
         trb(:,:,:,jpTr) = tsb(:,:,:,jp_tem)
         ! Added heat variables is initialized to 0.
         trn(:,:,:,jpTa) = 0.
         tra(:,:,:,jpTa) = 0.
         trb(:,:,:,jpTa) = 0.
      ENDIF 
      !
   END SUBROUTINE trc_ini_fafmip

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No fafmip model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_fafmip             ! Empty routine
   END SUBROUTINE trc_ini_fafmip
#endif

   !!======================================================================
END MODULE trcini_fafmip
