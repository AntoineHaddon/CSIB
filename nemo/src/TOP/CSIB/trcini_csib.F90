MODULE trcini_csib
   !!======================================================================
   !!                         ***  MODULE trcini_csib  ***
   !! TOP :   initialisation of the CSIB tracers
   !!======================================================================
   !! History :        !  2007  (C. Ethe, G. Madec) Original code
   !!                  !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_ini_csib   : CSIB model initialisation
   !!----------------------------------------------------------------------
   USE par_kind   !: access wp kind
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE par_csib
   USE trcnam_csib     ! csib SMS namelist
   USE trcsms_csib

   USE dom_oce, ONLY: glamt, gphit               ! latitude/longitude for funky initiation
   USE ice , ONLY: a_i
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_csib   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcini_my_trc.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_csib( Kmm )
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_csib  ***  
      !!
      !! ** Purpose :   initialization for CSIB model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kmm     ! time level indices
      INTEGER  ::   ji, jj             ! dummy loop indices
      !
      CALL trc_nam_csib
      !
      !                       ! Allocate CSIB arrays
      IF( trc_sms_csib_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_csib: unable to allocate CSIB arrays' )

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_csib: passive tracer unit vector'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      
      IF( .NOT. ln_rsttr ) THEN
         icedia(:,:,:)=0._wp
         ! icedia(:,:,:) = tr(:,:,1,jrdia,Kmm)
         WHERE( gphit(:,:) > 85._wp )   ;   icedia(:,:,3)=1._wp
         ELSEWHERE                     ;   icedia(:,:,3)=0._wp
         END WHERE
         icedia_gca(:,:,:) = icedia(:,:,:) * a_i(:,:,:)
      ENDIF

      flushrate(:,:,:)=0._wp
      flushdia(:,:,:)=0._wp
      
      !
   END SUBROUTINE trc_ini_csib

   !!======================================================================
END MODULE trcini_csib
