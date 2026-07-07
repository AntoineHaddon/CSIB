MODULE trcwri_dms
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!     trc_wri_dms   :  outputs of concentration fields for OCEAN DMS
   !!======================================================================
#if defined key_top && defined key_xios
   !!----------------------------------------------------------------------
   !! History :      ! 2026 (T. Sou, A. Haddon) Ocean DMS 
   !!----------------------------------------------------------------------
   USE par_trc         ! passive tracers common variables
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   USE dom_oce,  ONLY : e3t_0                !: t- vert. scale factor [m]

   USE trcsms_dms     ! zdmsflx
   USE par_dms

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_dms 

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_wri_dms( Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_dms  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  :: Kmm   ! time level indices
      CHARACTER (len=20)   :: cltra
      INTEGER              :: jn,jp
      REAL(wp)             :: zfact
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------

        CALL iom_put( "dmsflux" , zdmsflx(:,:))
        CALL iom_put( "smsdmspd", smsdmspd(:,:)*e3t_0(:,:,1) )
        CALL iom_put( "smsdms", smsdms(:,:)*e3t_0(:,:,1) )

        CALL iom_put( "sfcdmspd", tr(:,:,1,jrdmspd,Kmm) )
        CALL iom_put( "sfcdms", tr(:,:,1,jrdms,Kmm) )

        CALL iom_put( "dms", tr(:,:,:,jrdms,Kmm) )
        CALL iom_put( "dmspd", tr(:,:,:,jrdmspd,Kmm) )

      !
   END SUBROUTINE trc_wri_dms

#else

CONTAINS

   SUBROUTINE trc_wri_dms
      !
   END SUBROUTINE trc_wri_dms

#endif

END MODULE trcwri_dms
