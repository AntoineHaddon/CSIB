MODULE trcrst_fafmip
   !!======================================================================
   !!                       ***  MODULE trcrst_fafmip  ***
   !! TOP :   create, write, read the restart files of fafmip tracer
   !!======================================================================
   !! History :   1.0  !  2010-01 (C. Ethe)  Original
   !!----------------------------------------------------------------------
#if defined key_fafmip
   !!----------------------------------------------------------------------
   !!   'key_fafmip'                                         FAFMIP tracers
   !!----------------------------------------------------------------------
   !!   trc_rst_read_fafmip   : read  restart file
   !!   trc_rst_wri_fafmip    : write restart file
   !!----------------------------------------------------------------------

   IMPLICIT NONE
   PRIVATE

   PUBLIC  trc_rst_read_fafmip   ! called by trcini.F90 module
   PUBLIC  trc_rst_wri_fafmip   ! called by trcini.F90 module

CONTAINS
   
   SUBROUTINE trc_rst_read_fafmip( knum ) 
     INTEGER, INTENT(in)  :: knum
     WRITE(*,*) 'trc_rst_read_fafmip: No specific variables to read on unit', knum
   END SUBROUTINE trc_rst_read_fafmip

   SUBROUTINE trc_rst_wri_fafmip( kt, kirst, knum )
     INTEGER, INTENT(in)  :: kt, kirst, knum
     WRITE(*,*) 'trc_rst_wri_fafmip: No specific variables to write on unit' ,knum, ' at time ', kt, kirst
   END SUBROUTINE trc_rst_wri_fafmip

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rst_read_fafmip( knum )
     INTEGER, INTENT(in)  :: knum
     WRITE(*,*) 'trc_rst_read_fafmip: You should not have seen this print! error?', knum
   END SUBROUTINE trc_rst_read_fafmip

   SUBROUTINE trc_rst_wri_fafmip( kt, kirst, knum )
     INTEGER, INTENT(in)  :: kt, kirst, knum
     WRITE(*,*) 'trc_rst_wri_fafmip: You should not have seen this print! error?', kt, kirst, knum
   END SUBROUTINE trc_rst_wri_fafmip
#endif

   !!======================================================================
END MODULE trcrst_fafmip
