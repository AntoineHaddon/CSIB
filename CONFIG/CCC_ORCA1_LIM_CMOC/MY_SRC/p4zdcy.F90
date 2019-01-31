MODULE p4zdcy
   !!======================================================================
   !!                         ***  MODULE p4zdcy  ***
   !! TOP :   Calculate CMOC Radioactive decay
   !!======================================================================
   !! History :   
   !!----------------------------------------------------------------------
#if defined key_pisces
   USE dom_oce, only : nyear_len 
   USE oce_trc                      !  shared variables between ocean and passive tracers 
   USE trc                          !  passive tracers common variables
   USE sms_pisces                   !  PISCES Source Minus Sink variables
   USE par_pisces, only : jpdrc
   IMPLICIT NONE
   PRIVATE

   REAL(wp), SAVE :: lambda_c14 

   PUBLIC p4z_dcy

   CONTAINS

   !> Calculate the tendency due to radioactive decay
   SUBROUTINE p4z_dcy( kt )
   INTEGER, INTENT(in   ) :: kt ! ocean time step

   INTEGER :: ji, jj, jk
   IF ( kt == nit000 ) THEN
      lambda_c14 = LOG(2.)/(5700.*365.*86400.) ! 14C with 5700yr half-life converted to seconds
      WRITE(numout,*) '   Radioactive decay parameters'
      WRITE(numout,*) '     Carbon-14 decay constant:        lambda_c14 = ', lambda_c14
   ENDIF

   DO jk = 1,jpk
      DO jj = 1,jpj
         DO ji = 1,jpi
            tra(ji,jj,jk,jpdrc) = tra(ji,jj,jk,jpdrc) - lambda_c14*trn(ji,jj,jk,jpdrc) 
         ENDDO
      ENDDO
   ENDDO

   END SUBROUTINE
#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_dcy                   ! Empty routine
      WRITE(*,*) 'p4z_dcy: You should not have seen this print! error?'
   END SUBROUTINE p4z_dcy
#endif 

   !!======================================================================
END MODULE  p4zdcy
