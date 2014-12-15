MODULE trcwri_pisces
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    PISCES :   Output of PISCES tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && key_pisces && defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                           PISCES model
   !!----------------------------------------------------------------------
   !! trc_wri_pisces   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_pisces 

CONTAINS

   SUBROUTINE trc_wri_pisces
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      CHARACTER (len=20)   :: cltra
      REAL(wp)             :: zrfact
      INTEGER              :: jn
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------

      DO jn = 1, jptra

         zrfact = 1.0e+6_wp ! <CMOC OR 03/14/2014> 
! <CMOC OR 06/30/2014> Trimming code, tracers (jpnh4) !          IF( jn == jpnh4 )                                                                     zrfact = 1.0e+6 / 7.6 ! IF( jn == jpno3 .OR. jn == jpnh4 ) zrfact = 1.0e+6 / 7.6 <CMOC OR 10/28/2013> 
! <CMOC OR 06/27/2014> Trimming code, tracers (jppo4) !           IF( jn == jppo4  )                                                                    zrfact = 1.0e+6 / 122.
         IF( jn == jpoxy  )                                                                    zrfact = 1.0e+6 / 106._wp * 138._wp  ! <CMOC OR 11/30/2013> convert back to uM of O2 with the Redfield ratio
         IF( jn == jpno3 .OR. jn == jpphy .OR. jn == jpzoo .OR. jn == jppoc )                  zrfact = 1.0e+6 / 106._wp * 16._wp   ! <CMOC OR 10/29/2013> change the chemical currency from carbon to nitrogen

         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         CALL iom_put( cltra, trn(:,:,:,jn) * zrfact )

      END DO
      !
   END SUBROUTINE trc_wri_pisces

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_pisces
CONTAINS
   SUBROUTINE trc_wri_pisces                     ! Empty routine  
   END SUBROUTINE trc_wri_pisces
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcwri_pisces.F90 3160 2011-11-20 14:27:18Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_pisces
