MODULE trcwri_pisces
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    PISCES :   Output of PISCES tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!           CMOC1  !  2013-15 (O. Riche) code edition for consistency
   !!                                        with CMOC adaptation
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

         ! Scale CMOC tracers
         IF( jn >= jp_pcs0 .OR. jn <= jp_pcs1  ) THEN
             zrfact = 1.0e+6_wp 
         ELSE ! for all other passive tracers
             zrfact = 1.0_wp 
         ENDIF

         ! Change the chemical currency from carbon to nitrogen for certain CMOC variables
         IF( jn == jpno3 .OR. jn == jpphy .OR. jn == jpzoo .OR. jn == jppoc )  zrfact = 1.0e+6 / 106._wp * 16._wp   

         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         CALL iom_put( cltra, trn(:,:,:,jn) * zrfact )
          WRITE(*,*) cltra, zrfact

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
