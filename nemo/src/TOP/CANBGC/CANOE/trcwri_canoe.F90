MODULE trcwri_canoe
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!     trc_wri_canoe   :  outputs of concentration fields
   !!======================================================================
#if defined key_top && defined key_xios
   !!----------------------------------------------------------------------
   !! History :      !  2007  (C. Ethe, G. Madec)  Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   USE par_trc          ! passive tracers common variables
   USE trc              ! passive tracers common variables 
   USE iom              ! I/O manager
   
   USE sms_top_canbgc      ! access src2d/3d_dta
   USE trc_closea_canbgc   ! bgc-specific closea mask

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_canoe 

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcwri_canoe.F90 10069 2018-08-28 14:12:24Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_wri_canoe( Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_canoe  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  :: Kmm   ! time level indices
      CHARACTER (len=20)   :: cltra
      INTEGER              :: jn
      REAL(wp)             :: zfact
      !---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      DO jn = 1, jp_bgc+jp_canoe
      cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
      zfact = 1._wp
      IF ( cltra == 'DIC'      ) zfact = 1.e06_wp
      IF ( cltra == 'Alkalini' ) zfact = 1.e06_wp
      CALL iom_put( cltra, tr(:,:,:,jn,Kmm)*zfact ) ! O Riche June 6th 2022, manual scaling here as xml file issue not solved yet
      END DO
      !
      ! Testing trcopt diagnostics
      !CALL iom_put( "tmask", tmask(:,:,:) )
      CALL iom_put( "closea", tmask_bgc_closea(:,:,:) )
      !
      ! ---------------------------------------

   END SUBROUTINE trc_wri_canoe

#else

CONTAINS

   SUBROUTINE trc_wri_canoe
      !
   END SUBROUTINE trc_wri_canoe

#endif

END MODULE trcwri_canoe
