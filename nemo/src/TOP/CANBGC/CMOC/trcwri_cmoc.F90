MODULE trcwri_cmoc
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!     trc_wri_cmoc   :  outputs of concentration fields
   !!======================================================================
#if defined key_top && defined key_xios
   !!----------------------------------------------------------------------
   !! History :      !  2007  (C. Ethe, G. Madec)  Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   USE par_trc         ! passive tracers common variables
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   
   USE sms_top_canbgc             ! access src2d/3d_dta
   USE trc_closea_canbgc          ! bgc-specific closea mask
   USE sms_cmoc, ONLY: ncrr_cmoc  !  CMOC specific parameters declaration   

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_cmoc 

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcwri_cmoc.F90 10069 2018-08-28 14:12:24Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_wri_cmoc( Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_cmoc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  :: Kmm   ! time level indices
      CHARACTER (len=20)   :: cltra
      INTEGER              :: jn
      REAL(wp)             :: zfact
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      DO jn = 1, jp_bgc+jp_cmoc
      cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
      zfact = 1._wp
      IF ( cltra == 'DIC'      ) zfact = 1.e06_wp
      IF ( cltra == 'Alkalini' ) zfact = 1.e06_wp
      IF ( cltra == 'O2' )       zfact = 1.e06_wp
      IF ( cltra == 'NO3')       zfact = 1.e06_wp * ncrr_cmoc
      IF ( cltra == 'POC')       zfact = 1.e06_wp
      IF ( cltra == 'PHY')       zfact = 1.e06_wp
      IF ( cltra == 'ZOO')       zfact = 1.e06_wp
      IF ( cltra == 'NCHL')      zfact = 1.e06_wp
      CALL iom_put( cltra, tr(:,:,:,jn,Kmm)*zfact ) ! O Riche June 6th 2022, manual scaling here as xml file issue not solved yet
      END DO
      !
      ! Testing trcopt diagnostics
      CALL iom_put( "surf_chla", src2d_dta(:,:,js2d_chla))
      !CALL iom_put( "tmask", tmask(:,:,:) )
      CALL iom_put( "closea", tmask_bgc_closea(:,:,:) )
      !
      ! ---------------------------------------

   END SUBROUTINE trc_wri_cmoc

#else

CONTAINS

   SUBROUTINE trc_wri_cmoc
      !
   END SUBROUTINE trc_wri_cmoc

#endif

END MODULE trcwri_cmoc
