MODULE trcopt_canbgc
   !!======================================================================
   !!                         ***  MODULE trcopt  ***
   !! TOP - PISCES : Compute the light availability in the water column
   !!======================================================================
   !! History :  1.0  !  2004     (O. Aumont) Original code
   !!            2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!            3.2  !  2009-04  (C. Ethe, G. Madec)  optimisation
   !!            3.4  !  2011-06  (O. Aumont, C. Ethe) Improve light availability of nano & diat
   !!----------------------------------------------------------------------
   !!   trc_opt       : light availability in the water column
   !!----------------------------------------------------------------------
   USE trc            ! tracer variables
   USE oce_trc        ! tracer-ocean share variables

   ! USE sms_pisces   ! Source Minus Sink of PISCES
   ! USE sms_bgcm     ! Source Minus Sink of BGCMs

   USE iom            ! I/O manager
   USE fldread        !  time interpolation
   USE prtctl     !  print control for debugging

   ! read external file
   USE sms_top_canbgc
   USE trcsrc_canbgc     ! access to surface chlorophyll array from external file
   
   USE trc_closea_canbgc ! bgc-specific closea mask

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_opt              ! called in trcsms_canoe
   PUBLIC   trc_opt_1band        ! called in trcsms_cmoc
   PUBLIC   trc_opt_alloc        !
   PUBLIC   trc_opt_init         !

   REAL(wp) ::   parlux          ! Fraction of shortwave as PAR
   REAL(wp) ::   xparsw          ! parlux/3
   REAL(wp) ::   xsi0r           ! 1. /rn_si0
   
   REAL(wp) ::   kw_cmoc, kchl_cmoc ! 1-band PAR parameters

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_par      ! structure of input par
   INTEGER , PARAMETER :: nbtimes = 366  !: maximum number of times record in a file
   INTEGER  :: ntimes_par                ! number of time steps in a file
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:  ) :: par_varsw      ! PAR fraction of shortwave
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: ekb, ekg, ekr  ! wavelength (Red-Green-Blue)

   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  nelnb           !: number of T-levels + 1 in the euphotic layer

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  par_1band

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etotb           !: par (photosynthetic available radiation)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  par_3bands      !: 

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcopt.F90 13331 2020-07-22 14:00:04Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   ! SUBROUTINE trc_opt( kt )
   SUBROUTINE trc_opt( kt, knt , Kbb, Kmm)
   ! O Riche Aug 16th 2022
   ! knt is for time splitting, not implemented 
   ! at least for now
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_opt  ***
      !!
      !! ** Purpose :   Compute the light availability in the water column
      !!              depending on the depth and the chlorophyll concentration
      !!
      !! ** Method  : - Use a tabulated function to match 
      !!                chla concentration with light att.
      !!                based on Morel et al 1981
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step 
      INTEGER, INTENT(in) ::    Kbb, Kmm  ! time level indices
      ! INTEGER, INTENT(in) ::   kt        ! ocean time step 
      ! O Riche Aug 16th 2022
      ! knt is for time splitting in PISCES
      ! see trcsms_pisces s/routine for reference, look for jnt and p4z_bio call
      ! no need right now, later maybe?
      !
      INTEGER  ::   ji, jj, jk
      INTEGER  ::   irgb
      REAL(wp) ::   zchl
      REAL(wp) ::   zc0 , zc1 , zc2, zc3, z1_dep
      REAL(wp), DIMENSION(jpi,jpj    ) :: zdepmoy, zetmp1, zetmp2
      REAL(wp), DIMENSION(jpi,jpj    ) :: zqsr100, zqsr_corr
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpar, ze0, ze1, ze2, ze3, zchl3d
      ! O Riche Sept 13th 2022
      ! add an intermediate/working array to track total chla
      ! regardless of the BGCM used.
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ztotchla
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_opt')
      
      IF( lwp ) WRITE(numout,*) 'trc_opt: PAR attenuation'
      IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~'
      IF( lwp ) WRITE(numout,*)
      IF( lwp ) CALL FLUSH(numout)

      ! O Riche Aug 16th 2022
      ! knt/time splitting 

      !     Initialisation of variables used to compute PAR
      !     -----------------------------------------------
      ze1(:,:,:) = 0._wp
      ze2(:,:,:) = 0._wp
      ze3(:,:,:) = 0._wp
      !
      !                                        !* attenuation coef. function of Chlorophyll and wavelength (Red-Green-Blue)
      !                                        !  --------------------------------------------------------
      
      ! O Riche aug 16th 2022
      ! This will be reactivated when chla is available
      ! as a tracer
      ! for now read surface chlorophyll external file and 
      ! apply an e-folding of 30 m.
      !  zchl3d(:,:,:) = tr(:,:,:,jrnch, Kbb) + tr(:,:,:,jrdch, Kbb)
      !  CALL trc_src2d( kt, js2d_chla  )
      ! O Riche Sept 13th 2022
      ! this assumes that chlorophyll can be max 2 sizes
      ! and netcdf variable names are either NCHL or DCHL
      ! and NCHL is common to both CMOC and CanOE.
      !
      ! IFs can later be replaced by cpp key activation statements
      !
      ! Failsafe case
      IF( .NOT. ln_cmoc .AND. .NOT. ln_canoe) THEN
        IF( lwp ) WRITE(numout,*) 'trc_opt: no CanBGC selected:'
        IF( lwp ) WRITE(numout,*) 'prescribed 30-m e-folding chla'
        IF( lwp ) WRITE(numout,*) 'as a place holder here.'
        IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        IF( lwp ) WRITE(numout,*)
        IF( lwp ) CALL FLUSH(numout)      
        CALL trc_src2d( kt, js2d_chla )
        DO jk = 1, jpkm1
          ztotchla(:,:,jk) = src2d_dta(:,:,js2d_chla)*exp(-gdept(:,:,jk,Kmm)/30.)
        ENDDO
      ENDIF
      !
      IF( ln_cmoc ) THEN
        IF( lwp ) WRITE(numout,*) 'trc_opt: CMOC selected:'
        IF( lwp ) WRITE(numout,*) 'using CMOC chla-a tracer'
        IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~'
        IF( lwp ) WRITE(numout,*)
        IF( lwp ) CALL FLUSH(numout)         
        IF( iom_use("NCHL") ) THEN
          IF( lwp ) WRITE(numout,*), 'trc_opt: NCHL detected by iom_use S/R.'
          IF( lwp ) WRITE(numout,*), 'trc_opt: ztotchla assigned current tr(:,:,:,jrnch, Kmm) values'
          IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'          
          IF( lwp ) WRITE(numout,*)
          ztotchla(:,:,:) = tr(:,:,:,jrnch, Kmm)        
      ENDIF
      ENDIF
      IF( ln_canoe ) THEN
        IF( lwp ) WRITE(numout,*) 'trc_opt: CanOE selected:'
        IF( lwp ) WRITE(numout,*) 'using both CanOE chla-a tracers'
        IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        IF( lwp ) WRITE(numout,*)
        IF( lwp ) CALL FLUSH(numout)       
        IF( iom_use("NCHL") ) THEN
          IF( lwp ) WRITE(numout,*), 'trc_opt: NCHL detected by iom_use S/R.'
          IF( lwp ) WRITE(numout,*), 'trc_opt: ztotchla assigned current tr(:,:,:,jrnch, Kbb) values'
          IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'          
          IF( lwp ) WRITE(numout,*)        
          ztotchla(:,:,:) = tr(:,:,:,jrnch, Kbb)
        ENDIF
        IF( iom_use("DCHL") ) THEN
          IF( lwp ) WRITE(numout,*), 'trc_opt: DCHL detected by iom_use S/R.'
          IF( lwp ) WRITE(numout,*), 'trc_opt: ztotchla assigned current tr(:,:,:,jrdch, Kbb) values'
          IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'          
          IF( lwp ) WRITE(numout,*)        
          ztotchla(:,:,:) = ztotchla(:,:,:) + tr(:,:,:,jrdch, Kbb)
        ENDIF
      ENDIF
      
      DO jk = 1, jpkm1   
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! O Riche Aug 16th 2022
               ! krgb table is indexed by chla values in mg Chl m^-3 same units as in 
               ! the surface chlorophyll file. It has 61 rows for values varying between 0.01 to 10.
               ! O Riche Aug 17th 2022
               ! chl-a in the file is already in mg Chla m^-3 (ranging between 0.01 and 1)
               ! zchl = src2d_dta(ji,jj,js2d_chla)*exp(-gdept(ji,jj,jk,Kmm)/30.)
               ! O Riche Sept 13th 2022
               ! use chla arrays instead of mockup array
               zchl = ztotchla(ji,jj,jk)
               zchl = zchl + rtrn
               zchl = zchl * tmask(ji,jj,jk)
               zchl = MIN(  10. , MAX( 0.05, zchl )  )
               irgb = NINT( 41 + 20.* LOG10( zchl ) + rtrn )
               !
               ekb(ji,jj,jk) = rkrgb(1,irgb) * e3t(ji,jj,jk,Kmm)
               ekg(ji,jj,jk) = rkrgb(2,irgb) * e3t(ji,jj,jk,Kmm)
               ekr(ji,jj,jk) = rkrgb(3,irgb) * e3t(ji,jj,jk,Kmm)

            END DO
         END DO
      END DO
      !                                        !* Photosynthetically Available Radiation (PAR)
      !                                        !  --------------------------------------
      IF( l_trcdm2dc ) THEN                    !  diurnal cycle
         !
         zqsr_corr(:,:) = qsr_mean(:,:) 
         !
         CALL trc_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100 ) 
         !
         zqsr_corr(:,:) = max(0.,qsr(:,:))
         !
         CALL trc_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3 ) 
         !
         DO jk = 1, nksr      
            etotb(:,:,jk) =  ze1(:,:,jk) + ze2(:,:,jk) + ze3(:,:,jk)
         END DO
         !
      ELSE
         !
         zqsr_corr(:,:) = max(0.,qsr(:,:))
         !
         CALL trc_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100  ) 
         !
         DO jk = 1, nksr      
            etotb (:,:,jk) = ze1(:,:,jk) + ze2(:,:,jk) + ze3(:,:,jk)
         END DO
         !
         etotb     (:,:,:) =  etotb(:,:,:) * tmask_bgc_closea(:,:,:)
      ENDIF
      ! 
      ! O Riche Aug 30th 2022
      ! 
      par_3bands(:,:,:) = etotb(:,:,:)
      ! b is to distinguish etotb from etot in p4zopt.F90
      !
      IF( ln_qsr_bio ) THEN                    !* heat flux accros w-level (used in the dynamics)
         !                                     !  ------------------------
         zqsr_corr(:,:) = max(0.,qsr(:,:))
         CALL trc_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3, pe0=ze0 )
         !
         etot3(:,:,1) =  zqsr_corr(:,:) * tmask_bgc_closea(:,:,1)
         DO jk = 2, nksr + 1
            etot3(:,:,jk) =  ( ze0(:,:,jk) + ze1(:,:,jk) + ze2(:,:,jk) + ze3(:,:,jk) ) * tmask_bgc_closea(:,:,jk)
         END DO
         !                                     !  ------------------------
      ENDIF
      !
      ! O Riche Aug 16th 2022, time splitting not implemented at least for now.
!      IF( lk_iomput .AND. knt == qnrdttrc ) THEN
!         CALL iom_put( "PAR3" , etot3(:,:,:) * tmask_bgc_closea(:,:,:) )  ! Photosynthetically Available Radiation (no band att.)
!         CALL iom_put( "PAR2BIO", par_3bands(:,:,:) * tmask_bgc_closea(:,:,:) ) ! PAR w/o the mxl averaging and w/ the diurnal cycle if any
!      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('trc_opt')
      !
   END SUBROUTINE trc_opt

   SUBROUTINE trc_opt_1band( kt , knt, Kmm )
   !
      INTEGER, INTENT(in)  :: kt, knt            ! ocean time step
      INTEGER, INTENT(in) ::    Kmm  ! time level indices
      INTEGER              :: ierr, ji, jj, jk
      REAL(wp)             :: zchl               ! temporary value of chla
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zparsw  ! PAR/SW ratio  
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zetot   ! temporary SW downwelling rad. array
                                                         ! use qsr, penetrative solar radiation
      ! O Riche Sept 13th 2022
      ! add an intermediate/working array to track total chla
      ! regardless of the BGCM used.
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ztotchla                                                         
      !
      IF( ln_timing )  CALL timing_start('trc_opt_1band')
      !
      !
      IF( lwp ) WRITE(numout,*) 'trc_opt_1band: PAR attenuation'
      IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      IF( lwp ) CALL FLUSH(numout)  
      !
      ALLOCATE( zetot(jpi,jpj,jpk), zparsw(jpi,jpj), ztotchla(jpi,jpj), STAT=ierr)
      IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'trc_opt_1band: unable to allocate zetot' )
      !
      ! O Riche Aug 16th 2022
      ! knt/time splitting 
      !
      zparsw(:,:) = parlux
      zetot(:,:,:) = 0._wp
      ! O Riche Sept 13th 2022
      ! this assumes that chlorophyll can be max 2 sizes
      ! and netcdf variable names are either NCHL or DCHL
      ! and NCHL is common to both CMOC and CanOE.
      !
      ! IFs can later be replaced by cpp key activation statements
      !
      ! Failsafe case
      IF ( .NOT. ln_cmoc .AND. .NOT. ln_canoe) THEN
        IF( lwp ) WRITE(numout,*) 'trc_opt_1band: no CanBGC selected:'
        IF( lwp ) WRITE(numout,*) 'prescribed 30-m e-folding chla'
        IF( lwp ) WRITE(numout,*) 'as a place holder here.'
        IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        IF( lwp ) WRITE(numout,*)
        IF( lwp ) CALL FLUSH(numout)       
        ! OR Jan 20th 2023
        ! Temporary changes to test PAR and PP starting off
        CALL trc_src2d( kt, js2d_chla )
        ztotchla(:,:) = src2d_dta(:,:,js2d_chla)
      ENDIF
      !
      IF ( ln_cmoc ) THEN
        IF( lwp ) WRITE(numout,*) 'trc_opt_1band: CMOC selected:'
        IF( lwp ) WRITE(numout,*) 'using CMOC chla-a tracer'
        IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        IF( lwp ) WRITE(numout,*)
        IF( lwp ) CALL FLUSH(numout)      
        IF( iom_use("NCHL") ) THEN
          IF( lwp ) WRITE(numout,*), 'trc_opt_1band: NCHL detected by iom_use S/R.'
          IF( lwp ) WRITE(numout,*), 'trc_opt_1band: ztotchla assigned current tr(:,:,1,jqnch, Kmm) values'
          IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'          
          IF( lwp ) WRITE(numout,*)
          ztotchla(:,:) = tr(:,:,1,jqnch, Kmm)  !!! OR Jan 23rd 2023 ! Only use the surface ztotchla values
          ! IF( .NOT. ln_rsttr .AND. kt <= nittrc000 + nn_dttrc) THEN
            ! ! OR Jan 20th 2023
            ! ! Temporary changes to test PAR and PP starting off
            ! CALL trc_src2d( kt, js2d_chla )
            ! !
            ! IF( lwp ) WRITE(numout,*), 'trc_opt_1band: could use surf_chla on 1st time step from rest.'
            ! IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
            ! ztotchla(:,:) = src2d_dta(:,:,js2d_chla)
          ! ENDIF
        ENDIF
      ENDIF
      IF( ln_canoe ) THEN
        IF( lwp ) WRITE(numout,*) 'trc_opt_1band: CanOE selected:'
        IF( lwp ) WRITE(numout,*) 'using both CanOE chla-a tracers'
        IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        IF( lwp ) WRITE(numout,*)
        IF( lwp ) CALL FLUSH(numout)
        IF( iom_use("NCHL") ) THEN
          IF( lwp ) WRITE(numout,*), 'trc_opt_1band: NCHL detected by iom_use S/R.'
          IF( lwp ) WRITE(numout,*), 'trc_opt_1band: ztotchla assigned current tr(:,:,1,jrnch, Kmm) values'
          IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'          
          IF( lwp ) WRITE(numout,*)
          ztotchla(:,:) = tr(:,:,1,jrnch, Kmm)  !!! OR Jan 23rd 2023 ! Only use the surface ztotchla values
        ENDIF
        IF( iom_use("DCHL") ) THEN
          IF( lwp ) WRITE(numout,*), 'trc_opt_1band: DCHL detected by iom_use S/R.'
          IF( lwp ) WRITE(numout,*), 'trc_opt_1band: ztotchla added current tr(:,:,1,jrdch, Kmm) values'
          IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'          
          IF( lwp ) WRITE(numout,*)
          ztotchla(:,:) = ztotchla(:,:) + tr(:,:,1,jrdch, Kmm)  !!! OR Jan 23rd 2023 ! Only use the surface ztotchla values   
        ENDIF
      ENDIF
      !
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            ! O Riche Aug 17th 2022
            ! chl-a in the file is already in mg Chla m^-3 (ranging between 0.01 and 1)
            ! This is temporary as zchl/1st line should be replaced by
            ! tr(ji,jj,jk,jqdch, Kmm) + tr(ji,jj,jk,jqnch, Kmm) once they are available.
            ! zchl = src2d_dta(ji,jj,js2d_chla)*exp(-gdept(ji,jj,jk,Kmm)/30._wp)
            ! O Riche Sept 13th 2022
            ! use chla arrays instead of mockup array
            zchl = ztotchla(ji,jj)  !!! OR Jan 23rd 2023 ! Only use the surface ztotchla values
            zchl = zchl + rtrn
            zchl = zchl * tmask(ji,jj,jk)
            zetot(ji,jj,jk) = max(0.,qsr(ji,jj)) * zparsw(ji,jj)     & 
            &               * exp ( - ( (kw_cmoc + kchl_cmoc * zchl * 1e6_wp) * gdept(ji,jj,jk,Kmm) ) ) 
            !        
          ENDDO
        ENDDO
      ENDDO
      !
      par_1band(:,:,:) = zetot(:,:,:)
      !
      DEALLOCATE(zetot, zparsw)
      !
!      IF( lk_iomput .AND.  knt == qnrdttrc ) THEN
!        CALL iom_put( "PARDM",                     tmask_bgc_closea(:,:,:) )  ! diagnostic : PAR with no diurnal cycle (mixed layer mean within the mxl)
!        CALL iom_put( "PAR"  ,                     tmask_bgc_closea(:,:,:) )  ! Photosynthetically Available Radiation (3-band att., mxl meam within the mxl)
!        CALL iom_put( "PAR3" ,                     tmask_bgc_closea(:,:,:) )  ! Photosynthetically Available Radiation (no band att.)

!        CALL iom_put("PAR2BIO", par_1band(:,:,:) * tmask_bgc_closea(:,:,:) ) ! PAR to use for CMOC (or CanOE)
!      ENDIF
      !
      IF( ln_timing )  CALL timing_stop('trc_opt_1band')      
         
   END SUBROUTINE trc_opt_1band

   SUBROUTINE trc_opt_par( kt,Kmm, pqsr, pe1, pe2, pe3, pe0, pqsr100 ) 
      !!----------------------------------------------------------------------
      !!                  ***  routine trc_opt_par  ***
      !!
      !! ** purpose :   compute PAR of each wavelength (Red-Green-Blue)
      !!                for a given shortwave radiation
      !!
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in)              ::   kt,Kmm                ! ocean time-step
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(in   )           ::   pqsr              ! shortwave
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout)           ::   pe1 , pe2 , pe3   ! PAR ( R-G-B)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout), OPTIONAL ::   pe0               !
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(  out), OPTIONAL ::   pqsr100           !
      !
      INTEGER    ::   ji, jj, jk     ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj) ::  zqsr   ! shortwave
      !!----------------------------------------------------------------------

      !  Real shortwave
      zqsr(:,:) = xparsw         * pqsr(:,:)
      
      !  Light at the euphotic depth 
      IF( PRESENT( pqsr100 ) )   pqsr100(:,:) = 0.01 * 3. * zqsr(:,:)

      IF( PRESENT( pe0 ) ) THEN     !  W-level
         !
         pe0(:,:,1) = pqsr(:,:) - 3. * zqsr(:,:)    !   ( 1 - 3 * alpha ) * q
         pe1(:,:,1) = zqsr(:,:)         
         pe2(:,:,1) = zqsr(:,:)
         pe3(:,:,1) = zqsr(:,:)
         !
         DO jk = 2, nksr + 1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  pe0(ji,jj,jk) = pe0(ji,jj,jk-1) * EXP( -e3t(ji,jj,jk-1,Kmm) * xsi0r )
                  pe1(ji,jj,jk) = pe1(ji,jj,jk-1) * EXP( -ekb  (ji,jj,jk-1 )        )
                  pe2(ji,jj,jk) = pe2(ji,jj,jk-1) * EXP( -ekg  (ji,jj,jk-1 )        )
                  pe3(ji,jj,jk) = pe3(ji,jj,jk-1) * EXP( -ekr  (ji,jj,jk-1 )        )
               END DO
              !
            END DO
            !
         END DO
        !
      ELSE   ! T- level
        !
        pe1(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekb(:,:,1) )
        pe2(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekg(:,:,1) )
        pe3(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekr(:,:,1) )
        !
        DO jk = 2, nksr      
           DO jj = 1, jpj
              DO ji = 1, jpi
                 pe1(ji,jj,jk) = pe1(ji,jj,jk-1) * EXP( -0.5 * ( ekb(ji,jj,jk-1) + ekb(ji,jj,jk) ) )
                 pe2(ji,jj,jk) = pe2(ji,jj,jk-1) * EXP( -0.5 * ( ekg(ji,jj,jk-1) + ekg(ji,jj,jk) ) )
                 pe3(ji,jj,jk) = pe3(ji,jj,jk-1) * EXP( -0.5 * ( ekr(ji,jj,jk-1) + ekr(ji,jj,jk) ) )
              END DO
           END DO
        END DO    
        !
      ENDIF
      ! 
   END SUBROUTINE trc_opt_par

   SUBROUTINE trc_opt_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_opt_init  ***
      !!
      !! ** Purpose :   Initialization of tabulated attenuation coef
      !!                and of the percentage of PAR in Shortwave
      !!
      !! ** Input   :   external ascii and netcdf files
      !!----------------------------------------------------------------------
      INTEGER :: numpar, ierr, ios     ! Local integer 
      !
      CHARACTER(len=100) ::  cn_dir   ! Root directory for location of ssr files
      TYPE(FLD_N)        ::   sn_par  ! informations about the fields to be read
      !
      NAMELIST/namtrc_opt/ sn_par, cn_dir, parlux,      &
      &                    kw_cmoc, kchl_cmoc                         ! 1-band PAR parameters  
      !!----------------------------------------------------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_opt_init : '
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*)
      ENDIF

      READ  ( numnat_ref, namtrc_opt, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtrc_opt in top reference namelist' )

      READ  ( numnat_cfg, namtrc_opt, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtrc_opt in top configuration namelist' )
      IF(lwm) WRITE ( numonpb, namtrc_opt )

      IF(lwp) THEN
         WRITE(numout,*) '   Namelist : namtrc_opt '
         WRITE(numout,*) '      Default value for the PAR fraction                    parlux    = ', parlux
         WRITE(numout,*) '      1-band PAR att. coefficient by seawater  (m^-1)       kw_cmoc   = ', kw_cmoc
         WRITE(numout,*) '      1-band PAR att. coeff. by chla (m^-1) (mgChl m^-3)^-1 kchl_cmoc = ', kchl_cmoc
         WRITE(numout,*)
      ENDIF
      !
      xparsw = parlux / 3.0
      xsi0r  = 1.e0 / rn_si0
      !
      ! Variable PAR at the surface of the ocean
      ! ----------------------------------------
      !
                         ekr       (:,:,:) = 0._wp
                         ekb       (:,:,:) = 0._wp
                         ekg       (:,:,:) = 0._wp
                         etotb     (:,:,:) = 0._wp
                          par_1band(:,:,:) = 0._wp
                         par_3bands(:,:,:) = 0._wp
                         
      IF( ln_qsr_bio )   etot3     (:,:,:) = 0._wp
      ! 
   END SUBROUTINE trc_opt_init

   INTEGER FUNCTION trc_opt_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_opt_alloc  ***
      !!----------------------------------------------------------------------
      !
      ALLOCATE( ekb(jpi,jpj,jpk),   ekr(jpi,jpj,jpk),          &
        &       ekg(jpi,jpj,jpk), nelnb(jpi,jpj),              &
        &     etotb(jpi,jpj,jpk),                              &
        &     par_1band(jpi,jpj,jpk),par_3bands(jpi,jpj,jpk),  &
        &     STAT= trc_opt_alloc  ) 
      !
      IF( trc_opt_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_opt_alloc : failed to allocate arrays.' )
      !
   END FUNCTION trc_opt_alloc

   !!======================================================================
   
END MODULE trcopt_canbgc
