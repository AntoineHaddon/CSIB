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
   USE prtctl_trc     !  print control for debugging

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

   LOGICAL  ::   ln_varpar       ! boolean for variable PAR fraction
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
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heupb           !: euphotic layer depth
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heup_01b        !: absolute euphotic layer depth

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  par_1band
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  emoyb           !: averaged PAR iver the mixed layer
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etot_ndcyb      !: PAR over 24h in case of diurnal cycle

   REAL(wp), PUBLIC, TARGET, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etotb           !: par (photosynthetic available radiation)
   REAL(wp), PUBLIC, POINTER            , SAVE, DIMENSION(:,:,:) ::  par_3bands      !: pointer/alias for etotb

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcopt.F90 13331 2020-07-22 14:00:04Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_opt( kt )
   ! SUBROUTINE trc_opt( kt, knt )
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
      ! INTEGER, INTENT(in) ::   kt, knt   ! ocean time step 
      INTEGER, INTENT(in) ::   kt        ! ocean time step 
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
      ! IF( knt == 1 .AND. ln_varpar )   CALL trc_opt_sbc( kt )
      IF( ln_varpar )   CALL trc_opt_sbc( kt )

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
      !  zchl3d(:,:,:) = trb(:,:,:,jqnch) + trb(:,:,:,jqdch)
      !  CALL trc_src2d( kt, js2d_chla  )
      ! O Riche Sept 13th 2022
      ! this assumes that chlorophyll can be max 2 sizes
      ! and netcdf variable names are either NCHL or DCHL
      ! and NCHL is common to both CMOC and CanOE.
      !
      ! IFs can later be replaced by cpp key activation statements
      !
      ! Failsafe case
      CALL trc_src2d( kt, js2d_chla )
      DO jk = 1, jpkm1
        ztotchla(:,:,jk) = src2d_dta(:,:,js2d_chla)*exp(-gdept_n(:,:,jk)/30.)
      ENDDO
      !
      IF( iom_use("NCHL") )  ztotchla(:,:,:) = trn(:,:,:,jqnch)
      ! IF( iom_use("DCHL") )  ztotchla(:,:,:) = ztotchla(:,:,:) + trn(:,:,:,jqdch)    
      !
      DO jk = 1, jpkm1   
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! O Riche Aug 16th 2022
               ! krgb table is indexed by chla values in mg Chl m^-3 same units as in 
               ! the surface chlorophyll file. It has 61 rows for values varying between 0.01 to 10.
               ! O Riche Aug 17th 2022
               ! chl-a in the file is already in mg Chla m^-3 (ranging between 0.01 and 1)
               ! zchl = src2d_dta(ji,jj,js2d_chla)*exp(-gdept_n(ji,jj,jk)/30.)
               ! O Riche Sept 13th 2022
               ! use chla arrays instead of mockup array
               zchl = ztotchla(ji,jj,jk)
               zchl = zchl + rtrn
               ! the exponential is an ad-hoc e-folding as mentioned above
               ! O Riche Aug 30th 2022
               ! trn(ji,jj,jk,jqdch)+trn(ji,jj,jk,jqnch) will replace src2d_dta(ji,jj,jk,js2d_chla)
               ! once they are available. 
               zchl = zchl * tmask(ji,jj,jk)
               zchl = MIN(  10. , MAX( 0.05, zchl )  )
               irgb = NINT( 41 + 20.* LOG10( zchl ) + rtrn )
               !
               ekb(ji,jj,jk) = rkrgb(1,irgb) * e3t_n(ji,jj,jk)
               ekg(ji,jj,jk) = rkrgb(2,irgb) * e3t_n(ji,jj,jk)
               ekr(ji,jj,jk) = rkrgb(3,irgb) * e3t_n(ji,jj,jk)

            END DO
         END DO
      END DO
      !                                        !* Photosynthetically Available Radiation (PAR)
      !                                        !  --------------------------------------
      IF( l_trcdm2dc ) THEN                    !  diurnal cycle
         !
         zqsr_corr(:,:) = qsr_mean(:,:) / ( 1.-fr_i(:,:) + rtrn )
         !
         CALL trc_opt_par( kt, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100 ) 
         !
         DO jk = 1, nksr      
            etot_ndcyb(:,:,jk) = ze1(:,:,jk) + ze2(:,:,jk) + ze3(:,:,jk)
         END DO
         !
         zqsr_corr(:,:) = qsr(:,:) / ( 1.-fr_i(:,:) + rtrn )
         !
         CALL trc_opt_par( kt, zqsr_corr, ze1, ze2, ze3 ) 
         !
         DO jk = 1, nksr      
            etotb(:,:,jk) =  ze1(:,:,jk) + ze2(:,:,jk) + ze3(:,:,jk)
         END DO
         !
      ELSE
         !
         zqsr_corr(:,:) = qsr(:,:) / ( 1.-fr_i(:,:) + rtrn )
         !
         CALL trc_opt_par( kt, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100  ) 
         !
         DO jk = 1, nksr      
            etotb (:,:,jk) = ze1(:,:,jk) + ze2(:,:,jk) + ze3(:,:,jk)
         END DO
         !
         etotb     (:,:,:) =  etotb(:,:,:) * tmask_bgc_closea(:,:,:)
         etot_ndcyb(:,:,:) =  etotb(:,:,:) 
      ENDIF
      ! 
      ! O Riche Aug 30th 2022
      ! 
      par_3bands(:,:,:) = etotb(:,:,:)
      ! b is to distinguish etotb from etot in p4zopt.F90
      !
      IF( ln_qsr_bio ) THEN                    !* heat flux accros w-level (used in the dynamics)
         !                                     !  ------------------------
         CALL trc_opt_par( kt, qsr, ze1, ze2, ze3, pe0=ze0 )
         !
         etot3(:,:,1) =  qsr(:,:) * tmask_bgc_closea(:,:,1)
         DO jk = 2, nksr + 1
            etot3(:,:,jk) =  ( ze0(:,:,jk) + ze1(:,:,jk) + ze2(:,:,jk) + ze3(:,:,jk) ) * tmask_bgc_closea(:,:,jk)
         END DO
         !                                     !  ------------------------
      ENDIF
      !                                        !* Euphotic depth and level
      DO jk = 2, nksr
         DO jj = 1, jpj
           DO ji = 1, jpi
              IF( etot_ndcyb(ji,jj,jk) * tmask(ji,jj,jk) >=  zqsr100(ji,jj) )  THEN
                 nelnb(ji,jj) = jk+1                    ! Euphotic level : 1rst T-level strictly below Euphotic layer
                 !                                      ! nb: ensure the compatibility with nmld_trc definition in trd_mld_trc_zint
                 heupb(ji,jj) = gdepw_n(ji,jj,jk+1)     ! Euphotic layer depth
              ENDIF
              IF( etot_ndcyb(ji,jj,jk) * tmask(ji,jj,jk) >= 0.50 )  THEN
                 heup_01b(ji,jj) = gdepw_n(ji,jj,jk+1)  ! Euphotic layer depth (light level definition)
              ENDIF
           END DO
        END DO
      END DO
      !
      !
      heupb   (:,:) = MIN( 300., heupb   (:,:) ) 
      heup_01b(:,:) = MIN( 300., heup_01b(:,:) ) 
      !                                        !* mean light over the mixed layer
      zdepmoy(:,:)   = 0.e0                    !  -------------------------------
      zetmp1 (:,:)   = 0.e0
      zetmp2 (:,:)   = 0.e0

      DO jk = 1, nksr
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( gdepw_n(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                  zetmp1 (ji,jj) = zetmp1 (ji,jj) + etotb     (ji,jj,jk) * e3t_n(ji,jj,jk) ! remineralisation (?OR Aug 30th 2022)
                  zetmp2 (ji,jj) = zetmp2 (ji,jj) + etot_ndcyb(ji,jj,jk) * e3t_n(ji,jj,jk) ! production
                  zdepmoy(ji,jj) = zdepmoy(ji,jj) +                        e3t_n(ji,jj,jk)
               ENDIF
            END DO
         END DO
      END DO
      !
      emoyb(:,:,:) = etotb(:,:,:)       ! remineralisation (?OR Aug 30th 2022)
      zpar(:,:,:)  = etot_ndcyb(:,:,:)  ! diagnostic : PAR with no diurnal cycle
      !
      DO jk = 1, nksr
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( gdepw_n(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                  z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
                  emoyb (ji,jj,jk) = zetmp1(ji,jj) * z1_dep
                  zpar (ji,jj,jk)  = zetmp2(ji,jj) * z1_dep
               ENDIF
            END DO
         END DO
      END DO
      !
      ! O Riche Aug 16th 2022, time splitting not implemented at least for now.
      IF( lk_iomput ) THEN ! .AND.  knt == nrdttrc ) THEN
         CALL iom_put( "Heup" , heupb(:,:  ) * tmask_bgc_closea(:,:,1) )  ! euphotic layer depth
         CALL iom_put( "PARDM", zpar(:,:,: ) * tmask_bgc_closea(:,:,:) )  ! diagnostic : PAR with no diurnal cycle (mixed layer mean within the mxl)
         CALL iom_put( "PAR"  , emoyb(:,:,:) * tmask_bgc_closea(:,:,:) )  ! Photosynthetically Available Radiation (3-band att., mxl meam within the mxl)
         CALL iom_put( "PAR3" , etot3(:,:,:) * tmask_bgc_closea(:,:,:) )  ! Photosynthetically Available Radiation (no band att.)
         CALL iom_put( "PAR2BIO", par_3bands(:,:,:) * tmask_bgc_closea(:,:,:) ) ! PAR w/o the mxl averaging and w/ the diurnal cycle if any
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('trc_opt')
      !
   END SUBROUTINE trc_opt

   SUBROUTINE trc_opt_1band( kt )
   !
      INTEGER, INTENT(in)  :: kt                 ! ocean time step
      INTEGER              :: ierr, ji, jj, jk
      REAL(wp)             :: zchl               ! temporary value of chla
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zparsw  ! PAR/SW ratio  
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zetot   ! temporary SW downwelling rad. array
                                                         ! use qsr, penetrative solar radiation
      ! O Riche Sept 13th 2022
      ! add an intermediate/working array to track total chla
      ! regardless of the BGCM used.
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztotchla                                                         
      !
      IF( ln_timing )  CALL timing_start('trc_opt_1band')
      !
      !
      IF( lwp ) WRITE(numout,*) 'trc_opt_1band: PAR attenuation'
      IF( lwp ) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      IF( lwp ) CALL FLUSH(numout)  
      !
      ALLOCATE( zetot(jpi,jpj,jpk), zparsw(jpi,jpj), ztotchla(jpi,jpj,jpk), STAT=ierr)
      IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'trc_opt_1band: unable to allocate zetot' )
      !
      ! O Riche Aug 16th 2022
      ! knt/time splitting 
      ! IF( knt == 1 .AND. ln_varpar )   CALL trc_opt_sbc( kt )
      IF( ln_varpar )   CALL trc_opt_sbc( kt )      !
      !
      IF( ln_varpar ) THEN  ;  zparsw(:,:) = par_varsw(:,:) * 3._wp ! (as it meant for the 3-band PAR)
      ELSE                  ;  zparsw(:,:) = .43_wp
      ENDIF
      zetot(:,:,:) = 0._wp
      ! O Riche Sept 13th 2022
      ! this assumes that chlorophyll can be max 2 sizes
      ! and netcdf variable names are either NCHL or DCHL
      ! and NCHL is common to both CMOC and CanOE.
      !
      ! IFs can later be replaced by cpp key activation statements
      !
      ! Failsafe case
      CALL trc_src2d( kt, js2d_chla )
      DO jk = 1, jpkm1
        ztotchla(:,:,jk) = src2d_dta(:,:,js2d_chla)*exp(-gdept_n(:,:,jk)/30.)
      ENDDO
      !
      IF( iom_use("NCHL") )  ztotchla(:,:,:) = trn(:,:,:,jqnch)
      ! IF( iom_use("DCHL") )  ztotchla(:,:,:) = ztotchla(:,:,:) + trn(:,:,:,jqdch)      
      !
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            ! O Riche Aug 17th 2022
            ! chl-a in the file is already in mg Chla m^-3 (ranging between 0.01 and 1)
            ! This is temporary as zchl/1st line should be replaced by
            ! trn(ji,jj,jk,jqdch) + trn(ji,jj,jk,jqnch) once they are available.
            ! zchl = src2d_dta(ji,jj,js2d_chla)*exp(-gdept_n(ji,jj,jk)/30._wp)
            ! O Riche Sept 13th 2022
            ! use chla arrays instead of mockup array
            zchl = ztotchla(ji,jj,jk)
            zchl = zchl + rtrn
            zchl = zchl * tmask(ji,jj,jk)
            zetot(ji,jj,jk) = qsr(ji,jj) * zparsw(ji,jj)     & 
            &               * exp ( - ( (kw_cmoc + kchl_cmoc * zchl * 1e6_wp) * gdept_n(ji,jj,jk) ) ) 
            !        
          ENDDO
        ENDDO
      ENDDO
      !
      par_1band(:,:,:) = zetot(:,:,:)
      IF( lk_iomput )  CALL iom_put("PAR2BIO", par_1band(:,:,:) * tmask_bgc_closea(:,:,:) ) ! PAR to use for CMOC (or CanOE)
      !
      IF( ln_timing )  CALL timing_stop('trc_opt_1band')      
         
   END SUBROUTINE trc_opt_1band

   SUBROUTINE trc_opt_par( kt, pqsr, pe1, pe2, pe3, pe0, pqsr100 ) 
      !!----------------------------------------------------------------------
      !!                  ***  routine trc_opt_par  ***
      !!
      !! ** purpose :   compute PAR of each wavelength (Red-Green-Blue)
      !!                for a given shortwave radiation
      !!
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in)              ::   kt                ! ocean time-step
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(in   )           ::   pqsr              ! shortwave
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout)           ::   pe1 , pe2 , pe3   ! PAR ( R-G-B)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout), OPTIONAL ::   pe0               !
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(  out), OPTIONAL ::   pqsr100           !
      !
      INTEGER    ::   ji, jj, jk     ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj) ::  zqsr   ! shortwave
      !!----------------------------------------------------------------------

      !  Real shortwave
      IF( ln_varpar ) THEN  ;  zqsr(:,:) = par_varsw(:,:) * pqsr(:,:)
      ELSE                  ;  zqsr(:,:) = xparsw         * pqsr(:,:)
      ENDIF
      
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
                  pe0(ji,jj,jk) = pe0(ji,jj,jk-1) * EXP( -e3t_n(ji,jj,jk-1) * xsi0r )
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


   SUBROUTINE trc_opt_sbc( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine trc_opt_sbc  ***
      !!
      !! ** purpose :   read and interpolate the variable PAR fraction
      !!                of shortwave radiation
      !!
      !! ** method  :   read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('trc_optsbc')
      !
      ! Compute par_varsw at nit000 or only if there is more than 1 time record in par coefficient file
      IF( ln_varpar ) THEN
         IF( kt == nit000 .OR. ( kt /= nit000 .AND. ntimes_par > 1 ) ) THEN
            CALL fld_read( kt, 1, sf_par )
            par_varsw(:,:) = ( sf_par(1)%fnow(:,:,1) ) / 3.0
         ENDIF
      ENDIF
      !
      IF( ln_timing )  CALL timing_stop('trc_optsbc')
      !
   END SUBROUTINE trc_opt_sbc


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
      NAMELIST/namtrc_opt/ sn_par, cn_dir, ln_varpar, parlux,      &
      &                    kw_cmoc, kchl_cmoc                      ! 1-band PAR parameters  
      !!----------------------------------------------------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_opt_init : '
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*)
      ENDIF

      ! REWIND( numnat_ref )
      ! READ  ( numnat_ref, namtrc_opt, IOSTAT = ios, ERR = 901)
! 901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtrc_opt in reference namelist' )

      REWIND( numnat_cfg )
      READ  ( numnat_cfg, namtrc_opt, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtrc_opt in configuration namelist' )
      IF(lwm) WRITE ( numonpb, namtrc_opt )

      IF(lwp) THEN
         WRITE(numout,*) '   Namelist : namtrc_opt '
         WRITE(numout,*) '      PAR as a variable fraction of SW     ln_varpar      = ', ln_varpar
         WRITE(numout,*) '      Default value for the PAR fraction   parlux         = ', parlux
         WRITE(numout,*)
      ENDIF
      !
      xparsw = parlux / 3.0
      xsi0r  = 1.e0 / rn_si0
      !
      ! Variable PAR at the surface of the ocean
      ! ----------------------------------------
      IF( ln_varpar ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   initialize variable par fraction (ln_varpar=T)'
         !
         ALLOCATE( par_varsw(jpi,jpj) )
         !
         ALLOCATE( sf_par(1), STAT=ierr )           !* allocate and fill sf_sst (forcing structure) with sn_sst
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'trc_opt_init: unable to allocate sf_par structure' )
         !
         CALL fld_fill( sf_par, (/ sn_par /), cn_dir, 'trc_opt_init', 'Variable PAR fraction ', 'namcanopt' )
                                   ALLOCATE( sf_par(1)%fnow(jpi,jpj,1)   )
         IF( sn_par%ln_tint )      ALLOCATE( sf_par(1)%fdta(jpi,jpj,1,2) )

         CALL iom_open (  TRIM( sn_par%clname ) , numpar )
         ntimes_par = iom_getszuld( numpar )   ! get number of record in file
      ENDIF
      !
                         ekr       (:,:,:) = 0._wp
                         ekb       (:,:,:) = 0._wp
                         ekg       (:,:,:) = 0._wp
                         etotb     (:,:,:) = 0._wp
                         etot_ndcyb(:,:,:) = 0._wp
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
        &     heupb(jpi,jpj),  heup_01b(jpi,jpj),              & 
        &     etotb(jpi,jpj,jpk),                              &
        &     etot_ndcyb(jpi,jpj,jpk), emoyb(jpi,jpj,jpk),     &              
        &  par_1band(jpi,jpj,jpk),par_3bands(jpi,jpj,jpk),   STAT= trc_opt_alloc  ) 
      !
      IF( trc_opt_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_opt_alloc : failed to allocate arrays.' )
      !
   END FUNCTION trc_opt_alloc

   !!======================================================================
   
END MODULE trcopt_canbgc
