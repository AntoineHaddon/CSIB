MODULE trcwri_csib
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!     trc_wri_csib   :  outputs of concentration fields
   !!======================================================================
#if defined key_top && defined key_xios
   !!----------------------------------------------------------------------
   !! History :      !  2025 (A. Haddon) Original code
   !!----------------------------------------------------------------------
   USE par_trc         ! passive tracers common variables
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager

   USE par_csib  
   USE trcsms_csib  !CSIB variables
   USE ice              ! ice variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_csib

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcwri_my_trc.F90 14239 2020-12-23 08:57:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_wri_csib( Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  :: Kmm   ! time level indices
      ! CHARACTER (len=20)   :: cltra
      INTEGER              :: jn
      REAL(wp), DIMENSION(jpi,jpj) ::   z1_at_i
      !!---------------------------------------------------------------------
 
      ! 1 / sea ice concentration ; to compute average over sea ice thickness categories
      WHERE( at_i(:,:) > 1.e-4_wp )   ;   z1_at_i(:,:) = 1._wp / at_i(:,:)
      ELSEWHERE                       ;   z1_at_i(:,:) = 0._wp
      END WHERE

      ! write the tracer concentrations in the file
      ! ---------------------------------------
      DO jn=1,jp_csib
         CALL iom_put( TRIM(icetrcnm(jn))//'_cat'        , icetra(:,:,:,jn) )
         CALL iom_put( TRIM(icetrcnm(jn))//'_gca'  , icetra_gca(:,:,:,jn) )

         ! average over sea ice thickness categories
         CALL iom_put( icetrcnm(jn)                , SUM( icetra_gca(:,:,:,jn) , dim=3 ) * z1_at_i(:,:)  )
      ENDDO

      ! ! ocean surface BGC
      CALL iom_put( 'PHY2c_os'         , tr(:,:,1,jrdia,Kmm) )
      CALL iom_put( 'NO3_os'           , tr(:,:,1,jqno3,Kmm) )
      CALL iom_put( 'NH4_os'           , tr(:,:,1,jrnh4,Kmm) )
      CALL iom_put( 'PHY2N_os'         , tr(:,:,1,jrdn,Kmm) )
      CALL iom_put( 'PHY2CHL_os'       , tr(:,:,1,jrdch,Kmm) )
      CALL iom_put( 'GOC_os'           , tr(:,:,1,jrgoc,Kmm) )
      CALL iom_put( 'DIC_os'           , tr(:,:,1,jqdic,Kmm)*1.e06_wp )

      ! ! sea ice ocean exchanges
      CALL iom_put( 'flushrate'        , flushrate(:,:,:) )
      CALL iom_put( 'bogup'            , bogup(:,:,:) )

      CALL iom_put( 'dh_bom_cat'       , dh_bom_cat(:,:,:) )
      CALL iom_put( 'dh_sum_cat'       , dh_sum_cat(:,:,:) )
      CALL iom_put( 'da_lam_cat'       , da_lam_cat(:,:,:) )
      CALL iom_put( 'dh_snw_sum_cat'   , dh_snw_sum_cat(:,:,:) )
      CALL iom_put( 'dh_mpdrn_cat'     , dh_mpdrn_cat(:,:,:) )
      CALL iom_put( 'dh_bog_cat'       , dh_bog_cat(:,:,:) )
      CALL iom_put( 'da_lag_cat'       , da_lag_cat(:,:,:) )

      CALL iom_put( 'flush_dia'        , flush_dia(:,:,:) )
      CALL iom_put( 'slough_dia'       , slough_dia(:,:,:) )
      CALL iom_put( 'lamloss_dia'      , lamloss_dia(:,:,:) )
      CALL iom_put( 'heatexp_dia'      , heatexp_dia(:,:,:) )
      CALL iom_put( 'dt_i'             , dt_i(:,:,:) )
      CALL iom_put( 'bogup_dia'        , bogup_dia(:,:,:) )
      CALL iom_put( 'lagup_dia'        , lagup_dia(:,:,:) )
      CALL iom_put( 'nxsicedia'        , nxsicedia(:,:,:) )
      CALL iom_put( 'cxsicedia'        , cxsicedia(:,:,:) )
      
      CALL iom_put( 'flush_no3'        , flush_no3(:,:,:) )
      CALL iom_put( 'slough_no3'       , slough_no3(:,:,:) )
      CALL iom_put( 'lamloss_no3'      , lamloss_no3(:,:,:) )
      CALL iom_put( 'lagup_no3'        , lagup_no3(:,:,:) )
      CALL iom_put( 'bogup_no3'        , bogup_no3(:,:,:) )
      CALL iom_put( 'moldif_no3'       , moldif_no3(:,:,:) )
      CALL iom_put( 'fric_vel'         , fric_vel(:,:) )

      CALL iom_put( 'flush_nh4'        , flush_nh4(:,:,:) )
      CALL iom_put( 'slough_nh4'       , slough_nh4(:,:,:) )
      CALL iom_put( 'lamloss_nh4'      , lamloss_nh4(:,:,:) )
      CALL iom_put( 'moldif_nh4'       , moldif_nh4(:,:,:) )
      CALL iom_put( 'lagup_nh4'        , lagup_nh4(:,:,:) )
      CALL iom_put( 'bogup_nh4'        , bogup_nh4(:,:,:) )
      

      ! ! BGC process
      CALL iom_put( 'phot_dia'         , phot_dia(:,:,:) )
      CALL iom_put( 'lim_PAR'          , lim_PAR(:,:,:) )
      CALL iom_put( 'lim_nut'          , lim_nut(:,:,:) )
      CALL iom_put( 'lim_ice'          , lim_ice(:,:,:) )
      CALL iom_put( 'qtr_ice_bot_cat'  , qtr_ice_bot(:,:,:) )
      CALL iom_put( 'diaupn'           , diaupn(:,:,:) )
      CALL iom_put( 'chlsyn'           , chlsyn(:,:,:) )
      CALL iom_put( 'mortlin_dia'      , mortlin_dia(:,:,:) )
      CALL iom_put( 'mortquad_dia'     , mortquad_dia(:,:,:) )
      CALL iom_put( 'remin_dia'        , remin_dia(:,:,:) )
      CALL iom_put( 'nitri'            , nitri(:,:,:) )
      
      ! diagnostics
      CALL iom_put( 'icenpp_cat'       , phot_dia(:,:,:) - etares*diaupn(:,:,:) )
      CALL iom_put( 'icenpp'           , SUM( (phot_dia(:,:,:) - etares*diaupn(:,:,:)) *a_i(:,:,:) , dim=3 ) * z1_at_i(:,:)  )
      CALL iom_put( 'qchidia'          , qchidia(:,:,:) )
      CALL iom_put( 'qnidia'           , qnidia(:,:,:) )
      CALL iom_put( 'qnidiamax'        , qnidiamax(:,:,:) )

      IF (ln_dmsice) THEN
         CALL iom_put( "icedmspdrls"      , flush_dmspd(:,:,:) + slough_dmspd(:,:,:) + lamloss_dmspd(:,:,:))
         CALL iom_put( "icedmsrls"        , flush_dms(:,:,:) + slough_dms(:,:,:) + lamloss_dms(:,:,:))
         CALL iom_put( "dmsp_exud"        , dmsp_exud(:,:,:) )
         CALL iom_put( "dmsp_lysis"       , dmsp_lysis(:,:,:) )
         CALL iom_put( "dms_phot"         , dms_phot(:,:,:) )
      ENDIF

   END SUBROUTINE trc_wri_csib

#else

CONTAINS

   SUBROUTINE trc_wri_csib
      !
   END SUBROUTINE trc_wri_csib

#endif

END MODULE trcwri_csib
