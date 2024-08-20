MODULE trcwri_csib
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!     trc_wri_csib   :  outputs of concentration fields
   !!======================================================================
#if defined key_top && defined key_xios
   !!----------------------------------------------------------------------
   !! History :      !  2007  (C. Ethe, G. Madec)  Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   USE par_trc         ! passive tracers common variables
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
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
      ! INTEGER              :: jn
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      CALL iom_put( 'icedia'           , icedia(:,:,:) )
      CALL iom_put( 'icedia_gca'       , icedia_gca(:,:,:) )

      CALL iom_put( 'iceno3'           , iceno3(:,:,:) )
      CALL iom_put( 'iceno3_gca'       , iceno3_gca(:,:,:) )
      
      CALL iom_put( 'icenh4'           , icenh4(:,:,:) )
      CALL iom_put( 'icenh4_gca'       , icenh4_gca(:,:,:) )

      ! ! ocean surface BGC
      CALL iom_put( 'PHY2c_os'         , tr(:,:,1,jrdia,Kmm) )
      CALL iom_put( 'NO3_os'           , tr(:,:,1,jqno3,Kmm) )
      CALL iom_put( 'NH4_os'           , tr(:,:,1,jrnh4,Kmm) )
      ! CALL iom_put( 'NO3_i-NO3_os-cat1'     , iceno3(:,:,1)-tr(:,:,1,jqno3,Kmm) )
      CALL iom_put( 'NO3_i-NO3_os-cat1'     , tr(:,:,2,jqno3,Kmm) )

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
      CALL iom_put( 'lamloss_dia'      , lamloss_dia(:,:,:) )
      CALL iom_put( 'bogup_dia'        , bogup_dia(:,:,:) )
      CALL iom_put( 'lagup_dia'        , lagup_dia(:,:,:) )
      
      CALL iom_put( 'flush_no3'        , flush_no3(:,:,:) )
      CALL iom_put( 'lamloss_no3'      , lamloss_no3(:,:,:) )
      CALL iom_put( 'moldif_no3'       , moldif_no3(:,:,:) )
      CALL iom_put( 'lagup_no3'        , lagup_no3(:,:,:) )

      CALL iom_put( 'flush_nh4'        , flush_nh4(:,:,:) )
      CALL iom_put( 'lamloss_nh4'      , lamloss_nh4(:,:,:) )
      CALL iom_put( 'moldif_nh4'       , moldif_nh4(:,:,:) )
      CALL iom_put( 'lagup_nh4'        , lagup_nh4(:,:,:) )
      
      CALL iom_put( 'fric_vel'         , fric_vel(:,:) )


      ! ! BGC process
      CALL iom_put( 'growth_dia'       , growth_dia(:,:,:) )
      CALL iom_put( 'lim_lig'          , lim_lig(:,:,:) )
      CALL iom_put( 'qtr_ice_bot_cat'  , qtr_ice_bot(:,:,:) )

   END SUBROUTINE trc_wri_csib

#else

CONTAINS

   SUBROUTINE trc_wri_csib
      !
   END SUBROUTINE trc_wri_csib

#endif

END MODULE trcwri_csib
