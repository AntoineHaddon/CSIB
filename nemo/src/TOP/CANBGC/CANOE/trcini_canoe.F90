MODULE trcini_canoe
   !!======================================================================
   !!                         ***  MODULE trcini_canoe  ***
   !! TOP :   initialisation of the CANOE tracers
   !!======================================================================
   !! History :        !  2007  (C. Ethe, G. Madec) Original code
   !!                  !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!                  !  2022  (O. Riche) NEMO4 integration
   !!----------------------------------------------------------------------
   !! trc_ini_canoe   : CANOE model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   ! CanOE modules
   USE par_canoe
   USE trcnam_canoe        ! canoe SMS namelist
   USE trcsms_canoe
   ! BGCM modules
   USE trc_closea_canbgc ! bgc closea mask
   USE trcflx_canbgc     ! air-sea gas exch.
   USE trcche_canbgc     ! carbon chemistry 
   USE trcsrc_canbgc     ! external sources/other data 
   USE sms_top_canbgc    ! access ext. source arrays declaration
   !
   USE trcopt_canbgc     ! PAR attenuation
   USE trcsink_canbgc    ! CANBGC particules sinking package
   !
   USE sms_canoe         ! set elemental parameters
   USE canoetemp         ! CanOE temperature dependencies module
   USE canoeprod         ! CanOE PP module
   USE canoenzd          ! CanOE grazing/mortality/remineralization module
   !
   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_canoe   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcini_canoe.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_canoe( Kmm )
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_canoe  ***  
      !!
      !! ** Purpose :   initialization for CANOE model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   Kmm  ! time level indices
      INTEGER  :: jn, jp
      CHARACTER(len = 20)  ::  cltra
      REAL(wp), SAVE ::   sco2   =  2.312e-3_wp
      REAL(wp), SAVE ::   alka0  =  2.426e-3_wp
      REAL(wp), SAVE ::   oxyg0  =177.6_wp
      REAL(wp), SAVE ::   no30   = 31.04_wp
      REAL(wp), SAVE ::   bioma0 =  1.e-2_wp
      ! !
      CALL trc_nam_canoe   
	    !                       ! Allocate CANOE arrays
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_canoe: passive tracer unit vector'
      IF(lwp) WRITE(numout,*) ' To check conservation : '
      IF(lwp) WRITE(numout,*) '   1 - No sea-ice model '
      IF(lwp) WRITE(numout,*) '   2 - No runoff ' 
      IF(lwp) WRITE(numout,*) '   3 - precipitation and evaporation equal to 1 : E=P=1 ' 
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      ! assign an index in trc array for each prognostic variable
      DO jn = 1,jp_bgc
       IF( lwp ) THEN
         WRITE(numout,*) ctrcnm(jn)
       ENDIF
       !
       cltra = TRIM( ctrcnm(jn) )
       IF( cltra == 'DIC'      )   jqdic = jn      !: dissolved inorganic carbon concentration
       IF( cltra == 'Alkalini' )   jqtal = jn      !: total alkalinity
       IF( cltra == 'O2'       )   jqoxy = jn      !: oxygen concentration
       IF( cltra == 'NO3'      )   jqno3 = jn      !: NO3 concentration
      ENDDO
      !
      ! assign an index in trc array for each extra CMOC prognostic variable
      DO jp = 1,jp_canoe
       jn = jp + jp_bgc
       write(numout,*) ctrcnm(jn)
       cltra = TRIM( ctrcnm(jn) )     
       IF( cltra == 'NH4'      )   jrnh4 = jn      !: Ammonium Concentration
       IF( cltra == 'Fer'      )   jrfer = jn      !: Dissolved Iron Concentration
       IF( cltra == 'CaCO3'    )   jrcal = jn      !: calcite particules
       IF( cltra == 'POC'      )   jrpoc = jn      !: small sized POC
       IF( cltra == 'GOC'      )   jrgoc = jn      !: large sized POC
       IF( cltra == 'PHYC'     )   jrphy = jn      !: small sized phyto C biomass
       IF( cltra == 'NCHL'     )   jrnch = jn      !: small sized phyto chl-a
       IF( cltra == 'PHYN'     )   jrnn  = jn      !: Nanophytoplankton N concentration 
       IF( cltra == 'PHYFe'    )   jrnfe = jn      !: Nanophytoplankton Fe concentration 
       IF( cltra == 'PHY2C'    )   jrdia = jn      !: large sized phyto C biomass
       IF( cltra == 'DCHL'     )   jrdch = jn      !: large sized phyto chl-a
       IF( cltra == 'PHY2N'    )   jrdn  = jn      !: Diatoms N concentration 
       IF( cltra == 'PHY2Fe'   )   jrdfe = jn      !: Diatoms Fe concentration 
       IF( cltra == 'ZOO'      )   jrzoo = jn      !: small sized zoo C biomass
       IF( cltra == 'ZOO2'     )   jrmes = jn      !: large sized zoo C biomass         
       !       
      END DO
      !
      ! closea mask for BGCM
      CALL trc_closea_init(read_var_flag=.true.)
      !
      IF( .NOT. ln_rsttr ) THEN
        !
        tr(:,:,:,jqdic,Kmm) = sco2  * tmask_bgc_closea(:,:,:)
        tr(:,:,:,jqtal,Kmm) = alka0 * tmask_bgc_closea(:,:,:)
        tr(:,:,:,jqoxy,Kmm) = oxyg0 * tmask_bgc_closea(:,:,:)
        tr(:,:,:,jrcal,Kmm) = bioma0              * tmask_bgc_closea(:,:,:)  ! calcite particules
        tr(:,:,:,jrpoc,Kmm) = bioma0              * tmask_bgc_closea(:,:,:)  !: small sized POC
        tr(:,:,:,jrphy,Kmm) = bioma0              * tmask_bgc_closea(:,:,:)  !: small sized phyto C biomass
        tr(:,:,:,jrnn ,Kmm) = bioma0 * 12./106.   * tmask_bgc_closea(:,:,:)  ! Nanophytoplankton N concentration 
        tr(:,:,:,jrnfe,Kmm) = bioma0 * 5.         * tmask_bgc_closea(:,:,:)  ! Nanophytoplankton Fe concentration 
        tr(:,:,:,jrnch,Kmm) = bioma0 * 12./55.    * tmask_bgc_closea(:,:,:)  !: small sized phyto chl-a
        tr(:,:,:,jrdia,Kmm) = bioma0              * tmask_bgc_closea(:,:,:)  ! large sized phyto C by
        tr(:,:,:,jrdn ,Kmm) = bioma0 * 12./106.   * tmask_bgc_closea(:,:,:)  ! Diatoms N concentration 
        tr(:,:,:,jrdfe,Kmm) = bioma0 * 5.         * tmask_bgc_closea(:,:,:)  ! Diatoms Fe concentration 
        tr(:,:,:,jrdch,Kmm) = bioma0 * 12./55.    * tmask_bgc_closea(:,:,:)  ! large sized phyto chl-a
        tr(:,:,:,jrzoo,Kmm) = bioma0              * tmask_bgc_closea(:,:,:)  !: small sized zoo C biomass
        tr(:,:,:,jrmes,Kmm) = bioma0              * tmask_bgc_closea(:,:,:)  ! large sized zoo C biomass
        tr(:,:,:,jrfer,Kmm) = 6.e2_wp             * tmask_bgc_closea(:,:,:)  ! Dissolved Iron Concentration
        tr(:,:,:,jrgoc,Kmm) = bioma0              * tmask_bgc_closea(:,:,:)  ! large sized POC
        tr(:,:,:,jqno3,Kmm) = no30  * tmask_bgc_closea(:,:,:)
        tr(:,:,:,jrnh4,Kmm) = bioma0 * 16./106.   * tmask_bgc_closea(:,:,:)  ! Ammonium Concentration        
        !
      ENDIF
      !     
      ! Test allocation of space for CanOE arrays before initialization
      CALL canoe_alloc ! allocate arrays space, see end of this module
      !      
      ! O Riche Aug 4th 2022
      ! Initialise external sources reading
      ! Check namelist_top_cfg for &trcsrc_dta section
      ! call trc_src_init
      ! prep reading external sources
      ! open the files
      CALL trc_src_init
      !
      !
      !
      ! call all the BGC initialization subroutines in TOP tier
      !
      CALL trc_che_init_2D( Kmm )
      CALL trc_che_init_3D( Kmm )
      !
      CALL trc_flx_init
      !
      ! Calendar day counter for chemistry calls
      qndayflxtr = 0
      !
      CALL trc_opt_init
      !
      CALL canoe_temp_init
      !
      CALL canoe_prod_init
      !
      CALL canoe_nzd_init
      !
      CALL canoe_sink_init
      !
      CALL trc_n2fx_init_canoe
      !
      ! Set elemental ratios
      ! ---------------------
      rr_c2n  =  6.625_wp
      rr_fe2n = 33._wp
      rr_n2c  =      1./rr_c2n
      rr_n2fe =      1./rr_fe2n
      rr_fe2c = rr_fe2n/rr_c2n
      rr_c2fe =      1./rr_fe2c     
      !
      ! !
   END SUBROUTINE trc_ini_canoe

   SUBROUTINE canoe_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE pisces_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of CANOE 
      !!----------------------------------------------------------------------
      !
      USE lib_mpp, only: mpp_sum         ! ierr sum over all processors      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      !ierr =        sms_canoe_alloc()          ! Start of CANOE-related alloc routines...
      ierr =        trc_opt_alloc()
      ierr = ierr + sms_top_alloc()
      ierr = ierr + trc_che_alloc()
      ierr = ierr + trc_flx_alloc()
      ierr = ierr + canoe_sink_alloc()
      ierr = ierr + trc_sms_canoe_alloc()      
      ierr = ierr + canoe_nzd_alloc()      
      !
      IF( lk_mpp    )   CALL mpp_sum( 'canoe_alloc', ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'canoe_alloc: unable to allocate CANOE arrays' )
      !
   END SUBROUTINE canoe_alloc

   !!======================================================================
END MODULE trcini_canoe
