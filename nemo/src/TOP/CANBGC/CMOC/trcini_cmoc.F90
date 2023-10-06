MODULE trcini_cmoc
   !!======================================================================
   !!                         ***  MODULE trcini_cmoc  ***
   !! TOP :   initialisation of the cmoc tracers
   !!======================================================================
   !! History :        !  2007  (C. Ethe, G. Madec) Original code
   !!                  !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!                  !  2022  (O. Riche) NEMO4 integration
   !!----------------------------------------------------------------------
   !! trc_ini_cmoc   : cmoc model initialisation
   !!----------------------------------------------------------------------
   USE par_trc             ! TOP parameters
   USE oce_trc
   USE trc
   ! CMOC modules
   USE par_cmoc
   USE trcnam_cmoc         ! cmoc SMS namelist
   USE trcsms_cmoc
   ! BGCM modules
   USE trc_closea_canbgc   ! bgc closea mask
   USE trcflx_canbgc       ! air-sea gas exch.
   USE trcche_canbgc       ! carbon chemistry 
   USE trcsrc_canbgc       ! external sources/other data 
   USE sms_top_canbgc      ! access ext. source arrays declaration
   !
   USE trcopt_canbgc       ! PAR attenuation
   USE trcsink_canbgc      ! CANBGC particules sinking package
   !
   USE cmocprod            ! CMOC PP module
   ! USE cmocmort            ! CMOC phyto mortality module
   ! USE cmocrem             ! CMOC carbon remineralization
   ! USE cmoczoo             ! CMOC zooplankton grazing
   USE cmocnzd             ! consolidated module containing remineralization to (N)itrate
                           ! (Z)ooplankton grazing, and (D)etritus for phytoplankton mortality
   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_cmoc   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcini_cmoc.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_cmoc( Kmm )
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_cmoc  ***  
      !!
      !! ** Purpose :   initialization for cmoc model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   Kmm  ! time level indices
      INTEGER  :: jn, jp
      CHARACTER(len = 20)  ::  cltra
      REAL(wp), SAVE ::   sco2   =  2.312e-3_wp
      REAL(wp), SAVE ::   alka0  =  2.426e-3_wp
      REAL(wp), SAVE ::   oxyg0  =  177.6e-6_wp
      REAL(wp), SAVE ::   no30   =   31.04e-6_wp * 6.625_wp ! match value in CCC*CMOC/trcini_pisces.F90
      REAL(wp), SAVE ::   poc0   =    1.0e-8_wp
      REAL(wp), SAVE ::   phy0   =    1.0e-8_wp
      REAL(wp), SAVE ::   zoo0   =    1.0e-8_wp
      REAL(wp), SAVE ::   nch0   =    1.0e-8_wp * 12. / 55. ! match value in CCC*CMOC/trcini_pisces.F90
      !
      ! Load namelists for shared parameters
      CALL trc_nam_cmoc
	    !                           ! Allocate cmoc arrays
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_cmoc: passive tracer unit vector'
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
      END DO
      !
      ! assign an index in trc array for each extra CMOC prognostic variable
      DO jp = 1,jp_cmoc
       jn = jp + jp_bgc
       write(numout,*) ctrcnm(jn)
       cltra = TRIM( ctrcnm(jn) ) !! OR Jan 19th 2023
       IF( cltra == 'POC'      )   jqpoc = jn      !: particulate organic carbon
       IF( cltra == 'PHY'      )   jqphy = jn      !: phyto carbon
       IF( cltra == 'ZOO'      )   jqzoo = jn      !: zooplankton carbon
       IF( cltra == 'NCHL'     )   jqnch = jn      !: phyto chl-a
       ! IF( cltra == 'DICnat'   )   jqdnt = jn      !: natural DIC   
       ! IF( cltra == 'DICabio'  )   jqdab = jn      !: abiotic DIC 
       ! IF( cltra == 'ALKabio'  )   jqaab = jn      !: abiotic alkalinity
       ! IF( cltra == 'O2abio'   )   jqoab = jn      !: abiotic oxygen
       ! IF( cltra == 'DI14C'    )   jqdrc = jn      !: abiotic DI14C
      END DO
      ! closea mask for BGCM
      CALL trc_closea_init(read_var_flag=.true.)
      !
      IF( .NOT. ln_rsttr ) THEN
        tr(:,:,:,jqdic,Kmm) = sco2  * tmask_bgc_closea(:,:,:)
        tr(:,:,:,jqtal,Kmm) = alka0 * tmask_bgc_closea(:,:,:)
        tr(:,:,:,jqoxy,Kmm) = oxyg0 * tmask_bgc_closea(:,:,:)
        tr(:,:,:,jqno3,Kmm) = no30  * tmask_bgc_closea(:,:,:)
        tr(:,:,:,jqpoc,Kmm) = poc0  * tmask_bgc_closea(:,:,:)
        tr(:,:,:,jqphy,Kmm) = phy0  * tmask_bgc_closea(:,:,:)
        tr(:,:,:,jqzoo,Kmm) = zoo0  * tmask_bgc_closea(:,:,:)
        tr(:,:,:,jqnch,Kmm) = nch0  * tmask_bgc_closea(:,:,:)
      ENDIF
      !
      CALL cmoc_alloc ! allocate arrays space, see end of this module
      !
      ! O Riche Aug 4th 2022
      ! Initialise external sources reading
      ! Check namelist_top_cfg for &trcsrc_dta section
      ! call trc_src_init
      ! prep reading external sources
      ! open the files
      CALL trc_src_init
      !
      CALL trc_n2fx_init_cmoc ! load N2 fixation parameters
      !
      ! call all the BGC initialization subroutines in TOP tier
      !
      CALL trc_che_init_2D ( Kmm )
      !
      CALL trc_flx_init
      !
      ! Calendar day counter for chemistry calls
      qndayflxtr = 0
      !
      CALL trc_opt_init    ! PAR 
      !
      CALL cmoc_prod_init  ! Primary production
      !
      CALL cmoc_mort_init  ! Phytoplankton mortality
      !
      CALL cmoc_rem_init   ! NO3 remineralization
      !
      CALL cmoc_zoo_init   ! Zooplankton grazing and mortality
      !
      CALL cmoc_sink_init  ! Particule sinking
      !
      ! !
   END SUBROUTINE trc_ini_cmoc

   SUBROUTINE cmoc_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE pisces_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of cmoc 
      !!----------------------------------------------------------------------
      !
      USE lib_mpp, only: mpp_sum         ! ierr sum over all processors      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      !ierr =        sms_cmoc_alloc()          ! Start of cmoc-related alloc routines...
      ierr =           trc_opt_alloc()
      ierr = ierr +    sms_top_alloc()
      ierr = ierr +    trc_che_alloc()
      ierr = ierr +    trc_flx_alloc()
      ierr = ierr +  cmoc_sink_alloc()
      ierr = ierr +  trc_sms_cmoc_alloc()
      !
      IF( lk_mpp    )   CALL mpp_sum( 'cmoc_alloc', ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'cmoc_alloc: unable to allocate cmoc arrays' )
      !
   END SUBROUTINE cmoc_alloc

   !!======================================================================
END MODULE trcini_cmoc
