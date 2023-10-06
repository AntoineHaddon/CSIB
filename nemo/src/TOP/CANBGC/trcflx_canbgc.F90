MODULE trcflx_canbgc
   !!======================================================================
   !!                         ***  MODULE trcflx  ***
   !! TOP :   COMPUTES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!======================================================================
   !! History :   OPA  !  1988     (E. Maier-Reimer)  Original code
   !!              -   !  1998     (O. Aumont)  addition
   !!              -   !  1999     (C. Le Quere)  modification
   !!   NEMO      1.0  !  2004     (O. Aumont)  modification
   !!              -   !  2006     (R. Gangsto)  modification
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J.Orr ) update O2 solubility constants
   !!                  !  2022-2023(J. Christian, O. Riche) NEMO4 integration
   !!----------------------------------------------------------------------
   !!   trc_flx       :   main code
   !!   trc_flx_init  :   read the namelist
   !!   trc_flx_alloc :   allocate array space in memory
   !!======================================================================
! Module calls section
 
   USE trc                       ! time step in seconds whether or not euler is activated
   USE oce_trc                   ! give access to active tracers, wndm, ts, and fr_i (wind @10m, T/S, and ice fraction)
                                 ! oce_trc calls common OCE and TOP indices, e.g. jpi,jpj dimensions
   USE sms_top_canbgc            ! contains all common variables to Canadian BGCMs
 
   USE dom_oce									  ! give access to domain grid and z-levels
                                               ! grid cell area and tmask
   USE par_trc                                   ! par_trc calls par_kind and par among others, wp defined
   USE in_out_manager							  ! in_out_manager grants access to ln_timing variable among others
   USE trc_closea_canbgc
 
   USE iom                       ! to access iom_put for diagnostics
 
   USE trcche_canbgc  			      ! Carbon chemistry module
   USE trc_closea_canbgc         ! bgc-specific closea mask
   USE trcsrc_canbgc

! General scope section
    IMPLICIT NONE
    PRIVATE
    
! Procedures section

    PUBLIC   trc_flx  
    PUBLIC   trc_flx_init  
    PUBLIC   trc_flx_alloc 

! Variables section
	! partial pressures and fluxes
    REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: patmo      !: atmospheric pressure at kt [N/m2] but later on initialized to 1 atm
    REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: oce_co2g   !: CO2 flux
    REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: satmco2g   !: atmospheric pco2 
    REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: oce_o2g    !: O2  flux
    REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: satmo2g    !: atmospheric po2 

   ! !!* Substitution
#  include "domzgr_substitute.h90"

CONTAINS

   SUBROUTINE trc_flx(kt, Kmm, Krhs)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_flx  ***
      !!
      !! ** Purpose :   COMPUTES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
      !!
      !! ** Method  : 
      !!              - for now compute a partial pressure difference from
      !!                available tracers
      !!              - pCO2(w) is approximate as a percent of DIC for now
      !!              - no solubility is computed for now
      !!              - testing reading temperature and wind from OCE
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   !
      INTEGER, INTENT(in) ::   Kmm, Krhs  ! time level indices
      !
      INTEGER  ::   ji, jj, jm
      REAL(wp) ::   ztc, ztc2, ztc3, ztc4, zws, zkgwan
      REAL(wp) ::   zfld, zflu, zfld16, zflu16, zfact
      REAL(wp) ::   zsch_o2, zsch_co2
      REAL(wp), DIMENSION(jpi,jpj) :: zkgco2, zkgo2, zo2flx, zco2flx

      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('trc_flx')
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_flx:  air-sea processes' 
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~'

      ! SURFACE CHEMISTRY (PCO2 AND [H+] IN
      !     SURFACE LAYER); THE RESULT OF THIS CALCULATION
      !     IS USED TO COMPUTE AIR-SEA FLUX OF CO2

      ! 1. compute gas exchange velocities
      ! -------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            ztc  = ts(ji,jj,1,jp_tem,Kmm)
            ztc2 = ztc * ztc
            ztc3 = ztc * ztc2 
            ztc4 = ztc * ztc3 
			! O Riche June 10th 2022
			! Equations here are contemporary to NEMO 3.4/3.6?
			! They are original NEMO PISCES equations still used in NEMO4.2 PISCES
			! https://forge.nemo-ocean.eu/nemo/nemo/-/blob/main/src/TOP/PISCES/P4Z/p4zflx.F90
			! Wanninkhof 1992 is in the end-of-article biblio list in the Aumont et al 2015
			! Csanady 2001/ Sarmiento and Gruber 2006
			! The relationship here is based on the notion that
			! the gas concentration in seawater gets toward equilibrium with
			! the gas concentration on the air side, here we compute the 
			! gas exchange rate scale
            ! Compute the Schmidt Number both O2 and CO2
            zsch_co2 = 2116.8 - 136.25 * ztc + 4.7353 * ztc2 - 0.092307 * ztc3 + 0.0007555 * ztc4
            zsch_o2  = 1920.4 - 135.6  * ztc + 5.2122 * ztc2 - 0.10939  * ztc3 + 0.00093777 * ztc4
            !  wind speed 
            zws  = wndm(ji,jj) * wndm(ji,jj)
            ! Compute the piston velocity for O2 and CO2
			! Formulation of piston velocity is close to Wanninkhof 1992 found in Table 3.3.2
			! with a factor 0.251 smaller than 0.31 of Sarmiento and Gruber 2006 also citing Wanninkhof 1992
			! it's possible it's a coefficient adjusted to work best for both CO2 and O2
            zkgwan = 0.251 * zws  
            zkgwan = zkgwan * xconv0 * ( 1.- fr_i(ji,jj) ) * tmask_bgc_closea(ji,jj,1)
            ! compute gas exchange for CO2 and O2
            zkgco2(ji,jj) = zkgwan * SQRT( 660./ zsch_co2 )
            zkgo2 (ji,jj) = zkgwan * SQRT( 660./ zsch_o2 )
         END DO
      END DO

      ! 2. compute partial pressure differences and fluxes
      ! -------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi
            !zh2co3(ji,jj) = tr(ji,jj,1,jqdic, Kmm)*0.01_wp           ! set up fraction (1%) of DIC as a proxy for [H2CO3*]sw and convert in M
            ! 1% seems to be a reasonable fraction based on Zeebe et al textbook Fig. 1.6.27
            ! to represent total dissolved CO2g in seawater?	  
            ! CO2(g) <=> H2CO3 Zeebe Eqs 1.1.1 and H2CO3 in equilibrium with DIC (HCO3 and CO3 2-)
            ! Compute CO2 flux for the sea and air
            ! zfld flux is based on Henry's law giving the equilibrium concentration to reach
            ! in seawater and the partial pressure of CO2 g (O2 g in air at sea surface level)
            ! partial pressure in air is converted to seawater concentration in mol L^-1

            zfld = satmco2g(ji,jj) * K0CO2(ji,jj) * tmask_bgc_closea(ji,jj,1) * zkgco2(ji,jj)     ! (mol/L) * (m/s)
            zflu = qh2co3(ji,jj) * tmask_bgc_closea(ji,jj,1) * zkgco2(ji,jj)                        ! (mol/L) * (m/s) 

            oce_co2g(ji,jj) = ( zfld - zflu ) * qfact * e1e2t(ji,jj) * tmask_bgc_closea(ji,jj,1) * 1000. ! convert L^-1 to m^-3
            zco2flx(ji,jj)  = ( zfld - zflu ) * tmask_bgc_closea(ji,jj,1)
            tr(ji,jj,1,jqdic, Krhs) = tr(ji,jj,1,jqdic, Krhs) + zco2flx(ji,jj) / e3t(ji,jj,1,Kmm)
   
            ! Compute O2 flux 
            zfld16 = satmo2g(ji,jj) * K0O2(ji,jj) * tmask_bgc_closea(ji,jj,1) * zkgo2(ji,jj)     ! (mol/L) * (m/s)
            zflu16 = tr(ji,jj,1,jqoxy, Kmm) * tmask_bgc_closea(ji,jj,1) * zkgo2(ji,jj)                  ! (mol/L) * (m/s)

            oce_o2g(ji,jj) = ( zfld16 - zflu16 ) * qfact * e1e2t(ji,jj) * tmask_bgc_closea(ji,jj,1) * 1000. ! convert L^-1 to m^-3
            zo2flx(ji,jj)  = ( zfld16 - zflu16 ) * tmask_bgc_closea(ji,jj,1)
            tr(ji,jj,1,jqoxy, Krhs) = tr(ji,jj,1,jqoxy, Krhs) + zo2flx(ji,jj) / e3t(ji,jj,1,Kmm)
         END DO
      END DO

      ! O Riche Sept 7th 2022
      ! Diagnostics
      ! gas exchange rates
      CALL iom_put("KgCO2", zkgco2(:,:) )
      CALL iom_put("KgO2" ,  zkgo2(:,:) )
      ! net flux to the ocean
      CALL iom_put("Cflx", zco2flx(:,:) / qfact / e1e2t(:,:) )
      CALL iom_put("Oflx",  zo2flx(:,:) / qfact / e1e2t(:,:) )
      ! partial pressures
      CALL iom_put("DpCO2", ( satmco2g(:,:) - qh2co3(:,:) / ( K0CO2(:,:) + rtrn ) )   * tmask_bgc_closea(:,:,1) )
      CALL iom_put("pCO2" ,                 ( qh2co3(:,:) / ( K0CO2(:,:) + rtrn ) )   * tmask_bgc_closea(:,:,1) )
      CALL iom_put("DpO2" , ( satmo2g(:,:) - tr(:,:,1,jqoxy, Kmm) / ( K0O2(:,:) + rtrn ) ) * tmask_bgc_closea(:,:,1) )
      CALL iom_put("pO2"  ,                ( tr(:,:,1,jqoxy, Kmm) / ( K0O2(:,:) + rtrn ) ) * tmask_bgc_closea(:,:,1) )
      ! Carbonate system
      ! other fields will be set to 0s by default (compilation setting)
      ! CALL iom_put("CO3",      )
      ! CALL iom_put("CO3sat",   )
      ! CALL iom_put("HCO3",     )

      !
      IF( ln_timing )  CALL timing_stop('trc_flx')
      ! 
   END SUBROUTINE trc_flx
   
   
   SUBROUTINE trc_flx_init

	   ! Initialize partial pressure difference fields
       oce_co2g(:,:)  = 0._wp                ! Initialization co2(g)
       oce_o2g(:,:)   = 0._wp                ! Initialization o2(g)
       patmo(:,:)     = 1.e0                 ! Initialize patmo if no reading from a file
       satmco2g(:,:)  = atmco2*patmo(:,:)     ! Initialization with constant value default for ln_co2int = F
       satmo2g(:,:)   = atcoxy*patmo(:,:)    ! Initialization of atm pO2

       IF(lwp) WRITE(numout,*)
       IF(lwp) WRITE(numout,*) ' trc_flx_init : initialization'
       IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~'

       ! Silic acid clim
       ! qasi3 is used in CO2 flux calculation using
       CALL trc_src3d( 1 , js3d_si)
       qasi3(:,:,1) = src3d_dta(:,:,1,js3d_si)
    
   END SUBROUTINE trc_flx_init
 

   INTEGER FUNCTION trc_flx_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_flx_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to CANOE
      ALLOCATE( oce_co2g(jpi,jpj), satmco2g(jpi,jpj), patmo(jpi,jpj),    &
      &         oce_o2g(jpi,jpj),  satmo2g(jpi,jpj), STAT=trc_flx_alloc )

      IF( trc_flx_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_flx_alloc : failed to allocate arrays' )
      !
   END FUNCTION trc_flx_alloc


   !!======================================================================
END MODULE  trcflx_canbgc
