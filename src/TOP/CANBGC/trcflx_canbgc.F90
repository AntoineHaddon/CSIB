MODULE trcflx_canbgc
   !!======================================================================
   !!                         ***  MODULE trcflx  ***
   !! TOP :   COMPUTES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!======================================================================
   !!----------------------------------------------------------------------
   !!   trc_flx       :   main code
   !!   trc_flx_init  :   read the namelist
   !!   trc_flx_alloc :   allocate array space in memory
   !!----------------------------------------------------------------------

! Module calls section
	! internal calls to 
	USE trc                       ! time step in seconds whether or not euler is activated
  USE oce_trc                   ! give access to active tracers, wndm, tsn, and fr_i (wind @10m, T/S, and ice fraction)
                                ! oce_trc calls common OCE and TOP indices, e.g. jpi,jpj dimensions
	USE sms_top_canbgc            ! contains all common variables to Canadian BGCMs
	USE dom_oce									  ! give access to domain grid and z-levels
                             	  ! grid cell area and tmask
	USE par_trc                   ! par_trc calls par_kind and par among others, wp defined
	USE in_out_manager					  ! in_out_manager grants access to ln_timing variable among others

  USE iom                       ! to access iom_put for diagnostics

  USE trcche_canbgc  			      ! Carbon chemistry module
  USE trc_closea_canbgc         ! bgc-specific closea mask

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

CONTAINS

   SUBROUTINE trc_flx(kt)
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
      !
      INTEGER  ::   ji, jj, jm
      REAL(wp) ::   ztc, ztc2, ztc3, ztc4, zws, zkgwan
      REAL(wp) ::   zfld, zflu, zfld16, zflu16, zfact
	    REAL(wp) ::   zph, zah2, zbot, zdic, zalk, zalka               ! carbon chemistry wrk variables
      REAL(wp) ::   zph2, zph3, zpo4, zsi, zpd, zp0, zp1, zp3        ! coefficients added to account for P and Si contribution to TA
      REAL(wp) ::   zsch_o2, zsch_co2
      REAL(wp), DIMENSION(jpi,jpj)     :: zkg0, zkgco2, zkgo2, zh2co3, zo2flx, zco2flx 
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zph0
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('trc_flx')
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_flx:  air-sea processes'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

      ! SURFACE CHEMISTRY (PCO2 AND [H+] IN
      !     SURFACE LAYER); THE RESULT OF THIS CALCULATION
      !     IS USED TO COMPUTE AIR-SEA FLUX OF CO2

      DO jm = 1, 10
	  !
         DO jj = 1, jpj
	  !
            DO ji = 1, jpi

               ! DUMMY VARIABLES FOR DIC, H+, AND BORATE
               zbot  = qborat(ji,jj,1)
               zfact = rhop(ji,jj,1) / 1000. + rtrn
               zdic  = trn(ji,jj,1,jqdic) / zfact
               zph   = MAX( qhi(ji,jj,1), 1.e-10 ) / zfact
               zalka = trn(ji,jj,1,jqtal) / zfact
               zph2  = zph*zph
               zph3  = zph*zph2
               ! zpo4 = (trn(ji,jj,1,jqno3)+trn(ji,jj,1,jqnh4)) / 16. *0.000001 / zfact
               zpo4  = (5+.5) / 16. *0.000001 / zfact  
               ! O Riche Aug 16th 2022
               ! 5. and .5 are placeholder for nitrate and ammonium trn arrays
               ! which eventually will be used 
               ! in the meantime we can do better and assigned the phosphate IC
               ! values instead of these
               zsi   = qasi3(ji,jj,1) * 0.000001 / zfact                        ! silica is a static array based on initialization file, not a carried tracer

               ! CALCULATE P AND Si ION CONCENTRATIONS AS PER ORR ET AL (BPG EQUATIONS 43-47)
               ! zp3 = H3PO4, zp1 = HPO4(2-), zp0 = PO4(3-): denominator is the same for all 3 equations
               zpd = 1./ ( zph3 + qakp13(ji,jj,1)*zph2 + qakp13(ji,jj,1)*qakp23(ji,jj,1)*zph + qakp13(ji,jj,1)*qakp23(ji,jj,1)*qakp33(ji,jj,1) )
               zp3 = zph3*zpo4 * zpd
               zp1 = zph*zpo4*qakp13(ji,jj,1)*qakp23(ji,jj,1) * zpd
               zp0 = zpo4*qakp13(ji,jj,1)*qakp23(ji,jj,1)*qakp33(ji,jj,1) * zpd
               zsi = zsi / (1. + zph / qaksi3(ji,jj,1))

               ! CALCULATE [ALK]([CO3--], [HCO3-])
               zalk  = zalka - (  qakw3(ji,jj,1) / zph - zph + zbot / ( 1.+ zph / qakb3(ji,jj,1) ) + 2.*zp0 + zp1 - zp3 + zsi )

               ! CALCULATE [H+] AND [H2CO3]
               zah2   = SQRT(  (zdic-zalk)*(zdic-zalk) + 4.* ( zalk * qak23(ji,jj,1)   &
                  &                                        / qak13(ji,jj,1) ) * ( 2.* zdic - zalk )  )
               zah2   = 0.5 * qak13(ji,jj,1) / zalk * ( ( zdic - zalk ) + zah2 )
               zh2co3(ji,jj) = ( 2.* zdic - zalk ) / ( 2.+ qak13(ji,jj,1) / zah2 ) * zfact
               qhi(ji,jj,1)   = zah2 * zfact
            END DO
         END DO
      END DO

	  !
	  ! ----------------
	  ! compute fluxes
	  ! ----------------
	   
      ! 1. compute gas exchange velocities
      ! -------------------------------------------
      DO jj = 1, jpj

         DO ji = 1, jpi
            ztc  = tsn(ji,jj,1,jp_tem)
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
              zkg0(ji,jj) = zkgwan
            zkgco2(ji,jj) = zkgwan * SQRT( 660./ zsch_co2 )
            zkgo2 (ji,jj) = zkgwan * SQRT( 660./ zsch_o2  )
         END DO
      END DO

      ! 2. compute partial pressure differences and fluxes
      ! -------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi
			! zh2co3(ji,jj) = trn(ji,jj,1,jqdic)*.01_wp           ! set up fraction (1%) of DIC as a proxy for [H2CO3*]sw and convert in M
			! 1% seems to be a reasonable fraction based on Zeebe et al textbook Fig. 1.6.27
			! to represent total dissolved CO2g in seawater?	  
			! CO2(g) <=> H2CO3 Zeebe Eqs 1.1.1 and H2CO3 in equilibrium with DIC (HCO3 and CO3 2-)
            ! Compute CO2 flux for the sea and air
			! zfld flux is based on Henry's law giving the equilibrium concentration to reach
			! in seawater and the partial pressure of CO2 g (O2 g in air at sea surface level)
			! partial pressure in air is converted to seawater concentration in mol L^-1


            zfld = satmco2g(ji,jj) * chemc(ji,jj,1) * tmask_bgc_closea(ji,jj,1) * zkgco2(ji,jj)     ! (mol/L) * (m/s)
            zflu = zh2co3(ji,jj) * tmask_bgc_closea(ji,jj,1) * zkgco2(ji,jj)                        ! (mol/L) * (m/s) 

            oce_co2g(ji,jj) = ( zfld - zflu ) * rfact * e1e2t(ji,jj) * tmask_bgc_closea(ji,jj,1) * 1000. ! convert L^-1 to m^-3
            zco2flx(ji,jj)  = ( zfld - zflu ) * tmask_bgc_closea(ji,jj,1)
            tra(ji,jj,1,jqdic) = tra(ji,jj,1,jqdic) + zco2flx(ji,jj) / e3t_n(ji,jj,1)
			
            ! Compute O2 flux 
            zfld16 = satmo2g(ji,jj) * chemc(ji,jj,2) * tmask_bgc_closea(ji,jj,1) * zkgo2(ji,jj)     ! (mol/L) * (m/s)
            zflu16 = trn(ji,jj,1,jqoxy) * tmask_bgc_closea(ji,jj,1) * zkgo2(ji,jj)                  ! (mol/L) * (m/s)

            oce_o2g(ji,jj) = ( zfld16 - zflu16 ) * rfact * e1e2t(ji,jj) * tmask_bgc_closea(ji,jj,1) * 1000. ! convert L^-1 to m^-3
            zo2flx(ji,jj)  = ( zfld16 - zflu16 ) * tmask_bgc_closea(ji,jj,1)
            tra(ji,jj,1,jqoxy) = tra(ji,jj,1,jqoxy) + zo2flx(ji,jj) / e3t_n(ji,jj,1)
         END DO
      END DO

      ! O Riche Sept 7th 2022
      ! Diagnostics
      ! gas exchange rates
      CALL iom_put("Kg"   ,   zkg0(:,:) )
      CALL iom_put("KgCO2", zkgco2(:,:) )
      CALL iom_put("KgO2" ,  zkgo2(:,:) )
      ! net flux to the ocean
      CALL iom_put("Cflx", zco2flx(:,:) / rfact / e1e2t(:,:) )
      CALL iom_put("Oflx",  zo2flx(:,:) / rfact / e1e2t(:,:) )
      ! partial pressures
      CALL iom_put("DpCO2", ( satmco2g(:,:) - zh2co3(:,:) / ( chemc(:,:,1) + rtrn ) ) * tmask_bgc_closea(:,:,1) )
      CALL iom_put("pCO2" ,                 ( zh2co3(:,:) / ( chemc(:,:,1) + rtrn ) ) * tmask_bgc_closea(:,:,1) )
      CALL iom_put("DpO2" , ( satmo2g(:,:) - trn(:,:,1,jqoxy) / ( chemc(:,:,2) + rtrn ) ) * tmask_bgc_closea(:,:,1) )
      CALL iom_put("pO2"  ,                ( trn(:,:,1,jqoxy) / ( chemc(:,:,2) + rtrn ) ) * tmask_bgc_closea(:,:,1) )
      ! Carbonate system
      zph0(:,:,:) = rtrn
      zph0(:,:,1) = qhi(:,:,1) + rtrn   ! [H+] is 2D for now so just set to epsilon if deeper than level 1 
      CALL iom_put("pH",  -1. * LOG10( MAX( zph0(:,:,:), rtrn ) ) * tmask(:,:,:))
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
       satmco2g(:,:)  = atco2*patmo(:,:)     ! Initialization with constant value default for ln_co2int = F
       satmo2g(:,:)   = atcoxy*patmo(:,:)    ! Initialization of atm pO2

       IF(lwp) WRITE(numout,*)
       IF(lwp) WRITE(numout,*) ' trc_flx_init : initialization'
       IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~'
    
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
