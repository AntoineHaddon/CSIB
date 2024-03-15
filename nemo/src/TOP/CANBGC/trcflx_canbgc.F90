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
   !!                  !  2024     (J. Christian, N. Lambert) NEMO4 CanCPL integration
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
 
   USE dom_oce                   ! give access to domain grid and z-levels
                                 ! grid cell area and tmask
   USE par_trc                   ! par_trc calls par_kind and par among others, wp defined
   USE in_out_manager            ! in_out_manager grants access to ln_timing variable among others
 
   USE iom                       ! to access iom_put for diagnostics
   USE fldread                   ! read input fields
 
   USE trcche_canbgc             ! Carbon chemistry module
   USE sbcapr                    ! Ocean dynamic atm. press
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
    REAL(wp), PUBLIC, ALLOCATABLE,       DIMENSION(:,:) :: patmo      !: atmospheric pressure at kt [atm]
    REAL(wp), PUBLIC, ALLOCATABLE,       DIMENSION(:,:) :: oce_co2g    !: total CO2  flux
    REAL(wp), PUBLIC, ALLOCATABLE,       DIMENSION(:,:) :: satmco2g   !: atmospheric pco2 ( in atm )
    REAL(wp), PUBLIC, ALLOCATABLE,       DIMENSION(:,:) :: oce_o2g    !: total O2  flux
    REAL(wp), PUBLIC, ALLOCATABLE,       DIMENSION(:,:) :: satmo2g    !: atmospheric po2 ( in atm )

        !!  Variables related to reading atmospheric CO2 time history nn_co2int = 1
    INTEGER                                   ::   nmaxrec, numco2   !
    REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::   atcco2h, years    !

    !! input fields for the patm and co2 (nn_patmint & nn_co2int = 2)
    TYPE(FLD), ALLOCATABLE,       DIMENSION(:)   ::   sf_atmco2 ! structure of co2 input fields
    TYPE(FLD), ALLOCATABLE,       DIMENSION(:)   ::   sf_patm   ! structure of patm input fields

      ! gas exchange parameters in namelist_top (nambgcext)
      REAL(wp)           :: nn_atmco2   !: Constant atmospheric [ppm co2]
      INTEGER            :: nn_co2int   ! flag to set the atmospheric pressure
                                        ! =0, co2 constant (nn_atmco2) 
                                        ! =1, read atm co2 conc. from an annual file with no spatialized fields (clname)
                                        ! =2, read atm co2 conc. from files with spatialized fields (sn_atmco2)
                                        ! =3, get the fields from the coupler (sn_rcv_co2 must be 'coupled')
      CHARACTER(len=34)  :: clname      ! filename of co2 conc. values (nn_co2int = 1)
      INTEGER            ::   nn_offset ! Offset model-data start year (default = 0) 

      REAL(wp)           :: nn_patm     ! Constant atmospheric pressure [N/m2]
      INTEGER            :: nn_patmint  ! flag to set the atmospheric pressure
                                        ! =0, air pressure constant (nn_patm) 
                                        ! =2, read patm from files with spatialized fields (sn_patm)
                                        ! =3, get the fields from the coupler (sn_rcv_mslp must be 'coupled')


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

      ! SURFACE CHEMISTRY (PCO2 AND [H+] IN
      !     SURFACE LAYER); THE RESULT OF THIS CALCULATION
      !     IS USED TO COMPUTE AIR-SEA FLUX OF CO2

      ! 1. fill the gas exchange arrays 
      ! -------------------------------------------
      SELECT  CASE (nn_patmint) ! atm. press.
      CASE (0) ! atm press constant value
         patmo(:,:)=nn_patm/101325.0
      CASE (1,3) ! atm press from ocean dyn. or coupler
         patmo(:,:)=apr/101325.0
      CASE (2) ! atm press from inputs files
         CALL fld_read( kt, 1, sf_patm )               !* input Patm provided at kt + 1/2
         patmo(:,:) = sf_patm(1)%fnow(:,:,1)/101325.0     ! atmospheric pressure
      ENDSELECT



      ! 2. compute gas exchange velocities
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

      ! 3. compute partial pressure differences and fluxes
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
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_flx_init  ***
      !!
      !! ** Purpose :   Initialization of atmospheric conditions
      !!
      !! ** Method  :   Read the nambgcext namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nambgcext
      !!----------------------------------------------------------------------
      INTEGER ::   jm, ierr, ios   ! Local integer 
      !!
      ! gas exchange parameters in namelist_top (nambgcext)
      CHARACTER(len=100) ::   cn_dir    ! Root directory for location of ssr files
      TYPE(FLD_N)        ::   sn_patm   ! informations about the fields to be read for patm
      TYPE(FLD_N)        ::   sn_atmco2 ! informations about the fields to be read atmco2
      REAL(wp)           ::   nn_no3_sf ! scaling factor for NO3 (for estimation of PO4)


      NAMELIST/nambgcext/ nn_atmco2, nn_co2int, clname, nn_offset, nn_patm, nn_patmint, sn_patm, sn_atmco2, cn_dir, nn_no3_sf
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_flx_init : initialization'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~'

      READ  ( numnat_ref, nambgcext, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambgcext in reference namelist' )
      READ  ( numnat_cfg, nambgcext, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nambgcext in configuration namelist' )
      IF(lwm) WRITE ( numond, nambgcext )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : nambgcext --- parameters for air-sea exchange'
         WRITE(numout,*) '      atm press. constant value, from ocean dynamics or sn file:   nn_patmint =', nn_patmint
         WRITE(numout,*) '      atm CO2 constant value, from annual file or sn file:   nn_co2int =', nn_co2int
      ENDIF
         
      ! ---------------   Atmospheric pressure inputs ---------------------
      SELECT  CASE (nn_patmint)
      CASE (0) ! atm press constant value
         IF(lwp) WRITE(numout,*) '         Constant atmospheric pressure: ', nn_patm, ' Pa'
      CASE (1) ! atm press from ocean dyn. 
         IF(lwp) WRITE(numout,*) '         Atmospheric pressure from ocean dynamics (sbcapr)'
         IF (.not.ln_apr_dyn)  CALL ctl_stop( 'trc_flx_init: Patm gradien need to be added on the ocean equations ', &
                                  &  'to use the option nn_patmint=1 (ln_apr_dyn  = .true.).' )
      CASE (2) ! atm press from inputs files
         IF(lwp) WRITE(numout,*) '         Atmospheric pressure from inputs files (sn_patm) '
         ALLOCATE( sf_patm(1), STAT=ierr )           !* allocate and fill sf_patm (forcing structure) with sn_patm
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_flx: unable to allocate sf_patm structure' )
         !
         CALL fld_fill( sf_patm, (/ sn_patm /), cn_dir, 'trc_flx', 'Atmospheric pressure ', 'nambgcext' )
                                ALLOCATE( sf_patm(1)%fnow(jpi,jpj,1)   )
         IF( sn_patm%ln_tint )  ALLOCATE( sf_patm(1)%fdta(jpi,jpj,1,2) )
      CASE (3) ! atm press from coupler
         IF(lwp) WRITE(numout,*) '         Atmospheric pressure from coupler (sn_rcv_mslp = coupled)'
         IF (.not.l_aprcpl)  CALL ctl_stop( 'trc_flx_init: sea level air pressure need to be coupled ', &
                                  &  'to use the option nn_patmint=3 (sn_rcv_mslp = coupled).' )
      ENDSELECT
      !  Correct nn_patmint if l_aprcpl
      IF (l_aprcpl.and.nn_patmint.ne.3) THEN
          CALL ctl_warn( 'trc_flx_init: sn_rcv_mslp = coupled but nn_patmint/=3', &
               &         '===> nn_patmint forced to 3' )
          nn_patmint=3
      ENDIF

      !
      ! ---------------   Co2 concentration inputs ---------------------
      SELECT  CASE (nn_co2int)
      CASE (0) ! atm co2 constant value
         IF(lwp) WRITE(numout,*) '         Constant CO2 concentration: ', nn_atmco2, ' ppm'
      CASE (1) ! atm co2 from annual file
         IF(lwp) WRITE(numout,*) '         CO2 concentration read from an annual file, clname  =', TRIM( clname )
         IF(lwp) WRITE(numout,*) '         Offset model-data start year              nn_offset =', nn_offset
         CALL ctl_opn( numco2, TRIM( clname) , 'OLD', 'FORMATTED', 'SEQUENTIAL', -1 , numout, lwp )
         jm = 0                      ! Count the number of record in co2 file
         DO
           READ(numco2,*,END=100) 
           jm = jm + 1
         END DO
 100     nmaxrec = jm - 1 
         ALLOCATE( years  (nmaxrec) )   ;   years  (:) = 0._wp
         ALLOCATE( atcco2h(nmaxrec) )   ;   atcco2h(:) = 0._wp
         !
         REWIND(numco2)
         DO jm = 1, nmaxrec          ! get  xCO2 data
            READ(numco2, *)  years(jm), atcco2h(jm)
            IF(lwp) WRITE(numout, '(f6.0,f7.2)')  years(jm), atcco2h(jm)
         END DO
         CLOSE(numco2)
      CASE (2) ! co2 from inputs files
         IF(lwp) WRITE(numout,*) '         CO2 concentration from inputs files (sn_atmco2) '
         ALLOCATE( sf_atmco2(1), STAT=ierr )           !* allocate and fill sf_atmco2 (forcing structure) with sn_atmco2
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_flx: unable to allocate sf_atmco2 structure' )
         !
         CALL fld_fill( sf_atmco2, (/ sn_atmco2 /), cn_dir, 'trc_flx', 'CO2 concentration ', 'nambgcext' )
                                  ALLOCATE( sf_atmco2(1)%fnow(jpi,jpj,1)   )
         IF( sn_atmco2%ln_tint )  ALLOCATE( sf_atmco2(1)%fdta(jpi,jpj,1,2) )
      CASE (3) ! CO2 conc. from coupler
         IF(lwp) WRITE(numout,*) '         CO2 concentration from coupler (sn_rcv_co2 = coupled)'
         IF (.not.l_co2cpl)  CALL ctl_stop( 'trc_flx_init: CO2 concentration need to be coupled', &
                                  &  'to use the option nn_co2int=3 (sn_rcv_co2 = coupled).' )
      ENDSELECT
      !  Correct nn_patm if l_aprcpl
      IF (l_co2cpl.and.nn_co2int.ne.3) THEN
          CALL ctl_warn( 'trc_flx_init: sn_rcv_co2 = coupled but nn_co2int/=3', &
               &         '===> nn_co2int forced to 3' )
          nn_co2int=3
      ENDIF



       ! Silic acid clim
       ! qasi3 is used in CO2 flux calculation using
       ! CALL trc_src3d( 1 , js3d_si)
       ! qasi3(:,:,1) = src3d_dta(:,:,1,js3d_si)

       no3_sf         = nn_no3_sf          !  Initialization of no3_sf
       oce_co2g(:,:)  = 0._wp                ! Initialization of Flux of CO2
       oce_o2g(:,:)  = 0._wp                ! Initialization of Flux of oxygen 
    
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
