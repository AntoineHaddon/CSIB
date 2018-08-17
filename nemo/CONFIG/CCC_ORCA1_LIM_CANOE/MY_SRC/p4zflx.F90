MODULE p4zflx
   !!======================================================================
   !!                         ***  MODULE p4zflx  ***
   !! TOP :   PISCES CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr) Include total atm P correction 
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_flx       :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!   p4z_flx_init  :   Read the namelist
   !!   p4z_patm      :   Read sfc atm pressure [atm] for each grid cell
   !!----------------------------------------------------------------------
   USE dom_oce, only  : nyear, nyear_len, nsec_year   
   USE oce_trc                      !  shared variables between ocean and passive tracers 
   USE trc                          !  passive tracers common variables
   USE sms_pisces                   !  PISCES Source Minus Sink variables
   USE p4zche                       !  Chemical model
   USE prtctl_trc                   !  print control for debugging
   USE trc_util
   USE iom                          !  I/O manager
   USE fldread                      !  read input fields
#if defined key_cpl_carbon_cycle
   USE sbc_oce, ONLY :  atm_co2     !  atmospheric pCO2               
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_flx  
   PUBLIC   p4z_flx_init  
   PUBLIC   p4z_flx_alloc  

   !                                      !!** Namelist  nampisext  **
   REAL(wp)          ::  atcco2    = 278._wp       !: pre-industrial atmospheric [co2] (ppm) 	
   LOGICAL           ::  ln_co2int = .FALSE.       !: flag to read in a file and interpolate atmospheric pco2 or not
   CHARACTER(len=34) ::  clname    = 'atcco2.txt'  !: filename of pco2 values
   INTEGER           ::  nn_offset = 0             !: Offset model-data start year (default = 0) 

   !!  Variables related to reading atmospheric CO2 time history    
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) :: atcco2h, atcco2h_years
   INTEGER  :: nmaxrec, numco2

   !                                         !!* nampisatm namelist (Atmospheric PRessure) *
   LOGICAL, PUBLIC ::   ln_presatm = .true.  !: ref. pressure: global mean Patm (F) or a constant (F)

   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:)  ::  patm      ! atmospheric pressure at kt                 [N/m2]
   TYPE(FLD), ALLOCATABLE,       DIMENSION(:)    ::  sf_patm   ! structure of input fields (file informations, fields read)


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: oce_co2   !: ocean carbon flux 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: satmco2   !: atmospheric pco2 

   REAL(wp) ::  t_oce_co2_flx               !: Total ocean carbon flux 
   REAL(wp) ::  t_atm_co2_flx               !: global mean of atmospheric pco2
   REAL(wp) ::  area                        !: ocean surface
   REAL(wp) ::  xconv  = 0.01_wp / 3600._wp !: coefficients for conversion 

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zflx.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_flx ( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx  ***
      !!
      !! ** Purpose :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
      !!
      !! ** Method  : 
      !!              - Include total atm P correction via Esbensen & Kushnir (1981) 
      !!              - Pressure correction NOT done for key_cpl_carbon_cycle
      !!              - Remove Wanninkhof chemical enhancement;
      !!              - Add option for time-interpolation of atcco2.txt  
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   !
      !
      INTEGER  ::   ji, jj, jm, iind, iindm1
      REAL(wp) ::   ztc, ztc2, ztc3, ztc4, zws, zkgwan
      REAL(wp) ::   zfld, zflu, zfld16, zflu16, zfact
      REAL(wp) ::   zph, zah2, zbot, zdic, zalk, zsch_o2, zalka, zsch_co2
      REAL(wp) ::   zph2, zph3, zpo4, zsi, zpd, zp0, zp1, zp3        ! coefficients added to account for P and Si contribution to TA
      REAL(wp) ::   zyr_dec, zdco2dt, current_yearfrac
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:) :: zkgco2, zkgo2, zh2co3, zoflx 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_flx')
      !
      CALL wrk_alloc( jpi, jpj, zkgco2, zkgo2, zh2co3, zoflx )
      !

      ! SURFACE CHEMISTRY (PCO2 AND [H+] IN
      !     SURFACE LAYER); THE RESULT OF THIS CALCULATION
      !     IS USED TO COMPUTE AIR-SEA FLUX OF CO2

      IF( kt /= nit000 ) CALL p4z_patm( kt )    ! Get sea-level pressure (E&K [1981] climatology) for use in flux calcs

      IF( ln_co2int ) THEN 
         ! Linear temporal interpolation  of atmospheric pco2.  atcco2.txt has annual values.
         ! Caveats: First column of .txt must be in years, decimal  years preferably. 
         ! For nn_offset, if your model year is iyy, nn_offset=(years(1)-iyy) 
         ! then the first atmospheric CO2 record read is at years(1)
         current_yearfrac = nyear + (nsec_year / ( nyear_len * 86400.))
         satmco2(:,:) = lin_interp( current_yearfrac + nn_offset, atcco2h_years, atcco2h )
      ENDIF

#if defined key_cpl_carbon_cycle
      satmco2(:,:) = atm_co2(:,:)
#endif

      DO jm = 1, 10
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi

               ! DUMMY VARIABLES FOR DIC, H+, AND BORATE
               zbot  = borat(ji,jj,1)
               zfact = rhop(ji,jj,1) / 1000. + rtrn
               zdic  = trn(ji,jj,1,jpdic) / zfact
               zph   = MAX( hi(ji,jj,1), 1.e-10 ) / zfact
               zalka = trn(ji,jj,1,jptal) / zfact
               zph2 = zph*zph
               zph3 = zph*zph2
               zpo4 = (trn(ji,jj,1,jpno3)+trn(ji,jj,1,jpnh4)) / 16. *0.000001 / zfact
               zsi = asi3(ji,jj,1) * 0.000001 / zfact                        ! silica is a static array based on initialization file, not a carried tracer

               ! CALCULATE P AND Si ION CONCENTRATIONS AS PER ORR ET AL (BPG EQUATIONS 43-47)
               ! zp3 = H3PO4, zp1 = HPO4(2-), zp0 = PO4(3-): denominator is the same for all 3 equations
               zpd = 1./ ( zph3 + akp13(ji,jj,1)*zph2 + akp13(ji,jj,1)*akp23(ji,jj,1)*zph + akp13(ji,jj,1)*akp23(ji,jj,1)*akp33(ji,jj,1) )
               zp3 = zph3*zpo4 * zpd
               zp1 = zph*zpo4*akp13(ji,jj,1)*akp23(ji,jj,1) * zpd
               zp0 = zpo4*akp13(ji,jj,1)*akp23(ji,jj,1)*akp33(ji,jj,1) * zpd
               zsi = zsi / (1. + zph / aksi3(ji,jj,1))

               ! CALCULATE [ALK]([CO3--], [HCO3-])
               zalk  = zalka - (  akw3(ji,jj,1) / zph - zph + zbot / ( 1.+ zph / akb3(ji,jj,1) ) + 2.*zp0 + zp1 - zp3 + zsi )

               ! CALCULATE [H+] AND [H2CO3]
               zah2   = SQRT(  (zdic-zalk)*(zdic-zalk) + 4.* ( zalk * ak23(ji,jj,1)   &
                  &                                        / ak13(ji,jj,1) ) * ( 2.* zdic - zalk )  )
               zah2   = 0.5 * ak13(ji,jj,1) / zalk * ( ( zdic - zalk ) + zah2 )
               zh2co3(ji,jj) = ( 2.* zdic - zalk ) / ( 2.+ ak13(ji,jj,1) / zah2 ) * zfact
               hi(ji,jj,1)   = zah2 * zfact
            END DO
         END DO
      END DO


      ! --------------
      ! COMPUTE FLUXES
      ! --------------

      ! FIRST COMPUTE GAS EXCHANGE COEFFICIENTS
      ! -------------------------------------------

!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
!            ztc  = MIN( 35., tsn(ji,jj,1,jp_tem) )
            ztc2 = ztc * ztc
            ztc3 = ztc * ztc2 
            ztc4 = ztc * ztc3 

            ! Compute the schmidt Number both O2 and CO2
            zsch_co2 = 2116.8 - 136.25 * ztc + 4.7353 * ztc2 - 0.092307 * ztc3 + 0.0007555 * ztc4
            zsch_o2  = 1920.4 - 135.6  * ztc + 5.2122 * ztc2 - 0.10939  * ztc3 + 0.00093777 * ztc4
            !  wind speed 
            zws  = wndm(ji,jj) * wndm(ji,jj)
            ! Compute the piston velocity for O2 and CO2
            zkgwan = 0.251 * zws  
            zkgwan = zkgwan * xconv * ( 1.- fr_i(ji,jj) ) * tmask(ji,jj,1)
            ! compute gas exchange for CO2 and O2
            zkgco2(ji,jj) = zkgwan * SQRT( 660./ zsch_co2 )
            zkgo2 (ji,jj) = zkgwan * SQRT( 660./ zsch_o2 )
         END DO
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Compute CO2 flux for the sea and air
            zfld = satmco2(ji,jj) * patm(ji,jj) * tmask(ji,jj,1) * chemc(ji,jj,1) * zkgco2(ji,jj)   ! (mol/L) * (m/s)
            zflu = zh2co3(ji,jj) * tmask(ji,jj,1) * zkgco2(ji,jj)                                   ! (mol/L) (m/s) ?
            oce_co2(ji,jj) = ( zfld - zflu ) * rfact * e1e2t(ji,jj) * tmask(ji,jj,1) * 1000.
            ! compute the trend
            tra(ji,jj,1,jpdic) = tra(ji,jj,1,jpdic) + ( zfld - zflu ) / fse3t(ji,jj,1)

            ! Compute O2 flux 
            zfld16 = atcox * patm(ji,jj) * chemc(ji,jj,2) * tmask(ji,jj,1) * zkgo2(ji,jj)          ! (mol/L) * (m/s)
            zflu16 = trn(ji,jj,1,jpoxy) * tmask(ji,jj,1) * zkgo2(ji,jj) * 1.e-6   ! convert to mol/L
            zoflx(ji,jj) = zfld16 - zflu16
            tra(ji,jj,1,jpoxy) = tra(ji,jj,1,jpoxy) + zoflx(ji,jj) / fse3t(ji,jj,1) * 1.e6   ! convert to mmol m^-3
         END DO
      END DO

      t_oce_co2_flx = t_oce_co2_flx + glob_sum( oce_co2(:,:) )            ! Cumulative Total Flux of Carbon
      IF( kt == nitend ) THEN
         t_atm_co2_flx = glob_sum( satmco2(:,:) * patm(:,:) * e1e2t(:,:) )            ! Total atmospheric pCO2
         !
         t_oce_co2_flx = (-1.) * t_oce_co2_flx  * 12. / 1.e15             ! Conversion in PgC ; negative for out of the ocean
         t_atm_co2_flx = t_atm_co2_flx  / area                            ! global mean of atmospheric pCO2
         !
         IF( lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) ' Global mean of atmospheric pCO2 (ppm) at it= ', kt, ' date= ', ndastp
            WRITE(numout,*) '------------------------------------------------------- :  ',t_atm_co2_flx
            WRITE(numout,*)
            WRITE(numout,*) ' Cumulative total Flux of Carbon out of the ocean (PgC) :'
            WRITE(numout,*) '-------------------------------------------------------  ',t_oce_co2_flx
         ENDIF
         !
      ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('flx ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      IF( ln_diatrc ) THEN
         IF( lk_iomput ) THEN
            CALL iom_put( "Cflx" , oce_co2(:,:) / e1e2t(:,:) / rfact ) 
            CALL iom_put( "Oflx" , zoflx(:,:) * 1000 * tmask(:,:,1)  )
            CALL iom_put( "Kg"   , zkgco2(:,:) * tmask(:,:,1) )
            CALL iom_put( "Dpco2", ( zh2co3(:,:) / ( chemc(:,:,1) + rtrn ) - satmco2(:,:) * patm(:,:) ) * tmask(:,:,1) )
            CALL iom_put( "pco2" , ( zh2co3(:,:) / ( chemc(:,:,1) + rtrn ) ) * tmask(:,:,1) )
            CALL iom_put( "Dpo2" , ( atcox * patm(:,:) - trn(:,:,1,jpoxy) / ( chemc(:,:,2) + rtrn ) ) * tmask(:,:,1) )
         ELSE
            trc2d(:,:,jp_pcs0_2d    ) = oce_co2(:,:) / e1e2t(:,:) / rfact 
            trc2d(:,:,jp_pcs0_2d + 1) = zoflx(:,:) * 1000 * tmask(:,:,1) 
            trc2d(:,:,jp_pcs0_2d + 2) = zkgco2(:,:) * tmask(:,:,1) 
            trc2d(:,:,jp_pcs0_2d + 3) = ( satmco2(:,:) * patm(:,:) - zh2co3(:,:) / ( chemc(:,:,1) + rtrn ) ) * tmask(:,:,1) 
         ENDIF
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zkgco2, zkgo2, zh2co3, zoflx )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_flx')
      !
   END SUBROUTINE p4z_flx


   SUBROUTINE p4z_flx_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_flx_init  ***
      !!
      !! ** Purpose :   Initialization of atmospheric conditions
      !!
      !! ** Method  :   Read the nampisext namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !! ** input   :   Namelist nampisext
      !!----------------------------------------------------------------------
      NAMELIST/nampisext/ln_co2int, atcco2, clname, nn_offset
      INTEGER :: jm
      !!----------------------------------------------------------------------
      !
      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampisext )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for air-sea exchange, nampisext'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Choice for reading in the atm pCO2 file or constant value, ln_co2int =', ln_co2int
         WRITE(numout,*) ' '
      ENDIF
      IF( .NOT.ln_co2int ) THEN
         IF(lwp) THEN                         ! control print
            WRITE(numout,*) '    Constant Atmospheric pCO2 value  atcco2    =', atcco2
            WRITE(numout,*) ' '
         ENDIF
         satmco2(:,:)  = atcco2      ! Initialisation of atmospheric pco2
      ELSE
         IF(lwp)  THEN
            WRITE(numout,*) '    Atmospheric pCO2 value from file       clname         =', TRIM( clname )
            WRITE(numout,*) '    Atmospheric pCO2 variable name in file clvarname      =', TRIM( clvarname )
            WRITE(numout,*) '    Offset model-data start year           nn_offset   =', nn_offset
            WRITE(numout,*) ' '
         ENDIF
         CALL chkerr(nf90_open( clname, NF90_NOWRITE, ncid ), 'trcini_cfc', 0)
         CALL read_var2d( ncid, clvarname, tmp2d )
         CALL chkerr(nf90_close( ncid )) 
         ! Set the time-varying atmospheric history from the read in data
         ALLOCATE(atcco2h(ntime))
         atcco2h(:) = tmp2d(:,1) ! Sector '1' corresponds to global average
         DEALLOCATE(tmp2d)
         ! Input file for OMIP6 is in Gregorian days since 1 January 0000, manually overwrite
         ! so that atcco2h_years is in yearfraction
         DO jm = 1,ntime
            atcco2h_years(jm) = (jm-1) + 0.5
         ENDDO
         
      ENDIF
      !
      area = glob_sum( e1e2t(:,:) )        ! interior global domain surface
      !
      oce_co2(:,:)  = 0._wp                ! Initialization of Flux of Carbon
      t_atm_co2_flx = 0._wp
      t_oce_co2_flx = 0._wp
      !
      CALL p4z_patm( nit000 )
      !
   END SUBROUTINE p4z_flx_init

   SUBROUTINE p4z_patm( kt )

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_atm  ***
      !!
      !! ** Purpose :   Read and interpolate the external atmospheric sea-levl pressure
      !! ** Method  :   Read the files and interpolate the appropriate variables
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step
      !
      INTEGER            ::  ierr
      CHARACTER(len=100) ::  cn_dir   ! Root directory for location of ssr files
      TYPE(FLD_N)        ::  sn_patm  ! informations about the fields to be read
      !!
      NAMELIST/nampisatm/ ln_presatm, sn_patm, cn_dir

      !                                         ! -------------------- !
      IF( kt == nit000 ) THEN                   ! First call kt=nittrc000 !
         !                                      ! -------------------- !
         !                                            !* set file information (default values)
         ! ... default values (NB: frequency positive => hours, negative => months)
         !            !   file   ! frequency !  variable  ! time intep !  clim   ! 'yearly' or ! weights  ! rotation !
         !            !   name   !  (hours)  !   name     !   (T/F)    !  (T/F)  !  'monthly'  ! filename ! pairs    !
         sn_patm = FLD_N( 'pres'  ,    24     ,  'patm'    ,  .false.   , .true.  ,   'yearly'  , ''       , ''       )
         cn_dir  = './'          ! directory in which the Patm data are 

         REWIND( numnatp )                             !* read in namlist nampisatm
         READ  ( numnatp, nampisatm ) 
         !
         !
         IF(lwp) THEN                                 !* control print
            WRITE(numout,*)
            WRITE(numout,*) '   Namelist nampisatm : Atmospheric Pressure as external forcing'
            WRITE(numout,*) '      constant atmopsheric pressure (F) or from a file (T)  ln_presatm = ', ln_presatm
            WRITE(numout,*)
         ENDIF
         !
         IF( ln_presatm ) THEN
            ALLOCATE( sf_patm(1), STAT=ierr )           !* allocate and fill sf_patm (forcing structure) with sn_patm
            IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_flx: unable to allocate sf_patm structure' )
            !
            CALL fld_fill( sf_patm, (/ sn_patm /), cn_dir, 'p4z_flx', 'Atmospheric pressure ', 'nampisatm' )
                                   ALLOCATE( sf_patm(1)%fnow(jpi,jpj,1)   )
            IF( sn_patm%ln_tint )  ALLOCATE( sf_patm(1)%fdta(jpi,jpj,1,2) )
         ENDIF
         !                                         
         IF( .NOT.ln_presatm )   patm(:,:) = 1.e0    ! Initialize patm if no reading from a file
         !
      ENDIF
      !
      IF( ln_presatm ) THEN
         CALL fld_read( kt, 1, sf_patm )               !* input Patm provided at kt + 1/2
         patm(:,:) = sf_patm(1)%fnow(:,:,1)                        ! atmospheric pressure
      ENDIF
      !
   END SUBROUTINE p4z_patm

   INTEGER FUNCTION p4z_flx_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( oce_co2(jpi,jpj), satmco2(jpi,jpj), patm(jpi,jpj), STAT=p4z_flx_alloc )
      !
      IF( p4z_flx_alloc /= 0 )   CALL ctl_warn('p4z_flx_alloc : failed to allocate arrays')
      !
   END FUNCTION p4z_flx_alloc

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_flx( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p4z_flx: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_flx
#endif 

   !!======================================================================
END MODULE  p4zflx
