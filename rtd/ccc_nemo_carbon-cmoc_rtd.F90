PROGRAM nemo_ocean_diag
IMPLICIT NONE 
      
! ======================================================================
!  Purpose: Run-time diagnostics for NEMO (ORCA2) 
!
! HISTORY:
! -------
!
! N. Swart    Nov     2015 1. Remove annual mean calculation and rewrite
!                             all code to operate on monthly data.
!                          2. Remove North Fold point from computions.
!                            (i.e. sum to jmt -1)
!                          3. Rewrite of functions and code to freefrom. 
!                          4. Improve time axis in output NetCDF, include
!                            noleap_days subroutine.
!
!  N. Swart   Jul 15  2014 Update to CMOC only variables (exclude PISCES vars).
!  N. Swart   May 07  2014 Update to standard CMOC/CanESM2 RTD variable set. Major style revision to F90.
!  N. Swart   May 02  2014 Made resolution independent.
!  N. Swart   Apr 30  2014 Switched output format from CCCma to NetCDF.
!  N. Swart   Apr 29  2014 Complete re-write to remove extraneous variables and to improve ordering.
!  N. Swart   Apr 28  2014 Stripped out the physical ocean, fixed bugs, converted to ORCA1 section.
!  O.Riche    Feb 02. 2014 PISCES version 
!  O.Riche    Jan 11, 2014 Expand the CMOC set
!  O.Riche    Dec 18, 2013 Add a CMOC variable to test diagnostic of BGC tracers
!  WGL Sept 25, 2013 Conversion for output in CCCma format timeseries
!  O. Saenko (Sept 25, 2013)
!  O. Saenko (May 22, 2013)
!
! USAGE
! -----
! 
! nemo_ocean_diag.exe YYYY, MM
!
! where the first command line arg, YYYY, is the RTD year, and MM is the FIRST month in this RTD sequence.
!
! INPUT FILES
! -----------
! NEMO_PISCES NetCDF files, with the names:
!
!    - orca_mesh_mask
!    - diat_t : monthly frequency (_1m_)
!    - ptrc_t : monthly frequency (_1m_)
!
! OUTPUT FILES
! ------------
! nemo_carbon_rtd.nc - NetCDF output file with monthly timeseries for carbon variables.
!
!
! ======================================================================
! to compile:
!
! 1. Use build-nemo-rtd
!
! 2. xlf90_r -o nemo_physical_rtd.exe \
!   ccc_nemo_physical_rtd.F90 uvic_netcdf.f `nf-config --fflags --flibs`
! ======================================================================
      INTEGER :: imt, jmt, km, lm, iou

!         establish the size of the grid from the input file.
          CALL openfile ("ptrc_t", iou)
          CALL getdimlen ('x', iou, imt)
          CALL getdimlen ('y', iou, jmt)
          CALL getdimlen ('deptht', iou, km)
          call getdimlen ('time_counter', iou, lm)

!         do the calculations and save the output netcdf   
          CALL calc (imt, jmt, km, lm)

END PROGRAM nemo_ocean_diag

SUBROUTINE calc (imt, jmt, km, lm)
!     Does the required calculations and saves the output to netcdf

      IMPLICIT NONE
      integer, parameter:: dp=kind(0.d0) ! double precision


! ======================================================================
!     Input data 
! ======================================================================
      INTEGER imt, jmt, km, lm
      INTEGER i, j, k, l

!     Grid-related arrays
      REAL, DIMENSION(imt, jmt)     :: lon2d, lat2d, e1t, e2t
      REAL, DIMENSION(km)           :: deptht
      REAL, DIMENSION(imt, jmt, km) :: e3t(imt,jmt,km) 

!     Monthly DIC, CaCO3, TA, PH, O2
      REAL, DIMENSION(imt, jmt, km, lm) :: dic, caco3, tal, ph, oxy

!     Monthly POC, GOC, DOC
      REAL, DIMENSION(imt, jmt, km, lm) :: poc, goc, doc

!     Monthly NO3, NH4, PO4, Si
      REAL, DIMENSION(imt, jmt, km, lm) :: no3, nh4, po4, si

!     Monthly PHY and Zoo
      REAL, DIMENSION(imt, jmt, km, lm) :: phy, phy2, zoo, zoo2

!     Monthly primary production: PPPHY, PPPHY2 
      REAL, DIMENSION(imt, jmt, km, lm) :: ppphy, ppphy2

!     Monthly export fluxes of C (EPC100), CaCO3 (EPCAL100)
      REAL, DIMENSION(imt, jmt, lm)     :: epc100, epcal100

!     Monthly surface fluxes of DIC, O2, N2, Fe
      REAL, DIMENSION(imt, jmt, lm)     :: cflux, oflux, nfix, irondep

! ======================================================================
!     Pre-computed data  
! ======================================================================

!     3D variables
      REAL, DIMENSION(imt, jmt, km) :: t_mask

! ======================================================================
!     Working arrays / variables  
! ======================================================================
      REAL :: dum, dvol, vol 
      REAL :: dicz, caco3z, talz, phz, oxyz, pocz, gocz, docz,no3z, nh4z
      REAL :: po4z, siz, phyz, phy2z, zooz, zoo2z, ppphyz, ppphy2z

!     total ocean carbon, nitrogen      
      REAL, DIMENSION(lm) :: toc, ton
 
! g_mask
      REAL, DIMENSION(imt, jmt) :: g_mask

! ======================================================================
!     Output data 
! ======================================================================
! (1) Global-mean profiles for 3D data:
      REAL, DIMENSION(km, lm) :: dic_z, caco3_z, tal_z, ph_z, oxy_z, poc_z
      REAL, DIMENSION(km, lm) :: goc_z, doc_z, no3_z, nh4_z, po4_z, si_z, phy_z, phy2_z
      REAL, DIMENSION(km, lm) :: zoo_z, zoo2_z, ppphy_z, ppphy2_z

! (2) Global-mean (volume weighted) or integral

!     DIC, CaCO3 TA, PH, O2
      REAL, DIMENSION(lm) :: dicvol, caco3vol, talvol, phvol, oxyvol, pocvol, gocvol
      REAL, DIMENSION(lm) :: docvol, no3vol, nh4vol, po4vol, sivol, phyvol, phy2vol
      REAL, DIMENSION(lm) :: zoovol, zoo2vol, ppphyvol, ppphy2vol
      REAL, DIMENSION(lm) :: epc100glo, epcal100glo, cglo, ofluxglo, nfixglo, irondepglo

!----------------
!  NetCDF-output specific
      integer id_time, id_z, iou, ntrec, ntrec2, iyear, imon
      integer days_elapsed
      logical exists, exists1, notopen
      real tyear, tdays_elapsed
      CHARACTER(len=32) :: year_arg_in, mon_arg_in

!----------------
!     input file stuff
      character fname05*100, fname06*100, fname07*100  
      integer year, iou4, iou5, iou6, recn, nrecon

         year =0
         ntrec=0
         iou4 =0
         iou5 =0 
         iou6 =0 
         recn =12.
         nrecon = int(recn + 0.001)
!---------------------------------------------------
!    Define NetCDF files   
!---------------------------------------------------
        fname05='orca_mesh_mask'
        fname06='ptrc_t'   
        fname07='diad_t'   
!---------------------------------------------------
!    Get grid/mask data   
!---------------------------------------------------
      print*,'Reading data on NEMO grid...'
      print*,''
      print*,'Grid size: imt, jmt, km:', imt, jmt, km

      CALL openfile (fname05,iou4)
      CALL getvara ('e1t', iou4, imt*jmt, (/1,1,1/), (/imt,jmt,1/),e1t , 1., 0.)
      CALL getvara ('e2t', iou4, imt*jmt, (/1,1,1/),  (/imt,jmt,1/),e2t , 1., 0.)
      CALL getvara ('e3t', iou4, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),e3t , 1., 0.)
      CALL getvara ('tmask', iou4, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/), t_mask , 1., 0.)
      CALL closefile (iou4)

!     get some more grid information
      CALL openfile (fname06,iou5)
      CALL getvara ('nav_lon', iou5, imt*jmt, (/1,1,1/), (/imt,jmt,1/), lon2d, 1., 0.)
      CALL getvara ('nav_lat', iou5, imt*jmt, (/1,1,1/), (/imt,jmt,1/), lat2d, 1., 0.)
      CALL getvara ('deptht', iou5, km, (/1/), (/km/), deptht, 1., 0.)
      CALL closefile (iou5)

!---------------------------------------------------
! Read in from NetCDF
!---------------------------------------------------

!     read prognostic variables from the prtc_t file        
      CALL openfile (fname06,iou5)

!     read diagnostic variables from the diad_t file, if it exists
      inquire (file=trim(fname07), exist=exists)                          
      if (exists) then         
          CALL openfile (fname07,iou6)
      endif

!  DIC, TA, O2 
      CALL getvara('DIC', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), dic, 1., 0.)                                 
      CALL getvara('Alkalini', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), tal, 1., 0.)   

!  PH moved below for reading with other diat_t input

      CALL getvara('O2', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), oxy, 1., 0.)   

! POC, GOC, DOC
      CALL getvara('POC', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), poc, 1., 0.)   

!  NO3, NH4, PO4, Si 
      CALL getvara('NO3', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), no3, 1., 0.)

!  PHY, PHY2, ZOO, ZOO2
      CALL getvara('PHY',  iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), phy, 1., 0.)  
      CALL getvara('ZOO' , iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), zoo, 1., 0.)    
      
!  Diagnostic variables
      if (exists) then 
!       3-D: PH, PPPHY, PPPHY2, EPC100,
          CALL getvara('PH',       iou6, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), ph, 1., 0.)   
          CALL getvara('PPPHY',    iou6, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), ppphy, 1., 0.)   

!       2-D :  EPCAL100, DIC flux, Oflux, Nfix, Irondep
          CALL getvara('EPC100',   iou6, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/),   epc100, 1., 0.)    
          CALL getvara('EPCAL100', iou6, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/), epcal100, 1., 0.)    
          CALL getvara('Cflx',     iou6, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/),    cflux, 1., 0.)    
          CALL getvara('Oflx',     iou6, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/),    oflux, 1., 0.)    
          CALL getvara('Nfix',     iou6, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/),     nfix, 1., 0.)    
      endif 

      CALL closeall ! close all open netcdf files

! ********** Do some basic calculations ************

!---------------------------------------------------
! (1) Global annual mean T(z) and S(z)  
!---------------------------------------------------

!  DIC, TA, PH, O2 
      dicvol(:)   = 0.0_dp
      caco3vol(:) = 0.0_dp
      talvol(:)   = 0.0_dp
      phvol(:)    = 0.0_dp
      oxyvol(:)   = 0.0_dp
!  POC, GOC, DOC
      pocvol(:)   = 0.0_dp
      gocvol(:)   = 0.0_dp
      docvol(:)   = 0.0_dp
!  NO3, NH4, PO4, Si 
      no3vol(:)   = 0.0_dp
      nh4vol(:)   = 0.0_dp  
      po4vol(:)   = 0.0_dp  
      sivol(:)    = 0.0_dp 
!  PHY, PHY2, ZOO, ZOO2
      phyvol(:)   = 0.0_dp 
      phy2vol(:)  = 0.0_dp 
      zoovol(:)   = 0.0_dp 
      zoo2vol(:)  = 0.0_dp 
!  PPPHY, PPPHY2      
      ppphyvol(:) = 0.0_dp
      ppphy2vol(:) = 0.0_dp
                       
      do l = 1, lm  
          vol = 0.0_dp
          do k = 1, km   
             g_mask(:, :)  = t_mask(:, :, k) 

    !        DIC, TA, PH, O2 
             CALL area_ave(e1t, e2t, e3t, g_mask, dic(:, :, k, l), imt, jmt, km, dicz,   dvol, k)  
             CALL area_ave(e1t, e2t, e3t, g_mask, tal(:, :, k, l), imt, jmt, km, talz,   dvol, k)  
             CALL area_ave(e1t, e2t, e3t, g_mask, oxy(:, :, k, l), imt, jmt, km, oxyz,   dvol, k)

    !        POC, GOC, DOC
             CALL area_ave(e1t, e2t, e3t, g_mask, poc(:, :, k, l), imt, jmt, km, pocz, dvol, k)

    !        NO3, NH4, PO4, Si 
             CALL area_ave(e1t, e2t, e3t, g_mask, no3(:, :, k, l), imt, jmt, km, no3z, dvol, k)

    !        PHY, PHY2, ZOO, ZOO2
             CALL area_ave(e1t, e2t, e3t, g_mask, phy(:, :, k, l),  imt, jmt, km, phyz,  dvol, k)  
             CALL area_ave(e1t, e2t, e3t, g_mask, zoo(:, :, k, l),  imt, jmt, km, zooz,  dvol, k)  

             if (exists) then
    !            PPPHY, PPPHY2      
                 CALL area_ave(e1t, e2t, e3t, g_mask, ph(:, :, k, l),     imt, jmt, km, phz,     dvol, k)  
                 CALL area_ave(e1t, e2t, e3t, g_mask, ppphy(:, :, k, l),  imt, jmt, km, ppphyz,  dvol, k) 
             endif 
    !================================================================
    !        Assign outputs
    !================================================================
    !        DIC, TA, PH, O2 
             dic_z(k, l)   = dicz     
             tal_z(k, l)   = talz     
             oxy_z(k, l)   = oxyz    

    !        POC, GOC, DOC
             poc_z(k, l)    = pocz       

    !        NO3, NH4, PO4, Si 
             no3_z(k, l)   = no3z  

    !        PHY, PHY2, ZOO, ZOO2
             phy_z(k, l)   = phyz
             zoo_z(k, l)   = zooz
         
             if (exists) then
                 ph_z(k, l)    = phz  
                 ppphy_z(k, l) = ppphyz
             endif 

    !        DIC, TA, PH, O2 
             print*, dicz, dvol
             dicvol(l)   = dicvol(l) + dicz*dvol  
             talvol(l) = talvol(l) + talz*dvol  
             oxyvol(l) = oxyvol(l) + oxyz*dvol  

    !        POC, GOC, DOC
             pocvol(l) = pocvol(l) + pocz*dvol  

    !        NO3, NH4, PO4, Si 
             no3vol(l) = no3vol(l) + no3z*dvol  

    !        PHY, PHY2, ZOO, ZOO2
             phyvol(l)    = phyvol(l)  + phyz*dvol
             zoovol(l)    = zoovol(l)  + zooz*dvol

             if (exists) then 
                 phvol(l)     =     phvol(l) + phz*dvol    
                 ppphyvol(l)  =  ppphyvol(l) + ppphyz*dvol  
             endif
                       
              vol = vol + dvol
          enddo  ! depth, k        

    !     compute toc and ton
          toc(l) = dicvol(l)  + pocvol(l) + phyvol(l) + zoovol(l)                                          
    !     convert from mmol C to Pg C      
          toc(l) = toc(l) * 12.0e-18
          print*,'toc', toc(l)
          ton(l) = no3vol(l) + 16./122.*(phyvol(l) + zoovol(l)  + pocvol(l))                                             
    !     convert to Pg      
          ton(l) = ton(l) * 14.007e-18

          if (vol.ne.0.) then         
    !        DIC, TA, PH, O2 
             dicvol(l) = dicvol(l)/vol 
             talvol(l) = talvol(l)/vol 
             oxyvol(l) = oxyvol(l)/vol  
    !        POC, GOC, DOC  
             pocvol(l) = pocvol(l)/vol  
             no3vol(l) = no3vol(l)/vol  

             if (exists) then   
                 phvol(l)  = phvol(l)/vol 
    !            convert to PgC/yr
                 ppphyvol(l)  = ppphyvol(l) * 12.e-15 * 86400 * 365
             endif
          endif


    !---------------------------------------------------
    ! (2) Global surface fields (fluxes, etc...)  
    !---------------------------------------------------
          do i=1,imt
              do j=1,jmt
                  g_mask(i,j) = t_mask(i,j,1)
              enddo
          enddo

          if (exists) then 
    !        EPC100, EPCAL100  
             CALL area_ave_flx(e1t, e2t, g_mask, epc100(:,:,l),   imt, jmt,   epc100glo(l), dum) 
             CALL area_ave_flx(e1t, e2t, g_mask, epcal100(:,:,l), imt, jmt, epcal100glo(l), dum) 

    !        Cflux, Oflux, Nfix, Irondep  
             CALL area_ave_flx (e1t, e2t, g_mask, cflux(:,:,l),   imt, jmt, cglo(l),       dum) 
             CALL area_ave_flx (e1t, e2t, g_mask, oflux(:,:,l),   imt, jmt, ofluxglo(l),   dum) 
             CALL area_ave_flx (e1t, e2t, g_mask, nfix(:,:,l),    imt, jmt, nfixglo(l),   dum) 

    !       <PISCES OR 01/15/2014> convert into PgC/yr 
             cglo(l)         = cglo(l)        * dum * 12.e-15 * 86400 * 365
             epc100glo(l)    = epc100glo(l)   * dum * 12.e-15 * 86400 * 365
             epcal100glo(l)  = epcal100glo(l) * dum * 12.e-15 * 86400 * 365

    !        convert to TgN/yr (assuming this is N not N2)
             nfixglo(l)      = nfixglo(l)  * dum * 14.007e-12 * 86400 * 365

    !        convert to mol O2/yr
             ofluxglo(l)     = ofluxglo(l) * dum *  86400 * 365
          endif
      enddo ! time loop


!---------------------------------------------------------
!     NETCDF RTD output: Time series information section
!---------------------------------------------------------
!     Read in the year which is the first command line argument
      CALL getarg(1, year_arg_in )
      read (year_arg_in,'(I10)') iyear
!     Read in the 1st month which is the second command line argument
      CALL getarg(2, mon_arg_in )
      read (mon_arg_in,'(I10)') imon 
      print*," --- "
      print*," iyear, imon are:", iyear, imon 
      print*," --- "

      iou = 0
      id_time = 0
      id_z = 0

!     If the output file does not exist, create it and define dims and vars
      inquire (file="nemo_carbon_rtd.nc", exist=exists1)
      if (.not. exists1) then
      print*,"output file not found...creating a new file..."
      !CALL flush(6)
      CALL opennew ("nemo_carbon_rtd.nc", iou)
      ntrec = 1
        CALL redef (iou)
!         basic grid specification
          CALL defdim ('time', iou, 0, id_time)
          CALL defdim ('depth', iou, km, id_z)
          CALL defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'F'   &
            &        , 'time', 'time', 'days since 0000-01-01 00:00:00')

          call putatttext (iou, 'time', 'calendar', '365_day')

          CALL defvar ('depth', iou, 1, (/id_z/), 0., 0., 'Y', 'F'               &
     &       , 'depth of the t grid', 'depth', 'm')

!         DIC
          CALL defvar ('DIC', iou, 1, (/id_time/), -1.e4                         & 
     &        , 1.e4,' ', 'F', 'Global mean Dissolved Inorganic Carbon'          &
     &        , 'DIC', 'mmol m^-3')

          CALL defvar ('DICz', iou, 2, (/id_z, id_time/), -1.e4                  &
     &        , 1.e4,' ', 'F', 'Dissolved Inorganic Carbon by level'             &
     &        , 'DICz', 'mmol m^-3')

!         CaCO3
!          CALL defvar ('CaCO3', iou, 1, (/id_time/), -1.e4                       &  
!     &        , 1.e4,' ', 'F', 'Global mean calcite Concentration'               & 
!     &        , 'CaCO3', 'mmol m^-3')

!          CALL defvar ('CaCO3z', iou, 2, (/id_z, id_time/), -1.e4                &
!     &        , 1.e4,' ', 'F', 'Calcite concentration by level'                  &
!     &        , 'CaCO3z', 'mmol m^-3')

!         ALK
          CALL defvar ('TAL', iou, 1, (/id_time/), -1.e4                         &
     &        , 1.e4,' ', 'F', 'Global mean Total Alkalinity'                    &
     &        , 'TAL', 'mmol m^-3')

          CALL defvar ('TALz', iou, 2, (/id_z, id_time/), -1.e4                  &
     &        , 1.e4,' ', 'F', 'Total Alkalinity by level'                       &  
     &        , 'TALz', 'mmol m^-3')

!         O2
          CALL defvar ('O2', iou, 1, (/id_time/), -1.e4                          &
     &        , 1.e4,' ', 'F', 'Global mean O2'                                  &
     &        , 'O2', 'uM')

          CALL defvar ('O2z', iou, 2, (/id_z, id_time/), -1.e4                   &
     &        , 1.e4,' ', 'F', 'O2 by level'                                     &
     &        , 'O2z', 'uM')

!         POC
          CALL defvar ('POC', iou, 1, (/id_time/), -1.e4                         &
     &        , 1.e4,' ', 'F'                                                    & 
     &        , 'Global mean small Particulate Organic Carbon'                   & 
     &        , 'POC', 'mmol m^-3')

          CALL defvar ('POCz', iou, 2, (/id_z, id_time/), -1.e4                  &
     &        , 1.e4,' ', 'F', 'Small POC by level'                              &
     &        , 'POCz', 'mmol m^-3')

!         NO3
          CALL defvar ('NO3', iou, 1, (/id_time/), -1.e4                         &
     &        , 1.e4,' ', 'F', 'Global mean NO3'                                 &  
     &        , 'NO3', 'uM')

          CALL defvar ('NO3z', iou, 2, (/id_z, id_time/), -1.e4                  &
     &        , 1.e4,' ', 'F', 'NO3 by level'                                    & 
     &        , 'NO3z', 'uM')

!         TC
          CALL defvar ('TC', iou, 1, (/id_time/), -1.e4                          &
     &        , 1.e4,' ', 'F', 'Total Ocean Carbon'                              & 
     &        , 'TC', 'Pg')
!         TN
          CALL defvar ('TN', iou, 1, (/id_time/), -1.e4                          & 
     &        , 1.e4,' ', 'F', 'Total Ocean (fixed) Nitrogen'                    & 
     &        , 'TN', 'Pg')

          if (exists) then
!             PH
              CALL defvar ('PH', iou, 1, (/id_time/), -1.e4                      &
     &            , 1.e4,' ', 'F', 'Global mean pH'                              & 
     &            , 'PH', '')

              CALL defvar ('PHz', iou, 2, (/id_z, id_time/), -1.e4               &
     &            , 1.e4,' ', 'F', 'pH by level'                                 & 
     &            , 'PHz', '')

!             PPPHY
              CALL defvar ('PPPHY', iou, 1, (/id_time/), -1.e4                   &
     &            , 1.e4,' ', 'F', 'Primary production of nanophyto'             & 
     &            , 'PPPHY', 'PgC/yr')

!             EPC100
              CALL defvar ('EPC100', iou, 1, (/id_time/), -1.e4                  &  
     &            , 1.e4,' ', 'F', 'Export of carbon particles at 100m'          &  
     &            , 'EPC100', 'PgC/yr')

!             EPCAL100
              CALL defvar ('EPCAL100', iou, 1, (/id_time/), -1.e4                &  
     &            , 1.e4,' ', 'F', 'Export of Calcite at 100m'                   &
     &            , 'EPCAL100', 'PgC/yr')

!             Cflux
              CALL defvar ('CFLX', iou, 1, (/id_time/), -1.e4                    &
     &            , 1.e4,' ', 'F', 'Global mean surface flux of DIC'             &
     &            , 'CFLX', 'PgC/yr')
!             Oflux
              CALL defvar ('OFLX', iou, 1, (/id_time/), -1.e4                    & 
     &           , 1.e4,' ', 'F', 'Global mean surface flux of oxygen'           &
     &           , 'OFLX', 'mol/yr')
!             Nfix
              CALL defvar ('NFIX', iou, 1, (/id_time/), -1.e4                    & 
     &            , 1.e4,' ', 'F', 'Nitrogen fixation at surface'                &
     &            , 'NFIX', 'TgN/yr')
          endif

          CALL enddef (iou)
!         define the depth axis
          CALL putvara ('depth', iou, km, (/1/), (/km/)                          & 
     &      , deptht, 1., 0.)
 
      else

!      if the file does exist then open it for writing at the next record
       print*,"output file found...opening existing file for appending"
       !CALL flush(6)
       CALL opennext ("nemo_carbon_rtd.nc", tyear, ntrec, iou)
      endif

!       append variables
!--------------------------------------------------------------------------
      do l = 1, lm
!       Convert the date into days since 01-01-0001        
        CALL noleap_days(iyear, imon+l-1, 1, days_elapsed)
        tdays_elapsed = float(days_elapsed)
!        print*, days_elapsed
!       time
        ntrec2 = ntrec + l - 1
        call putvars ('time', iou, ntrec2, tdays_elapsed, 1., 0.)

!       DIC
        CALL putvars ('DIC', iou, ntrec2, dicvol(l), 1., 0.)
        CALL putvara ('DICz', iou, km, (/1, ntrec2/), (/km, 1/), dic_z(:,l), 1., 0.) 

!       ALK
        CALL putvars ('TAL', iou, ntrec2, talvol(l), 1., 0.)
        CALL putvara ('TALz', iou, km, (/1, ntrec2/), (/km, 1/), tal_z(:,l), 1., 0.)

!       OXY
        CALL putvars ('O2', iou, ntrec2, oxyvol(l), 1., 0.)
        CALL putvara ('O2z', iou, km, (/1, ntrec2/), (/km, 1/), oxy_z(:,l), 1., 0.)

!       POC
        CALL putvars ('POC', iou, ntrec2, pocvol(l), 1., 0.)
        CALL putvara ('POCz', iou, km, (/1, ntrec2/), (/km, 1/), poc_z(:,l), 1., 0.)

!       NO3
        CALL putvars ('NO3', iou, ntrec2, no3vol(l), 1., 0.)
        CALL putvara ('NO3z', iou, km, (/1, ntrec2/), (/km, 1/), no3_z(:,l), 1., 0.)

!       Total C
        CALL putvars ('TC', iou, ntrec2, toc(l), 1., 0.)

!       Total N
        CALL putvars ('TN', iou, ntrec2, ton(l), 1., 0.)

!       Diagnostic variables
        if (exists) then 
!           PH
            CALL putvars ('PH', iou, ntrec2, phvol(l), 1., 0.)
            CALL putvara ('PHz', iou, km, (/1, ntrec2/), (/km, 1/), ph_z(:,l), 1., 0.)

!           PHY
            CALL putvars ('PPPHY', iou, ntrec2, ppphyvol(l), 1., 0.)

!           EPC100
            CALL putvars ('EPC100', iou, ntrec2, epc100glo(l), 1., 0.)

!           EPCAL100
            CALL putvars ('EPCAL100', iou, ntrec2, epcal100glo(l), 1., 0.)

!           Cflux
            CALL putvars ('CFLX', iou, ntrec2, cglo(l), 1., 0.)

!           Oflux
            CALL putvars ('OFLX', iou, ntrec2, ofluxglo(l), 1., 0.)

!           Nfix
            CALL putvars ('NFIX', iou, ntrec2, nfixglo(l), 1., 0.)

        endif
      enddo
      print*, 'closing netcdf'
      !CALL flush(6)

      CALL closefile (iou)
      
END SUBROUTINE calc
!=========================================================
! Area averaging of 3d field over the selected regions 
!=========================================================
SUBROUTINE area_ave (e1,e2,e3, mask,a, imt,jmt,km,a_mean,ss,kk)
      implicit none
      integer imt, jmt, km, i, j, kk
      real e1(imt,jmt),e2(imt,jmt), e3(imt,jmt,km) 
      real a(imt,jmt), mask(imt,jmt) 
      real a_mean, ss, s1, vol

          s1=0.
          ss=0.
          do i=1,imt-2  ! not to double count the cyclic boundary
              do j=1,jmt-1
                  if (mask(i,j).gt.0.5) then  ! mask the region of interst
                      vol = e1(i,j)*e2(i,j)*e3(i,j,kk)
                      ss=ss+vol 
                      s1=s1+a(i,j)*vol  
                  endif
              enddo
          enddo

          a_mean = 0.
          if (ss.ne.0.) then 
              a_mean =s1/ss
          endif

      return
END SUBROUTINE area_ave
!=========================================================
! Area averaging of 2d field over the selected regions 
!=========================================================
SUBROUTINE area_ave_flx(e1,e2, mask, a, imt, jmt,a_mean,ss)
      implicit none
      integer imt, jmt, i, j
      real e1(imt,jmt),e2(imt,jmt) 
      real a(imt,jmt), mask(imt,jmt) 
      real a_mean, ss, s1, arc

          s1=0.
          ss=0.
          do i=1,imt-2  ! not to double count the cyclic boundary
              do j=1,jmt-1
                  if (mask(i,j).gt.0.5) then  ! mask the region of interst
                      arc = e1(i,j)*e2(i,j)
                      ss=ss+arc 
                      s1=s1+a(i,j)*arc  
                  endif
              enddo
          enddo

          a_mean = 0.
          if (ss.ne.0.) then 
              a_mean =s1/ss
          endif

      return
END SUBROUTINE area_ave_flx

SUBROUTINE noleap_days(year, mon, day, days_elapsed)
!    Given a year, mon, day, returns the number of days elapsed
!    since 01-01-0000 (using a noleap/365_day calendar)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: year, mon, day
    INTEGER, INTENT(OUT) :: days_elapsed
    INTEGER              :: doy, dpy, mmon, myear

    ! Adjust if mon > 12. This is a potential here. No correction
    ! for days being off. I'm assuming day will mostly be 1 anyway.

    IF (mon > 12) THEN
        myear = year + mon/12
        mmon = MOD(mon, 12)
    ELSE
       myear = year
       mmon = mon
    ENDIF

    ! In a given (noleap) year, compute the day of year    
    ! http://www.davidgsimpson.com/software/greg2doy_f90.txt
    !
    doy = ((275*mmon)/9) - ((mmon+9)/6) + day - 30

    ! Compute the number of days in the preceeding years
    dpy = 365 * year

    ! Tally for the final result
    days_elapsed = dpy + doy
    return
END SUBROUTINE noleap_days

