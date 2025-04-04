PROGRAM nemo_ocean_diag
IMPLICIT NONE 
      
! ======================================================================
!  Purpose: Run-time diagnostics for NEMO (ORCA2) 
!
! HISTORY:
! -------
! N. Lambert  July   2023   Include the time variation of e3t
!
! O. Riche    Jan    2023   change PH to pH when reading the variable from the o/p file.
!                           read e3t from grid_t instead of from mesh_mask (not present).
!
! O Riche     Jun    2022   Bare-bones version to test global tracer mass cons.
!                           in NEMO4
! N. Swart    Dec    2015   Abstract all calculations to ccc_nemo_rtd_utils
!                           module, which is shared between all rtd.
!
! N. Swart    Nov    2015 1. Remove annual mean calculation and rewrite
!                            all code to operate on monthly data.
!                         2. Remove North Fold point from computions.
!                            (i.e. sum to jmt -1)
!                         3. Rewrite of functions and code to freefrom. 
!                         4. Improve time axis in output NetCDF, include
!                            noleap_days subroutine.
!
!  J. Christian Apr 07 2016 Update to CanOE variable set.
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
          CALL getdimlen ('time_counter', iou, lm) 
      
!         do the calculations and save the output netcdf   
          CALL calc (imt, jmt, km, lm)

END PROGRAM nemo_ocean_diag

SUBROUTINE calc (imt, jmt, km, lm)
!     Does the required calculations and saves the output to netcdf
      USE ccc_nemo_rtd_utils, only: area_ave, area_ave_flx, noleap_days

      IMPLICIT NONE
      integer, parameter:: dp=kind(0.d0) ! double precision



! ======================================================================
!     Input data 
! ======================================================================
      INTEGER imt, jmt, km, lm
      INTEGER i, j, k, l

!     Grid-related arrays
      REAL, DIMENSION(:, :), ALLOCATABLE    :: lon2d, lat2d, e1t, e2t
      REAL, DIMENSION(:, :, :), ALLOCATABLE :: t_mask
      REAL, DIMENSION(:), ALLOCATABLE       :: deptht

!     Monthly DIC, CaCO3, TA, PH, O2
      REAL, DIMENSION(:, :, :, :), ALLOCATABLE :: dic, caco3, tal, oxy

!     Monthly POC, GOC
      REAL, DIMENSION(:, :, :, :), ALLOCATABLE :: e3t, poc, goc

!     Monthly NO3, NH4, dFe
      REAL, DIMENSION(:, :, :, :), ALLOCATABLE :: no3, nh4, dfe

!     Monthly PHY and Zoo
      REAL, DIMENSION(:, :, :, :), ALLOCATABLE :: phy, phy2, phyn, phy2n, zoo, zoo2

!     Monthly primary production: PPPHY, PPPHY2 
      REAL, DIMENSION(:, :, :, :), ALLOCATABLE :: ppphy, ppphy2, nfix, irondep

!     Monthly export fluxes of C (EPC100), CaCO3 (EPCAL100)
      REAL, DIMENSION(:, :, :), ALLOCATABLE :: epc100, epcal100

!     Monthly surface fluxes of DIC, O2
      REAL, DIMENSION(:, :, :), ALLOCATABLE :: cflux, oflux, ph

! ======================================================================
!     Working arrays / variables  
! ======================================================================
      REAL :: dum, dvol, vol 
      REAL :: dicz, caco3z, talz, oxyz, pocz, gocz, no3z, nh4z
      REAL :: dfez, phyz, phy2z, phynz, phy2nz, zooz, zoo2z 
      REAL :: ppphyz, ppphy2z, nfixz, irondepz
      REAL :: test_var

!     total ocean carbon, nitrogen      
      REAL, DIMENSION(:), ALLOCATABLE :: toc, ton
 
! g_mask
      REAL, DIMENSION(:, :), ALLOCATABLE :: g_mask

! ======================================================================
!     Output data 
! ======================================================================
! (1) Global-mean profiles for 3D data:

      REAL, DIMENSION(:, :), ALLOCATABLE :: dic_z, caco3_z, tal_z, oxy_z, poc_z
      REAL, DIMENSION(:, :), ALLOCATABLE :: goc_z, no3_z, nh4_z, dfe_z, phy_z, phy2_z, phyn_z, phy2n_z
      REAL, DIMENSION(:, :), ALLOCATABLE :: zoo_z, zoo2_z, ppphy_z, ppphy2_z, nfix_z, irondep_z

! (2) Global-mean (volume weighted) or integral

!     DIC, CaCO3 TA, PH, O2
      REAL, DIMENSION(:), ALLOCATABLE :: dicvol, caco3vol, talvol, oxyvol, pocvol, gocvol
      REAL, DIMENSION(:), ALLOCATABLE :: no3vol, nh4vol, dfevol, phyvol, phy2vol, phynvol, phy2nvol
      REAL, DIMENSION(:), ALLOCATABLE :: zoovol, zoo2vol, ppphyvol, ppphy2vol, nfixvol, irondepvol
      REAL, DIMENSION(:), ALLOCATABLE :: epc100glo, epcal100glo, cglo, ofluxglo, phglo 

!----------------
!  NetCDF-output specific
      integer id_time, id_z, iou, ntrec, ntrec2, iyear, imon, days_elapsed
      logical exists, exists1, notopen
      real tyear, tdays_elapsed
      CHARACTER(len=32) :: year_arg_in, mon_arg_in
      integer, dimension(7) :: ierr


!----------------
!     input file stuff
      character fname05*100, fname06*100, fname07*100, fname08*100
      integer year, iou4, iou5, iou6, iou7, recn, nrecon

         year =0
         ntrec=0
         iou4 =0
         iou5 =0 
         iou6 =0 
         iou7 =0 
         recn =12.
         nrecon = int(recn + 0.001)
!---------------------------------------------------
!    Define NetCDF files   
!---------------------------------------------------
        fname05='orca_mesh_mask'
        fname06='ptrc_t'   
        fname07='diad_t'   
        fname08='grid_t' 
!---------------------------------------------------
!    Get grid/mask data   
!---------------------------------------------------
      ALLOCATE( lon2d(imt,jmt), lat2d(imt,jmt), e1t(imt,jmt), e2t(imt,jmt),      &
         &      g_mask(imt,jmt), STAT=ierr(1) )
      ALLOCATE( e3t(imt,jmt,km,lm), t_mask(imt,jmt,km), STAT=ierr(2) )
      ALLOCATE( deptht(km), STAT=ierr(3) )
      ALLOCATE( dic(imt,jmt,km,lm), caco3(imt,jmt,km,lm),tal (imt,jmt,km,lm),    &
         &      oxy(imt,jmt,km,lm), poc(imt,jmt,km,lm),                          &
         &      goc(imt,jmt,km,lm), no3(imt,jmt,km,lm), nh4(imt,jmt,km,lm),      &
         &      dfe(imt,jmt,km,lm), phy(imt,jmt,km,lm), phy2(imt,jmt,km,lm),     &
         &      phyn(imt,jmt,km,lm), phy2n(imt,jmt,km,lm), zoo(imt,jmt,km,lm),   &
         &      zoo2(imt,jmt,km,lm), ppphy(imt,jmt,km,lm), ppphy2(imt,jmt,km,lm),&
         &      nfix(imt,jmt,km,lm), irondep(imt,jmt,km,lm), STAT=ierr(4) )
      ALLOCATE( epc100(imt,jmt,lm), epcal100(imt,jmt,lm), cflux(imt,jmt,lm),     &
         &      oflux(imt,jmt,lm), ph(imt,jmt,lm),    STAT=ierr(5) ) 
      ALLOCATE( dic_z(km, lm), caco3_z(km, lm), tal_z(km, lm),                   &
         &      oxy_z(km, lm), poc_z(km, lm), goc_z(km, lm), no3_z(km, lm),      &
         &      nh4_z(km, lm), dfe_z(km, lm), phy_z(km, lm), phy2_z(km, lm),     &
         &      phyn_z(km, lm), phy2n_z(km, lm), zoo_z(km, lm), zoo2_z(km, lm),  &
         &      ppphy_z(km, lm), ppphy2_z(km, lm), nfix_z(km, lm),               &
         &      irondep_z(km, lm), STAT=ierr(6) )
      ALLOCATE( toc(lm), ton(lm), dicvol(lm), caco3vol(lm), talvol(lm),          &
         &      phglo(lm), oxyvol(lm), pocvol(lm), gocvol(lm), no3vol(lm),       &
         &      nh4vol(lm), dfevol(lm), phyvol(lm), phy2vol(lm), phynvol(lm),    &
         &      phy2nvol(lm), zoovol(lm), zoo2vol(lm), ppphyvol(lm),             &
         &      ppphy2vol(lm), nfixvol(lm), irondepvol(lm), epc100glo(lm),       &
         &      epcal100glo(lm), cglo(lm), ofluxglo(lm), STAT=ierr(7) )

         IF (MAXVAL(ierr) /=0) THEN
           STOP 'Memory allocation error in carbon RTD'
         ENDIF

      print*,'Reading data on NEMO grid...'
      print*,''
      print*,'Grid size: imt, jmt, km:', imt, jmt, km

      CALL openfile (fname05,iou4)
      CALL getvara ('e1t', iou4, imt*jmt, (/1,1,1/), (/imt,jmt,1/),e1t , 1., 0.)
      CALL getvara ('e2t', iou4, imt*jmt, (/1,1,1/),  (/imt,jmt,1/),e2t , 1., 0.)
      CALL getvara ('tmask', iou4, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/), t_mask , 1., 0.)
      CALL closefile (iou4)

!     get some more grid information
      CALL openfile (fname06,iou5)
      CALL getvara ('nav_lon', iou5, imt*jmt, (/1,1,1/), (/imt,jmt,1/), lon2d, 1., 0.)
      CALL getvara ('nav_lat', iou5, imt*jmt, (/1,1,1/), (/imt,jmt,1/), lat2d, 1., 0.)
      CALL getvara ('deptht', iou5, km, (/1/), (/km/), deptht, 1., 0.)
      CALL closefile (iou5)
      CALL openfile(fname08,iou7) 
      CALL getvara ('thkcello', iou7, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/),e3t , 1., 0.)
      CALL closefile(iou7)  

!---------------------------------------------------
! Read from NetCDF
!---------------------------------------------------

!     read prognostic variables from the prtc_t file        
      CALL openfile (fname06,iou5)

!     read diagnostic variables from the diad_t file, if it exists
      inquire (file=trim(fname07), exist=exists)
! O Riche June 6th 2022
! no diad_t here for this test and for now
!	  exists=.false.
      if (exists) then         
          CALL openfile (fname07,iou6)
      endif

!        DIC, TA, O2 
         CALL getvara('DIC',      iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/),   dic, 1., 0.)                                 
         CALL getvara('CaCO3',    iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), caco3, 1., 0.)  
         CALL getvara('Alkalini',     iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/),   tal, 1., 0.)   

!        PH moved below for reading with other diat_t input

         CALL getvara('O2', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), oxy, 1., 0.)   

!        POC, GOC
         CALL getvara('POC', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), poc, 1., 0.)   
         CALL getvara('GOC', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), goc, 1., 0.)   

!        NO3, NH4, dFe
         CALL getvara('NO3', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), no3, 1., 0.)
         CALL getvara('NH4', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), nh4, 1., 0.)    
         CALL getvara('Fer', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), dfe, 1., 0.)    

!        PHY, PHY2, ZOO, ZOO2
         CALL getvara('PHYC',  iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/),  phy, 1., 0.)  
         CALL getvara('PHY2C', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), phy2, 1., 0.)    
         CALL getvara('PHYN',  iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), phyn, 1., 0.)  
         CALL getvara('PHY2N', iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/),phy2n, 1., 0.)    
         CALL getvara('ZOO' ,  iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/),  zoo, 1., 0.)    
         CALL getvara('ZOO2',  iou5, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), zoo2, 1., 0.)    
      
!        Diagnostic variables
         if (exists) then 
!            3-D: PPPHY, PPPHY2, EPC100, Nfix, Irondep
             CALL getvara('PPPHY',    iou6, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/),   ppphy, 1., 0.)   
             CALL getvara('PPPHY2',   iou6, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/),  ppphy2, 1., 0.)   
             CALL getvara('Nfix',     iou6, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/),    nfix, 1., 0.)   
             CALL getvara('Irondep',  iou6, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), irondep, 1., 0.)   
!            2-D: PH, EPCAL100, DIC flux, Oflux
             CALL getvara('EPC100',   iou6, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/),   epc100, 1., 0.)    
             CALL getvara('EPCALC100',iou6, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/), epcal100, 1., 0.)    
             CALL getvara('Cflx',     iou6, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/),    cflux, 1., 0.)    
             CALL getvara('Oflx',     iou6, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/),    oflux, 1., 0.)    
             CALL getvara('pH',       iou6, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/),       ph, 1., 0.)   
         endif 
      CALL closeall ! close all open netcdf files

! ********** Do some basic calculations ************

!---------------------------------------------------
! (1) Global annual mean T(z) and S(z)  
!---------------------------------------------------

! DIC, TA, PH, O2 
      dicvol(:)   = 0.0_dp
      caco3vol(:) = 0.0_dp
      talvol(:)   = 0.0_dp
      phglo(:)    = 0.0_dp
      oxyvol(:)   = 0.0_dp
! POC, GOC
      pocvol(:)   = 0.0_dp
      gocvol(:)   = 0.0_dp
! NO3, NH4
      no3vol(:)   = 0.0_dp
      nh4vol(:)   = 0.0_dp  
      dfevol(:)   = 0.0_dp  
! PHY, PHY2, ZOO, ZOO2
      phyvol(:)   = 0.0_dp 
      phy2vol(:)  = 0.0_dp 
      phynvol(:)  = 0.0_dp 
      phy2nvol(:) = 0.0_dp 
      zoovol(:)   = 0.0_dp 
      zoo2vol(:)  = 0.0_dp 
! PPPHY, PPPHY2      
      ppphyvol(:) = 0.0_dp
      ppphy2vol(:) = 0.0_dp
      nfixvol(:) = 0.0_dp
      irondepvol(:) = 0.0_dp

      do l = 1, lm                             
          vol =0.0_dp
          do k = 1, km   
             g_mask(:, :)  = t_mask(:, :, k) 

    !        DIC, TA, PH, O2 
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, dic(:, :, k, l),   imt, jmt, km, dicz,   dvol, k)  
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, caco3(:, :, k, l), imt, jmt, km, caco3z, dvol, k)  
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, tal(:, :, k, l),     imt, jmt, km, talz,   dvol, k)  
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, oxy(:, :, k, l),     imt, jmt, km, oxyz,   dvol, k)

    !        POC, GOC
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, poc(:, :, k, l), imt, jmt, km, pocz, dvol, k)
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, goc(:, :, k, l), imt, jmt, km, gocz, dvol, k)

    !        NO3, NH4
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, no3(:, :, k, l), imt, jmt, km, no3z, dvol, k)
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, nh4(:, :, k, l), imt, jmt, km, nh4z, dvol, k)
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, dfe(:, :, k, l), imt, jmt, km, dfez, dvol, k)  

    !        PHY, PHY2, ZOO, ZOO2
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, phy(:, :, k, l),   imt, jmt, km, phyz,   dvol, k)  
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, phy2(:, :, k, l),  imt, jmt, km, phy2z,  dvol, k)  
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, phyn(:, :, k, l),  imt, jmt, km, phynz,  dvol, k)  
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, phy2n(:, :, k, l), imt, jmt, km, phy2nz, dvol, k)  
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, zoo(:, :, k, l),   imt, jmt, km, zooz,   dvol, k)  
             CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, zoo2(:, :, k, l),  imt, jmt, km, zoo2z,  dvol, k) 

             if (exists) then
    !            PPPHY, PPPHY2      
                 !CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, ph(:, :, k, l),      imt, jmt, km, phz,      dvol, k)  
                 CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, ppphy(:, :, k, l),   imt, jmt, km, ppphyz,   dvol, k) 
                 CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, ppphy2(:, :, k, l),  imt, jmt, km, ppphy2z,  dvol, k) 
                 CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, nfix(:, :, k, l),    imt, jmt, km, nfixz,    dvol, k) 
                 CALL area_ave(e1t, e2t, e3t(:,:,:,l), g_mask, irondep(:, :, k, l), imt, jmt, km, irondepz, dvol, k) 
             endif 
    !================================================================
    !        Assign outputs
    !================================================================
    !        DIC, TA, PH, O2 
             dic_z(k, l)   = dicz     
             caco3_z(k, l) = caco3z     
             tal_z(k, l)   = talz     
             oxy_z(k, l)   = oxyz    

    !        POC, GOC
             poc_z(k, l)    = pocz       
             goc_z(k, l)    = gocz       

    !        NO3, NH4
             no3_z(k, l)   = no3z  
             nh4_z(k, l)   = nh4z   
             dfe_z(k, l)   = dfez  

    !        PHY, PHY2, ZOO, ZOO2
             phy_z(k, l)   = phyz
             phy2_z(k, l)  = phy2z
             phyn_z(k, l)   = phynz
             phy2n_z(k, l)  = phy2nz
             zoo_z(k, l)   = zooz
             zoo2_z(k, l)  = zoo2z
         
             if (exists) then
                 !ph_z(k, l)    = phz  
                 ppphy_z(k, l) = ppphyz
                 ppphy2_z(k, l)= ppphy2z
                 nfix_z(k, l)= nfixz
                 irondep_z(k, l)= irondepz
             endif 

    !        DIC, TA, PH, O2 
             dicvol(l)   = dicvol(l)  + dicz*dvol  
             caco3vol(l)  = caco3vol(l)  + caco3z*dvol  
             talvol(l)  = talvol(l)  + talz*dvol  
             oxyvol(l)  = oxyvol(l)  + oxyz*dvol  

    !        POC, GOC
             pocvol(l)  = pocvol(l)  + pocz*dvol  
             gocvol(l)  = gocvol(l)  + gocz*dvol  

    !        NO3, NH4, dFe
             no3vol(l)  = no3vol(l)  + no3z*dvol  
             nh4vol(l)  = nh4vol(l)  + nh4z*dvol  
             dfevol(l)  = dfevol(l)  + dfez*dvol  

    !        PHY, PHY2, ZOO, ZOO2
             phyvol(l)     = phyvol(l)   + phyz*dvol
             phy2vol(l)    = phy2vol(l)  + phy2z*dvol 
             phynvol(l)    = phynvol(l)  + phynz*dvol
             phy2nvol(l)   = phy2nvol(l) + phy2nz*dvol 
             zoovol(l)     = zoovol(l)   + zooz*dvol
             zoo2vol(l)    = zoo2vol(l)  + zoo2z*dvol 

             if (exists) then 
                 ppphyvol(l)   =  ppphyvol(l)  + ppphyz*dvol  
                 ppphy2vol(l)  = ppphy2vol(l)  + ppphy2z*dvol  
                 nfixvol(l)  = nfixvol(l)  + nfixz*dvol  
                 irondepvol(l)  = irondepvol(l)  + irondepz*dvol  
             endif
                       
              vol = vol + dvol
          enddo  ! depth, k        

    !     compute toc and ton
          toc(l)  = dicvol(l)  + caco3vol(l)  + pocvol(l)  + gocvol(l)                          &
         &      + phyvol(l)  + phy2vol(l)  + zoovol(l)  + zoo2vol(l) 
    !     convert from mmol C to Pg C      
          toc(l)  = toc(l)  * 12.0e-18
          ton(l)  = no3vol(l) + nh4vol(l) + phynvol(l) + phy2nvol(l)  +                                         &
         &      (zoovol(l)  + zoo2vol(l)  + pocvol(l)  + gocvol(l)) * 16./106.
    !     convert to Pg      
          ton(l)  = ton(l)  * 14.007e-18

          if (vol.ne.0.) then         
    !        DIC, TA, PH, O2 
             dicvol(l)  = dicvol(l) /vol 
             caco3vol(l)  = caco3vol(l) /vol 
             talvol(l)  = talvol(l) /vol 
             oxyvol(l)  = oxyvol(l) /vol  
    !        POC, GOC  
             pocvol(l)  = pocvol(l) /vol  
             gocvol(l)  = gocvol(l) /vol  
    !        NO3, NH4
             no3vol(l)  = no3vol(l) /vol  
             nh4vol(l)  = nh4vol(l) /vol 
             dfevol(l)  = dfevol(l) /vol

             if (exists) then   
    !            convert to PgC/y
                 ppphyvol(l)   = ppphyvol(l)  * 12.e-15 * 86400. * 365.
                 ppphy2vol(l)  = ppphy2vol(l)  * 12.e-15 * 86400. * 365.
    !            convert to TgN/y
                 nfixvol(l)       = nfixvol(l)   * 14.007e-12 * 86400. * 365.
    !            convert to mol Fe/y 
                 irondepvol(l)    = irondepvol(l)   * 86400. * 365.

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
             CALL area_ave_flx (e1t, e2t, g_mask, epc100(:,:,l),   imt, jmt,   epc100glo(l), dum) 
             CALL area_ave_flx (e1t, e2t, g_mask, epcal100(:,:,l), imt, jmt, epcal100glo(l), dum) 

    !        Cflux, Oflux, Nfix, Irondep  
             CALL area_ave_flx (e1t, e2t, g_mask, cflux(:,:,l),   imt, jmt, cglo(l),       dum) 
             CALL area_ave_flx (e1t, e2t, g_mask, oflux(:,:,l),   imt, jmt, ofluxglo(l),   dum) 
             CALL area_ave_flx (e1t, e2t, g_mask, ph(:,:,l),   imt, jmt, phglo(l),   dum) 

    !       <PISCES OR 01/15/2014> convert into PgC/y 
             cglo(l)          = cglo(l)         * dum * 12.e-15 * 86400. * 365.
             epc100glo(l)     = epc100glo(l)    * dum * 12.e-15 * 86400. * 365.
             epcal100glo(l)   = epcal100glo(l)  * dum * 12.e-15 * 86400. * 365.

    !        convert to mol O2/yr
             ofluxglo(l)      = ofluxglo(l)  * dum *  86400. * 365.

          endif
      enddo !main time loop

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
          CALL defvar ('CaCO3', iou, 1, (/id_time/), -1.e4                       &  
     &        , 1.e4,' ', 'F', 'Global mean calcite Concentration'               & 
     &        , 'CaCO3', 'mmol m^-3')

          CALL defvar ('CaCO3z', iou, 2, (/id_z, id_time/), -1.e4                &
     &        , 1.e4,' ', 'F', 'Calcite concentration by level'                  &
     &        , 'CaCO3z', 'mmol m^-3')

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

!         GOC
          CALL defvar ('GOC', iou, 1, (/id_time/), -1.e4                         &
     &        , 1.e4,' ', 'F'                                                    &
     &        , 'Global mean large Particulate Organic Carbon'                   &
     &        , 'GOC', 'mmol m^-3')

          CALL defvar ('GOCz', iou, 2, (/id_z, id_time/), -1.e4                  & 
     &        , 1.e4,' ', 'F', 'Large POC by level'                              &
     &        , 'GOCz', 'mmol m^-3')

!         NO3
          CALL defvar ('NO3', iou, 1, (/id_time/), -1.e4                         &
     &        , 1.e4,' ', 'F', 'Global mean NO3'                                 &  
     &        , 'NO3', 'uM')

          CALL defvar ('NO3z', iou, 2, (/id_z, id_time/), -1.e4                  &
     &        , 1.e4,' ', 'F', 'NO3 by level'                                    & 
     &        , 'NO3z', 'uM')

!         NH4
          CALL defvar ('NH4', iou, 1, (/id_time/), -1.e4                         &
     &        , 1.e4,' ', 'F', 'Global mean NH4'                                 &  
     &        , 'NH4', 'uM')

          CALL defvar ('NH4z', iou, 2, (/id_z, id_time/), -1.e4                  &
     &        , 1.e4,' ', 'F', 'NH4 by level'                                    &
     &        , 'NH4z', 'uM')

!         dFe
          CALL defvar ('dFe', iou, 1, (/id_time/), -1.e4                         &
     &        , 1.e4,' ', 'F', 'Global mean dFe'                                 & 
     &        , 'dFe', 'uM')

          CALL defvar ('dFez', iou, 2, (/id_z, id_time/), -1.e4                  &
     &        , 1.e4,' ', 'F', 'dFe by level'                                    &
     &        , 'dFez', 'uM')

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

!              CALL defvar ('PHz', iou, 2, (/id_z, id_time/), -1.e4               &
!     &            , 1.e4,' ', 'F', 'pH by level'                                 & 
!     &            , 'PHz', '')

!             PPPHY
              CALL defvar ('PPPHY', iou, 1, (/id_time/), -1.e4                   &
     &            , 1.e4,' ', 'F', 'Primary production of nanophyto'             & 
     &            , 'PPPHY', 'PgC/yr')


!             PPPHY2
              CALL defvar ('PPPHY2', iou, 1, (/id_time/), -1.e4                  &  
     &            , 1.e4,' ', 'F', 'Primary production of Diatoms'               & 
     &            , 'PPPHY2', 'PgC/yr')

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
     &            , 1.e4,' ', 'F', 'Nitrogen fixation'                           &
     &            , 'NFIX', 'TgN/yr')
!             Irondep
              CALL defvar ('Irondep', iou, 1, (/id_time/), -1.e4                 & 
     &           , 1.e4,' ', 'F', 'Iron deposition'                              & 
     &           , 'Irondep', 'mol/yr')
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

!       time
        ntrec2 = ntrec + l - 1
        call putvars ('time', iou, ntrec2, tdays_elapsed, 1., 0.)

!       DIC
        CALL putvars ('DIC', iou, ntrec2, dicvol(l), 1., 0.)
        CALL putvara ('DICz', iou, km, (/1, ntrec2/), (/km, 1/), dic_z(:, l), 1., 0.) 

!       CaCO3
        CALL putvars ('CaCO3', iou, ntrec2, caco3vol(l), 1., 0.)
        CALL putvara ('CaCO3z', iou, km, (/1, ntrec2/), (/km, 1/), caco3_z(:, l), 1., 0.) 

!       ALK
        CALL putvars ('TAL', iou, ntrec2, talvol(l), 1., 0.)
        CALL putvara ('TALz', iou, km, (/1, ntrec2/), (/km, 1/), tal_z(:, l), 1., 0.)

!       OXY
        CALL putvars ('O2', iou, ntrec2, oxyvol(l), 1., 0.)
        CALL putvara ('O2z', iou, km, (/1, ntrec2/), (/km, 1/), oxy_z(:, l), 1., 0.)

!       POC
        CALL putvars ('POC', iou, ntrec2, pocvol(l), 1., 0.)
        CALL putvara ('POCz', iou, km, (/1, ntrec2/), (/km, 1/), poc_z(:, l), 1., 0.)

!       GOC
        CALL putvars ('GOC', iou, ntrec2, gocvol(l), 1., 0.)
        CALL putvara ('GOCz', iou, km, (/1, ntrec2/), (/km, 1/), goc_z(:, l), 1., 0.)

!       NO3
        CALL putvars ('NO3', iou, ntrec2, no3vol(l), 1., 0.)
        CALL putvara ('NO3z', iou, km, (/1, ntrec2/), (/km, 1/), no3_z(:, l), 1., 0.)

!       NH4
        CALL putvars ('NH4', iou, ntrec2, nh4vol(l), 1., 0.)
        CALL putvara ('NH4z', iou, km, (/1, ntrec2/), (/km, 1/), nh4_z(:, l), 1., 0.)

!       dFe
        CALL putvars ('dFe', iou, ntrec2, dfevol(l), 1., 0.)
        CALL putvara ('dFez', iou, km, (/1, ntrec2/), (/km, 1/), dfe_z(:, l), 1., 0.)

!       Total C
        CALL putvars ('TC', iou, ntrec2, toc(l), 1., 0.)

!       Total N
        CALL putvars ('TN', iou, ntrec2, ton(l), 1., 0.)

!       Diagnostic variables
        if (exists) then 
!           PH
            CALL putvars ('PH', iou, ntrec2, phglo(l), 1., 0.)
            !CALL putvara ('PHz', iou, km, (/1, ntrec2/), (/km, 1/), ph_z(:, l), 1., 0.)

!           PHY
            CALL putvars ('PPPHY', iou, ntrec2, ppphyvol(l), 1., 0.)

!           PHY2
            CALL putvars ('PPPHY2', iou, ntrec2, ppphy2vol(l), 1., 0.)

!           EPC100
            CALL putvars ('EPC100', iou, ntrec2, epc100glo(l), 1., 0.)

!           EPCAL100
            CALL putvars ('EPCAL100', iou, ntrec2, epcal100glo(l), 1., 0.)

!           Cflux
            CALL putvars ('CFLX', iou, ntrec2, cglo(l), 1., 0.)

!           Oflux
            CALL putvars ('OFLX', iou, ntrec2, ofluxglo(l), 1., 0.)

!           Nfix
            CALL putvars ('NFIX', iou, ntrec2, nfixvol(l), 1., 0.)

!           Irondep
            CALL putvars ('Irondep', iou, ntrec2, irondepvol(l), 1., 0.)
        endif
      enddo ! time loop

      print*, 'closing netcdf'
      !CALL flush(6)

      CALL closefile (iou)
      
END SUBROUTINE calc

