PROGRAM nemo_ocean_diag
IMPLICIT NONE 
      
! ======================================================================
!  Purpose: Run-time diagnostics for NEMO (ORCA2) 
!  O. Riche   Dec 07  2015 Fix ton and toc; assume all tracers are in nitrogen currency
!                          and convert accordingly with classic Redfield ratios (ll 448 and 452)
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
! ======================================================================
! to compile loCALLy:
! to compile:
!   $F77 -o nemo_ocean_diag nemo_ocean_diag.f uvic_netcdf.f $LINK -L/home/rls/wrk/AR5/CMOR/lib -lnetcdf -I/home/rls/wrk/AR5/CMOR/include/
!   
! UPDATE - 2013Feb20 - MB
!   gfortran -o nemo_ocean_diag nemo_ocean_diag.f uvic_netcdf.f -I/usr/local/netcdf-4.2.1_GF/include -L/usr/local/netcdf-4.2.1_GF/lib/ -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
! OR
!   pgf90 -m64 -O -D_LARGE_FILES -DpgiFortran -o nemo_ocean_diag nemo_ocean_diag.f uvic_netcdf.f -I/usr/local/netcdf-4.2.1_PG/include -L/usr/local/netcdf-4.2.1_PG/lib/ -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz

! ======================================================================
      INTEGER :: imt, jmt, km, ll, iou

!         establish the size of the grid from the input file.
          CALL openfile ("ptrc_t", iou)
          CALL getdimlen ('x', iou, imt)
          CALL getdimlen ('y', iou, jmt)
          CALL getdimlen ('deptht', iou, km)
      
          ll = 1

!         do the calculations and save the output netcdf   
          CALL calc (imt, jmt, km, ll)

END PROGRAM nemo_ocean_diag

SUBROUTINE calc (imt, jmt, km, ll)
!     Does the required calculations and saves the output to netcdf

      IMPLICIT NONE

! ======================================================================
!     Input data 
! ======================================================================
      INTEGER imt, jmt, km, ll
      INTEGER i, j, k, mon

!     Grid-related arrays
      REAL, DIMENSION(imt, jmt)     :: lon2d, lat2d, e1t, e2t
      REAL, DIMENSION(km)           :: deptht
      REAL, DIMENSION(imt, jmt, km) :: e3t(imt,jmt,km) 

!     Monthly DIC, CaCO3, TA, PH, O2
      REAL, DIMENSION(imt, jmt, km) :: dic, caco3, tal, ph, oxy

!     Monthly POC, GOC, DOC
      REAL, DIMENSION(imt, jmt, km) :: poc, goc, doc

!     Monthly NO3, NH4, PO4, Si
      REAL, DIMENSION(imt, jmt, km) :: no3, nh4, po4, si

!     Monthly PHY and Zoo
      REAL, DIMENSION(imt, jmt, km) :: phy, phy2, zoo, zoo2

!     Monthly primary production: PPPHY, PPPHY2 
      REAL, DIMENSION(imt, jmt, km) :: ppphy, ppphy2

!     Monthly export fluxes of C (EPC100), CaCO3 (EPCAL100)
      REAL, DIMENSION(imt, jmt)     :: epc100, epcal100

!     Monthly surface fluxes of DIC, O2, N2, Fe
      REAL, DIMENSION(imt, jmt)     :: cflux, oflux, nfix, irondep

! ======================================================================
!     Pre-computed data  
! ======================================================================

!     3D variables
      REAL, DIMENSION(imt, jmt, km) :: t_mask, dic_ann, caco3_ann
      REAL, DIMENSION(imt, jmt, km) :: tal_ann, ph_ann, oxy_ann, poc_ann, goc_ann, doc_ann, no3_ann
      REAL, DIMENSION(imt, jmt, km) :: nh4_ann, po4_ann, si_ann, phy_ann, phy2_ann, zoo_ann
      REAL, DIMENSION(imt, jmt, km) :: zoo2_ann, ppphy_ann, ppphy2_ann  

!    2D variables
      REAL, DIMENSION(imt, jmt) ::  epc100_ann, epcal100_ann
     REAL, DIMENSION(imt, jmt) ::   cflux_ann, oflux_ann, nfix_ann, irondep_ann

!     weight for annual mean calculations, DY, 2/OCT/2013
      INTEGER DIM(12)
      DATA DIM/31,28,31,30,31,30,31,31,30,31,30,31/
! ======================================================================
!     Working arrays / variables  
! ======================================================================
      REAL :: dum, dvol, vol 
      REAL :: dicz, caco3z, talz, phz, oxyz, pocz, gocz, docz,no3z, nh4z
      REAL :: po4z, siz, phyz, phy2z, zooz, zoo2z, ppphyz, ppphy2z

!     total ocean carbon, nitrogen      
      REAL, DIMENSION(imt, jmt, km) :: toc, ton
 
! g_mask
      REAL, DIMENSION(imt, jmt) :: g_mask

! ======================================================================
!     Output data 
! ======================================================================
! (1) Global-mean profiles for 3D data:

      REAL, DIMENSION(km) :: dic_z, caco3_z, tal_z, ph_z, oxy_z, poc_z
      REAL, DIMENSION(km) :: goc_z, doc_z, no3_z, nh4_z, po4_z, si_z, phy_z, phy2_z
      REAL, DIMENSION(km) :: zoo_z, zoo2_z, ppphy_z, ppphy2_z

! (2) Global-mean (volume weighted) or integral

!     DIC, CaCO3 TA, PH, O2
      REAL :: dicvol, caco3vol, talvol, phvol, oxyvol, pocvol, gocvol
      REAL :: docvol, no3vol, nh4vol, po4vol, sivol, phyvol, phy2vol
      REAL :: zoovol, zoo2vol, ppphyvol, ppphy2vol

      REAL epc100glo, epcal100glo, cglo, ofluxglo, nfixglo, irondepglo

!----------------
!  NetCDF-output specific
      integer id_time, id_z, iou, ntrec, iyear
      logical exists, exists1, notopen
      real tyear
      CHARACTER(len=32) :: year_arg_in

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
      CALL getvara ('tmask', iou4, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),t_mask , 1., 0.)
      CALL closefile (iou4)

!     get some more grid information
      CALL openfile (fname06,iou5)
      CALL getvara ('nav_lon', iou5, imt*jmt, (/1,1,1/), (/imt,jmt,1/), lon2d, 1., 0.)
      CALL getvara ('nav_lat', iou5, imt*jmt, (/1,1,1/), (/imt,jmt,1/), lat2d, 1., 0.)
      CALL getvara ('deptht', iou5, km, (/1/), (/km/), deptht, 1., 0.)
      CALL closefile (iou5)

!---------------------------------------------------
!    Set to zero arrays for annual accumulation   
!---------------------------------------------------

!     3D fields
!     DIC, TA, PH, O2
      dic_ann(:, :, :)       = 0.
      caco3_ann(:, :, :)     = 0.
      tal_ann(:, :, :)       = 0.
      ph_ann(:, :, :)        = 0.    
      oxy_ann(:, :, :)       = 0.
!     POC, GOC, DOC
      poc_ann(:, :, :)       = 0.
      goc_ann(:, :, :)       = 0.
      doc_ann(:, :, :)       = 0.
!     NO3, NH4, PO4, Si 
      no3_ann(:, :, :)       = 0. 
      nh4_ann(:, :, :)       = 0.
      po4_ann(:, :, :)       = 0. 
      si_ann(:, :, :)        = 0.
!     PHY, PHY2, ZOO, ZOO2
      phy_ann(:, :, :)       = 0.    
      phy2_ann(:, :, :)      = 0
      zoo_ann(:, :, :)       = 0.
      zoo2_ann(:, :, :)      = 0.
!     PPPHY, PPPHY2
      ppphy_ann(:, :, :)     = 0.    
      ppphy2_ann(:, :, :)    = 0.    

!     2D fields
!     EPC100, EPCAL100
      epc100_ann(:, :)       = 0.    
      epcal100_ann(:, :)     = 0.    
!     DIC flux 
      cflux_ann(:, :)        = 0.         
      oflux_ann(:, :)        = 0.         
      nfix_ann(:, :)         = 0.         
      irondep_ann(:, :)      = 0.         

!---------------------------------------------------
! reading time-dependent (monthly) data and compute 
! annual-mean (the latter does not take into consideration 
! the different number of days in each month; has to be 
! corrected at some point...)
!---------------------------------------------------

!     read prognostic variables from the prtc_t file        
      CALL openfile (fname06,iou5)

!     read diagnostic variables from the diad_t file, if it exists
      inquire (file=trim(fname07), exist=exists)                          
      if (exists) then         
          CALL openfile (fname07,iou6)
      endif

      do mon=1, nrecon !  12         
!        DIC, TA, O2 
         CALL getvara ('DIC',      iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/),   dic, 1., 0.)                                 
!         CALL getvara ('CaCO3',    iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), caco3, 1., 0.)  
         CALL getvara ('Alkalini', iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/),   tal, 1., 0.)   

!        PH moved below for reading with other diat_t input

         CALL getvara ('O2', iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), oxy, 1., 0.)   

!        POC, GOC, DOC
         CALL getvara ('POC', iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), poc, 1., 0.)   
!         CALL getvara ('GOC', iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), goc, 1., 0.)   
!         CALL getvara ('DOC', iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), doc, 1., 0.)   

!        NO3, NH4, PO4, Si 
         CALL getvara ('NO3', iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), no3, 1., 0.)
!         CALL getvara ('NH4', iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), nh4, 1., 0.)    
!         CALL getvara ('PO4', iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), po4, 1., 0.)    
!         CALL getvara ('Si' , iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/),  si, 1., 0.)      

!        PHY, PHY2, ZOO, ZOO2
         CALL getvara ('PHY',  iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/),  phy, 1., 0.)  
!         CALL getvara ('PHY2', iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), phy2, 1., 0.)    
         CALL getvara ('ZOO' , iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/),  zoo, 1., 0.)    
!         CALL getvara ('ZOO2', iou5, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), zoo2, 1., 0.)    
      
!        Diagnostic variables
         if (exists) then 
!            3-D: PH, PPPHY, PPPHY2, EPC100,
             CALL getvara ('PH',       iou6, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/),     ph, 1., 0.)   
             CALL getvara ('PPPHY',    iou6, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/),  ppphy, 1., 0.)   
!             CALL getvara ('PPPHY2',   iou6, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), ppphy2, 1., 0.)   

!            2-D :  EPCAL100, DIC flux, Oflux, Nfix, Irondep
             CALL getvara ('EPC100',   iou6, imt*jmt, (/1,1,mon/), (/imt,jmt,ll/),   epc100, 1., 0.)    
!             CALL getvara ('EPCAL100', iou6, imt*jmt, (/1,1,mon/), (/imt,jmt,ll/), epcal100, 1., 0.)    
             CALL getvara ('Cflx',     iou6, imt*jmt, (/1,1,mon/), (/imt,jmt,ll/),    cflux, 1., 0.)    
             CALL getvara ('Oflx',     iou6, imt*jmt, (/1,1,mon/), (/imt,jmt,ll/),    oflux, 1., 0.)    
             CALL getvara ('Nfix',     iou6, imt*jmt, (/1,1,mon/), (/imt,jmt,ll/),     nfix, 1., 0.)    
!             CALL getvara ('Irondep',  iou6, imt*jmt, (/1,1,mon/), (/imt,jmt,ll/),  irondep, 1., 0.)    
         endif 

!        calculate the annual mean from the monthly data.
!         3D data
!         DIC, TA, PH, O2 
          dic_ann   = dic_ann  + dic*DIM(mon)/365.
!          caco3_ann = caco3_ann+ caco3*DIM(mon)/365.
          tal_ann   = tal_ann  + tal*DIM(mon)/365. 
          oxy_ann   = oxy_ann  + oxy*DIM(mon)/365.

!         POC, GOC, DOC
          poc_ann   = poc_ann  + poc*DIM(mon)/365.
!          goc_ann   = goc_ann  + goc*DIM(mon)/365.
!          doc_ann   = doc_ann  + doc*DIM(mon)/365.

!         NO3, NH4, PO4, Si 
          no3_ann   = no3_ann + no3*DIM(mon)/365. 
!          nh4_ann   = nh4_ann + nh4*DIM(mon)/365. 
!          po4_ann   = po4_ann + po4*DIM(mon)/365.  
!          si_ann    =  si_ann + si*DIM(mon)/365. 

!         PHY, PHY2, ZOO, ZOO2
          phy_ann   = phy_ann  + phy*DIM(mon)/365.  
!          phy2_ann  = phy2_ann + phy2*DIM(mon)/365. 
          zoo_ann   = zoo_ann  + zoo*DIM(mon)/365. 
!          zoo2_ann  = zoo2_ann + zoo2*DIM(mon)/365. 

!         Diagnostic variables
          if (exists) then   
              ph_ann    = ph_ann   + ph*DIM(mon)/365. 
              ppphy_ann = ppphy_ann  + ppphy*DIM(mon)/365.
!              ppphy2_ann= ppphy2_ann + ppphy2*DIM(mon)/365.

!             2D data: 
              epc100_ann    = epc100_ann   + epc100*DIM(mon)/365. 
!              epcal100_ann  = epcal100_ann + epcal100*DIM(mon)/365.

              cflux_ann    = cflux_ann   + cflux*DIM(mon)/365. 
              oflux_ann    = oflux_ann   + oflux*DIM(mon)/365. 
              nfix_ann     = nfix_ann    + nfix*DIM(mon)/365. 
!             irondep_ann  = irondep_ann + irondep*DIM(mon)/365. 
          endif 
      enddo  ! time  (month) 
      CALL closeall ! close all open netcdf files

! ********** Do some basic calculations ************

!---------------------------------------------------
! (1) Global annual mean T(z) and S(z)  
!---------------------------------------------------
          vol =0.

!         DIC, TA, PH, O2 
          dicvol   = 0.
          caco3vol = 0.
          talvol   = 0.
          phvol    = 0.
          oxyvol   = 0.
!         POC, GOC, DOC
          pocvol   = 0.
          gocvol   = 0.
          docvol   = 0.
!         NO3, NH4, PO4, Si 
          no3vol   = 0.
          nh4vol   = 0.  
          po4vol   = 0.  
          sivol    = 0. 
!         PHY, PHY2, ZOO, ZOO2
          phyvol   = 0. 
          phy2vol  = 0. 
          zoovol   = 0. 
          zoo2vol  = 0. 
!         PPPHY, PPPHY2      
          ppphyvol = 0.
          ppphy2vol = 0.
                         
      do k=1, km   
         g_mask(:, :)  = t_mask(:, :, k) 

!        DIC, TA, PH, O2 
         CALL area_ave (e1t, e2t, e3t, g_mask, dic_ann(: , : , k),   imt, jmt, km, dicz,   dvol, k)  
!         CALL area_ave (e1t, e2t, e3t, g_mask, caco3_ann(: , : , k), imt, jmt, km, caco3z, dvol, k)  
         CALL area_ave (e1t, e2t, e3t, g_mask, tal_ann(:, :, k),     imt, jmt, km, talz,   dvol, k)  
         CALL area_ave (e1t, e2t, e3t, g_mask, oxy_ann(:, :, k),     imt, jmt, km, oxyz,   dvol, k)

!        POC, GOC, DOC
         CALL area_ave (e1t, e2t, e3t, g_mask, poc_ann(:, :, k), imt, jmt, km, pocz, dvol, k)
!         CALL area_ave (e1t, e2t, e3t, g_mask, goc_ann(:, :, k), imt, jmt, km, gocz, dvol, k)
!         CALL area_ave (e1t, e2t, e3t, g_mask, doc_ann(:, :, k), imt, jmt, km, docz, dvol, k)

!        NO3, NH4, PO4, Si 
         CALL area_ave (e1t, e2t, e3t, g_mask, no3_ann(:, :, k), imt, jmt, km, no3z, dvol, k)
!         CALL area_ave (e1t, e2t, e3t, g_mask, nh4_ann(:, :, k), imt, jmt, km, nh4z, dvol, k)
!         CALL area_ave (e1t, e2t, e3t, g_mask, po4_ann(:, :, k), imt, jmt, km, po4z, dvol, k)  
!         CALL area_ave (e1t, e2t, e3t, g_mask, si_ann(:, :, k),  imt, jmt, km, siz,  dvol, k)

!        PHY, PHY2, ZOO, ZOO2
         CALL area_ave (e1t, e2t, e3t, g_mask, phy_ann(:, :, k),  imt, jmt, km, phyz,  dvol, k)  
!         CALL area_ave (e1t, e2t, e3t, g_mask, phy2_ann(:, :, k), imt, jmt, km, phy2z, dvol, k)  
         CALL area_ave (e1t, e2t, e3t, g_mask, zoo_ann(:, :, k),  imt, jmt, km, zooz,  dvol, k)  
!         CALL area_ave (e1t, e2t, e3t, g_mask, zoo2_ann(:, :, k), imt, jmt, km, zoo2z, dvol, k) 

         if (exists) then
!            PPPHY, PPPHY2      
             CALL area_ave (e1t, e2t, e3t, g_mask, ph_ann(:, :, k),     imt, jmt, km, phz,     dvol, k)  
             CALL area_ave (e1t, e2t, e3t, g_mask, ppphy_ann(:, :, k),  imt, jmt, km, ppphyz,  dvol, k) 
 !            CALL area_ave (e1t, e2t, e3t, g_mask, ppphy2_ann(:, :, k), imt, jmt, km, ppphy2z, dvol, k) 
         endif 
!================================================================
!        Assign outputs
!================================================================
!        DIC, TA, PH, O2 
         dic_z(k)   = dicz     
!         caco3_z(k) = caco3z     
         tal_z(k)   = talz     
         oxy_z(k)   = oxyz    

!        POC, GOC, DOC
         poc_z(k)    = pocz       
!         goc_z(k)    = gocz       
!         doc_z(k)    = docz       

!        NO3, NH4, PO4, Si 
         no3_z(k)   = no3z  
!         nh4_z(k)   = nh4z   
!         po4_z(k)   = po4z  
!         si_z(k)    = siz       

!        PHY, PHY2, ZOO, ZOO2
         phy_z(k)   = phyz
!         phy2_z(k)  = phy2z
         zoo_z(k)   = zooz
!         zoo2_z(k)  = zoo2z
     
         if (exists) then
             ph_z(k)    = phz  
             ppphy_z(k) = ppphyz
!             ppphy2_z(k)= ppphy2z
         endif 

!        DIC, TA, PH, O2 
         dicvol   = dicvol + dicz*dvol  
!         caco3vol = caco3vol + caco3z*dvol  
         talvol = talvol + talz*dvol  
         oxyvol = oxyvol + oxyz*dvol  

!        POC, GOC, DOC
         pocvol = pocvol + pocz*dvol  
!         gocvol = gocvol + gocz*dvol  
!         docvol = docvol + docz*dvol  

!        NO3, NH4, PO4, Si 
         no3vol = no3vol + no3z*dvol  
!         nh4vol = nh4vol + nh4z*dvol  
!         po4vol = po4vol + po4z*dvol  
!         sivol  = sivol  + siz*dvol  


!        PHY, PHY2, ZOO, ZOO2
         phyvol    = phyvol  + phyz*dvol
!         phy2vol   = phy2vol + phy2z*dvol 
         zoovol    = zoovol  + zooz*dvol
!         zoo2vol   = zoo2vol + zoo2z*dvol 

         if (exists) then 
             phvol     =     phvol + phz*dvol    
             ppphyvol  =  ppphyvol + ppphyz*dvol  
!             ppphy2vol = ppphy2vol + ppphy2z*dvol  
         endif
                   
          vol = vol + dvol
      enddo  ! depth, k        

!     compute toc and ton
      toc = dicvol  + 106./16. * ( pocvol + phyvol + zoovol )                                                       ! &
!    &               + gocvol + docvol + caco3vol + phy2vol  + zoo2vol
!     convert from mmol C to Pg C      
      toc = toc * 12.0e-18
      ton = no3vol + phyvol + zoovol  + pocvol                                                                      ! &
!    &              + nh4vol + 16./122.*( phy2vol + zoo2vol + gocvol + docvol )   
!     convert to Pg      
      ton = ton * 14.007e-18

      if (vol.ne.0.) then         
!        DIC, TA, PH, O2 
         dicvol = dicvol/vol 
!         caco3vol = caco3vol/vol 
         talvol = talvol/vol 
         oxyvol = oxyvol/vol  
!        POC, GOC, DOC  
         pocvol = pocvol/vol  
!         gocvol = gocvol/vol  
!         docvol = docvol/vol  
!        NO3, NH4, PO4, Si 
         no3vol = no3vol/vol  
!         nh4vol = nh4vol/vol 
!         po4vol = po4vol/vol
!         sivol  = sivol/vol 

         if (exists) then   
             phvol  = phvol/vol 
!            convert to PgC/yr
             ppphyvol  = ppphyvol * 12.e-15 * 86400 * 365
!             ppphy2vol = ppphy2vol * 12.e-15 * 86400 * 365
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
         CALL area_ave_flx (e1t, e2t, g_mask, epc100_ann,   imt, jmt,   epc100glo, dum) 
!         CALL area_ave_flx (e1t, e2t, g_mask, epcal100_ann, imt, jmt, epcal100glo, dum) 

!        Cflux, Oflux, Nfix, Irondep  
         CALL area_ave_flx (e1t, e2t, g_mask, cflux_ann,   imt, jmt, cglo,       dum) 
         CALL area_ave_flx (e1t, e2t, g_mask, oflux_ann,   imt, jmt, ofluxglo,   dum) 
         CALL area_ave_flx (e1t, e2t, g_mask, nfix_ann,    imt, jmt, nfixglo ,   dum) 
!         CALL area_ave_flx (e1t, e2t, g_mask, irondep_ann, imt, jmt, irondepglo, dum) 

!       <PISCES OR 01/15/2014> convert into PgC/yr 
         cglo         = cglo        * dum * 12.e-15 * 86400 * 365
         epc100glo    = epc100glo   * dum * 12.e-15 * 86400 * 365
!         epcal100glo  = epcal100glo * dum * 12.e-15 * 86400 * 365

!        convert to TgN/yr (assuming this is N not N2)
         nfixglo      = nfixglo  * dum * 14.007e-12 * 86400 * 365

!        convert to mol O2/yr
         ofluxglo     = ofluxglo * dum *  86400 * 365

!        convert to mol Fe/yr 
!         irondepglo   = irondepglo  * dum * 86400 * 365
      endif

!     Read in the year which is the first command line argument
      CALL getarg(1, year_arg_in )
      read (year_arg_in,'(I10)') iyear

      tyear = float( iyear )
      print*," --- "
      print*," Tyear is:", tyear
      print*," --- "

!---------------------------------------------------------
!     NETCDF RTD output: Time series information section
!---------------------------------------------------------
      iou = 0
      id_time = 0
      id_z = 0

!     If the output file does not exist, create it and define dims and vars
      inquire (file="nemo_carbon_rtd.nc", exist=exists1)
      if (.not. exists1) then
      print*,"output file not found...creating a new file..."
      CALL flush(6)
      CALL opennew ("nemo_carbon_rtd.nc", iou)
      ntrec = 1
        CALL redef (iou)
!         basic grid specification
          CALL defdim ('time', iou, 0, id_time)
          CALL defdim ('depth', iou, km, id_z)
          CALL defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'D'             &
     &        , 'time', 'time', 'common_year since 1-1-0 00:00:0.0')
          CALL putatttext (iou, 'time', 'calendar', 'noleap')

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

!         GOC
!          CALL defvar ('GOC', iou, 1, (/id_time/), -1.e4                         &
!     &        , 1.e4,' ', 'F'                                                    &
!     &        , 'Global mean large Particulate Organic Carbon'                   &
!     &        , 'GOC', 'mmol m^-3')

!          CALL defvar ('GOCz', iou, 2, (/id_z, id_time/), -1.e4                  & 
!     &        , 1.e4,' ', 'F', 'Large GOC by level'                              &
!     &        , 'GOCz', 'mmol m^-3')

!         DOC
!          CALL defvar ('DOC', iou, 1, (/id_time/), -1.e4                         &   
!     &        , 1.e4,' ', 'F'                                                    &
!     &        , 'Global mean Dissolved Organic Carbon'                           &  
!     &        , 'DOC', 'mmol m^-3')

!          CALL defvar ('DOCz', iou, 2, (/id_z, id_time/), -1.e4                  &
!     &        , 1.e4,' ', 'F', 'DOC by level'                                    &
!     &        , 'DOCz', 'mmol m^-3')


!         NO3
          CALL defvar ('NO3', iou, 1, (/id_time/), -1.e4                         &
     &        , 1.e4,' ', 'F', 'Global mean NO3'                                 &  
     &        , 'NO3', 'uM')

          CALL defvar ('NO3z', iou, 2, (/id_z, id_time/), -1.e4                  &
     &        , 1.e4,' ', 'F', 'NO3 by level'                                    & 
     &        , 'NO3z', 'uM')

!         NH4
!          CALL defvar ('NH4', iou, 1, (/id_time/), -1.e4                         &
!     &        , 1.e4,' ', 'F', 'Global mean NH4'                                 &  
!     &        , 'NH4', 'uM')

!          CALL defvar ('NH4z', iou, 2, (/id_z, id_time/), -1.e4                  &
!     &        , 1.e4,' ', 'F', 'NH4 by level'                                    &
!     &        , 'NH4z', 'uM')

!         PO4
!          CALL defvar ('PO4', iou, 1, (/id_time/), -1.e4                         &
!     &        , 1.e4,' ', 'F', 'Global mean PO4'                                 & 
!     &        , 'PO4', 'uM')

!          CALL defvar ('PO4z', iou, 2, (/id_z, id_time/), -1.e4                  &
!     &        , 1.e4,' ', 'F', 'PO4 by level'                                    &
!     &        , 'PO4z', 'uM')

!         SI
!          CALL defvar ('SI', iou, 1, (/id_time/), -1.e4                          &  
!     &        , 1.e4,' ', 'F', 'Global mean Si'                                  & 
!     &        , 'SI', 'uM') 

!          CALL defvar ('SIz', iou, 2, (/id_z, id_time/), -1.e4                   &
!     &        , 1.e4,' ', 'F', 'Si by level'                                     &
!     &        , 'SIz', 'uM')

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


!             PPPHY2
!              CALL defvar ('PPPHY2', iou, 1, (/id_time/), -1.e4                  &  
!     &            , 1.e4,' ', 'F', 'Primary production of Diatoms'               & 
!     &            , 'PPPHY2', 'PgC/yr')

!             EPC100
              CALL defvar ('EPC100', iou, 1, (/id_time/), -1.e4                  &  
     &            , 1.e4,' ', 'F', 'Export of carbon particles at 100m'          &  
     &            , 'EPC100', 'PgC/yr')

!             EPCAL100
!              CALL defvar ('EPCAL100', iou, 1, (/id_time/), -1.e4                &  
!     &            , 1.e4,' ', 'F', 'Export of Calcite at 100m'                   &
!     &            , 'EPCAL100', 'PgC/yr')

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
!             Irondep
!              CALL defvar ('Irondep', iou, 1, (/id_time/), -1.e4                 & 
!     &           , 1.e4,' ', 'F', 'Iron deposition'                              & 
!     &           , 'Irondep', 'mol/yr')
          endif

          CALL enddef (iou)
!         define the depth axis
          CALL putvara ('depth', iou, km, (/1/), (/km/)                          & 
     &      , deptht, 1., 0.)
 
      else

!      if the file does exist then open it for writing at the next record
       print*,"output file found...opening existing file for appending"
       CALL flush(6)
       CALL opennext ("nemo_carbon_rtd.nc", tyear, ntrec, iou)
      endif

!       append variables
!--------------------------------------------------------------------------
!       time
        CALL putvars ('time', iou, ntrec, tyear, 1., 0.)

!       DIC
        CALL putvars ('DIC', iou, ntrec, dicvol, 1., 0.)
        CALL putvara ('DICz', iou, km, (/1, ntrec/), (/km, 1/), dic_z, 1., 0.) 

!       CaCO3
!        CALL putvars ('CaCO3', iou, ntrec, caco3vol, 1., 0.)
!        CALL putvara ('CaCO3z', iou, km, (/1, ntrec/), (/km, 1/), caco3_z, 1., 0.) 

!       ALK
        CALL putvars ('TAL', iou, ntrec, talvol, 1., 0.)
        CALL putvara ('TALz', iou, km, (/1, ntrec/), (/km, 1/), tal_z, 1., 0.)

!       OXY
        CALL putvars ('O2', iou, ntrec, oxyvol, 1., 0.)
        CALL putvara ('O2z', iou, km, (/1, ntrec/), (/km, 1/), oxy_z, 1., 0.)

!       POC
        CALL putvars ('POC', iou, ntrec, pocvol, 1., 0.)
        CALL putvara ('POCz', iou, km, (/1, ntrec/), (/km, 1/), poc_z, 1., 0.)
!       GOC
!        CALL putvars ('GOC', iou, ntrec, gocvol, 1., 0.)
!        CALL putvara ('GOCz', iou, km, (/1, ntrec/), (/km, 1/), goc_z, 1., 0.)
!       DOC
!        CALL putvars ('DOC', iou, ntrec, docvol, 1., 0.)
!        CALL putvara ('DOCz', iou, km, (/1, ntrec/), (/km, 1/), doc_z, 1., 0.)

!       NO3
        CALL putvars ('NO3', iou, ntrec, no3vol, 1., 0.)
        CALL putvara ('NO3z', iou, km, (/1, ntrec/), (/km, 1/), no3_z, 1., 0.)

!       NH4
!        CALL putvars ('NH4', iou, ntrec, nh4vol, 1., 0.)
!        CALL putvara ('NH4z', iou, km, (/1, ntrec/), (/km, 1/), nh4_z, 1., 0.)

!       PO4
!        CALL putvars ('PO4', iou, ntrec, po4vol, 1., 0.)
!        CALL putvara ('PO4z', iou, km, (/1, ntrec/), (/km, 1/), po4_z, 1., 0.)

!       SI
!        CALL putvars ('SI', iou, ntrec, sivol, 1., 0.)
!        CALL putvara ('SIz', iou, km, (/1, ntrec/), (/km, 1/), si_z, 1., 0.)

!       Total C
        CALL putvars ('TC', iou, ntrec, toc, 1., 0.)

!       Total N
        CALL putvars ('TN', iou, ntrec, ton, 1., 0.)

!       Diagnostic variables
        if (exists) then 
!           PH
            CALL putvars ('PH', iou, ntrec, phvol, 1., 0.)
            CALL putvara ('PHz', iou, km, (/1, ntrec/), (/km, 1/), ph_z, 1., 0.)

!           PHY
            CALL putvars ('PPPHY', iou, ntrec, ppphyvol, 1., 0.)

!           PHY2
!            CALL putvars ('PPPHY2', iou, ntrec, ppphy2vol, 1., 0.)

!           EPC100
            CALL putvars ('EPC100', iou, ntrec, epc100glo, 1., 0.)

!           EPCAL100
!            CALL putvars ('EPCAL100', iou, ntrec, epcal100glo, 1., 0.)

!           Cflux
            CALL putvars ('CFLX', iou, ntrec, cglo, 1., 0.)

!           Oflux
            CALL putvars ('OFLX', iou, ntrec, ofluxglo, 1., 0.)

!           Nfix
            CALL putvars ('NFIX', iou, ntrec, nfixglo, 1., 0.)

!           Irondep
!            CALL putvars ('Irondep', iou, ntrec, irondepglo, 1., 0.)
        endif
      print*, 'closing netcdf'
      CALL flush(6)

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
              do j=1,jmt
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
SUBROUTINE area_ave_flx (e1,e2, mask, a, imt, jmt,a_mean,ss)
      implicit none
      integer imt, jmt, i, j
      real e1(imt,jmt),e2(imt,jmt) 
      real a(imt,jmt), mask(imt,jmt) 
      real a_mean, ss, s1, arc

          s1=0.
          ss=0.
          do i=1,imt-2  ! not to double count the cyclic boundary
              do j=1,jmt
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
