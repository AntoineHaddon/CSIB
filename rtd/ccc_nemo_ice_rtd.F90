! Compute CCCma RTDs from NEMO-CICEi/LIM timeseries files.
!
PROGRAM nemo_ice_diag 
! ======================================================================
! PROGRAM nemo_ice_diag
! ---------------------------------------------------
!
! Description:
!     Reads the 'original' NEMO icemod timeseries files (defined in 
!     iodef), and computes various 1D run time diagnostics for 
!     plottting at CCCma. As long at the required fields are output
!     this will work for LIM or CICE.
!  
!
! Method:
!     Area averaging based on NEMO physical RTD.
!     NetCDF writing with the aid of the UVic library.
!
! Input files:
!     icemod - the NEMO monthly frequency (_1m_) icemod.nc file in 
!               NetCDF
!     
!     orca_mesh_mask - the mesh mask file defining the ORCA grid,
!              specifically the mask_t and e1t, e2t fields.
!
! Output files:
!     nemo_ice_rtd.nc - the 1-D RTDs (noleap time axis).
! 
! History:
!
! Version    Date         Comment
! -------    ----         -------
!            Dec 15/15   Abstract all calculations to ccc_nemo_rtd_utils
!                        module, which is shared between all rtd.
!
!            OCT 29/15    Initial implementation based on physical RTD
!                         <Neil Swart>
!
! ======================================================================
      IMPLICIT NONE 
      INTEGER :: imt, jmt, lm, iou

!         establish the size of the grid from the input file.
          call openfile ("icemod", iou) 
          call getdimlen ('x', iou, imt) 
          call getdimlen ('y', iou, jmt) 
          call getdimlen ('time_counter', iou, lm)

!         do the calculations and save the output netcdf   
          call calc (imt, jmt, lm)

END PROGRAM nemo_ice_diag

SUBROUTINE calc (imt, jmt, lm)
!     Does the required calculations and saves the output to netcdf
      USE ccc_nemo_rtd_utils, only: area_ave, area_ave_flx, noleap_days
      IMPLICIT NONE 
      integer, parameter:: dp=kind(0.d0) ! double precision
      INTEGER :: i, j, k, imt, jmt, lm, cur_mon
      INTEGER :: iou1, iou2

! ======================================================================
!     Input data 
! ======================================================================
!     Grid-related arrays
      REAL, DIMENSION(imt, jmt)     :: lon2d, lat2d
      REAL, DIMENSION(imt, jmt)     :: e1t, e2t, t_mask
!     Monthly ice fraction, thickness, surface temp, i-vel, j-vel
      REAL, DIMENSION(imt, jmt, lm) :: soicecov, iicethic, &
          & iicetemp, iicevelu, iicevelv
!     Monthly fluxes: uflx, vflx, ocean heat flux at base, sublimation 
      REAL, DIMENSION(imt, jmt, lm) :: iicestru, iicestrv, &
          & qt_ice_oce, sublim_over_sea_ice, qtr_ice_bot  
!     Monthly snow fields: snow thickness, snow precip, snow precip 
!                          over ice
      REAL, DIMENSION(imt, jmt, lm) :: isnowthi, isnowpre, hflx_snow_ai_cea, hflx_err, & 
          &  snow_over_sea_ice, aicesflx, aicenflx, iicesflx, iicenflx, iicetflx

! ======================================================================
!     Output data 
! ======================================================================
!     Integral quantities: ice area, extent, volume
      REAL, DIMENSION(lm)      :: area_nh, area_sh, extent_nh      &
          & , extent_sh, vol_nh, vol_sh, vol_sno_nh, vol_sno_sh

!     Average quantities: ice thickness
      REAL, DIMENSION(lm)      :: thick_nh, thick_sh
      REAL, DIMENSION(lm)      :: surf_temp_nh, surf_temp_sh
      REAL, DIMENSION(lm)      :: ivelu_nh, ivelu_sh, ivelv_nh, ivelv_sh
      REAL, DIMENSION(lm)      :: itauu_nh, itauu_sh, itauv_nh, itauv_sh
      REAL, DIMENSION(lm)      :: isnowthick_nh, isnowthick_sh
      REAL, DIMENSION(lm)      :: iohflx_nh, iohflx_sh
      REAL, DIMENSION(lm)      :: test_calc
      REAL, DIMENSION(lm)      :: aicesflx_ave, aicenflx_ave,hflx_snow_ai_cea_ave,hflx_err_ave
      REAL, DIMENSION(lm)      :: iicesflx_ave, iicenflx_ave,iicetflx_ave

! ======================================================================
!     Working Arrays
! ======================================================================

      REAL, DIMENSION(imt, jmt)       :: nh_mask, sh_mask, tarea
      REAL, DIMENSION(imt, jmt, lm)   :: extent_tf, thick_over_ice, snow_thick_over_ice
      REAL, DIMENSION(imt, jmt, lm)   :: ipres_mask
      REAL, DIMENSION(imt, jmt)   :: ipres_mask_nh, ipres_mask_sh
      REAL                        :: ss                      ! area

!----------------
!  NetCDF-output specific
      INTEGER :: id_time, id_z, iou, ntrec, ntrec2, iyear, imon
      INTEGER :: days_elapsed
      LOGICAL :: exists, exists1, notopen
      REAL    :: tyear, tmon, norec, tdays_elapsed
      CHARACTER :: fname01*100,fname02*100
      CHARACTER(len=32) :: year_arg_in, mon_arg_in

!----------------
           print*,'Reading data on NEMO grid...'

         iou1 =0
         iou2 =0

!---------------------------------------------------
!    Define NetCDF files   
!---------------------------------------------------
        fname01='icemod'
        fname02='orca_mesh_mask'
!---------------------------------------------------
!    Open the defined NetCDF files   
!---------------------------------------------------
      call openfile (fname01,iou1)
      call openfile (fname02,iou2)

!---------------------------------------------------
!    Get grid/mask data   
!---------------------------------------------------
      call getvara ('nav_lon', iou1, imt*jmt, (/1,1,1/), (/imt,jmt,1/) &
          &         , lon2d, 1., 0.)
      call getvara ('nav_lat', iou1, imt*jmt, (/1,1,1/), (/imt,jmt,1/) &
          &         , lat2d, 1., 0.)
      call getvara ('e1t', iou2, imt*jmt, (/1,1,1/), (/imt,jmt,1/)     &
          &         , e1t, 1., 0.)
      call getvara ('e2t', iou2, imt*jmt, (/1,1,1/), (/imt,jmt,1/)     &
          &         , e2t , 1., 0.)
      call getvara ('tmask', iou2, imt*jmt, (/1,1,1,1/)             &
          &         , (/imt,jmt,1,1/),t_mask , 1., 0.)

! compute total grid cell area for valid ocean points   
      tarea(:, :)   = e1t(:, :)*e2t(:, :)*t_mask(:, :)
!---------------------------------------------------
!    Construct masks for NH and SH sub-regions   
!---------------------------------------------------
      nh_mask(:,:) = 0.0_dp
      sh_mask(:,:) = 0.0_dp
      
      where (lat2d .gt. 0.0_dp) ! logical array.
          nh_mask = 1.0_dp
      elsewhere
          sh_mask = 1.0_dp
      endwhere

      nh_mask = nh_mask * t_mask
      sh_mask = sh_mask * t_mask

      print*,'NH ocean area (x10**6 km**2): ', SUM(nh_mask*tarea)/1e12
      print*,'SH ocean area (x10**6 km**2): ', SUM(sh_mask*tarea)/1e12
      print*,'Ocean area (x10**6 km**2): ',    SUM(tarea)/1e12

!---------------------------------------------------
!    Load ice data from NetCDF
!---------------------------------------------------
! Ice Fraction
      call getvara ('siconc', iou1, imt*jmt*lm, (/1,1,1/)            &
          & , (/imt,jmt,lm/), soicecov, 1., 0.)
      soicecov=soicecov/100 ! from % to factor
! Ice thickness (cell average)
      call getvara ('sithick', iou1, imt*jmt*lm, (/1,1,1/)            &
          & , (/imt,jmt,lm/), iicethic, 1., 0.)
! Ice surface temperature (cell average)
      call getvara ('iicetemp', iou1, imt*jmt*lm           &
          & ,(/1,1,1/), (/imt,jmt,lm/), iicetemp, 1., 0.)
! Ice velocity along i-axis at I-point (ice presence average)
      call getvara ('siu', iou1, imt*jmt*lm                       &
          & ,(/1,1,1/), (/imt,jmt,lm/), iicevelu, 1., 0.)
! Ice velocity along j-axis at I-point (ice presence average)
      call getvara ('siv', iou1, imt*jmt*lm                       &
          & ,(/1,1,1/), (/imt,jmt,lm/), iicevelv, 1., 0.)
! Wind stress along i-axis over the ice at i-point
      call getvara ('sistrxdtop', iou1, imt*jmt*lm                       &
          & ,(/1,1,1/), (/imt,jmt,lm/), iicestru, 1., 0.)
! Wind stress along j-axis over the ice at i-point
      call getvara ('sistrydtop', iou1, imt*jmt*lm                       &
          & ,(/1,1,1/), (/imt,jmt,lm/), iicestrv, 1., 0.)
! Oceanic total heat flux at ice base 
     call getvara ('qt_ice_oce', iou1, imt*jmt*lm                       &
         & ,(/1,1,1/), (/imt,jmt,lm/), qt_ice_oce, 1., 0.)
! Snow thickness (cell average)
      call getvara ('sisnthick', iou1, imt*jmt*lm                       &
          & ,(/1,1,1/), (/imt,jmt,lm/), isnowthi, 1., 0.)
! Solar heat flux over ice
      call getvara ('qsr_ice', iou1, imt*jmt*lm                       &
          & ,(/1,1,1/), (/imt,jmt,lm/), aicesflx, 1., 0.)
! Non Solar heat flux over ice
      call getvara ('qns_ice', iou1, imt*jmt*lm                       &
          & ,(/1,1,1/), (/imt,jmt,lm/), aicenflx, 1., 0.)
! Solar heat flux under the ice 
      call getvara ('qtr_ice_bot', iou1, imt*jmt*lm                       &
          & ,(/1,1,1/), (/imt,jmt,lm/), iicesflx, 1., 0.)
! Heat flux from the snow precipitation of ice
      call getvara ('hflx_snow_ai_cea', iou1, imt*jmt*lm,                   &
          & (/1,1,1/), (/imt,jmt,lm/), hflx_snow_ai_cea, 1., 0.)
! Heat flux error
      call getvara ('hfxerr', iou1, imt*jmt*lm,                   &
          & (/1,1,1/), (/imt,jmt,lm/), hflx_err, 1., 0.)

! Non solar heat flux under the ice (total - solar )
      iicenflx= qt_ice_oce - iicesflx
! Non solar radiation in ice 
      aicenflx = aicenflx - hflx_snow_ai_cea

! Hold these for now.
! Sublimation over sea-ice (cell average)
!      call getvara ('sublim_over_sea_ice', iou1, imt*jmt*lm            &
!          & ,(/1,1,1/), (/imt,jmt,lm/), sublim_over_sea_ice, 1., 0.)
! Snow precipitation
!      call getvara ('snowpre', iou1, imt*jmt*lm                       &
!          & ,(/1,1,1/), (/imt,jmt,lm/), isnowpre, 1., 0.)
! Snow over sea-ice (cell average) [kg/m2/s]
!      call getvara ('snow_over_sea_ice', iou1, imt*jmt*lm              &
!          & ,(/1,1,1/), (/imt,jmt,lm/), snow_over_sea_ice, 1., 0.)

!---------------------------------------------------
!    Pre-calculate boolean mask for extent and thicknes over ice
!---------------------------------------------------

! find where ice fraction > 0.15 for extent calculation
      extent_tf(:,:,:) = 0.0_dp
      where (soicecov .gt. 0.15_dp) ! logical array.
          extent_tf = 1.0_dp
      elsewhere
          extent_tf = 0.0_dp
      end where

! Thickness is given as a grid-cell ave. Convert to thickness over ice.
      thick_over_ice(:,:,:) = 0.0_dp
      where (soicecov .gt. 0.0_dp) ! logical array.
          thick_over_ice = iicethic / soicecov
      endwhere
! Repeat for snow Thickness is given as a grid-cell ave. Convert to thickness over ice.
      snow_thick_over_ice(:,:,:) = 0.0_dp
      where (soicecov .gt. 0.0_dp) ! logical array.
          snow_thick_over_ice = isnowthi / soicecov
      endwhere

! Now make a mask so we don't include 0s (ice free) in mean calculations
      ipres_mask(:,:,:) = 0.0_dp
      where (soicecov .gt. 0.0_dp) ! logical array.
          ipres_mask = 1.0_dp
      endwhere
          
!---------------------------------------------------
!    Main loop over all months
!---------------------------------------------------
      do cur_mon = 1, lm
!
!   AREA INTEGRALS for each hemisphere
!
! calculate sea-ice area in each hemisphere
!          area_nh(cur_mon) = SUM(soicecov(1:imt-2, :, cur_mon)*e1t(1:imt-2,:)*e2t(1:imt-2,:), MASK=nh_mask(1:imt-2,:) > 0.5)
!          area_nh(cur_mon) = SUM(soicecov(1:imt-2, :, cur_mon)*e1t(1:imt-2,:)*e2t(1:imt-2,:)*nh_mask(1:imt-2,:))
          area_nh(cur_mon) = SUM(soicecov(:,:,cur_mon)*nh_mask*tarea)
! calculate sea-ice area differently to test the calc
!          CALL area_ave_flx(e1t, e2t, nh_mask                           &
!              &            , soicecov(:, :, cur_mon)                   & 
!              &            , imt, jmt, test_calc(cur_mon), ss)
!          print*, area_nh(cur_mon), test_calc(cur_mon)*ss    
!          print*, (area_nh(cur_mon) - test_calc(cur_mon)*ss) / area_nh(cur_mon)*100
     
          area_sh(cur_mon) = SUM(soicecov(:, :, cur_mon)*sh_mask*tarea)
! calculate sea-ice extent in each hemisphere
          extent_nh(cur_mon) = SUM(extent_tf(:, :, cur_mon)*nh_mask*tarea)
          extent_sh(cur_mon) = SUM(extent_tf(:, :, cur_mon)*sh_mask*tarea)
! calculate sea-ice volume in each hemisphere
          vol_nh(cur_mon) = SUM(iicethic(:, :, cur_mon)*nh_mask*tarea)
          vol_sh(cur_mon) = SUM(iicethic(:, :, cur_mon)*sh_mask*tarea)
! calculate snow volume in each hemisphere
          vol_sno_nh(cur_mon) = SUM(isnowthi(:, :, cur_mon)*nh_mask*tarea)
          vol_sno_sh(cur_mon) = SUM(isnowthi(:, :, cur_mon)*sh_mask*tarea)
!
!     AREA AVERAGES for each hemisphere
!

! Construct masks to make sure averages only where sea ice is present
          ipres_mask_nh = ipres_mask(:, :, cur_mon) * nh_mask * t_mask
          ipres_mask_sh = ipres_mask(:, :, cur_mon) * sh_mask * t_mask

! calculate sea-ice thickness 
          CALL area_ave_flx(e1t, e2t, ipres_mask_nh*soicecov(:,:,cur_mon)                     &
              &            , thick_over_ice(:, :, cur_mon)             & 
              &            , imt, jmt, thick_nh(cur_mon), ss)
          CALL area_ave_flx(e1t, e2t, ipres_mask_sh*soicecov(:,:,cur_mon)   &
              &            , thick_over_ice(:, :, cur_mon)             & 
              &            , imt, jmt, thick_sh(cur_mon), ss)
! Snow thickness
          CALL area_ave_flx(e1t, e2t, ipres_mask_nh*soicecov(:,:,cur_mon),                         &
              & snow_thick_over_ice(:, :, cur_mon), imt, jmt                      &
              &            , isnowthick_nh(cur_mon), ss)
          CALL area_ave_flx(e1t, e2t, ipres_mask_sh*soicecov(:,:,cur_mon),                         &
              & snow_thick_over_ice(:, :, cur_mon), imt, jmt,                     &
              & isnowthick_sh(cur_mon), ss)
! calculate surface temp 
          CALL area_ave_flx(e1t, e2t, ipres_mask_nh*soicecov(:,:,cur_mon),                    &
              &  iicetemp(:, :, cur_mon), imt, jmt                      &
              &            , surf_temp_nh(cur_mon), ss)
          CALL area_ave_flx(e1t, e2t, ipres_mask_sh*soicecov(:,:,cur_mon),                    &
              & iicetemp(:, :, cur_mon), imt, jmt,                      &
              & surf_temp_sh(cur_mon), ss)
! calculate ice u velocities 
          CALL area_ave_flx(e1t, e2t, ipres_mask_nh*soicecov(:,:,cur_mon),                    &
              &  iicevelu(:, :, cur_mon), imt, jmt                      &
              &            , ivelu_nh(cur_mon), ss)
          CALL area_ave_flx(e1t, e2t, ipres_mask_sh,                    &
              & iicevelu(:, :, cur_mon), imt, jmt,                      &
              & ivelu_sh(cur_mon), ss)
! calculate ice v velocities 
          CALL area_ave_flx(e1t, e2t, ipres_mask_nh*soicecov(:,:,cur_mon),                    &
              &  iicevelv(:, :, cur_mon), imt, jmt                      &
              &            , ivelv_nh(cur_mon), ss)
          CALL area_ave_flx(e1t, e2t, ipres_mask_sh,                    &
              & iicevelv(:, :, cur_mon), imt, jmt,                      &
              & ivelv_sh(cur_mon), ss)
! calculate ice u wind stress
          CALL area_ave_flx(e1t, e2t, ipres_mask_nh*soicecov(:,:,cur_mon),                    &
              &  iicestru(:, :, cur_mon), imt, jmt                      &
              &            , itauu_nh(cur_mon), ss)
          CALL area_ave_flx(e1t, e2t, ipres_mask_sh*soicecov(:,:,cur_mon),                    &
              & iicestru(:, :, cur_mon), imt, jmt,                      &
              & itauu_sh(cur_mon), ss)
! calculate ice v wind stress
          CALL area_ave_flx(e1t, e2t, ipres_mask_nh*soicecov(:,:,cur_mon),                    &
              &  iicestrv(:, :, cur_mon), imt, jmt                      &
              &            , itauv_nh(cur_mon), ss)
          CALL area_ave_flx(e1t, e2t, ipres_mask_sh*soicecov(:,:,cur_mon),                    &
              & iicestrv(:, :, cur_mon), imt, jmt,                      &
              & itauv_sh(cur_mon), ss)
! calculate oceanic heat flux at ice base  (remove for now, variable not in SI3 , NL)
         CALL area_ave_flx(e1t, e2t, ipres_mask_nh*soicecov(:,:,cur_mon),                    &
             &  qt_ice_oce(:, :, cur_mon), imt, jmt                      &
             &            , iohflx_nh(cur_mon), ss)
         CALL area_ave_flx(e1t, e2t, ipres_mask_sh*soicecov(:,:,cur_mon),                    &
             & qt_ice_oce(:, :, cur_mon), imt, jmt,                      &
             & iohflx_sh(cur_mon), ss)
!  Solar and non solar heat fluxes from atmosphere
          CALL area_ave_flx(e1t, e2t, t_mask,                    &
              & aicenflx(:, :, cur_mon), imt, jmt,                      &
              & aicenflx_ave(cur_mon), ss)
          CALL area_ave_flx(e1t, e2t, t_mask,                    &
              & hflx_snow_ai_cea(:, :, cur_mon), imt, jmt,                      &
              & hflx_snow_ai_cea_ave(cur_mon), ss)
          CALL area_ave_flx(e1t, e2t, t_mask,                    &
              & hflx_err(:, :, cur_mon), imt, jmt,                      &
              & hflx_err_ave(cur_mon), ss)
          CALL area_ave_flx(e1t, e2t, t_mask,                    &
              & aicesflx(:, :, cur_mon), imt, jmt,                      &
              & aicesflx_ave(cur_mon), ss)
!  Solar and non solar heat fluxes from ice to ocean 
          CALL area_ave_flx(e1t, e2t, t_mask,                    &
              & qt_ice_oce(:, :, cur_mon), imt, jmt,                      &
              & iicetflx_ave(cur_mon), ss)
          CALL area_ave_flx(e1t, e2t, t_mask,                    &
              & iicesflx(:, :, cur_mon), imt, jmt,                      &
              & iicesflx_ave(cur_mon), ss)
          iicenflx_ave(cur_mon) = iicetflx_ave(cur_mon) - iicesflx_ave(cur_mon)
      enddo 

!---------------------------------------------------
!   Main outputs 
!--------------------------------------------------

      print*,'-------------------------------------'
      print*,'    Sea ice area (x10^6 km^2)        '
      print*,'-------------------------------------'
      print*,'NH     ', area_nh/1e12
      print*,'SH     ', area_sh/1e12
      print*,'-------------------------------------'
      print*,'    Sea ice extent (x10^6 km^2)      '
      print*,'-------------------------------------'
      print*,'NH     ', extent_nh/1e12
      print*,'SH     ', extent_sh/1e12
      print*,'-------------------------------------'
      print*,'    Sea ice volume (x10^3 km^3)      '
      print*,'-------------------------------------'
      print*,'NH     ', vol_nh/1e12
      print*,'SH     ', vol_sh/1e12
      print*,'-------------------------------------'
      print*,'    Average sea ice thickness (m)    '
      print*,'-------------------------------------'
      print*,'NH     ', thick_nh
      print*,'SH     ', thick_sh
      print*,'-------------------------------------'
      print*,'    Average sea ice surface temp (C) '
      print*,'-------------------------------------'
      print*,'NH     ', surf_temp_nh
      print*,'SH     ', surf_temp_sh
      print*,'-------------------------------------'
      print*,'    Average ice u-velocities (m/s) '
      print*,'-------------------------------------'
      print*,'NH     ', ivelu_nh
      print*,'SH     ', ivelu_sh
      print*,'-------------------------------------'
      print*,'    Average ice v-velocities (m/s) '
      print*,'-------------------------------------'
      print*,'NH     ', ivelv_nh
      print*,'SH     ', ivelv_sh
      print*,'-------------------------------------'
      print*,'    Average u windstress on ice (Pa) '
      print*,'-------------------------------------'
      print*,'NH     ', itauu_nh
      print*,'SH     ', itauu_sh
      print*,'-------------------------------------'
      print*,'    Average v windstress on ice (Pa) '
      print*,'-------------------------------------'
      print*,'NH     ', itauv_nh
      print*,'SH     ', itauv_sh
      print*,'-------------------------------------'
      print*,'   Oceanic heat flux at base (W/m^2) '
      print*,'-------------------------------------'
      print*,'NH     ', iohflx_nh
      print*,'SH     ', iohflx_sh
      print*,'-------------------------------------'
      print*,'    Average snow thickness on ice(m) '
      print*,'-------------------------------------'
      print*,'NH     ', isnowthick_nh
      print*,'SH     ', isnowthick_sh
      print*,'-------------------------------------'
      print*,'   ATM-ICE Solar/nonsolar flux  (W/m^2) '
      print*,'-------------------------------------'
      print*,'S      ', aicesflx_ave
      print*,'NS     ', aicenflx_ave
      print*,'NET    ', aicenflx_ave + aicesflx_ave
      print*,'SNOW   ', hflx_snow_ai_cea_ave
      print*,'Error  ', hflx_err_ave
      print*,'-------------------------------------'
      print*,'   ICE-OCE Solar/nonsolar flux  (W/m^2) '
      print*,'-------------------------------------'
      print*,'S      ', iicesflx_ave
      print*,'NS     ', iicenflx_ave
      print*,'NET    ', iicenflx_ave + iicesflx_ave
     

!---------------------------------------------------------
!     RTD output: Time series information section
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

!     NETCDF output
!     If the output file does not exist, create it and define dims and vars
      inquire (file="nemo_ice_rtd.nc", exist=exists1)

      if (.not. exists1) then 
          print*,"output file not found...creating a new file..."
          call opennew ("nemo_ice_rtd.nc", iou) 
          ntrec = 1
          call redef (iou)

!         basic grid specification
          call defdim ('time', iou, 0, id_time)
          call defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'F'   &    
     &        , 'time', 'time', 'days since 0000-01-01 00:00:00')

          call putatttext (iou, 'time', 'calendar', '365_day')

!   Area NH
          call defvar ('Area_NH', iou, 1, (/id_time/), 0               &    
     &        , 1.e15,' ', 'D', 'Northern Hemisphere sea-ice area'     &    
     &        , 'sea_ice_area', 'm^2')     
!   Area SH
          call defvar ('Area_SH', iou, 1, (/id_time/), 0               &    
     &        , 1.e15,' ', 'D', 'Southern Hemisphere sea-ice area'     &    
     &        , 'sea_ice_area', 'm^2')     
!   Extent NH
          call defvar ('Extent_NH', iou, 1, (/id_time/), 0              &    
     &        , 1.e15,' ', 'D', 'Northern Hemisphere sea-ice extent'    &    
     &        , 'sea_ice_extent', 'm^2')     
!   Extent SH
          call defvar ('Extent_SH', iou, 1, (/id_time/), 0              &    
     &        , 1.e15,' ', 'D', 'Southern Hemisphere sea-ice extent'    &    
     &        , 'sea_ice_extent', 'm^2')     
!   Volume NH
          call defvar ('Volume_NH', iou, 1, (/id_time/), 0              &    
     &        , 1.e15,' ', 'D', 'Northern Hemisphere sea-ice volume'    &    
     &        , 'sea_ice_volume', 'm^3')     
!   Volume SH
          call defvar ('Volume_SH', iou, 1, (/id_time/), 0              &    
     &        , 1.e15,' ', 'D', 'Southern Hemisphere sea-ice volume'    &    
     &        , 'sea_ice_volume', 'm^3')     
!   Volume NH snow
          call defvar ('Volume_snow_NH', iou, 1, (/id_time/), 0              &    
     &        , 1.e15,' ', 'D', 'Northern Hemisphere snow volume on ice'    &    
     &        , '', 'm^3')     
!   Volume SH snow
          call defvar ('Volume_snow_SH', iou, 1, (/id_time/), 0              &    
     &        , 1.e15,' ', 'D', 'Southern Hemisphere sea-ice volume on ice'    &    
     &        , '', 'm^3')     
!   Thickness NH
          call defvar ('Thickness_NH', iou, 1, (/id_time/), 0           &    
     &        , 1.e15,' ', 'D', 'Northern Hemisphere sea-ice thickness' &    
     &        , 'sea_ice_thickness', 'm')     
!   Thickness SH
          call defvar ('Thickness_SH', iou, 1, (/id_time/), 0           &    
     &        , 1.e15,' ', 'D', 'Southern Hemisphere sea-ice thickness' &    
     &        , 'sea_ice_thickness', 'm')     
!   Snow thickness NH
          call defvar ('Snow_thickness_NH', iou, 1, (/id_time/), 0      &    
     &        , 1.e15,' ', 'D'                                          &
     &        , 'Northern Hemisphere snow thickness over sea ice'       &    
     &        , 'surface_snow_thickness_where_sea_ice', 'm')     
!   Snow thickness SH
          call defvar ('Snow_thickness_SH', iou, 1, (/id_time/), 0      &    
     &        , 1.e15,' ', 'D'                                          &
     &        , 'Southern Hemisphere snow thickness over sea ice'       &    
     &        , 'surface_snow_thickness_where_sea_ice', 'm')     

!   Solar heat flux
          call defvar ('aicesflx', iou, 1, (/id_time/), -1e15           &    
     &        , 1.e15,' ', 'D', 'Atmosphere-Ice solar flux'             &    
     &        , '', 'W m^{-2}')     
!   Non Solar heat flux
          call defvar ('aicenflx', iou, 1, (/id_time/), -1e15           &    
     &        , 1.e15,' ', 'D', 'Atmosphere-Ice non solar flux'         &    
     &        , '', 'W m^{-2}')     
!   Solar heat flux ice-oce
          call defvar ('iicesflx', iou, 1, (/id_time/), -1e15           &    
     &        , 1.e15,' ', 'D', 'Ice-ocean solar flux'             &    
     &        , '', 'W m^{-2}')     
!   Non Solar heat flux ice-oce
          call defvar ('iicenflx', iou, 1, (/id_time/), -1e15           &    
     &        , 1.e15,' ', 'D', 'Ice-ocean non solar flux'         &    
     &        , '', 'W m^{-2}')     

          call enddef (iou)

      else
!      if the file does exist then open it for writing at the next record
       print*,"output file found...opening existing file for appending"
       call opennext ("nemo_ice_rtd.nc", tyear, ntrec, iou)
      endif

!       append variables
!--------------------------------------------------------------------------
      do cur_mon = 1, lm 

!       Adjust the time-step    
        norec = float(ntrec) + float(cur_mon) - 1.0_dp
        ntrec2 = ntrec + cur_mon - 1
      
!       Convert the date into days since 01-01-0001        
        CALL noleap_days(iyear, imon+cur_mon-1, 1, days_elapsed)
        tdays_elapsed = float(days_elapsed)  
!        print*, days_elapsed
!       time
        call putvars ('time', iou, ntrec2, tdays_elapsed, 1., 0.)

!       Area
        call putvars ('Area_NH', iou, ntrec2, area_nh(cur_mon), 1., 0.)
        call putvars ('Area_SH', iou, ntrec2, area_sh(cur_mon), 1., 0.)
!       Extent
        call putvars ('Extent_NH', iou, ntrec2, extent_nh(cur_mon), 1., 0.)
        call putvars ('Extent_SH', iou, ntrec2, extent_sh(cur_mon), 1., 0.)
!       Volume
        call putvars ('Volume_NH', iou, ntrec2, vol_nh(cur_mon), 1., 0.)
        call putvars ('Volume_SH', iou, ntrec2, vol_sh(cur_mon), 1., 0.)
        call putvars ('Volume_snow_NH', iou, ntrec2, vol_sno_nh(cur_mon), 1., 0.)
        call putvars ('Volume_snow_SH', iou, ntrec2, vol_sno_sh(cur_mon), 1., 0.)
!       Thickness
        call putvars ('Thickness_NH', iou, ntrec2, thick_nh(cur_mon), 1., 0.)
        call putvars ('Thickness_SH', iou, ntrec2, thick_sh(cur_mon), 1., 0.)
!       Snow Thickness
        call putvars ('Snow_thickness_NH', iou, ntrec2, isnowthick_nh(cur_mon), 1., 0.)
        call putvars ('Snow_thickness_SH', iou, ntrec2, isnowthick_sh(cur_mon), 1., 0.)
!       Heat fluxes
        call putvars ('aicenflx', iou, ntrec2, aicenflx_ave(cur_mon), 1., 0.)
        call putvars ('aicesflx', iou, ntrec2, aicesflx_ave(cur_mon), 1., 0.)
        call putvars ('iicenflx', iou, ntrec2, iicenflx_ave(cur_mon), 1., 0.)
        call putvars ('iicesflx', iou, ntrec2, iicesflx_ave(cur_mon), 1., 0.)
      enddo
        
      print*, 'closing netcdf'
      call closeall

END SUBROUTINE calc

