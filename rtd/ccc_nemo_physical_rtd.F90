PROGRAM nemo_ocean_diag
! ======================================================================
! PROGRAM nemo_ocean_diag
! ---------------------------------------------------
! JUN 24/14 - N. Swart    1. Make resolution independent (handles ORCA1,
!                            ORCA2, ORCA0.25)
!                         2. Calculate Indo. Pass. transport indirectly
!                            as the transport btwn Africa and Aus. 
!                         3. Updated, partially, to F90                           
!
! May 30/14 - N. Swart    1. Remove writes to CCCma format, add writes to 
!                            netCDF.
!                         2. Modify annual mean calculation of DY so that 
!                            it is applied to all variables (not just T,S).
!
! OCT 18/13 - D. YANG (1. Added IMPLICIT NONE and related declarations;
!                      2. Revised calculations for global means of T,S & SSH:
!                         a. accounted for volume from SSH; b. masked row 
!                         (:,jmt) from area calculation; c. added weighted 
!                         annual mean; d. generated new mesh_mask file not 
!                         including Caspian Sea. 
!                         Note that c and d have slight impact on T & S
!                         related fields.)
! May 22/13 - O. Saenko
! 
! AUTHOR  - O. Saenko
!
! PURPOSE - Run-time diagnostics for NEMO (ORCA1)
!
! INPUT FILE(S)...
!
! OUTPUT FILE...
! ======================================================================
! to compile locally:
! to compile:
!   $F77 -o nemo_ocean_diag nemo_ocean_diag.f uvic_netcdf.f $LINK -L/home/rls/wrk/AR5/CMOR/lib -lnetcdf -I/home/rls/wrk/AR5/CMOR/include/
!   
! UPDATE - 2013Feb20 - MB
!   gfortran -o nemo_ocean_diag nemo_ocean_diag.f uvic_netcdf.f -I/usr/local/netcdf-4.2.1_GF/include -L/usr/local/netcdf-4.2.1_GF/lib/ -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz
! OR
!   pgf90 -m64 -O -D_LARGE_FILES -DpgiFortran -o nemo_ocean_diag nemo_ocean_diag.f uvic_netcdf.f -I/usr/local/netcdf-4.2.1_PG/include -L/usr/local/netcdf-4.2.1_PG/lib/ -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz

! to compile on backend (spica/hadar)
!$F77 -o nemo_diag_global_mean nemo_diag_global_mean.f uvic_netcdf.f -I /users/tor/acrn/rls/local/aix64/netcdf3.6.3/include -L /users/tor/acrn/rls/local/aix64/netcdf3.6.3/lib -lnetcdf $LINK

! ======================================================================
      IMPLICIT NONE
      INTEGER :: imt, jmt, km, ll, iou

!         establish the size of the grid from the input file.
          call openfile ("grid_t", iou)
          call getdimlen ('x', iou, imt)
          call getdimlen ('y', iou, jmt)
          call getdimlen ('deptht', iou, km)

          ll = 1

!         do the calculations and save the output netcdf   
          call calc (imt, jmt, km, ll)

END PROGRAM nemo_ocean_diag

SUBROUTINE calc (imt, jmt, km, ll)
!     Does the required calculations and saves the output to netcdf
      IMPLICIT NONE
      INTEGER :: i, j, k, imt, jmt, km, ll, year, mon, nrecon
      INTEGER :: iou0, iou1, iou2, iou3, iou4, iou5, iou11
      INTEGER :: j_20N, j_20S, j_eq, k60, k500, k2000, i_DP, j_DP_S 
      INTEGER :: j_DP_N, i_IN_E1, i_IN_W1, i_IN_E2, i_IN_W2
      INTEGER :: i_AN_E, i_AN_W, i_AS_E, i_AS_W, i_PN_E, i_PN_W
      REAL    :: recn, cp, tz, sz, w_meanx, w_meany, area1, area2
      REAL    :: area, arc, arcn, arcs
! ======================================================================
!     Input data 
! ======================================================================
!     Grid-related arrays
      REAL, DIMENSION(imt, jmt)     :: lon2d, lat2d
      REAL, DIMENSION(imt, jmt)     :: e1t, e2t, e1u, e2u, e1v, e2v
      REAL, DIMENSION(imt, jmt, km) :: e3t, e3u, e3v
      REAL, DIMENSION(km)           :: depthw, deptht
!     Monthly T,S,u,v,w, eddy-induced u,v,w
      REAL, DIMENSION(imt, jmt, km) :: theta, salt, u, v, w
      REAL, DIMENSION(imt, jmt, km) :: gmu, gmv, gmw
!     Monthly fluxes of heat, water, and momentum, sea surface height, 
!     and mixed layer depth (MLD) 
      REAL, DIMENSION(imt, jmt)     :: hflux, wflux, tau_x, tau_y, ssh
      REAL, DIMENSION(imt, jmt)     :: mld10
! ======================================================================
!     Pre-computed data  
! ======================================================================
!     weight for annual mean calculations, DY, 2/OCT/2013
      INTEGER DIM(12)
      DATA DIM/31,28,31,30,31,30,31,31,30,31,30,31/
!     t,u,v masks 
      REAL, DIMENSION(imt, jmt, km) :: t_mask, u_mask, v_mask 
!     Annual T,S
      REAL, DIMENSION(imt, jmt, km) ::  theta_ann, salt_ann
!     Annual u, v, w
      REAL, DIMENSION(imt, jmt, km) ::  u_ann, v_ann, w_ann 
!     Annual eddy-induced u, v, w
      REAL, DIMENSION(imt, jmt, km) :: gmu_ann, gmv_ann, gmw_ann 
!     Annual fluxes of heat, water, and momentum
      REAL, DIMENSION(imt, jmt) ::  hflux_ann, wflux_ann
      REAL, DIMENSION(imt, jmt) ::  tau_x_ann,tau_y_ann
!     Annual sea surface height  
      REAL, DIMENSION(imt, jmt) ::  ssh_ann
!     Winter and summer MLD
      REAL, DIMENSION(imt, jmt) ::  mld10_win, mld10_sum
!     Wind energy input 
      REAL, DIMENSION(imt, jmt) ::  wind_x, wind_y
!     Global and annual mean meridional overturning circulation (MOC)
!     (Note: meaningful values are only south of 20N)
      REAL, DIMENSION(jmt, km) :: over_psi, over_psi_eddy
!     Masks for some regions 
      REAL, DIMENSION(imt, jmt) :: trop_up_mask ! for tropical Pacific upwelling 
      REAL, DIMENSION(imt, jmt) :: g_mask, g_mask1, g_mask2
      REAL, DIMENSION(imt, jmt) :: nino3_mask, nino34_mask, nino4_mask
! ======================================================================
!     Working arrays / variables  
! ======================================================================
      REAL, DIMENSION(imt, jmt) :: arr2d1, arr2d2, tarea, zarea_ssh
      REAL                      :: dum, dvol, volssh, volt
      REAL                      :: area_tot, vol0, zztmp, vol
! ======================================================================
!     Output data 
! ======================================================================
! (1) Global-mean profiles of T(z) and S(z) (C, g/kg)
      REAL, DIMENSION(km) :: theta_z, salt_z
!     Global-mean T and S  (C, g/kg)) 
      REAL                :: tvol, svol 
! (2) Global surface fields 
      REAL ::  hglo    ! heat flux (W/m2) 
      REAL ::  wglo    ! water flux (1.e+7 kg/m2/s) 
      REAL ::  sshglo  ! sea level (cm) 
! (3) Energetics
      REAL :: wind_work_glb ! net wind energy input to the ocean (TW)
      REAL :: wind_work_so  ! wind energy input south of 40S     (TW)
! (4) Tropical Pacific dynamics/therodynamics 
      REAL :: trp_up  ! Upwelling across 60m (Sv), 150E - 75W,  2S - 2N
      REAL :: t_nino3 ! Nino3   SST,  150W - 90W,  5S - 5N 
      REAL :: t_nino34! Nino3.4 SST,  170W - 120W, 5S - 5N 
      REAL :: t_nino4 ! Nino4   SST,  160E - 150W, 5S - 5N 
      REAL :: euc_max ! Max speed of EUC    (m/s)
! (5) Transports through key passages 
      REAL :: dp_tran  ! Drake Passage (Sv)
      REAL :: pi_tran  ! Indonesian Passage (Sv)
      REAL :: be_tran  ! Net transport across 20N in Atlantic
      REAL :: be_tran2 ! Net transport across 20N in Pacific
                    ! (can be used as proxies to Bering Strait tran.)                  
! (6) Mixed layer depth ( for values > 200m) 
      REAL :: win_mld      ! mean February MLD (m)     
      REAL :: sum_mld      ! mean August MLD (m)
      REAL :: win_mld_max  ! max. February MLD (m)  
      REAL :: sum_mld_max  ! max. August MLD (m)      
      REAL :: win_area     ! area of February MLD (1.e+14 m2)
      REAL :: sum_area     ! area of February MLD (1.e+14 m2)
! (7) Meridional overturning circulation (MOC)   
!     Maximum upper ocean MOC at 20N and 20S (Sv)
      REAL :: over_max_20N, over_max_20S
!     Minimum lower ocean MOC at 20N and 20S (Sv) 
      REAL :: over_min_20N, over_min_20S
!     Upper Southern Ocean MOC, net and eddy-induced (Sv, south of 40S) 
      REAL :: over_max_SO_net, over_min_SO_eddy  
! (8) Heat transport  (PW)   
!     across 20N, global ocean and Atlantic 
      REAL :: h_tran_20N, h_tran_20NA
!     across 20S, global ocean and Atlantic
      REAL :: h_tran_20S, h_tran_20SA
!----------------
!  NetCDF-output specific
      INTEGER :: id_time, id_z, iou, ntrec, iyear
      LOGICAL :: exists, exists1, notopen
      REAL    :: tyear
      CHARACTER :: fname01*100,fname02*100, fname03*100
      CHARACTER :: fname04*100, fname05*100
      CHARACTER(len=32) :: year_arg_in

!----------------
           print*,'Reading data on NEMO grid...'

         year =0
         ntrec=0
         iou0 =0
         iou1 =0
         iou2 =0
         iou3 =0
         iou4 =0
         recn =12.
         nrecon = int(recn + 0.001)
         Cp   = 4.2e+6  ! J/m3/K

! ======================================================================
!  Determine resolution and set parameters
! ======================================================================
!     ORCA2
      if ( imt == 182 ) then
        print *, "Using ORCA2 configuration"
        j_20N   =  92; j_20S   = 56; j_eq    = 74
        k60     =   7; k500    = 20; k2000   = 24
        i_DP    = 107; j_DP_S  = 17; j_DP_N  = 34
        i_IN_E1 =   1; i_IN_W1 = 23
        i_IN_E2 = 159; i_IN_W2 = imt-2
        i_AN_E  =  92; i_AN_W  = 135
        i_AS_E  = 121; i_AS_W  = 149
        i_PN_E  =  16; i_PN_W  = 91
!     ORCA1
      else if ( imt == 362 ) then
        print *, "Using ORCA1 configuration"
        j_20N   = 182; j_20S   = 112; j_eq    = 147
        k60     =   8; k500    =  20; k2000   =  29
        i_DP    = 221; j_DP_S  =  41; j_DP_N  =  66
        i_IN_E1 =   1; i_IN_W1 =  49
        i_IN_E2 = 322; i_IN_W2 = imt-2
        i_AN_E  = 191; i_AN_W  = 274
        i_AS_E  = 247; i_AS_W  = 302
        i_PN_E  =  34; i_PN_W  = 185
!     ORCA0.25
      else if ( imt == 1442 ) then
        print *, "Using ORCA0.25 configuration"
        j_20N   =  581; j_20S  =  419; j_eq    = 499
        k60     =    8; k500    =  20; k2000   =  29
        i_DP    =  880; j_DP_S  = 132; j_DP_N  = 238
        i_IN_E1 =    1; i_IN_W1 = 194
        i_IN_E2 = 1284; i_IN_W2 = imt-2
        i_AN_E  =  798; i_AN_W  = 1090
        i_AS_E  =  985; i_AS_W  = 1205
        i_PN_E  =  127; i_PN_W  =  735
      else
        print *, "Dont recognize the configuration."
        print *, "Only ORCA2, ORCA1 and ORCA0.25 compatible"
        stop
      endif

!---------------------------------------------------
!    Define NetCDF files   
!---------------------------------------------------
        fname01='grid_t'
        fname02='grid_u'
        fname03='grid_v'
        fname04='grid_w'
        fname05='orca_mesh_mask'

!---------------------------------------------------
!    Open the defined NetCDF files   
!---------------------------------------------------
! open files with ocean physical variables
      call openfile (fname01,iou0)
      call openfile (fname02,iou1)
      call openfile (fname03,iou2)
      call openfile (fname04,iou3)
! open file with mask/grid info
      call openfile (fname05,iou4)

!---------------------------------------------------
!    Get grid/mask data   
!---------------------------------------------------
      call getvara ('nav_lon', iou0, imt*jmt, (/1,1,1/), (/imt,jmt,1/), lon2d, 1., 0.)
      call getvara ('nav_lat', iou0, imt*jmt, (/1,1,1/), (/imt,jmt,1/), lat2d, 1., 0.)
      call getvara ('deptht', iou0, km, (/1/), (/km/), deptht, 1., 0.)
      call getvara ('depthw', iou3, km, (/1/), (/km/), depthw, 1., 0.)
      call getvara ('e1t', iou4, imt*jmt, (/1,1,1/), (/imt,jmt,1/),e1t , 1., 0.)
      call getvara ('e2t', iou4, imt*jmt, (/1,1,1/), (/imt,jmt,1/),e2t , 1., 0.)
      call getvara ('e3t', iou4, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),e3t , 1., 0.)
      call getvara ('e1v', iou4, imt*jmt, (/1,1,1/), (/imt,jmt,1/),e1v , 1., 0.)
      call getvara ('e2v', iou4, imt*jmt, (/1,1,1/), (/imt,jmt,1/),e2v , 1., 0.)
      call getvara ('e3v', iou4, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),e3v , 1., 0.)
      call getvara ('e1u', iou4, imt*jmt, (/1,1,1/), (/imt,jmt,1/),e1u , 1., 0.)
      call getvara ('e2u', iou4, imt*jmt, (/1,1,1/), (/imt,jmt,1/),e2u , 1., 0.)
      call getvara ('e3u', iou4, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),e3u , 1., 0.)
      call getvara ('tmask', iou4, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),t_mask , 1., 0.)
      call getvara ('umask', iou4, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),u_mask , 1., 0.)
      call getvara ('vmask', iou4, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),v_mask , 1., 0.)
!---------------------------------------------------
!    Construct masks for some sub-regions   
!    (will be done only once at some point...)
!---------------------------------------------------
      do i = 1, imt
          do j = 1, jmt
! Mask for tropical Pacific upwelling
              trop_up_mask(i,j) = 0. 
              if (lon2d(i,10).ge.150..or.lon2d(i,10).le.-75.) then 
                  if (lat2d(10,j).gt.-2..and.lat2d(10,j).lt.2.) then 
                      if (t_mask(i,j,8).eq.1.) then ! ~ 60m
                          trop_up_mask(i,j) = 1. 
                      endif
                  endif           
              endif   

! Masks for Nino
             nino3_mask(i,j) = 0. 
             nino34_mask(i,j)= 0. 
             nino4_mask(i,j) = 0.
             if (t_mask(i,j,1).eq.1.) then ! sst 
               if (lat2d(10,j).ge.-5..and.lat2d(10,j).le.5.) then 
! Nino3
                 if (lon2d(i,10).ge.-150..and.lon2d(i,10).le.-90.) then
                   nino3_mask(i,j) =  1. 
                 endif
! Nino3.4
                 if (lon2d(i,10).ge.-170..and.lon2d(i,10).le.-120.) then
                   nino34_mask(i,j) = 1. 
                 endif
! Nino4
                 if (lon2d(i,10).ge.160..or.lon2d(i,10).le.-150.) then 
                   nino4_mask(i,j) = 1. 
                 endif
               endif           
             endif  
          enddo
      enddo
      print*, SUM( nino3_mask(:,:) )
!---------------------------------------------------
!    Set to zero arrays for annual accumulation   
!---------------------------------------------------
      hflux_ann(:,:)    = 0.
      wflux_ann(:,:)    = 0.
      tau_x_ann(:,:)    = 0. 
      tau_y_ann(:,:)    = 0. 
      ssh_ann(:,:)      = 0.  
      mld10_win(:,:)    = 0. 
      mld10_sum(:,:)    = 0.
      wind_x(:,:)       = 0.
      wind_y(:,:)       = 0.          

      theta_ann(:,:,:)  = 0.
      salt_ann(:,:,:)   = 0.
      u_ann(:,:,:)      = 0.
      v_ann(:,:,:)      = 0.
      w_ann(:,:,:)      = 0.
      gmu_ann(:,:,:)    = 0.
      gmv_ann(:,:,:)    = 0.
      gmw_ann(:,:,:)    = 0.
!---------------------------------------------------
! reading time-dependent (monthly) data and compute 
! annual-mean (the latter does not take into consideration 
! the different number of days in each month; has to be 
! corrected at some point...)
!---------------------------------------------------
        
      do mon = 1, nrecon !  12  
! temperature
          call getvara ('votemper', iou0, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), theta, 1., 0.)
! salinity
          call getvara ('vosaline', iou0, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), salt, 1., 0.)
! net heat flux
          call getvara ('sohefldo', iou0, imt*jmt, (/1,1,mon/), (/imt,jmt,ll/), hflux, 1., 0.)
! net water flux
          call getvara ('sowaflup', iou0, imt*jmt, (/1,1,mon/), (/imt,jmt,ll/), wflux, 1., 0.)
! u-velocity 
          call getvara ('vozocrtx', iou1, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), u, 1., 0.)
! v-velocity 
          call getvara ('vomecrty', iou2, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), v, 1., 0.)
! w-velocity 
          call getvara ('vovecrtz', iou3, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), w, 1., 0.)
! EI u-velocity 
          call getvara ('vozoeivu', iou1, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), gmu, 1., 0.)
! EI v-velocity 
          call getvara ('vomeeivv', iou2, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), gmv, 1., 0.)
! EI w-velocity 
          call getvara ('voveeivw', iou3, imt*jmt*km, (/1,1,1,mon/), (/imt,jmt,km,ll/), gmw, 1., 0.)
! Wind Stress along i-axis
          call getvara ('sozotaux', iou1, imt*jmt, (/1,1,mon/), (/imt,jmt,ll/), tau_x, 1., 0.)
! Wind Stress along j-axis
          call getvara ('sometauy', iou2, imt*jmt, (/1,1,mon/), (/imt,jmt,ll/), tau_y, 1., 0.)
! Sea surface height
          call getvara ('sossheig', iou0, imt*jmt, (/1,1,mon/), (/imt,jmt,ll/), ssh, 1., 0.)
! Mixed Layer Depth 0.01 ref.10m
          call getvara ('somxl010', iou0, imt*jmt, (/1,1,mon/), (/imt,jmt,ll/), mld10, 1., 0.)

          print*,'month = ',mon ! ,theta(120,70,1)
! 2D data
          hflux_ann(:,:) = hflux_ann(:,:) + hflux(:,:)*DIM(mon)/365.
          wflux_ann(:,:) = wflux_ann(:,:) + wflux(:,:)*DIM(mon)/365.
          tau_x_ann(:,:) = tau_x_ann(:,:) + tau_x(:,:)*DIM(mon)/365.
          tau_y_ann(:,:) = tau_y_ann(:,:) + tau_y(:,:)*DIM(mon)/365.
          ssh_ann(:,:)   = ssh_ann(:,:)   + ssh(:,:)   * DIM(mon) / 365.
! MLD
          if (mon.eq.2) then 
              mld10_win(:,:)    =  mld10(:,:) 
          endif
 
          if (mon.eq.8) then 
              mld10_sum(:,:)    =  mld10(:,:) 
          endif 
! Wind enery input
          wind_x(:,:) = wind_x(:,:) + tau_x(:,:)*u(:,:,1)*DIM(mon)/365.
          wind_y(:,:) = wind_y(:,:) + tau_y(:,:)*v(:,:,1)*DIM(mon)/365.

! 3D data
         theta_ann(:,:,:) = theta_ann(:,:,:)+ theta(:,:,:)*DIM(mon)/365.
         salt_ann(:,:,:)  = salt_ann(:,:,:) + salt(:,:,:)*DIM(mon)/365.
         u_ann(:,:,:)     = u_ann(:,:,:)    + u(:,:,:)*DIM(mon)/365.
         v_ann(:,:,:)     = v_ann(:,:,:)    + v(:,:,:)*DIM(mon)/365.
         w_ann(:,:,:)     = w_ann(:,:,:)    + w(:,:,:)*DIM(mon)/365.
         gmu_ann(:,:,:)   = gmu_ann(:,:,:)  + gmu(:,:,:)*DIM(mon)/365.
         gmv_ann(:,:,:)   = gmv_ann(:,:,:)  + gmv(:,:,:)*DIM(mon)/365.
         gmw_ann(:,:,:)   = gmw_ann(:,:,:)  + gmw(:,:,:)*DIM(mon)/365.
      enddo  ! time  (month) 

! ********** Do some basic calculations ************

! start DY, 11/OCT/2013
!---------------------------------------------------
! (1) Global annual mean T and S and SSH
!---------------------------------------------------
      volssh   = 0.
      area_tot = 0.
! ---------------------------- total area    
      tarea(:, :)   = e1t(:, :)*e2t(:, :)*t_mask(:, :, 1)
      tarea(:, jmt) = 0. ! not to count the wrap row added to the northmost.
      tarea(1, :)   = 0. ! not to count the 2 wrap columns for the cyclic boundary
      tarea(imt, :) = 0. ! sshglo is not identical when using area(imt-1:imt,:)=0.
! ---------------------------- volume due to ssh   
      zarea_ssh(:, :) = tarea(:, :)*ssh_ann(:, :)
      volssh = SUM( zarea_ssh(:, :) )
      area_tot = SUM( tarea(:, :) )
! ---------------------------- Global mean ssh (cm)
      sshglo = volssh/area_tot
      sshglo = sshglo*1.e+2 ! cm
! ---------------------------- Global volume (not count ssh)
      vol0 = 0.
      do k = 1, km
          vol0 = vol0 + SUM( tarea(:, :)*t_mask(:, :, k)*e3t(:, :, k) )
      enddo
! ---------------------------- Total volume
      vol = vol0 + volssh
! ---------------------------- Global mean temperature (C) & salinity (psu)
      tvol = 0.
      svol = 0.
      tvol = SUM( zarea_ssh(:, :)*theta_ann(:, :, 1) )
      svol=SUM( zarea_ssh(:, :)*salt_ann(:, :, 1) )
      do k = 1, km
          do j = 1, jmt
              do i = 1, imt
                  zztmp = tarea(i,j)*e3t(i,j,k)
                  tvol = tvol + zztmp*theta_ann(i, j, k)
                  svol = svol + zztmp*salt_ann(i, j, k)
              enddo
          enddo
      enddo
      if (vol.ne.0.) then
        tvol = tvol / vol ! C 
        svol = svol / vol ! psu 
      endif
! end DY, 11/OCT/2013

!---------------------------------------------------
! (1) Global annual mean T(z) and S(z)
!---------------------------------------------------
      do k = 1, km
          call area_ave (e1t, e2t, e3t, t_mask(:, :, k), theta_ann(:, :, k), imt, jmt, km, tz, dvol, k)
          call area_ave (e1t, e2t, e3t, t_mask(:, :, k), salt_ann(:, :, k), imt, jmt, km, sz, dvol, k)
          theta_z(k) = tz
          salt_z(k) = sz
      enddo

!---------------------------------------------------
! (2) Global surface fields (fluxes, SSH, etc...)  
!---------------------------------------------------
      g_mask(:,:) = t_mask(:,:,1)
      call area_ave_flx (e1t, e2t, g_mask, hflux_ann, imt, jmt, hglo, dum)
      call area_ave_flx (e1t, e2t, g_mask, wflux_ann, imt, jmt, wglo, dum)
!DY      call area_ave_flx (e1t,e2t,g_mask,ssh_ann,sshglo,dum)
      wglo   = wglo*1.e+7   !  1.e-7 kg/m2/s
!DY      sshglo = sshglo*1.e+2 ! cm         

!---------------------------------------------------
! (3) Energetics
!---------------------------------------------------
!  Global Ocean 
      g_mask1(:, :) = u_mask(:, :, 1)
      g_mask2(:, :) = v_mask(:, :, 1)

      call area_ave_flx (e1u, e2u, g_mask1, wind_x, imt, jmt, w_meanx, area1)
      call area_ave_flx (e1v, e2v, g_mask2, wind_y, imt, jmt, w_meany, area2)
      wind_work_glb = (w_meanx*area1 + w_meany*area2) *1.e-12 ! TW 

!   Southern Ocean         
      do i = 1, imt
          do j = 1, jmt
              g_mask1(i,j) = 0.
              g_mask2(i,j) = 0.
              if (lat2d(10,j).le.-40.) then ! south of 40S 
                  g_mask1(i,j) = u_mask(i,j,1)
                  g_mask2(i,j) = v_mask(i,j,1)
              endif
          enddo
      enddo

      call area_ave_flx (e1u, e2u, g_mask1, wind_x, imt, jmt, w_meanx, area1)
      call area_ave_flx (e1v, e2v, g_mask2, wind_y, imt, jmt, w_meany, area2)
      wind_work_so = (w_meanx*area1 + w_meany*area2) *1.e-12 ! TW 
!---------------------------------------------------
! (4) Tropical Pacific dynamics/therodynamics 
!---------------------------------------------------
      arr2d1(:, :)  = theta_ann(:, :, 1)
      arr2d2(:, :)  = w_ann(:, :, k60)
      call area_ave_flx (e1t, e2t, trop_up_mask, arr2d2, imt, jmt, trp_up, area)
      trp_up = trp_up*area*1.e-6 ! Sv
 
      call area_ave (e1t, e2t, e3t, nino3_mask, arr2d1, imt, jmt, km,t_nino3, dvol, 1)  
      call area_ave (e1t, e2t, e3t, nino34_mask, arr2d1, imt, jmt, km, t_nino34, dvol, 1) 
      call area_ave (e1t, e2t, e3t, nino4_mask, arr2d1, imt, jmt, km, t_nino4, dvol, 1) 

      euc_max = 0. 
      do i = 1, imt-2
           if (lon2d(i,10).ge.150..or.lon2d(i,10).le.-75.) then
               do k = 1, k500 
                   if (u_ann(i, j_eq, k).gt. euc_max) then
                       euc_max = u_ann(i, j_eq, k)
                   endif
               enddo 
           endif
      enddo            
!---------------------------------------------------
! (5) Transports through key passages 
!---------------------------------------------------
! Drake Passage
      dp_tran = 0.
      do j = j_DP_S, j_DP_N
          do k=1,km
              arc = e2u(i_DP, j)*e3u(i_DP, j, k)*u_mask(i_DP, j, k)
              dp_tran = dp_tran + u_ann(i_DP, j, k)*arc
          enddo
      enddo
      dp_tran = dp_tran * 1.e-6 ! Sv         
! Indonesian Passage : calculated as transport btwn Africa and Aus. at 20S.
      pi_tran  = 0. 
        
!       Do two loops, because we cross the edge of the domain
!       From index 1 to Australia         
      do i = i_IN_E1, i_IN_W1
          do k = 1, km
              arc = e1v(i, j_20S)*e3v(i, j_20S ,k )*v_mask(i, j_20S, k)
              pi_tran = pi_tran + v_ann(i, j_20S, k)*arc
          enddo
      enddo
!       From Africa to the edge of the domain 
!      (excluding the periodic boundary )         
      do i = i_IN_E2, i_IN_W2
          do k = 1, km
              arc = e1v(i, j_20S)*e3v(i, j_20S ,k )*v_mask(i, j_20S, k)
              pi_tran = pi_tran + v_ann(i, j_20S, k)*arc
          enddo
      enddo

      pi_tran = pi_tran * 1.e-6 ! Sv        

! Net barotropic transports across 20N in Atlantic and Pacific
!   (can be used as proxies to Baring Passage transport)
      be_tran  = 0. 
      be_tran2 = 0. 
      do k = 1, km
! Atlantic, 20N
          do i = i_AN_E, i_AN_W
              arc = e1v(i, j_20N)*e3v(i, j_20N, k)*v_mask(i, j_20N, k)
              be_tran = be_tran + v_ann(i, j_20N, k)*arc
         enddo
! Pacific , 20N
          do i = i_PN_E, i_PN_W
              arc = e1v(i, j_20N)*e3v(i, j_20N, k)*v_mask(i, j_20N, k)
              be_tran2 = be_tran2 + v_ann(i, j_20N, k)*arc
          enddo
      enddo
      be_tran  = be_tran  * 1.e-6 ! Sv       
      be_tran2 = be_tran2 * 1.e-6 ! Sv  

!---------------------------------------------------
! (6) Mixed layer depth  
!---------------------------------------------------
      do i =1, imt
          do j =1, jmt
              g_mask1(i,j) = t_mask(i,j,1)
              g_mask2(i,j) = t_mask(i,j,1)
              if (mld10_win(i,j).lt.200.) then  ! mask for NH MLD > 200m
                  g_mask1(i,j) = 0.
              endif  
              if (mld10_sum(i,j).lt.200.) then  ! mask for SH MLD > 200m
                  g_mask2(i,j) = 0.
              endif              
          enddo
      enddo

      win_mld_max = 0.
      sum_mld_max = 0.

      do i=1,imt
          do j=1,jmt
              if (t_mask(i,j,1).gt.0.5) then 
! Winter max MLD
                  if (mld10_win(i,j).gt.win_mld_max) then
                      win_mld_max=mld10_win(i,j)
                  endif
! Summer max MLD 
                  if (mld10_sum(i,j).gt.sum_mld_max) then
                      sum_mld_max=mld10_sum(i,j)
                  endif
              endif       
          enddo
      enddo

      call area_ave_flx (e1t, e2t, g_mask1, mld10_win, imt, jmt, win_mld, win_area)
      call area_ave_flx (e1t, e2t, g_mask2, mld10_sum, imt, jmt, sum_mld, sum_area)          
      win_area = win_area*1.e-14 ! 1.e+14 m2
      sum_area = sum_area*1.e-14 ! 1.e+14 m2

!---------------------------------------------------
! (7) Meridional overturning circulation (MOC) 
!---------------------------------------------------

! Constract net, or residual velocities
      v_ann(:, :, :) = v_ann(:, :, :) + gmv_ann(:, :, :) 
      w_ann(:, :, :) = w_ann(:, :, :) + gmw_ann(:, :, :) 

      call moc(e1v, e3v, v_ann, imt, jmt, km, over_psi)
      call moc(e1v, e3v, gmv_ann, imt, jmt, km, over_psi_eddy)
! NADW
      over_max_20N = 0. 
      over_max_20S = 0.
      do k = k500, km ! below ~ 500 m 
          if (over_psi(j_20N,k).gt.over_max_20N) then  ! 20 N 
              over_max_20N = over_psi(j_20N,k) 
          endif
!
          if (over_psi(j_20S,k).gt.over_max_20S) then  ! 20 S 
              over_max_20S = over_psi(j_20S,k) 
         endif
      enddo    
! AABW 
      over_min_20N = 0. 
      over_min_20S = 0.
      do k = k2000, km ! below ~ 2000 m 
          if (over_psi(j_20N,k).lt.over_min_20N) then  ! 20 N 
              over_min_20N = over_psi(j_20N,k) 
          endif

          if (over_psi(j_20S,k).lt.over_min_20S) then  ! 20 S 
              over_min_20S = over_psi(j_20S,k) 
          endif
      enddo         
! Upper Southern Ocean MOC
      over_max_SO_net  = 0. 
      over_min_SO_eddy = 0.
      do k = 1, km 
          do j =1, jmt   
              if (lat2d(10,j).le.-40.) then ! south of 40S 
                  if (over_psi(j,k).gt.over_max_SO_net) then        ! Net  
                      over_max_SO_net = over_psi(j,k) 
                  endif
!
                  if (over_psi_eddy(j,k).lt.over_min_SO_eddy) then  ! Eddy-induced  
                      over_min_SO_eddy = over_psi_eddy(j,k) 
                  endif
              endif
          enddo 
      enddo    
!---------------------------------------------------
! (8) Heat transport (PW) 
!---------------------------------------------------
      h_tran_20N  = 0. 
      h_tran_20S  = 0.
      h_tran_20NA = 0. 
      h_tran_20SA = 0.

      do k = 1, km
          do i = 1, imt-2
! Global ocean at 20N 
              arcn = e1v(i, j_20N)*e3v(i, j_20N, k)*v_mask(i, j_20N, k)
              arcn = arcn*theta_ann(i, j_20N, k)*t_mask(i, j_20N, k)
              h_tran_20N = h_tran_20N + v_ann(i,j_20N,k)*arcn
! Global ocean at 20S 
              arcs = e1v(i, j_20S)*e3v(i, j_20S, k)*v_mask(i, j_20S, k)
              arcs = arcs*theta_ann(i, j_20S, k)*t_mask(i, j_20S, k)
              h_tran_20S = h_tran_20S + v_ann(i, j_20S, k)*arcs
! Atlantic at 20N
              if (i.ge.i_AN_E.and.i.le.i_AN_W) then
                  h_tran_20NA = h_tran_20NA + v_ann(i, j_20N, k)*arcn
              endif
! Atlantic at 20S
              if (i.ge.i_AS_E.and.i.le.i_AS_W) then
                  h_tran_20SA = h_tran_20SA + v_ann(i, j_20S, k)*arcs
              endif
          enddo
      enddo
      h_tran_20N  = Cp*h_tran_20N  *1.e-15 ! PW       
      h_tran_20NA = Cp*h_tran_20NA *1.e-15 ! PW  
      h_tran_20S  = Cp*h_tran_20S  *1.e-15 ! PW       
      h_tran_20SA = Cp*h_tran_20SA *1.e-15 ! PW  

!---------------------------------------------------
!   Main outputs 
!--------------------------------------------------

      print*,'-------------------------------------'
      print*,'    Temperature (C)                  '
      print*,'-------------------------------------'      
      print*,'SST     ', theta_z(1)
      print*,'T(220m) ', theta_z(15)
      print*,'T(450m) ', theta_z(19)
      print*,'T(850m) ', theta_z(23)
      print*,'T(1300m)', theta_z(26)
      print*,'T(2500m)', theta_z(32)
      print*,'Global T', tvol 

      print*,'-------------------------------------'
      print*,'    Salinity (g/kg)                  '
      print*,'-------------------------------------'     
      print*,'SSS     ', salt_z(1)
      print*,'S(220m) ', salt_z(15)
      print*,'S(450m) ', salt_z(19)
      print*,'S(850m) ', salt_z(23)
      print*,'S(1300m)', salt_z(26)
      print*,'S(2500m)', salt_z(32)
      print*,'Global S', svol 
      print*,'-------------------------------------'
      print*,'  Surface fields (fluxes, SSH, etc)  '
      print*,'-------------------------------------' 
      print*,'Heat  (W/m2)           ', hglo      
      print*,'Water (kg/m2/s)*1.e-7  ', wglo
      print*,'Sea surface height (cm)', sshglo

      print*,'------------------------------------'
      print*,'    Upper MOC (Sv)                  '
      print*,'------------------------------------'      
      print*,'Max at 20N ', over_max_20N   
      print*,'Max at 20S ', over_max_20S   

      print*,'------------------------------------'
      print*,'    Lower MOC (Sv)                  '
      print*,'------------------------------------'      
      print*,'Min at 20N ', over_min_20N   
      print*,'Min at 20S ', over_min_20S   
      print*,'------------------------------------'
      print*,'    Upper Southern Ocean MOC (Sv)   '
      print*,'------------------------------------'      
      print*,'Max net    ', over_max_SO_net   
      print*,'Min eddy   ', over_min_SO_eddy  
      print*,'------------------------------------'
      print*,'    Tropical Pacific                '
      print*,'------------------------------------'   
      print*,'Upwelling across 60m (Sv)',  trp_up
      print*,'Nino3   SST           (C)',  t_nino3
      print*,'Nino3.4 SST           (C)',  t_nino34
      print*,'Nino4   SST           (C)',  t_nino4
      print*,'Max speed of EUC    (m/s)',  euc_max
      print*,'------------------------------------'
      print*,'    Transports through key passages '
      print*,'------------------------------------'   
      print*,'Drake Passage              (Sv)',dp_tran
      print*,'Indonesian Passage         (Sv)',pi_tran
      print*,'Net across 20N in Atlantic (Sv)',be_tran
      print*,'Net across 20N in Pacific  (Sv)',be_tran2
      print*,'------------------------------------'
      print*,'    Deep (>200m) mixed layer        '  
      print*,'------------------------------------'
      print*,' Winter MLD, mean       (m) ', win_mld
      print*,' Summer MLD, mean       (m)' , sum_mld
      print*,' Winter MLD area (m2)*1.e+14', win_area
      print*,' Summer MLD area (m2)*1.e+14', sum_area
      print*,' Winter MLD, max        (m)' , win_mld_max
      print*,' Summer MLD, max        (m)' , sum_mld_max
      print*,'------------------------------------'
      print*,'    Wind energy input (TW)          '  
      print*,'------------------------------------'
      print*,' Global Ocean ', wind_work_glb
      print*,' South of 40S ', wind_work_so
      print*,'------------------------------------'
      print*,'    Heat transport  (PW)            '  
      print*,'------------------------------------'
      print*,' Global Ocean   at 20N ', h_tran_20N 
      print*,' Atlantic Ocean at 20N ', h_tran_20NA
      print*,' Global Ocean   at 20S ', h_tran_20S 
      print*,' Atlantic Ocean at 20S ', h_tran_20SA 

          print*,''             
          print*,'Done'

      call closeall

!---------------------------------------------------------
!     RTD output: Time series information section
!---------------------------------------------------------
!     Read in the year which is the first command line argument
      CALL getarg(1, year_arg_in )
      read (year_arg_in,'(I10)') iyear
      
      tyear = float(iyear)
      print*," --- "
      print*," Tyear is:", tyear
      print*," --- "

      iou = 0
      id_time = 0
      id_z = 0

!     NETCDF output
!     If the output file does not exist, create it and define dims and vars
      inquire (file="nemo_physical_rtd.nc", exist=exists1)

      if (.not. exists1) then
          print*,"output file not found...creating a new file..."
          call flush(6)
          call opennew ("nemo_physical_rtd.nc", iou)
          ntrec = 1
          call redef (iou)

!         basic grid specification
          call defdim ('time', iou, 0, id_time)
          call defdim ('depth', iou, km, id_z)
          call defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'D'           &
     &        , 'time', 'time', 'common_year since 1-1-0 00:00:0.0')

          call putatttext (iou, 'time', 'calendar', 'noleap')

          call defvar ('depth', iou, 1, (/id_z/), 0., 0., 'Y', 'F'             &
     &       , 'depth of the t grid', 'depth', 'm')

!         Temp
          call defvar ('T', iou, 1, (/id_time/), -1.e4                         &
     &        , 1.e4,' ', 'F', 'Global mean ocean temperature'                 &
     &        , 'T', 'k')                                                      

          call defvar ('Tz', iou, 2, (/id_z, id_time/), -1.e4                  &
     &        , 1.e4,' ', 'F', 'Temperature by level'                          &
     &        , 'Tz', 'k')
!         Salt
          call defvar ('S', iou, 1, (/id_time/), -1.e4                         &
     &        , 1.e4,' ', 'F', 'Global mean ocean salinity'                    &
     &        , 'S', 'psu')

          call defvar ('Sz', iou, 2, (/id_z, id_time/), -1.e4                  &
     &        , 1.e4,' ', 'F', 'Salinity by level'                             &
     &        , 'Sz', 'psu')
!-----------------------------------------------------------------
!        
!         Drake Passage transport
          call defvar ('dp_tran', iou, 1, (/id_time/), -1.e4                   &
     &        , 1.e4,' ', 'F', 'Drake Passage transport'                       & 
     &        , 'dp_tran', 'Sv')
!         Indo Passage transport
          call defvar ('pi_tran', iou, 1, (/id_time/), -1.e4                   &
     &        , 1.e4,' ', 'F', 'Indonesian Passage transport'                  &
     &        , 'pi_tran', 'Sv')
!          Transport across 20N Atlantic
          call defvar ('atl_tran_20n', iou, 1, (/id_time/), -1.e4              &
     &        , 1.e4,' ', 'F', 'Atlantic transport across 20N'                 &
     &        , 'atl_tran_20n', 'Sv')
!          Transport across 20N Pacific
          call defvar ('pac_tran_20n', iou, 1, (/id_time/), -1.e4              &
     &        , 1.e4,' ', 'F', 'Pacific transport across 20N'                  &
     &        , 'pac_tran_20n', 'Sv')
!-----------------------------------------------------------------
!        
!          MOC UPPER MAX at 20N
          call defvar ('mocmax20n', iou, 1, (/id_time/), -1.e4                 &
     &        , 1.e4,' ', 'F', 'MOC Max. at 20N'                               &
     &        , 'mocmax20n', 'Sv')
!          MOC UPPER MAX at 20S
          call defvar ('mocmax20s', iou, 1, (/id_time/), -1.e4                 &
     &        , 1.e4,' ', 'F', 'MOC Max. at 20S'                               & 
     &        , 'mocmax20s', 'Sv')
!          MOC LOWER MIN at 20N
          call defvar ('mocmin20n', iou, 1, (/id_time/), -1.e4                 &
     &        , 1.e4,' ', 'F', 'MOC Min. at 20N'                               &
     &        , 'mocmin20n', 'Sv')
!          MOC LOWER MIN at 20S
          call defvar ('mocmin20s', iou, 1, (/id_time/), -1.e4                 &
     &        , 1.e4,' ', 'F', 'MOC Min. at 20S'                               &
     &        , 'mocmin20s', 'Sv')
!          MOC UPPER SOU OCN MAX NET (SV)
          call defvar ('moc_max_SO_net', iou, 1, (/id_time/), -1.e4            &
     &        , 1.e4,' ', 'F', 'S.O. MOC Max. NET'                             &  
     &        , 'moc_max_SO_net', 'Sv')
!          MOC UPPER SOU OCN MIN EDDY (SV)
          call defvar ('moc_min_SO_eddy', iou, 1, (/id_time/), -1.e4           &
     &        , 1.e4,' ', 'F', 'S.O. MOC Min. Eddy'                            &
     &        , 'moc_min_SO_eddy', 'Sv')
!-----------------------------------------------------------------
!        
!          HEAT TRANS GLOBAL AT 20N (PW)
          call defvar ('h_tran_20N', iou, 1, (/id_time/), -1.e4                &
     &        , 1.e4,' ', 'F', 'Global heat transport at 20N'                  & 
     &        , 'h_tran_20N', 'PW')
!          HEAT TRANS GLOBAL AT 20S (PW)
          call defvar ('h_tran_20S', iou, 1, (/id_time/), -1.e4                &
     &        , 1.e4,' ', 'F', 'Global heat transport at 20S'                  &
     &        , 'h_tran_20S', 'PW')
!          Atlantic HEAT TRANS AT 20N (PW)
          call defvar ('h_tran_20NA', iou, 1, (/id_time/), -1.e4               &
     &        , 1.e4,' ', 'F', 'Atlantic heat transport at 20N'                &  
     &        , 'h_tran_20NA', 'PW')
!          Atlantic HEAT TRANS AT 20S (PW)
          call defvar ('h_tran_20SA', iou, 1, (/id_time/), -1.e4               & 
     &        , 1.e4,' ', 'F', 'Atlantic heat transport at 20S'                &
     &        , 'h_tran_20SA', 'PW')
!-----------------------------------------------------------------
!        
!       MEAN HEAT FLUX SURFACE (W/M^2)
          call defvar ('hglo', iou, 1, (/id_time/), -1.e4                      &
     &        , 1.e4,' ', 'F', 'Global mean heat flux at the surface'          &
     &        , 'hglo', 'W m^-^2')
!        MEAN FRESHWATER FLUX (KG/M^2/S)
          call defvar ('wglo', iou, 1, (/id_time/), -1.e4                      &
     &        , 1.e4,' ', 'F'                                                  &
     &        , 'Global mean freshwater flux at the surface'                   & 
     &        , 'wglo', 'kg m ^-^2 s^-^1')
!        SEA SURFACE HEIGHT (M)  
          call defvar ('sshglo', iou, 1, (/id_time/), -1.e4                    & 
     &        , 1.e4,' ', 'F', 'Global mean sea surface height'                &
     &        , 'sshglo', 'm')
!-----------------------------------------------------------------
!        
!        WIND ENERGY INPUT GLOBAL (TW)
          call defvar ('wind_work_glb', iou, 1, (/id_time/), -1.e4             &
     &        , 1.e4,' ', 'F', 'Global wind energy input '                     &
     &        , 'wind_work_glb', 'TW')                             
!        WIND ENERGY INPUT SOU OF 40S (TW)
          call defvar ('wind_work_so', iou, 1, (/id_time/), -1.e4              & 
     &        , 1.e4,' ', 'F', 'S.O. wind energy input s/o 40S '               & 
     &        , 'wind_work_so', 'TW')
!-----------------------------------------------------------------
!        
!        MLD MEAN WINTER (M)
          call defvar ('feb_mld_mean', iou, 1, (/id_time/), -1.e4              &
     &        , 1.e4,' ', 'F', 'February mean mixed layer depth '              &
     &        , 'feb_mld_mean', 'm')
!        MLD MEAN SUMMER
          call defvar ('aug_mld_mean', iou, 1, (/id_time/), -1.e4              &
     &        , 1.e4,' ', 'F', 'August mean mixed layer depth '                &
     &        , 'aug_mld', 'm')
!        MLD MAX WINTER (M) 
          call defvar ('feb_mld_max', iou, 1, (/id_time/), -1.e4               &
     &        , 1.e4,' ', 'F', 'February max. mixed layer depth '              &  
     &        , 'feb_mld_max', 'm')
!        MLD MAX SUMMER (M) 
          call defvar ('aug_mld_max', iou, 1, (/id_time/), -1.e4               &
     &        , 1.e4,' ', 'F', 'August max. mixed layer depth '                &
     &        , 'aug_mld_max', 'm')
!        MLD MEAN WINTER AREA>200M (M^2)
          call defvar ('feb_mld_area', iou, 1, (/id_time/), -1.e4              &
     &        , 1.e4,' ', 'F'                                                  & 
     &         , 'February area of mixed layer depth > 200 m'                  &
     &        , 'feb_mld_area', '10^{14} m^2')                                 
!        MLD MEAN SUMMER AREA>200M (M^2)
          call defvar ('aug_mld_area', iou, 1, (/id_time/), -1.e4              & 
     &        , 1.e4,' ', 'F'                                                  &
     &        , 'August area of mixed layer depth > 200 m'                     &
     &        , 'aug_mld_area', '10^{14} m^2')
!-----------------------------------------------------------------
!        
!        NINO3 SST (C) 
          call defvar ('nino3', iou, 1, (/id_time/), -1.e4                     & 
     &        , 1.e4,' ', 'F', 'NINO3 SST'                                     & 
     &        , 'nino3', '^oC')
!        NINO3.4 SST (C) 
          call defvar ('nino34', iou, 1, (/id_time/), -1.e4                    &
     &        , 1.e4,' ', 'F', 'NINO3.4 SST'                                   & 
     &        , 'nino34', '^oC')
!        NINO4 SST (C) 
          call defvar ('nino4', iou, 1, (/id_time/), -1.e4                     & 
     &        , 1.e4,' ', 'F', 'NINO4 SST'                                     &
     &        , 'nino4', '^oC')
!        UPWELL ACROSS 60M TROP PAC (SV)
          call defvar ('trop_upwell', iou, 1, (/id_time/), -1.e4               & 
     &        , 1.e4,' ', 'F', 'Topical Pacific upwelling across 60m'          &
     &        , 'trop_upwell', 'Sv')
!         MAX SPEED OF EUC (M/S) 
          call defvar ('euc_max', iou, 1, (/id_time/), -1.e4                   &
     &        , 1.e4,' ', 'F', 'Max. speed of the EUC'                         &
     &        , 'euc_max', 'm s^-^1')
!-----------------------------------------------------------------

          call enddef (iou)
!         define the depth axis
          call putvara ('depth', iou, km, (/1/), (/km/)                        &
     &      , deptht, 1., 0.)

      else
!      if the file does exist then open it for writing at the next record
       print*,"output file found...opening existing file for appending"
       call flush(6)
       call opennext ("nemo_physical_rtd.nc", tyear, ntrec, iou)
      endif

!       append variables
!--------------------------------------------------------------------------
!       time
        call putvars ('time', iou, ntrec, tyear, 1., 0.)

!       Temp
        call putvars ('T', iou, ntrec, tvol, 1., 0.)
        call putvara ('Tz', iou, km, (/1, ntrec/), (/km, 1/), theta_z, 1., 0.)
!       Salt
        call putvars ('S', iou, ntrec, svol, 1., 0.)
        call putvara ('Sz', iou, km, (/1, ntrec/), (/km, 1/), salt_z, 1., 0.)

!       Drake Passage transport
        call putvars ('dp_tran', iou, ntrec, dp_tran, 1., 0.)
!       Indo Passage transport
        call putvars ('pi_tran', iou, ntrec, pi_tran, 1., 0.)
!       Transport across 20N Atlantic
        call putvars ('atl_tran_20n', iou, ntrec, be_tran, 1., 0.)
!       Transport across 20N Pacific
        call putvars ('pac_tran_20n', iou, ntrec, be_tran2, 1., 0.)
!--------------------------------------------------------------------------
!       MOC UPPER MAX at 20N
        call putvars ('mocmax20n', iou, ntrec, over_max_20N, 1., 0.)
!       MOC UPPER MAX at 20S
        call putvars ('mocmax20s', iou, ntrec, over_max_20S, 1., 0.)
!       MOC LOWER MIN at 20N
        call putvars ('mocmin20n', iou, ntrec, over_min_20N, 1., 0.)
!       MOC LOWER MIN at 20S
        call putvars ('mocmin20s', iou, ntrec, over_min_20S, 1., 0.)
!       MOC UPPER SOU OCN MAX NET (SV)
        call putvars ('moc_max_SO_net', iou, ntrec, over_max_SO_net, 1., 0.)
!       MOC UPPER SOU OCN MIN EDDY (SV)
        call putvars ('moc_min_SO_eddy', iou, ntrec, over_min_SO_eddy, 1., 0.)
!--------------------------------------------------------------------------
!       HEAT TRANS GLOBAL AT 20N (PW)
        call putvars ('h_tran_20N', iou, ntrec, h_tran_20N, 1., 0.)
!       HEAT TRANS GLOBAL AT 20S (PW)
        call putvars ('h_tran_20S', iou, ntrec, h_tran_20S, 1., 0.)
!       Atlantic HEAT TRANS AT 20N (PW)
        call putvars ('h_tran_20NA', iou, ntrec, h_tran_20NA, 1., 0.)
!       Atlantic HEAT TRANS AT 20S (PW)
        call putvars ('h_tran_20SA', iou, ntrec, h_tran_20SA, 1., 0.)
!--------------------------------------------------------------------------
!       MEAN HEAT FLUX SURFACE (W/M^2)
        call putvars ('hglo', iou, ntrec, hglo, 1., 0.)
!       MEAN FRESHWATER FLUX (KG/M^2/S)
        call putvars ('wglo', iou, ntrec, wglo, 1., 0.)
!       SEA SURFACE HEIGHT (M)
        call putvars ('sshglo', iou, ntrec, sshglo, 1., 0.)
!--------------------------------------------------------------------------
!        WIND ENERGY INPUT GLOBAL (TW)
        call putvars ('wind_work_glb', iou, ntrec, wind_work_glb, 1., 0.)
!        WIND ENERGY INPUT SOU OF 40S (TW)
        call putvars ('wind_work_so', iou, ntrec, wind_work_so, 1., 0.)
!--------------------------------------------------------------------------
!         MLD MEAN WINTER (M)
        call putvars ('feb_mld_mean', iou, ntrec, win_mld, 1., 0.)
!         MLD MEAN SUMMER
        call putvars ('aug_mld_mean', iou, ntrec, sum_mld, 1., 0.)
!         MLD MAX WINTER (M)
        call putvars ('feb_mld_max', iou, ntrec, win_mld_max, 1., 0.)
!         MLD MAX SUMMER (M)
        call putvars ('aug_mld_max', iou, ntrec, sum_mld_max, 1., 0.)
!         MLD MEAN WINTER AREA>200M (M^2)
        call putvars ('feb_mld_area', iou, ntrec, win_area, 1., 0.)
!         MLD MEAN WINTER AREA>200M (M^2)
        call putvars ('aug_mld_area', iou, ntrec, sum_area, 1., 0.)
!--------------------------------------------------------------------------
!        NINO3 SST (C) 
        call putvars ('nino3', iou, ntrec, t_nino3, 1., 0.)
!        NINO3.4 SST (C) 
        call putvars ('nino34', iou, ntrec, t_nino34, 1., 0.)
!        NINO4 SST (C) 
        call putvars ('nino4', iou, ntrec, t_nino4, 1., 0.)
!        UPWELL ACROSS 60M TROP PAC (SV)
        call putvars ('trop_upwell', iou, ntrec, trp_up, 1., 0.)
!        MAX SPEED OF EUC (M/S) 
        call putvars ('euc_max', iou, ntrec, euc_max, 1., 0.)
!--------------------------------------------------------------------------
        print*, 'closing netcdf'
        call flush(6)
        call closefile (iou)
END SUBROUTINE  calc

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

!=========================================================
! Global meridional overturning (Sv) (valid only south of 20N)  
!=========================================================
SUBROUTINE moc(e1v, e3v, v, imt, jmt, km, over_psi)
      implicit none
      integer imt, jmt, km, i, j, k
      REAL, DIMENSION(imt, jmt) :: e1v 
      REAL, DIMENSION(imt, jmt, km) :: e3v, v
      REAL, DIMENSION(jmt, km) :: over_tran, over_psi
      REAL s


      do j = 1, jmt
          do k = km, 1, -1
              s=0.
              do i = 1, imt - 2  ! not to double count the cyclic boundary
                  s = s + v(i, j, k)*e1v(i, j)*e3v(i, j, k)   
              enddo
              over_tran(j, k) = s 
              over_psi(j, k)  = 0.
          enddo     
      enddo

      do j = 1, jmt
          do k = km, 1, -1
              if (k.eq.km) then 
                  over_psi(j,k) = -over_tran(j,k)
              else
                  over_psi(j,k) = over_psi(j, k + 1) - over_tran(j,k)
              endif
          enddo
      enddo

      do j = 1, jmt
          do k = 1, km
              over_psi(j, k)  =  over_psi(j, k)*1.e-6 ! to Sv                      
          enddo
      enddo

      return
END SUBROUTINE moc

