PROGRAM nemo_ocean_diag
! ======================================================================
! PROGRAM nemo_ocean_diag
! ---------------------------------------------------
!
! HISTORY
!--------
!
! Nov 5/15  - N. Swart    1. Remove annual mean calculation and rewrite
!                            all code to operate on monthly data.
!                         2. Remove North Fold point from computions.
!                            (i.e. sum to jmt -1)
!                         3. Rewrite of functions and code to freefrom. 
!                         4. Improve time axis in output NetCDF, include
!                            noleap_days subroutine.
!
! JUN 24/14 - N. Swart    1. Make resolution independent (handles ORCA1,
!                            ORCA2, ORCA0.25)
!                         2. Calculate Indo. Pass. transport indirectly
!                            as the transport btwn Africa and Aus. 
!                         3. Transition to freeform F90 (incomplete).                          
!
! May 30/14 - N. Swart    1. Major cleanup and reorganization of code.
!                         2. Remove writes to CCCma format, and replaces 
!                            with writes to netCDF.
!                         3. Modify annual mean calculation of DY so that 
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
!
! May 22/13 - O. Saenko   Original code.
! 
!
! PURPOSE - Run-time diagnostics for NEMO (ORCA1)
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
!    - grid_t : monthly frequency (_1m_)
!    - grid_u : monthly frequency (_1m_)
!    - grid_v : monthly frequency (_1m_)
!    - grid_w : monthly frequency (_1m_)
!
! OUTPUT FILES
! ------------
! nemo_physical_rtd.nc - NetCDF output file with monthly timeseries for physics variables.
!
! ======================================================================
! to compile:
!
! 1. Use build-nemo-rtd
!
! 2. xlf90_r -o nemo_physical_rtd.exe \
!   ccc_nemo_physical_rtd.F90 uvic_netcdf.f `nf-config --fflags --flibs`
! ======================================================================
      IMPLICIT NONE
      INTEGER :: imt, jmt, km, lm, iou

!         establish the size of the grid from the input file.
          call openfile ("grid_t", iou)
          call getdimlen ('x', iou, imt)
          call getdimlen ('y', iou, jmt)
          call getdimlen ('deptht', iou, km)
          call getdimlen ('time_counter', iou, lm)

!         do the calculations and save the output netcdf   
          call calc (imt, jmt, km, lm)

END PROGRAM nemo_ocean_diag

SUBROUTINE calc (imt, jmt, km, lm)
!     Does the required calculations and saves the output to netcdf
      IMPLICIT NONE
      integer, parameter:: dp=kind(0.d0) ! double precision
      INTEGER :: i, j, k, l, imt, jmt, km, lm, year, mon, nrecon
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
      REAL, DIMENSION(imt, jmt)         :: lon2d, lat2d
      REAL, DIMENSION(imt, jmt)         :: e1t, e2t, e1u, e2u, e1v, e2v
      REAL, DIMENSION(imt, jmt, km)     :: e3t, e3u, e3v
      REAL, DIMENSION(km)               :: depthw, deptht
!     Monthly T,S,u,v,w, eddy-induced u,v,w
      REAL, DIMENSION(imt, jmt, km, lm) :: theta, salt, u, v, w
      REAL, DIMENSION(imt, jmt, km, lm) :: gmu, gmv, gmw
!     Monthly fluxes of heat, water, and momentum, sea surface height, 
!     and mixed layer depth (MLD) 
      REAL, DIMENSION(imt, jmt, lm)     :: hflux, wflux, tau_x, tau_y
      REAL, DIMENSION(imt, jmt, lm)     :: mld10, ssh
! ======================================================================
!     Pre-computed data  
! ======================================================================
!     t,u,v masks 
      REAL, DIMENSION(imt, jmt, km)     :: t_mask, u_mask, v_mask 
!     Winter and summer MLD
      REAL, DIMENSION(imt, jmt, lm) ::  mld10_win, mld10_sum
!     Wind energy input 
      REAL, DIMENSION(imt, jmt, lm) ::  wind_x, wind_y
!     Global meridional overturning circulation (MOC)
!     (Note: meaningful values are only south of 20N)
      REAL, DIMENSION(jmt, km, lm) :: over_psi, over_psi_eddy
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
      REAL, DIMENSION(km, lm) :: theta_z, salt_z
!     Global-mean T and S  (C, g/kg)) 
      REAL, DIMENSION(lm)     :: tvol, svol 
! (2) Global surface fields 
      REAL, DIMENSION(lm) :: hglo    ! heat flux (W/m2) 
      REAL, DIMENSION(lm) :: wglo    ! water flux (1.e+7 kg/m2/s) 
      REAL, DIMENSION(lm) :: sshglo  ! sea level (cm) 
! (3) Energetics
      REAL, DIMENSION(lm) :: wind_work_glb ! net wind energy input to the ocean (TW)
      REAL, DIMENSION(lm) :: wind_work_so  ! wind energy input south of 40S     (TW)
! (4) Tropical Pacific dynamics/therodynamics 
      REAL, DIMENSION(lm) :: trp_up  ! Upwelling across 60m (Sv), 150E - 75W,  2S - 2N
      REAL, DIMENSION(lm) :: t_nino3 ! Nino3   SST,  150W - 90W,  5S - 5N 
      REAL, DIMENSION(lm) :: t_nino34! Nino3.4 SST,  170W - 120W, 5S - 5N 
      REAL, DIMENSION(lm) :: t_nino4 ! Nino4   SST,  160E - 150W, 5S - 5N 
      REAL, DIMENSION(lm) :: euc_max ! Max speed of EUC    (m/s)
! (5) Transports through key passages 
      REAL, DIMENSION(lm) :: dp_tran  ! Drake Passage (Sv)
      REAL, DIMENSION(lm) :: pi_tran  ! Indonesian Passage (Sv)
      REAL, DIMENSION(lm) :: be_tran  ! Net transport across 20N in Atlantic
      REAL, DIMENSION(lm) :: be_tran2 ! Net transport across 20N in Pacific
                    ! (can be used as proxies to Bering Strait tran.)                  
! (6) Mixed layer depth ( for values > 200m) 
      REAL, DIMENSION(lm) :: win_mld      ! mean February MLD (m)     
      REAL, DIMENSION(lm) :: sum_mld      ! mean August MLD (m)
      REAL, DIMENSION(lm) :: win_mld_max  ! max. February MLD (m)  
      REAL, DIMENSION(lm) :: sum_mld_max  ! max. August MLD (m)      
      REAL, DIMENSION(lm) :: win_area     ! area of February MLD (1.e+14 m2)
      REAL, DIMENSION(lm) :: sum_area     ! area of February MLD (1.e+14 m2)
! (7) Meridional overturning circulation (MOC)   
!     Maximum upper ocean MOC at 20N and 20S (Sv)
      REAL, DIMENSION(lm) :: over_max_20N, over_max_20S
!     Minimum lower ocean MOC at 20N and 20S (Sv) 
      REAL, DIMENSION(lm) :: over_min_20N, over_min_20S
!     Upper Southern Ocean MOC, net and eddy-induced (Sv, south of 40S) 
      REAL, DIMENSION(lm) :: over_max_SO_net, over_min_SO_eddy  
! (8) Heat transport  (PW)   
!     across 20N, global ocean and Atlantic 
      REAL, DIMENSION(lm) :: h_tran_20N, h_tran_20NA
!     across 20S, global ocean and Atlantic
      REAL, DIMENSION(lm) :: h_tran_20S, h_tran_20SA
!----------------
!  NetCDF-output specific
      INTEGER :: id_time, id_z, iou, ntrec, ntrec2, iyear, imon
      INTEGER :: days_elapsed
      LOGICAL :: exists, exists1, notopen
      REAL    :: tyear, tdays_elapsed
      CHARACTER :: fname01*100,fname02*100, fname03*100
      CHARACTER :: fname04*100, fname05*100
      CHARACTER(len=32) :: year_arg_in, mon_arg_in

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

!---------------------------------------------------
! Read in the monthly data from NetCDF
!---------------------------------------------------
        
! temperature
      CALL getvara ('votemper', iou0, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), theta, 1., 0.)
! salinity
      CALL getvara ('vosaline', iou0, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), salt, 1., 0.)
! net heat flux
      CALL getvara ('sohefldo', iou0, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/), hflux, 1., 0.)
! net water flux
      CALL getvara ('sowaflup', iou0, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/), wflux, 1., 0.)
! u-velocity 
      CALL getvara ('vozocrtx', iou1, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), u, 1., 0.)
! v-velocity 
      CALL getvara ('vomecrty', iou2, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), v, 1., 0.)
! w-velocity 
      CALL getvara ('vovecrtz', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), w, 1., 0.)
! EI u-velocity 
      CALL getvara ('vozoeivu', iou1, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), gmu, 1., 0.)
! EI v-velocity 
      CALL getvara ('vomeeivv', iou2, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), gmv, 1., 0.)
! EI w-velocity 
      CALL getvara ('voveeivw', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), gmw, 1., 0.)
! Wind Stress along i-axis
      CALL getvara ('sozotaux', iou1, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/), tau_x, 1., 0.)
! Wind Stress along j-axis
      CALL getvara ('sometauy', iou2, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/), tau_y, 1., 0.)
! Sea surface height
      CALL getvara ('sossheig', iou0, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/), ssh, 1., 0.)
! Mixed Layer Depth 0.01 ref.10m
       CALL getvara ('somxl010', iou0, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/), mld10, 1., 0.)

! ********** Do calculations ************
      tvol(:)          = 0.0_dp
      svol(:)          = 0.0_dp
      euc_max(:)       = 0.0_dp 
      dp_tran(:)       = 0.0_dp
      pi_tran(:)       = 0.0_dp 
      be_tran(:)       = 0.0_dp
      be_tran2(:)      = 0.0_dp 
      win_mld(:)       = 0.0_dp
      sum_mld(:)       = 0.0_dp
      win_mld_max(:)   = 0.0_dp
      sum_mld_max(:)   = 0.0_dp
      win_area(:)      = 0.0_dp
      sum_area(:)      = 0.0_dp
      over_max_20N(:)  = 0.0_dp 
      over_max_20S(:)  = 0.0_dp
      over_min_20N     = 0.0_dp 
      over_min_20S     = 0.0_dp
      over_max_SO_net  = 0.0_dp 
      over_min_SO_eddy = 0.0_dp
      h_tran_20N       = 0.0_dp
      h_tran_20S       = 0.0_dp
      h_tran_20NA      = 0.0_dp
      h_tran_20SA      = 0.0_dp

    ! ---------------------------- total area    
      tarea(:, :)   = e1t(:, :)*e2t(:, :)*t_mask(:, :, 1)
      tarea(:, jmt) = 0. ! not to count the wrap row added to the northmost.
      tarea(1, :)   = 0. ! not to count the 2 wrap columns for the cyclic boundary
      tarea(imt, :) = 0. ! sshglo is not identical when using area(imt-1:imt,:)=0.
      area_tot = 0.
      area_tot = SUM(tarea(:, :))
    ! ---------------------------- Global volume (not counting ssh)
          vol0 = 0.
          do k = 1, km
              vol0 = vol0 + SUM( tarea(:, :)*t_mask(:, :, k)*e3t(:, :, k) )
          enddo
! Wind enery input
          wind_x =  tau_x(:,:,:)*u(:,:,1,:)
          wind_y =  tau_y(:,:,:)*v(:,:,1,:)
! MLD
! Seasonal separation doesnt make sense with monthly RTDs, but define both (as the same),
! in order to maintain backwards compatibility in the plotting routine.
          mld10_win(:,:,:)    =  mld10(:,:,:)
          mld10_sum(:,:,:)    =  mld10(:,:,:)

    ! Main compute loop over all months
      do l = 1, lm 

    ! start DY, 11/OCT/2013
    !---------------------------------------------------
    ! (1) Global annual mean T and S and SSH
    !---------------------------------------------------
          volssh   = 0.
    ! ---------------------------- volume due to ssh   
          zarea_ssh(:, :) = tarea(:, :)*ssh(:, :, l)
          volssh = SUM( zarea_ssh(:, :) )
    ! ---------------------------- Global mean ssh (cm)
          sshglo(l) = (volssh/area_tot)*1.0e2
    ! ---------------------------- Total volume
          vol = vol0 + volssh
    ! ---------------------------- Global mean temperature (C) & salinity (psu)
          tvol(l) = SUM(zarea_ssh(:, :)*theta(:, :, 1, l))
          svol(l) = SUM(zarea_ssh(:, :)* salt(:, :, 1, l))
          do k = 1, km
              do j = 1, jmt
                  do i = 1, imt
                      zztmp = tarea(i,j)*e3t(i,j,k)
                      tvol(l) = tvol(l) + zztmp*theta(i, j, k, l)
                      svol(l) = svol(l) + zztmp*salt(i, j, k, l)
                  enddo
              enddo
          enddo
          if (vol.ne.0.) then
            tvol(l) = tvol(l) / vol ! C 
            svol(l) = svol(l) / vol ! psu 
          endif
    ! end DY, 11/OCT/2013

    !---------------------------------------------------
    ! (1) Global mean T(z) and S(z)
    !---------------------------------------------------
          do k = 1, km
              call area_ave(e1t, e2t, e3t, t_mask(:, :, k)             &
                &           , theta(:, :, k, l), imt, jmt              &
                &           , km, theta_z(k, l), dvol, k)
              call area_ave(e1t, e2t, e3t, t_mask(:, :, k)             &
                &           , salt(:, :, k, l), imt, jmt               &
                &           , km, salt_z(k, l), dvol, k)
          enddo

    !---------------------------------------------------
    ! (2) Global surface fields (fluxes, SSH, etc...)  
    !---------------------------------------------------
          g_mask(:,:) = t_mask(:,:,1)
          call area_ave_flx(e1t, e2t, g_mask, hflux(:, :, l), imt      &
            &               , jmt, hglo(l), dum)
          call area_ave_flx(e1t, e2t, g_mask, wflux(:, :, l), imt      &
            &              , jmt, wglo(l), dum)
    !DY      call area_ave_flx (e1t,e2t,g_mask,ssh_ann,sshglo,dum)
          wglo(l)   = wglo(l)*1.e+7   !  1.e-7 kg/m2/s
    !DY      sshglo = sshglo*1.e+2 ! cm         

    !---------------------------------------------------
    ! (3) Energetics
    !---------------------------------------------------
    !  Global Ocean 
          g_mask1(:, :) = u_mask(:, :, 1)
          g_mask2(:, :) = v_mask(:, :, 1)

          call area_ave_flx(e1u, e2u, g_mask1, wind_x(:, :, l), imt, jmt &
            &               , w_meanx, area1)
          call area_ave_flx(e1v, e2v, g_mask2, wind_y(:, :, l), imt      &
            &               , jmt, w_meany, area2)
          wind_work_glb(l) = (w_meanx*area1 + w_meany*area2) *1.e-12 ! TW 

    !   Southern Ocean         
          do i = 1, imt
              do j = 1, jmt
                  g_mask1(i,j) = 0.
                  g_mask2(i,j) = 0.
                  if (lat2d(i, j).le.-40.) then ! south of 40S 
                      g_mask1(i, j) = u_mask(i, j, 1)
                      g_mask2(i, j) = v_mask(i, j, 1)
                  endif
              enddo
          enddo

          call area_ave_flx(e1u, e2u, g_mask1, wind_x(:,:,l), imt      & 
            &               , jmt, w_meanx, area1)
          call area_ave_flx(e1v, e2v, g_mask2, wind_y(:,:,l), imt      &
            &               , jmt, w_meany, area2)
          wind_work_so(l) = (w_meanx*area1 + w_meany*area2) *1.e-12 ! TW 

    !---------------------------------------------------
    ! (4) Tropical Pacific dynamics/therodynamics 
    !---------------------------------------------------
          arr2d1(:, :)  = theta(:, :, 1, l)
          arr2d2(:, :)  = w(:, :, k60, l)
          call area_ave_flx(e1t, e2t, trop_up_mask, arr2d2, imt, jmt   &
            &               , trp_up(l), area)
          trp_up(l) = trp_up(l)*area*1.e-6 ! Sv
     
          call area_ave(e1t, e2t, e3t, nino3_mask, arr2d1, imt, jmt    &
            &           , km,t_nino3(l), dvol, 1)  
          call area_ave(e1t, e2t, e3t, nino34_mask, arr2d1, imt, jmt   &
            &           , km, t_nino34(l), dvol, 1) 
          call area_ave(e1t, e2t, e3t, nino4_mask, arr2d1, imt, jmt    &
            &           , km, t_nino4(l), dvol, 1) 

          do i = 1, imt-2
               if (lon2d(i,10).ge.150..or.lon2d(i,10).le.-75.) then
                   do k = 1, k500 
                       if (u(i, j_eq, k, l).gt. euc_max(l)) then
                           euc_max(l) = u(i, j_eq, k, l)
                       endif
                   enddo 
               endif
          enddo            
    !---------------------------------------------------
    ! (5) Transports through key passages 
    !---------------------------------------------------
    ! Drake Passage
          do j = j_DP_S, j_DP_N
              do k = 1, km
                  arc = e2u(i_DP, j)*e3u(i_DP, j, k)*u_mask(i_DP, j, k)
                  dp_tran(l) = dp_tran(l) + u(i_DP, j, k, l)*arc
              enddo
          enddo
          dp_tran(l) = dp_tran(l) * 1.e-6 ! Sv         

    ! Indonesian Passage : calculated as transport btwn Africa and Aus. at 20S.
    !       Do two loops, because we cross the edge of the domain
    !       From index 1 to Australia         
          do i = i_IN_E1, i_IN_W1
              do k = 1, km
                  arc = e1v(i, j_20S)*e3v(i, j_20S ,k )*v_mask(i, j_20S, k)
                  pi_tran(l) = pi_tran(l) + v(i, j_20S, k, l)*arc
              enddo
          enddo
    !       From Africa to the edge of the domain 
    !      (excluding the periodic boundary )         
          do i = i_IN_E2, i_IN_W2
              do k = 1, km
                  arc = e1v(i, j_20S)*e3v(i, j_20S ,k )*v_mask(i, j_20S, k)
                  pi_tran(l) = pi_tran(l) + v(i, j_20S, k, l)*arc
              enddo
          enddo

          pi_tran(l) = pi_tran(l) * 1.e-6 ! Sv        

    ! Net barotropic transports across 20N in Atlantic and Pacific
    !   (can be used as proxies to Baring Passage transport)
          do k = 1, km
    ! Atlantic, 20N
              do i = i_AN_E, i_AN_W
                  arc = e1v(i, j_20N)*e3v(i, j_20N, k)*v_mask(i, j_20N, k)
                  be_tran(l) = be_tran(l) + v(i, j_20N, k, l)*arc
             enddo
    ! Pacific , 20N
              do i = i_PN_E, i_PN_W
                  arc = e1v(i, j_20N)*e3v(i, j_20N, k)*v_mask(i, j_20N, k)
                  be_tran2(l) = be_tran2(l) + v(i, j_20N, k, l)*arc
              enddo
          enddo
          be_tran(l)  = be_tran(l)  * 1.e-6 ! Sv       
          be_tran2(l) = be_tran2(l) * 1.e-6 ! Sv  

    !---------------------------------------------------
    ! (6) Mixed layer depth  
    !---------------------------------------------------
          do i =1, imt
              do j =1, jmt
                  g_mask1(i,j) = t_mask(i,j,1)
                  g_mask2(i,j) = t_mask(i,j,1)
                  if (mld10_win(i,j,l).lt.200.) then  ! mask for NH MLD > 200m
                      g_mask1(i,j) = 0.
                  endif  
                  if (mld10_sum(i,j,l).lt.200.) then  ! mask for SH MLD > 200m
                      g_mask2(i,j) = 0.
                  endif              
              enddo
          enddo

          do i=1,imt
              do j=1,jmt
                  if (t_mask(i,j,1).gt.0.5) then 
    ! Winter max MLD
                      if (mld10_win(i, j, l).gt.win_mld_max(l)) then
                          win_mld_max(l) = mld10_win(i, j, l)
                      endif
    ! Summer max MLD 
                      if (mld10_sum(i, j, l).gt.sum_mld_max(l)) then
                          sum_mld_max(l) = mld10_sum(i, j, l)
                      endif
                  endif       
              enddo
          enddo
          call area_ave_flx(e1t, e2t, g_mask1, mld10(:, :, l), imt &
            &               , jmt, win_mld(l), win_area(l))
          call area_ave_flx(e1t, e2t, g_mask2, mld10(:, :, l), imt &
            &               , jmt, sum_mld(l), sum_area(l))          
          win_area(l) = win_area(l)*1.e-14 ! 1.e+14 m2
          sum_area(l) = sum_area(l)*1.e-14 ! 1.e+14 m2

    !---------------------------------------------------
    ! (7) Meridional overturning circulation (MOC) 
    !---------------------------------------------------

    ! Constract net, or residual velocities
          v(:, :, :, l) = v(:, :, :, l) + gmv(:, :, :, l) 
          w(:, :, :, l) = w(:, :, :, l) + gmw(:, :, :, l) 

          call moc(e1v, e3v, v(:, :, :, l), imt, jmt, km               & 
            &     , over_psi(:, :, l))
          call moc(e1v, e3v, gmv(:,:,:,l), imt, jmt, km                &
            &     , over_psi_eddy(:, :, l))
    ! NADW
          do k = k500, km ! below ~ 500 m 
              if (over_psi(j_20N, k, l).gt.over_max_20N(l)) then   
                  over_max_20N(l) = over_psi(j_20N, k, l) 
              endif
    !
              if (over_psi(j_20S, k, l).gt.over_max_20S(l)) then   
                  over_max_20S(l) = over_psi(j_20S, k, l) 
             endif
          enddo    
    ! AABW 
          do k = k2000, km ! below ~ 2000 m 
              if (over_psi(j_20N, k, l).lt.over_min_20N(l)) then   
                  over_min_20N(l) = over_psi(j_20N, k, l) 
              endif

              if (over_psi(j_20S, k, l).lt.over_min_20S(l)) then  
                  over_min_20S(l) = over_psi(j_20S, k, l) 
              endif
          enddo         
    ! Upper Southern Ocean MOC
          do k = 1, km 
              do j =1, jmt   
                  if (lat2d(10,j).le.-40.) then ! south of 40S 
                      if (over_psi(j, k, l).gt.over_max_SO_net(l)) then      
                          over_max_SO_net(l) = over_psi(j, k, l) 
                      endif
    !
                      if (over_psi_eddy(j, k, l).lt.over_min_SO_eddy(l)) then    
                          over_min_SO_eddy(l) = over_psi_eddy(j, k, l) 
                      endif
                  endif
              enddo 
          enddo    
    !---------------------------------------------------
    ! (8) Heat transport (PW) 
    !---------------------------------------------------
          do k = 1, km
              do i = 1, imt-2
    ! Global ocean at 20N 
                  arcn = e1v(i, j_20N)*e3v(i, j_20N, k)*v_mask(i, j_20N, k)
                  arcn = arcn*theta(i, j_20N, k, l)*t_mask(i, j_20N, k)
                  h_tran_20N(l) = h_tran_20N(l) + v(i, j_20N, k, l)*arcn
    ! Global ocean at 20S 
                  arcs = e1v(i, j_20S)*e3v(i, j_20S, k)*v_mask(i, j_20S, k)
                  arcs = arcs*theta(i, j_20S, k, l)*t_mask(i, j_20S, k)
                  h_tran_20S(l) = h_tran_20S(l) + v(i, j_20S, k, l)*arcs
    ! Atlantic at 20N
                  if (i.ge.i_AN_E.and.i.le.i_AN_W) then
                      h_tran_20NA(l) = h_tran_20NA(l) + v(i, j_20N, k, l)*arcn
                  endif
    ! Atlantic at 20S
                  if (i.ge.i_AS_E.and.i.le.i_AS_W) then
                      h_tran_20SA(l) = h_tran_20SA(l) + v(i, j_20S, k, l)*arcs
                  endif
              enddo
          enddo 
          h_tran_20N(l)  = Cp*h_tran_20N(l)  *1.e-15 ! PW       
          h_tran_20NA(l) = Cp*h_tran_20NA(l) *1.e-15 ! PW  
          h_tran_20S(l)  = Cp*h_tran_20S(l)  *1.e-15 ! PW       
          h_tran_20SA(l) = Cp*h_tran_20SA(l) *1.e-15 ! PW  

      enddo ! Main computing loop
!---------------------------------------------------
!   Main outputs 
!--------------------------------------------------

!      print*,'-------------------------------------'
!      print*,'    Temperature (C)                  '
!      print*,'-------------------------------------'      
!      print*,'SST     ', theta_z(1)
!      print*,'T(220m) ', theta_z(15)
!      print*,'T(450m) ', theta_z(19)
!      print*,'T(850m) ', theta_z(23)
!      print*,'T(1300m)', theta_z(26)
!      print*,'T(2500m)', theta_z(32)
!      print*,'Global T', tvol 
!
!      print*,'-------------------------------------'
!      print*,'    Salinity (g/kg)                  '
!      print*,'-------------------------------------'     
!      print*,'SSS     ', salt_z(1)
!      print*,'S(220m) ', salt_z(15)
!      print*,'S(450m) ', salt_z(19)
!      print*,'S(850m) ', salt_z(23)
!      print*,'S(1300m)', salt_z(26)
!      print*,'S(2500m)', salt_z(32)
!      print*,'Global S', svol 
!      print*,'-------------------------------------'
!      print*,'  Surface fields (fluxes, SSH, etc)  '
!      print*,'-------------------------------------' 
!      print*,'Heat  (W/m2)           ', hglo      
!      print*,'Water (kg/m2/s)*1.e-7  ', wglo
!      print*,'Sea surface height (cm)', sshglo
!
!      print*,'------------------------------------'
!      print*,'    Upper MOC (Sv)                  '
!      print*,'------------------------------------'      
!      print*,'Max at 20N ', over_max_20N   
!      print*,'Max at 20S ', over_max_20S   
!
!      print*,'------------------------------------'
!      print*,'    Lower MOC (Sv)                  '
!      print*,'------------------------------------'      
!      print*,'Min at 20N ', over_min_20N   
!      print*,'Min at 20S ', over_min_20S   
!      print*,'------------------------------------'
!      print*,'    Upper Southern Ocean MOC (Sv)   '
!      print*,'------------------------------------'      
!      print*,'Max net    ', over_max_SO_net   
!      print*,'Min eddy   ', over_min_SO_eddy  
!      print*,'------------------------------------'
!      print*,'    Tropical Pacific                '
!      print*,'------------------------------------'   
!      print*,'Upwelling across 60m (Sv)',  trp_up
!      print*,'Nino3   SST           (C)',  t_nino3
!      print*,'Nino3.4 SST           (C)',  t_nino34
!      print*,'Nino4   SST           (C)',  t_nino4
!      print*,'Max speed of EUC    (m/s)',  euc_max
!      print*,'------------------------------------'
!      print*,'    Transports through key passages '
!      print*,'------------------------------------'   
!      print*,'Drake Passage              (Sv)',dp_tran
!      print*,'Indonesian Passage         (Sv)',pi_tran
!      print*,'Net across 20N in Atlantic (Sv)',be_tran
!      print*,'Net across 20N in Pacific  (Sv)',be_tran2
!      print*,'------------------------------------'
!      print*,'    Deep (>200m) mixed layer        '  
!      print*,'------------------------------------'
!      print*,' Winter MLD, mean       (m) ', win_mld
!      print*,' Summer MLD, mean       (m)' , sum_mld
!      print*,' Winter MLD area (m2)*1.e+14', win_area
!      print*,' Summer MLD area (m2)*1.e+14', sum_area
!      print*,' Winter MLD, max        (m)' , win_mld_max
!      print*,' Summer MLD, max        (m)' , sum_mld_max
!      print*,'------------------------------------'
!      print*,'    Wind energy input (TW)          '  
!      print*,'------------------------------------'
!      print*,' Global Ocean ', wind_work_glb
!      print*,' South of 40S ', wind_work_so
!      print*,'------------------------------------'
!      print*,'    Heat transport  (PW)            '  
!      print*,'------------------------------------'
!      print*,' Global Ocean   at 20N ', h_tran_20N 
!      print*,' Atlantic Ocean at 20N ', h_tran_20NA
!      print*,' Global Ocean   at 20S ', h_tran_20S 
!      print*,' Atlantic Ocean at 20S ', h_tran_20SA 
!
!          print*,''             
!          print*,'Done'

      call closeall

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
      id_z = 0

!     NETCDF output
!     If the output file does not exist, create it and define dims and vars
      inquire (file="nemo_physical_rtd.nc", exist=exists1)

      if (.not. exists1) then
          print*,"output file not found...creating a new file..."
          !call flush(6)
          call opennew ("nemo_physical_rtd.nc", iou)
          ntrec = 1
          call redef (iou)

!         basic grid specification
          call defdim ('time', iou, 0, id_time)
          call defdim ('depth', iou, km, id_z)
          call defvar ('time', iou, 1, (/id_time/), 0., 0., 'T', 'F'   &
     &        , 'time', 'time', 'days since 0000-01-01 00:00:00')

          call putatttext (iou, 'time', 'calendar', '365_day')

          call defvar ('depth', iou, 1, (/id_z/), 0., 0., 'Y', 'F'     &
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
       !call flush(6)
       call opennext ("nemo_physical_rtd.nc", tyear, ntrec, iou)
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

!       Temp
        call putvars ('T', iou, ntrec2, tvol(l), 1., 0.)
        call putvara ('Tz', iou, km, (/1, ntrec2/), (/km, 1/), theta_z(:, l), 1., 0.)
!       Salt
        call putvars ('S', iou, ntrec2, svol, 1., 0.)
        call putvara ('Sz', iou, km, (/1, ntrec2/), (/km, 1/), salt_z(:, l), 1., 0.)

!       Drake Passage transport
        call putvars ('dp_tran', iou, ntrec2, dp_tran(l), 1., 0.)
!       Indo Passage transport
        call putvars ('pi_tran', iou, ntrec2, pi_tran(l), 1., 0.)
!       Transport across 20N Atlantic
        call putvars ('atl_tran_20n', iou, ntrec2, be_tran(l), 1., 0.)
!       Transport across 20N Pacific
        call putvars ('pac_tran_20n', iou, ntrec2, be_tran2(l), 1., 0.)
!--------------------------------------------------------------------------
!       MOC UPPER MAX at 20N
        call putvars ('mocmax20n', iou, ntrec2, over_max_20N(l), 1., 0.)
!       MOC UPPER MAX at 20S
        call putvars ('mocmax20s', iou, ntrec2, over_max_20S(l), 1., 0.)
!       MOC LOWER MIN at 20N
        call putvars ('mocmin20n', iou, ntrec2, over_min_20N(l), 1., 0.)
!       MOC LOWER MIN at 20S
        call putvars ('mocmin20s', iou, ntrec2, over_min_20S(l), 1., 0.)
!       MOC UPPER SOU OCN MAX NET (SV)
        call putvars ('moc_max_SO_net', iou, ntrec2, over_max_SO_net(l), 1., 0.)
!       MOC UPPER SOU OCN MIN EDDY (SV)
        call putvars ('moc_min_SO_eddy', iou, ntrec2, over_min_SO_eddy(l), 1., 0.)
!--------------------------------------------------------------------------
!       HEAT TRANS GLOBAL AT 20N (PW)
        call putvars ('h_tran_20N', iou, ntrec2, h_tran_20N(l), 1., 0.)
!       HEAT TRANS GLOBAL AT 20S (PW)
        call putvars ('h_tran_20S', iou, ntrec2, h_tran_20S(l), 1., 0.)
!       Atlantic HEAT TRANS AT 20N (PW)
        call putvars ('h_tran_20NA', iou, ntrec2, h_tran_20NA(l), 1., 0.)
!       Atlantic HEAT TRANS AT 20S (PW)
        call putvars ('h_tran_20SA', iou, ntrec2, h_tran_20SA(l), 1., 0.)
!--------------------------------------------------------------------------
!       MEAN HEAT FLUX SURFACE (W/M^2)
        call putvars ('hglo', iou, ntrec2, hglo(l), 1., 0.)
!       MEAN FRESHWATER FLUX (KG/M^2/S)
        call putvars ('wglo', iou, ntrec2, wglo(l), 1., 0.)
!       SEA SURFACE HEIGHT (M)
        call putvars ('sshglo', iou, ntrec2, sshglo(l), 1., 0.)
!--------------------------------------------------------------------------
!        WIND ENERGY INPUT GLOBAL (TW)
        call putvars ('wind_work_glb', iou, ntrec2, wind_work_glb(l), 1., 0.)
!        WIND ENERGY INPUT SOU OF 40S (TW)
        call putvars ('wind_work_so', iou, ntrec2, wind_work_so(l), 1., 0.)
!--------------------------------------------------------------------------
!         MLD MEAN WINTER (M)
        call putvars ('feb_mld_mean', iou, ntrec2, win_mld(l), 1., 0.)
!         MLD MEAN SUMMER
        call putvars ('aug_mld_mean', iou, ntrec2, sum_mld(l), 1., 0.)
!         MLD MAX WINTER (M)
        call putvars ('feb_mld_max', iou, ntrec2, win_mld_max(l), 1., 0.)
!         MLD MAX SUMMER (M)
        call putvars ('aug_mld_max', iou, ntrec2, sum_mld_max(l), 1., 0.)
!         MLD MEAN WINTER AREA>200M (M^2)
        call putvars ('feb_mld_area', iou, ntrec2, win_area(l), 1., 0.)
!         MLD MEAN WINTER AREA>200M (M^2)
        call putvars ('aug_mld_area', iou, ntrec2, sum_area(l), 1., 0.)
!--------------------------------------------------------------------------
!        NINO3 SST (C) 
        call putvars ('nino3', iou, ntrec2, t_nino3(l), 1., 0.)
!        NINO3.4 SST (C) 
        call putvars ('nino34', iou, ntrec2, t_nino34(l), 1., 0.)
!        NINO4 SST (C) 
        call putvars ('nino4', iou, ntrec2, t_nino4(l), 1., 0.)
!        UPWELL ACROSS 60M TROP PAC (SV)
        call putvars ('trop_upwell', iou, ntrec2, trp_up(l), 1., 0.)
!        MAX SPEED OF EUC (M/S) 
        call putvars ('euc_max', iou, ntrec2, euc_max(l), 1., 0.)
!--------------------------------------------------------------------------
      enddo ! netcdf writing time loop
        print*, 'closing netcdf'
        !call flush(6)
        call closefile (iou)
END SUBROUTINE  calc

!=========================================================
! Area averaging of 3d field over the selected regions 
!=========================================================
SUBROUTINE area_ave(e1,e2,e3, mask, a, imt,jmt,km,a_mean,ss,kk)
      implicit none
      integer imt, jmt, km, i, j, kk
      real e1(imt,jmt),e2(imt,jmt), e3(imt,jmt,km)
      real a(imt,jmt), mask(imt,jmt)
      real a_mean, ss, s1, vol

          s1=0.
          ss=0.
          do i=1,imt-2  ! not to double count the cyclic boundary
              do j=1,jmt-1 !avoid north fold
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
              do j=1,jmt-1 ! north fold
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


      do j = 1, jmt-1 ! north fold
          do k = km, 1, -1
              s=0.
              do i = 1, imt - 2  ! not to double count the cyclic boundary
                  s = s + v(i, j, k)*e1v(i, j)*e3v(i, j, k)   
              enddo
              over_tran(j, k) = s 
              over_psi(j, k)  = 0.
          enddo     
      enddo

      do j = 1, jmt-1 !north fold
          do k = km, 1, -1
              if (k.eq.km) then 
                  over_psi(j,k) = -over_tran(j,k)
              else
                  over_psi(j,k) = over_psi(j, k + 1) - over_tran(j,k)
              endif
          enddo
      enddo

      do j = 1, jmt-1 !north fold
          do k = 1, km
              over_psi(j, k)  =  over_psi(j, k)*1.e-6 ! to Sv                      
          enddo
      enddo

      return
END SUBROUTINE moc

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

