PROGRAM nemo_diag
   
   !!===============================================================
   !!                 ***    PROGRAM nemo_diag    ***              
   !!                  CMIP6 nemo offline diagnostics
   !!===============================================================
   !! 2018-04 (D. Yang): Original code
   !! 2020-02 (D. Yang): Revise tbnds to bnds in defdim to avoid crash
   !!                    when later "cdo mergetime" in canesm_nemo_diag.sh
   !!                    makes the unwanted change from tbnds to bnds.
   !! 2020-03 (D. Yang): Add getdimnm
   !!---------------------------------------------------------------
   !!
   !!---------------------------------------------------------------
   !! INPUT FIELDS
   !! 2D: e2u, e1v   
   !! 3D: e3u, e3v, e3t, umask, vmask, tmask
   !! 2D: ssh
   !! 3D: tnp, tnc, snp, snc 
   !! 4D: u, v
   !! 
   !! OUTPUT FIELDS
   !! 2D: mfo
   !! 2D: msftbarot
   !! 3D: opottemptend & osalttend
   !!
   !! INPUT FILES
   !! grid_u
   !! grid_v
   !! grid_t
   !! tnc, snc, tnp, snp
   !! orca_mesh_mask
   !!
   !! OUTPUT FILES
   !! mfo.nc
   !! msftbarot.nc
   !! tstend.nc
   !!---------------------------------------------------------------
   USE nemo_diag_glovars     ! global variable declarations
   USE nemo_diag_cal         ! diagnostics calculations
 
   IMPLICIT NONE

   CHARACTER(len=100) :: fname01, fname02, fname03, fname04
   CHARACTER(len=100) :: fname05, fname06, fname07, fname08
   CHARACTER(len=100) :: axis, standard_name, units, calendar, title
   CHARACTER(len=100) :: long_name, time_origin, bounds
   CHARACTER(len=100) :: bndsnm
   INTEGER   :: iou, iou1, iou2, iou3, iou4, iou5, iou6, iou7, iou8
   INTEGER   :: ntrec, id_time, id_tbnds, id_l, id_s, id_x, id_y, id_z
   INTEGER   :: ntbnds
   INTEGER   :: i, l
   LOGICAL   :: exists
   !INTEGER   :: strlen
   INTEGER, DIMENSION(10)            :: ierr
   REAL, DIMENSION(:), ALLOCATABLE   :: time
   REAL, DIMENSION(:), ALLOCATABLE   :: ytime
   REAL, DIMENSION(:), ALLOCATABLE   :: deptht
   REAL, DIMENSION(nline)            :: line
   REAL, DIMENSION(:), ALLOCATABLE   :: x, y
   REAL, DIMENSION(:,:), ALLOCATABLE :: nav_lon_u, nav_lat_u
   REAL, DIMENSION(:,:), ALLOCATABLE :: nav_lon_t, nav_lat_t
   REAL, DIMENSION(:,:), ALLOCATABLE :: time_bnds
   !!-------------------------------------
   !! Establish grid size from input files.
   !!-------------------------------------
   CALL openfile  ("grid_t", iou)
   CALL getdimlen ('x', iou, imt)
   CALL getdimlen ('y', iou, jmt)
   call getdimlen ('deptht', iou, km)
   CALL getdimlen ('time_counter', iou, lm)      
   ! dimid = 1: x; 2: y; 3: deptht; 4: time_counter; 5: tbnds or bnds
   CALL getdimnm  (bndsnm, iou, 5, ntbnds)

   ly = lm / 12

   !!----------------
   !! Allocate Arrays
   !!----------------
   ALLOCATE( e2u(imt,jmt), e1v(imt,jmt), STAT=ierr(1) )
   ALLOCATE( e3u(imt,jmt,km), e3v(imt,jmt,km), e3t(imt,jmt,km),    &
             umask(imt,jmt,km), vmask(imt,jmt,km), tmask(imt,jmt,km), STAT=ierr(2) )
   ALLOCATE( time(lm), ytime(ly), deptht(km), x(imt), y(jmt), STAT=ierr(3) )
   ALLOCATE( nav_lon_u(imt,jmt), nav_lat_u(imt,jmt), STAT=ierr(4) )
   ALLOCATE( nav_lon_t(imt,jmt), nav_lat_t(imt,jmt), STAT=ierr(5) )
   ALLOCATE( ssh(imt,jmt,lm), STAT=ierr(6) )
   ALLOCATE( u(imt,jmt,km,lm), v(imt,jmt,km,lm), STAT=ierr(7) )
   ALLOCATE( tnc(imt,jmt,km), tnp(imt,jmt,km), STAT=ierr(8) )
   ALLOCATE( snc(imt,jmt,km), snp(imt,jmt,km), STAT=ierr(9) )
   ALLOCATE( time_bnds(ntbnds,lm), STAT=ierr(10) )
      
   IF (MAXVAL(ierr) /=0) THEN
      STOP 'Memory allocation error in cmip6_nemo_offl'
   ENDIF
      
   iou1 =0
   iou2 =0
   iou3 =0
   iou4 =0
   iou5 =0
   iou6 =0
   iou7 =0
   iou8 =0

   !!--------------------
   !! Define NetCDF files   
   !!--------------------
   print*,'Reading data on NEMO grid...'
   fname01='orca_mesh_mask'
   fname02='grid_u'
   fname03='grid_v'
   fname04='grid_t'
   fname05='tnp.nc'
   fname06='snp.nc'
   fname07='tnc.nc'
   fname08='snc.nc'
   
   !!------------------------------
   !! Open the defined NetCDF files   
   !!------------------------------
   ! mask/grid info
   CALL openfile (fname01,iou1)
   ! ocean physical variables
   CALL openfile (fname02,iou2)
   CALL openfile (fname03,iou3)
   CALL openfile (fname04,iou4)
   CALL openfile (fname05,iou5)
   CALL openfile (fname06,iou6)
   CALL openfile (fname07,iou7)
   CALL openfile (fname08,iou8)

   !!-------------------
   !! Get grid/mask data   
   !!-------------------
   CALL getvara ('e2u', iou1, imt*jmt, (/1,1,1/), (/imt,jmt,1/),e2u , 1., 0.)
   CALL getvara ('e1v', iou1, imt*jmt, (/1,1,1/), (/imt,jmt,1/),e1v , 1., 0.)
   CALL getvara ('e3u', iou1, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),e3u , 1., 0.)
   CALL getvara ('e3v', iou1, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),e3v , 1., 0.)
   CALL getvara ('e3t', iou1, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),e3t , 1., 0.)
   CALL getvara ('umask', iou1, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),umask , 1., 0.)
   CALL getvara ('vmask', iou1, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),vmask , 1., 0.)
   CALL getvara ('tmask', iou1, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),tmask , 1., 0.)
   ! WRITE(*,*) 'umask'
   ! umask(144,45,1) = 1
   ! WRITE(*,*) umask(144,45,1)
   ! WRITE(*,*) 'e3u'
   ! e3u(144,45,1)=6.19056967037
   ! WRITE(*,*) e3u(144,45,1)

   !!-------------------------------------
   !! Read in the monthly data from NetCDF
   !!-------------------------------------
   ! time_counter
   CALL getvara ('time_counter', iou2, lm, (/1/), (/lm/), time, 1., 0.)
   !WRITE(*,*) 'time'
   !WRITE(*,*) time
   ! time_counter attribute
   CALL getatttext (iou2, 'time_counter', 'axis', axis)
   CALL getatttext (iou2, 'time_counter', 'standard_name', standard_name)
   CALL getatttext (iou2, 'time_counter', 'units', units)
   CALL getatttext (iou2, 'time_counter', 'calendar', calendar)
   CALL getatttext (iou2, 'time_counter', 'title', title)
   CALL getatttext (iou2, 'time_counter', 'long_name', long_name)
   CALL getatttext (iou2, 'time_counter', 'time_origin', time_origin)
   CALL getatttext (iou2, 'time_counter', 'bounds', bounds)
   ! time_counter_bnds
   CALL getvara ('time_counter_bnds', iou2, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
   !WRITE(*,*) 'time_bnds'
   !WRITE(*,*) time_bnds
   ! nav_lon on grid_U
   CALL getvara ('nav_lon', iou2, imt*jmt, (/1,1/), (/imt,jmt/), nav_lon_u, 1., 0.)
   ! nav_lat on grid_U
   CALL getvara ('nav_lat', iou2, imt*jmt, (/1,1/), (/imt,jmt/), nav_lat_u, 1., 0.)
   ! u-velocity
   CALL getvara ('vozocrtx', iou2, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), u, 1., 0.)
   ! WRITE(*,*) 'vozocrtx'
   ! u(144,45,1,1)=0.0234638 m/s
   ! WRITE(*,*) u(144,45,1,1)
   ! v-velocity 
   CALL getvara ('vomecrty', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), v, 1., 0.)
   ! WRITE(*,*) 'vomecrty'
   ! v(144,45,1,1)=-0.0440831 m/s     
   ! WRITE(*,*) v(144,45,1,1) 
   ! nav_lon on grid_T
   CALL getvara ('nav_lon', iou4, imt*jmt, (/1,1/), (/imt,jmt/), nav_lon_t, 1., 0.)
   !WRITE(*,*) 'nav_lon_t'
   !WRITE(*,*) nav_lon_t(144,45)
   ! nav_lat on grid_T
   CALL getvara ('nav_lat', iou4, imt*jmt, (/1,1/), (/imt,jmt/), nav_lat_t, 1., 0.)
   ! deptht 
   CALL getvara ('deptht', iou4, km, (/1/), (/km/), deptht, 1., 0.)
   ! ssh
   CALL getvara ('sossheig', iou4, imt*jmt*lm, (/1,1,1/), (/imt,jmt,lm/), ssh, 1., 0.)
   ! WRITE(*,*) 'sossheig'
   ! ssh(144,45,1,1)=-1.77589 m
   ! WRITE(*,*) ssh(144,45,1)
   ! tn from last time step of previous year
   CALL getvara ('tn', iou5, imt*jmt*km*1, (/1,1,1,1/), (/imt,jmt,km,1/), tnp, 1., 0.)
   !WRITE(*,*) 'tnp'
   ! tnp(144,45,1)=-0.0207554732464
   !WRITE(*,*) tnp(144,45,1)
   ! sn from last time step of previous year
   CALL getvara ('sn', iou6, imt*jmt*km*1, (/1,1,1,1/), (/imt,jmt,km,1/), snp, 1., 0.)
   !WRITE(*,*) 'snp'
   ! snp(144,45,1)=34.0309282726
   !WRITE(*,*) snp(144,45,1)
   ! tn from last time step of current year
   CALL getvara ('tn', iou7, imt*jmt*km*1, (/1,1,1,1/), (/imt,jmt,km,1/), tnc, 1., 0.)
   !WRITE(*,*) 'tnc'
   ! tnc(144,45,1)=2.14057999538
   !WRITE(*,*) tnc(144,45,1)
   ! sn from last time step of current year
   CALL getvara ('sn', iou8, imt*jmt*km*1, (/1,1,1,1/), (/imt,jmt,km,1/), snc, 1., 0.)
   !WRITE(*,*) 'snc'
   ! snc(144,45,1)=34.1197096452
   !WRITE(*,*) snc(144,45,1)
   print*, '-------------------'
   print*, 'Input data read OK!'
   print*, '-------------------'

   CALL closeall

   !!---------------------------------------------------------
   !! Computations of mfo, msftbarot, opottemptend & osalttend
   !!---------------------------------------------------------
   CALL cmip6_mfo
   CALL cmip6_msftbarot
   DO l = 1, ly
      CALL cmip6_tstend(l)
   END DO

   !!-----------------------------------------------------------------
   !! Output mfo, msftbarot, opottemptend & osalttend in NetCDF format
   !!-----------------------------------------------------------------
   iou = 0
   ntrec = 0
   id_time = 0
   id_tbnds = 0
   id_l = 0
   id_s = 0
   id_x = 0
   id_y = 0
   id_z = 0
   ytime = 15768000
   !strlen = 31
      
   ! NETCDF output mfo.nc
   DO i = 1, nline
      line(i) = i
   ENDDO

   ! If the output file does not exist, abort
   INQUIRE (file="mfo.nc", exist=exists)

   IF (.not. exists) THEN
      print*,"output file mfo.nc not found...creating a new file..."
      !CALL flush(6)
      CALL opennew ("mfo.nc", iou)
      ntrec = 1
      CALL redef (iou)

      ! basic grid specification
      CALL defdim ('time_counter', iou, 0, id_time)
      CALL defdim ('bnds', iou, ntbnds, id_tbnds)
      CALL defdim ('line', iou, nline, id_l)
      !CALL defdim ('strlen', iou, 31, id_s)
      CALL defvar ('time_counter', iou, 1, (/id_time/), 0., 0., 'T', 'D'   &
                   , long_name, standard_name, units)
      CALL putatttext (iou, 'time_counter', 'calendar', calendar)
      CALL putatttext (iou, 'time_counter', 'title', title)
      CALL putatttext (iou, 'time_counter', 'time_origin', time_origin)
      CALL putatttext (iou, 'time_counter', 'bounds', bounds)
      CALL defvar ('time_counter_bnds', iou, 2, (/id_tbnds, id_time/), 0., 0., ' ', 'D' &
                   , '', '', '')
      !CALL defvar ('passage', iou, 2, (/id_s, id_l/), 1, 1, 'Y', 'Tstrlen'     &
     !&            , 'ocean passage', 'region', ' ')
      CALL defvar ('line', iou, 1, (/id_l/), 1, 16, ' ', 'I'     &
     &             , 'sections', 'sections', ' ')
      CALL defvar ('mfo', iou, 2, (/id_l, id_time/), 0., 0., ' ', 'F',  &
                   'Sea Water Transport', 'sea_water_transport_across_line', 'kg/s')
      CALL enddef (iou)
      ! define the section axis
      !CALL putvara ('passage', iou, nline, (/1/), (/nline/)           &
     !&              , line, 1., 0.)
      CALL putvara ('line', iou, nline, (/1/), (/nline/)           &
     &              , line, 1., 0.)

      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      ! WRITE(*,*) 'time_bnds'
      ! WRITE(*,*) time_bnds
      ! mfo    
      CALL putvara ('mfo', iou, nline*lm, (/1,1/), (/nline, lm/), mfo(:,:), 1., 0.)
      !ENDDO
      print*, '------------------'
      print*, 'mfo.nc written OK!'
      print*, '------------------'
      call closefile (iou)
   ELSE 
      print*, 'mfo.nc already exists'
   ENDIF

   ! NETCDF output msftbarot.nc
   ! If the output file does not exist, abort
   INQUIRE (file="msftbarot.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file msftbarot.nc not found...creating a new file..."
      !call flush(6)
      CALL opennew ("msftbarot.nc", iou)
      ntrec = 1
      CALL redef (iou)

      ! basic grid specification
      CALL defdim ('time_counter', iou, 0, id_time)
      CALL defdim ('bnds', iou, ntbnds, id_tbnds)
      CALL defdim ('x', iou, imt, id_x)
      CALL defdim ('y', iou, jmt, id_y)
      CALL defvar ('time_counter', iou, 1, (/id_time/), 0., 0., 'T', 'D'   &
                   , long_name, standard_name, units)
      CALL putatttext (iou, 'time_counter', 'calendar', calendar)
      CALL putatttext (iou, 'time_counter', 'title', title)
      CALL putatttext (iou, 'time_counter', 'time_origin', time_origin)
      CALL putatttext (iou, 'time_counter', 'bounds', bounds)
      CALL defvar ('time_counter_bnds', iou, 2, (/id_tbnds, id_time/), 0., 0., ' ', 'D' &
                   , '', '', '')
      !CALL defvar ('nav_lon', iou, 2, (/id_x, id_y/), 0., 0., ' ', 'F', &
      !             'Longitude', 'longitude', 'degrees_east')
      !CALL putatttext (iou, 'nav_lon', 'valid_min', '-180.f')
      !CALL putatttext (iou, 'nav_lon', 'valid_max', '179.9999f')
      !CALL putatttext (iou, 'nav_lon', 'nav_model', 'Default grid')
      !CALL defvar ('nav_lat', iou, 2, (/id_x, id_y/), 0., 0., ' ', 'F', &
      !             'Latitude', 'latitude', 'degrees_north')
      !CALL putatttext (iou, 'nav_lat', 'valid_min', '-78.3935f')
      !CALL putatttext (iou, 'nav_lat', 'valid_max', '89.78919f')
      !CALL putatttext (iou, 'nav_lat', 'nav_model', 'Default grid')
      ! msftbarot
      CALL defvar ('msftbarot', iou, 3, (/id_x, id_y, id_time/), 0., 0., ' ', 'F', &
                   'Ocean Barotropic Mass Streamfunction', 'ocean_barotropic_mass_streamfunction', 'kg/s')
      CALL putatttext (iou, 'msftbarot', 'comment', 'written on u grid')
      CALL putatttext (iou, 'msftbarot', 'coordinates', 'time_counter nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_u(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_u(:,:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('msftbarot', iou, imt*jmt*lm, (/1,1,1/), (/imt, jmt, lm/), &
                     msftbarot(:,:,:), 1., 0.)
      print*, '------------------------'
      print*, 'msftbarot.nc written OK!'
      print*, '------------------------'

      CALL closefile (iou)
   ELSE
      print*, 'msftbarot.nc already exists'
   ENDIF

   ! NETCDF output tstend.nc
   ! If the output file does not exist, abort
   INQUIRE (file="tstend.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file tstend.nc not found...creating a new file..."
      CALL opennew ("tstend.nc", iou)
      ntrec = 1
      CALL redef (iou)

      ! basic grid specification
      CALL defdim ('time_counter', iou, 0, id_time)
      CALL defdim ('bnds', iou, ntbnds, id_tbnds)
      CALL defdim ('x', iou, imt, id_x)
      CALL defdim ('y', iou, jmt, id_y)
      CALL defdim ('deptht', iou, km, id_z)
      CALL defvar ('time_counter', iou, 1, (/id_time/), 0., 0., 'T', 'D'   &
                   , long_name, standard_name, units)
      CALL putatttext (iou, 'time_counter', 'calendar', calendar)
      CALL putatttext (iou, 'time_counter', 'title', title)
      CALL putatttext (iou, 'time_counter', 'time_origin', time_origin)
      CALL putatttext (iou, 'time_counter', 'bounds', bounds)
      CALL defvar ('time_counter_bnds', iou, 2, (/id_tbnds, id_time/), 0., 0., ' ', 'D' &
             , '', '', '')
      !CALL defvar ('nav_lon', iou, 2, (/id_x, id_y/), 0., 0, ' ', 'F', &
      !             'Longitude', 'longitude', 'degrees_east')
      !CALL putatttext (iou, 'nav_lon', 'valid_min', '-179.9942f')
      !CALL putatttext (iou, 'nav_lon', 'valid_max', '179.9903f')
      !CALL putatttext (iou, 'nav_lon', 'nav_model', 'Default grid')
      !CALL defvar ('nav_lat', iou, 2, (/id_x, id_y/), 0., 0., ' ', 'F', &
      !             'Latitude', 'latitude', 'degrees_north')
      !CALL putatttext (iou, 'nav_lat', 'valid_min', '-78.3935f')
      !CALL putatttext (iou, 'nav_lat', 'valid_max', '89.74177f')
      !CALL putatttext (iou, 'nav_lat', 'nav_model', 'Default grid')
      CALL defvar ('deptht', iou, 1, id_z, 0., 0., ' ', 'F', &
                   'Vertical T levels', 'model_level_number', 'm')
      CALL putatttext (iou, 'deptht', 'axis', 'Z')
      CALL putatttext (iou, 'deptht', 'positive', 'down')
      CALL putatttext (iou, 'deptht', 'valid_min', '3.046773f')
      CALL putatttext (iou, 'deptht', 'valid_max', '5875.141f')
      CALL putatttext (iou, 'deptht', 'title', 'deptht')
      ! opottemptend
      CALL defvar ('opottemptend', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'opottemptend', 'tendency_of_sea_water_potential_temperature_expressed_as_heat_content', 'Wm-2')
      CALL putatttext (iou, 'opottemptend', 'coordinates', 'time_counter deptht nav_lat nav_lon')
      ! osalttend
      CALL defvar ('osalttend', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'osalttend', 'tendency_of_sea_water_salinity_expressed_as_salt_content', 'kg m-2 s-1 ')
      CALL putatttext (iou, 'osalttend', 'coordinates', 'time_counter deptht nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_t(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_t(:,:), 1., 0.)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, ly, (/1/), (/ly/), ytime, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*ly, (/1,1/), (/ntbnds,ly/), time_bnds, 1., 0.)
      CALL putvara ('opottemptend', iou, imt*jmt*km*ly, (/1,1,1,1/), (/imt, jmt, km, ly/), opottemptend(:,:,:,:), 1., 0.)
      CALL putvara ('osalttend', iou, imt*jmt*km*ly, (/1,1,1,1/), (/imt, jmt, km, ly/), osalttend(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'tstend.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'tstend.nc already exists'
   ENDIF

END PROGRAM nemo_diag
