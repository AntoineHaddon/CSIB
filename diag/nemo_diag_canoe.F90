PROGRAM nemo_diag_canoe
   
   !!===============================================================
   !!                 ***    PROGRAM nemo_diag_canoe    ***              
   !!                  CMIP6 nemo offline diagnostics
   !!===============================================================
   !! 2018-04 (D. Yang): Original code
   !! 2019-01 (J. Christian): biogeochemistry (CanOE) version
   !!---------------------------------------------------------------
   !!
   !!---------------------------------------------------------------
   !! INPUT FIELDS
   !! 3D: e3t, tmask
   !! 
   !! OUTPUT FIELDS
   !! 3D: [CO3--]sat, [CO3--], pH, [O2]sat
   !!
   !! INPUT FILES
   !! grid_t
   !! ptrc_t
   !! orca_mesh_mask
   !!
   !! OUTPUT FILES
   !! CO3sata.nc, CO3sata.nc, CO3.nc pH.nc, O2sat.nc
   !!---------------------------------------------------------------
   USE nemo_diag_glovars_canoe     ! global variable declarations
   USE nemo_diag_cal_canoe         ! diagnostics calculations
 
   IMPLICIT NONE

   CHARACTER(len=100) :: fname01, fname02, fname03, fname04
   CHARACTER(len=100) :: axis, standard_name, units, calendar, title
   CHARACTER(len=100) :: long_name, time_origin, bounds
   INTEGER   :: iou, iou1, iou2, iou3, iou4
   INTEGER   :: ntrec, id_time, id_tbnds, id_l, id_s, id_x, id_y, id_z
   INTEGER   :: ntbnds
   INTEGER   :: i, l
   LOGICAL   :: exists
   !INTEGER   :: strlen
   INTEGER, DIMENSION(10)            :: ierr
   REAL, DIMENSION(:), ALLOCATABLE   :: time
   REAL, DIMENSION(:), ALLOCATABLE   :: ytime
   REAL, DIMENSION(:), ALLOCATABLE   :: x, y
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
   CALL getdimlen ('tbnds', iou, ntbnds)
   ly = lm / 12

   !!----------------
   !! Allocate Arrays
   !!----------------
   ALLOCATE( e3t(imt,jmt,km), tmask(imt,jmt,km), time_bnds(ntbnds,lm), STAT=ierr(1) )
   ALLOCATE( time(lm), ytime(ly), deptht(km), x(imt), y(jmt), STAT=ierr(2) )
   ALLOCATE( nav_lon_t(imt,jmt), nav_lat_t(imt,jmt), STAT=ierr(3) )
   ALLOCATE( borat(imt,jmt,km,lm), ak13(imt,jmt,km,lm), ak23(imt,jmt,km,lm), akb3(imt,jmt,km,lm), &
     &       akw3(imt,jmt,km,lm), akp13(imt,jmt,km,lm), akp23(imt,jmt,km,lm), akp33(imt,jmt,km,lm), &
     &       aksi3(imt,jmt,km,lm), asi3(imt,jmt,km,lm), STAT=ierr(4) )
   ALLOCATE( TT(imt,jmt,km,lm), SS(imt,jmt,km,lm), CC(imt,jmt,km,lm), AA(imt,jmt,km,lm), &
     &       NO3(imt,jmt,km,lm), NH4(imt,jmt,km,lm), STAT=ierr(5) )
 
   IF (MAXVAL(ierr) /=0) THEN
      STOP 'Memory allocation error in cmip6_nemo_offl'
   ENDIF
      
   iou1 =0
   iou2 =0
   iou3 =0
   iou4 =0

   !!--------------------
   !! Define NetCDF files   
   !!--------------------
   print*,'Reading data on NEMO grid...'
   fname01='orca_mesh_mask'
   fname02='grid_t'
   fname03='ptrc_t'
   fname04='si.nc'
   
   !!------------------------------
   !! Open the defined NetCDF files   
   !!------------------------------
   ! mask/grid info
   CALL openfile (fname01,iou1)
   ! ocean physical variables
   CALL openfile (fname02,iou2)
   CALL openfile (fname03,iou3)
   CALL openfile (fname04,iou4)

   !!-------------------
   !! Get grid/mask data   
   !!-------------------
   CALL getvara ('e3t', iou1, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),e3t , 1., 0.)
   CALL getvara ('tmask', iou1, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),tmask , 1., 0.)

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
   ! nav_lon on grid_T
   CALL getvara ('nav_lon', iou2, imt*jmt, (/1,1/), (/imt,jmt/), nav_lon_t, 1., 0.)
   ! nav_lat on grid_T
   CALL getvara ('nav_lat', iou2, imt*jmt, (/1,1/), (/imt,jmt/), nav_lat_t, 1., 0.)
   ! deptht 
   CALL getvara ('deptht', iou2, km, (/1/), (/km/), deptht, 1., 0.)
   ! temperature
   CALL getvara ('votemper', iou2, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), TT, 1., 0.)
   ! salinity
   CALL getvara ('vosaline', iou2, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), SS, 1., 0.)
   ! DIC
   CALL getvara ('DIC', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), CC, 1., 0.)
   ! alkalinity
   CALL getvara ('TAlk', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), AA, 1., 0.)
   ! Nitrate
   CALL getvara ('NO3', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), NO3, 1., 0.)
   ! Ammonium
   CALL getvara ('NH4', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), NH4, 1., 0.)
   ! Silicate
   CALL getvara ('Si', iou4, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), asi3, 1., 0.)
   print*, '-------------------'
   print*, 'Input data read OK!'
   print*, '-------------------'

   CALL closeall

   !!---------------------------------------------------------
   !! Computations of solubility product, O2 solubility, pH, carbonate ion
   !!---------------------------------------------------------
   CALL cmip6_co3sat
   CALL cmip6_o2sol
   CALL cmip6_cchem

   !!-----------------------------------------------------------------
   !! Output  in NetCDF format
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
      
   ! If the output file does not exist, abort
   INQUIRE (file="CO3sata.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file CO3sata.nc not found...creating a new file..."
      CALL opennew ("CO3sata.nc", iou)
      ntrec = 1
      CALL redef (iou)

      ! basic grid specification
      CALL defdim ('time_counter', iou, 0, id_time)
      CALL defdim ('tbnds', iou, ntbnds, id_tbnds)
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
      CALL defvar ('deptht', iou, 1, id_z, 0., 0., ' ', 'F', &
                   'Vertical T levels', 'model_level_number', 'm')
      CALL putatttext (iou, 'deptht', 'axis', 'Z')
      CALL putatttext (iou, 'deptht', 'positive', 'down')
      CALL putatttext (iou, 'deptht', 'valid_min', '3.046773f')
      CALL putatttext (iou, 'deptht', 'valid_max', '5875.141f')
      CALL putatttext (iou, 'deptht', 'title', 'deptht')
      ! CO3_sat (aragonite)
      CALL defvar ('CO3_sata', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'CO3_sata', '[CO3--] at Aragonite Saturation', 'mol m-3')
      CALL putatttext (iou, 'CO3_sata', 'coordinates', 'time_counter deptht nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_t(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_t(:,:), 1., 0.)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('CO3_sata', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), CO3_sata(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'CO3sata.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'CO3sata.nc already exists'
   ENDIF

   ! If the output file does not exist, abort
   INQUIRE (file="CO3satc.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file CO3satc.nc not found...creating a new file..."
      CALL opennew ("CO3satc.nc", iou)
      ntrec = 1
      CALL redef (iou)

      ! basic grid specification
      CALL defdim ('time_counter', iou, 0, id_time)
      CALL defdim ('tbnds', iou, ntbnds, id_tbnds)
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
      CALL defvar ('deptht', iou, 1, id_z, 0., 0., ' ', 'F', &
                   'Vertical T levels', 'model_level_number', 'm')
      CALL putatttext (iou, 'deptht', 'axis', 'Z')
      CALL putatttext (iou, 'deptht', 'positive', 'down')
      CALL putatttext (iou, 'deptht', 'valid_min', '3.046773f')
      CALL putatttext (iou, 'deptht', 'valid_max', '5875.141f')
      CALL putatttext (iou, 'deptht', 'title', 'deptht')
      ! CO3_sat (calcite)
      CALL defvar ('CO3_satc', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'CO3_satc', '[CO3--] at Calcite Saturation', 'mol m-3')
      CALL putatttext (iou, 'CO3_satc', 'coordinates', 'time_counter deptht nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_t(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_t(:,:), 1., 0.)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('CO3_satc', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), CO3_satc(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'CO3satc.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'CO3satc.nc already exists'
   ENDIF

   ! If the output file does not exist, abort
   INQUIRE (file="o2sol.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file o2sol.nc not found...creating a new file..."
      CALL opennew ("o2sol.nc", iou)
      ntrec = 1
      CALL redef (iou)

      ! basic grid specification
      CALL defdim ('time_counter', iou, 0, id_time)
      CALL defdim ('tbnds', iou, ntbnds, id_tbnds)
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
      CALL defvar ('deptht', iou, 1, id_z, 0., 0., ' ', 'F', &
                   'Vertical T levels', 'model_level_number', 'm')
      CALL putatttext (iou, 'deptht', 'axis', 'Z')
      CALL putatttext (iou, 'deptht', 'positive', 'down')
      CALL putatttext (iou, 'deptht', 'valid_min', '3.046773f')
      CALL putatttext (iou, 'deptht', 'valid_max', '5875.141f')
      CALL putatttext (iou, 'deptht', 'title', 'deptht')
      ! o2sol
      CALL defvar ('o2sol', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'o2sol', 'Oxygen concentration at saturation', 'mol m-3')
      CALL putatttext (iou, 'o2sol', 'coordinates', 'time_counter deptht nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_t(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_t(:,:), 1., 0.)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('o2sol', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), o2sol(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'o2sol.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'o2sol.nc already exists'
   ENDIF

   ! If the output file does not exist, abort
   INQUIRE (file="CO3.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file CO3.nc not found...creating a new file..."
      CALL opennew ("CO3.nc", iou)
      ntrec = 1
      CALL redef (iou)

      ! basic grid specification
      CALL defdim ('time_counter', iou, 0, id_time)
      CALL defdim ('tbnds', iou, ntbnds, id_tbnds)
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
      CALL defvar ('deptht', iou, 1, id_z, 0., 0., ' ', 'F', &
                   'Vertical T levels', 'model_level_number', 'm')
      CALL putatttext (iou, 'deptht', 'axis', 'Z')
      CALL putatttext (iou, 'deptht', 'positive', 'down')
      CALL putatttext (iou, 'deptht', 'valid_min', '3.046773f')
      CALL putatttext (iou, 'deptht', 'valid_max', '5875.141f')
      CALL putatttext (iou, 'deptht', 'title', 'deptht')
      ! CO3
      CALL defvar ('CO3', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'CO3', 'Carbonate ion concentration', 'mol m-3')
      CALL putatttext (iou, 'CO3', 'coordinates', 'time_counter deptht nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_t(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_t(:,:), 1., 0.)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('CO3', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), CO3(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'CO3.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'CO3.nc already exists'
   ENDIF

   ! If the output file does not exist, abort
   INQUIRE (file="pH.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file pH.nc not found...creating a new file..."
      CALL opennew ("pH3D.nc", iou)
      ntrec = 1
      CALL redef (iou)

      ! basic grid specification
      CALL defdim ('time_counter', iou, 0, id_time)
      CALL defdim ('tbnds', iou, ntbnds, id_tbnds)
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
      CALL defvar ('deptht', iou, 1, id_z, 0., 0., ' ', 'F', &
                   'Vertical T levels', 'model_level_number', 'm')
      CALL putatttext (iou, 'deptht', 'axis', 'Z')
      CALL putatttext (iou, 'deptht', 'positive', 'down')
      CALL putatttext (iou, 'deptht', 'valid_min', '3.046773f')
      CALL putatttext (iou, 'deptht', 'valid_max', '5875.141f')
      CALL putatttext (iou, 'deptht', 'title', 'deptht')
      ! pH
      CALL defvar ('pH3D', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'pH3D', 'pH3D', ' ')
      CALL putatttext (iou, 'pH3D', 'coordinates', 'time_counter deptht nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_t(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_t(:,:), 1., 0.)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('pH3D', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), pH(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'pH3D.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'pH3D.nc already exists'
   ENDIF

END PROGRAM nemo_diag_canoe
