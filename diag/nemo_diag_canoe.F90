PROGRAM nemo_diag_canoe
   
   !!===============================================================
   !!                 ***    PROGRAM nemo_diag_canoe    ***              
   !!                  CMIP6 nemo offline diagnostics
   !!===============================================================
   !! 2018-04 (D. Yang): Original code
   !! 2019-01 (J. Christian): biogeochemistry (CanOE) version
   !! 2020-02 (D. Yang): Revise tbnds to bnds in defdim to avoid crash
   !!                    when later "cdo mergetime" in canesm_nemo_bgc_diag.sh
   !!                    makes the unwanted change from tbnds to bnds.
   !! 2020-03 (D. Yang): Add getdimnm to accommodate both tbnds & bnds
   !!---------------------------------------------------------------
   !!
   !!---------------------------------------------------------------
   !! INPUT FIELDS
   !! 3D: e3t, tmask
   !! 
   !! OUTPUT FIELDS
   !! 3D: [CO3--]sat, [CO3--], pH, [O2]sat, Omega_A, Omega_C
   !! 2D: calcite and aragonite saturation depth, minimum [O2], depth of minimum [O2]
   !!
   !! INPUT FILES
   !! grid_t
   !! ptrc_t
   !! orca_mesh_mask
   !!
   !! OUTPUT FILES
   !! CO3sata.nc, CO3sata.nc, CO3.nc pH3D.nc, O2sat.nc, Omega_C.nc, Omega_A.nc
   !!---------------------------------------------------------------
   USE nemo_diag_glovars_canoe     ! global variable declarations
   USE nemo_diag_cal_canoe         ! diagnostics calculations
 
   IMPLICIT NONE

   CHARACTER(len=100) :: fname01, fname02, fname03, fname04
   CHARACTER(len=100) :: axis, standard_name, units, calendar, title
   CHARACTER(len=100) :: long_name, time_origin, bounds
   CHARACTER(len=100) :: dimnm
   INTEGER   :: iou, iou1, iou2, iou3, iou4
   INTEGER   :: ntrec, id_time, id_tbnds, id_l, id_s, id_x, id_y, id_z
   INTEGER   :: ntbnds, ndim, ntdim
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
   ! number of total dimensions (x, y, deptht, time_counter & tbnds or bnds)
   ntdim = 5
   CALL openfile  ("grid_t", iou)
   CALL getdimlen ('x', iou, imt)
   CALL getdimlen ('y', iou, jmt)
   call getdimlen ('deptht', iou, km)
   CALL getdimlen ('time_counter', iou, lm)      
   ! get grid size for tbnds/bnds
   ntbnds=2
   Do i = 1, ntdim
      CALL getdimnm  (dimnm, iou, i, ndim)
      print*, 'DIM',i,':',dimnm, 'length:', ndim
      IF (dimnm .eq. 'tbnds' .or. dimnm .eq. 'bnds' .or. dimnm .eq. 'axis_nbounds' ) ntbnds = ndim
   END DO

   ly = lm / 12

   !!----------------
   !! Allocate Arrays
   !!----------------
   ierr=0
   ALLOCATE( e3t(imt,jmt,km,lm), tmask(imt,jmt,km), time_bnds(ntbnds,lm), STAT=ierr(1) )
   ALLOCATE( time(lm), ytime(ly), deptht(km), x(imt), y(jmt), STAT=ierr(2) )
   ALLOCATE( nav_lon_t(imt,jmt), nav_lat_t(imt,jmt), STAT=ierr(3) )
   ALLOCATE( TT(imt,jmt,km,lm), SS(imt,jmt,km,lm), CC(imt,jmt,km,lm), AA(imt,jmt,km,lm), &
     &       NO3(imt,jmt,km,lm), NH4(imt,jmt,km,lm), O2(imt,jmt,km,lm), asi3(imt,jmt,km,lm), STAT=ierr(4) )
   ALLOCATE( K_sp_cal(imt,jmt,km,lm), K_sp_arag(imt,jmt,km,lm), Omega_C(imt,jmt,km,lm), Omega_A(imt,jmt,km,lm), STAT=ierr(5) )
   ALLOCATE( prhop(imt,jmt,km,lm), pH(imt,jmt,km,lm), CO3(imt,jmt,km,lm), co3_satc(imt,jmt,km,lm), &
     &       co3_sata(imt,jmt,km,lm), o2sol(imt,jmt,km,lm), STAT=ierr(6) )
   ALLOCATE( zsat_c(imt,jmt,lm), zsat_a(imt,jmt,lm), o2min(imt,jmt,lm), zo2min(imt,jmt,lm), STAT=ierr(7) )
 
   IF (MAXVAL(ierr) /=0) THEN
      print*,'ierr=',ierr
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
   !CALL getvara ('e3t', iou1, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),e3t , 1., 0.)
   CALL getvara ('tmask', iou1, imt*jmt*km, (/1,1,1,1/), (/imt,jmt,km,1/),tmask , 1., 0.)
   ! Mask out Caspian in CCCma ORCA1 grid
   IF ( (imt == 362) .AND. (jmt == 292) ) THEN
     tmask( 332:344, 203:235, 1:km ) = 0.
   ELSEIF ( (imt == 362) .AND. (jmt == 332) ) THEN ! eORCA1 grid (no jstart)
     tmask( 332:344, 243:275, 1:km ) = 0.
   ELSE
     stop "NEMO BGC diagnostic deck expects ORCA R1 grid"
   ENDIF

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
   CALL getatttext (iou2, 'time_counter', 'long_name', long_name)
   CALL getatttext (iou2, 'time_counter', 'time_origin', time_origin)
   CALL getatttext (iou2, 'time_counter', 'bounds', bounds)
   ! time_counter_bnds
   CALL getvara ('time_counter_bounds', iou2, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
   ! nav_lon on grid_T
   CALL getvara ('nav_lon', iou2, imt*jmt, (/1,1/), (/imt,jmt/), nav_lon_t, 1., 0.)
   ! nav_lat on grid_T
   CALL getvara ('nav_lat', iou2, imt*jmt, (/1,1/), (/imt,jmt/), nav_lat_t, 1., 0.)
   ! deptht 
   CALL getvara ('deptht', iou2, km, (/1/), (/km/), deptht, 1., 0.)
   CALL getvara ('thkcello', iou2, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/),e3t , 1., 0.)
   ! temperature
   CALL getvara ('thetao', iou2, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), TT, 1., 0.)
   ! salinity
   CALL getvara ('so', iou2, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), SS, 1., 0.)
   ! DIC
   CALL getvara ('DIC', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), CC, 1., 0.)
   ! alkalinity
   CALL getvara ('Alkalini', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), AA, 1., 0.)
   ! Nitrate
   CALL getvara ('NO3', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), NO3, 1., 0.)
   ! Ammonium
   CALL getvara ('NH4', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), NH4, 1., 0.)
   ! Oxygen
   CALL getvara ('O2', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), O2, 1., 0.)
   ! Silicate
   CALL getvara ('Si', iou4, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), asi3, 1., 0.)
   print*, '-------------------'
   print*, 'Input data read OK!'
   print*, '-------------------'

   CALL closeall

   !!---------------------------------------------------------
   !! Computations of solubility product, O2 solubility, pH, carbonate ion
   !!---------------------------------------------------------
   NO3=NO3+NH4      ! nitrate is used to estimate PO4 for carbon chemistry: in CanOE should include NH4
   CALL density
   CALL cmip6_co3sat(K_sp_arag,K_sp_cal)
   CALL cmip6_o2sol
   CALL cmip6_zo2min
   CALL cmip6_cchem(CC,AA,K_sp_cal,K_sp_arag,CO3,pH,Omega_C,Omega_A)
   CALL saturation_depth(Omega_C,Omega_A)

   DEALLOCATE( TT, SS, CC, AA, NO3, NH4, O2, asi3 )

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
      CALL defdim ('bnds', iou, ntbnds, id_tbnds)
      CALL defdim ('x', iou, imt, id_x)
      CALL defdim ('y', iou, jmt, id_y)
      CALL defdim ('deptht', iou, km, id_z)
      CALL defvar ('time_counter', iou, 1, (/id_time/), 0., 0., 'T', 'D'   &
                   , long_name, standard_name, units)
      CALL putatttext (iou, 'time_counter', 'calendar', calendar)
      CALL putatttext (iou, 'time_counter', 'time_origin', time_origin)
      CALL putatttext (iou, 'time_counter', 'bounds', bounds)
      CALL putatttext (iou, 'time_counter', 'axis', axis)
      CALL putatttext (iou, 'time_counter', 'standard_name', standard_name)
      CALL putatttext (iou, 'time_counter', 'units', units)
      CALL putatttext (iou, 'time_counter', 'long_name', long_name)
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
      CALL defvar ('CO3sata', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'CO3sata', '[CO3--] at Aragonite Saturation', 'mol m-3')
      CALL putatttext (iou, 'CO3sata', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_t(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_t(:,:), 1., 0.)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('CO3sata', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), CO3_sata(:,:,:,:), 1., 0.)
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
      CALL defdim ('bnds', iou, ntbnds, id_tbnds)
      CALL defdim ('x', iou, imt, id_x)
      CALL defdim ('y', iou, jmt, id_y)
      CALL defdim ('deptht', iou, km, id_z)
      CALL defvar ('time_counter', iou, 1, (/id_time/), 0., 0., 'T', 'D'   &
                   , long_name, standard_name, units)
      CALL putatttext (iou, 'time_counter', 'calendar', calendar)
      CALL putatttext (iou, 'time_counter', 'axis', axis)
      CALL putatttext (iou, 'time_counter', 'standard_name', standard_name)
      CALL putatttext (iou, 'time_counter', 'units', units)
      CALL putatttext (iou, 'time_counter', 'long_name', long_name)
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
      CALL defvar ('CO3satc', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'CO3satc', '[CO3--] at Calcite Saturation', 'mol m-3')
      CALL putatttext (iou, 'CO3satc', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_t(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_t(:,:), 1., 0.)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('CO3satc', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), CO3_satc(:,:,:,:), 1., 0.)
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
      CALL defdim ('bnds', iou, ntbnds, id_tbnds)
      CALL defdim ('x', iou, imt, id_x)
      CALL defdim ('y', iou, jmt, id_y)
      CALL defdim ('deptht', iou, km, id_z)
      CALL defvar ('time_counter', iou, 1, (/id_time/), 0., 0., 'T', 'D'   &
                   , long_name, standard_name, units)
      CALL putatttext (iou, 'time_counter', 'calendar', calendar)
      CALL putatttext (iou, 'time_counter', 'axis', axis)
      CALL putatttext (iou, 'time_counter', 'standard_name', standard_name)
      CALL putatttext (iou, 'time_counter', 'units', units)
      CALL putatttext (iou, 'time_counter', 'long_name', long_name)
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
      CALL putatttext (iou, 'o2sol', 'coordinates', 'nav_lat nav_lon')
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
      CALL defdim ('bnds', iou, ntbnds, id_tbnds)
      CALL defdim ('x', iou, imt, id_x)
      CALL defdim ('y', iou, jmt, id_y)
      CALL defdim ('deptht', iou, km, id_z)
      CALL defvar ('time_counter', iou, 1, (/id_time/), 0., 0., 'T', 'D'   &
                   , long_name, standard_name, units)
      CALL putatttext (iou, 'time_counter', 'calendar', calendar)
      CALL putatttext (iou, 'time_counter', 'axis', axis)
      CALL putatttext (iou, 'time_counter', 'standard_name', standard_name)
      CALL putatttext (iou, 'time_counter', 'units', units)
      CALL putatttext (iou, 'time_counter', 'long_name', long_name)
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
      CALL putatttext (iou, 'CO3', 'coordinates', 'nav_lat nav_lon')
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
   INQUIRE (file="pH3D.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file pH3D.nc not found...creating a new file..."
      CALL opennew ("pH3D.nc", iou)
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
      CALL putatttext (iou, 'time_counter', 'axis', axis)
      CALL putatttext (iou, 'time_counter', 'standard_name', standard_name)
      CALL putatttext (iou, 'time_counter', 'units', units)
      CALL putatttext (iou, 'time_counter', 'long_name', long_name)
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
      CALL putatttext (iou, 'pH3D', 'coordinates', 'nav_lat nav_lon')
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

   ! If the output file does not exist, abort
   INQUIRE (file="Omega_C.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file Omega_C.nc not found...creating a new file..."
      CALL opennew ("Omega_C.nc", iou)
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
      CALL putatttext (iou, 'time_counter', 'axis', axis)
      CALL putatttext (iou, 'time_counter', 'standard_name', standard_name)
      CALL putatttext (iou, 'time_counter', 'units', units)
      CALL putatttext (iou, 'time_counter', 'long_name', long_name)
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
      ! Omega_C
      CALL defvar ('Omega_C', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'Omega_C', 'Calcite saturation state', '1')
      CALL putatttext (iou, 'Omega_C', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('Omega_C', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), Omega_C(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'Omega_C.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'Omega_C.nc already exists'
   ENDIF

   ! If the output file does not exist, abort
   INQUIRE (file="Omega_A.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file Omega_A.nc not found...creating a new file..."
      CALL opennew ("Omega_A.nc", iou)
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
      CALL putatttext (iou, 'time_counter', 'axis', axis)
      CALL putatttext (iou, 'time_counter', 'standard_name', standard_name)
      CALL putatttext (iou, 'time_counter', 'units', units)
      CALL putatttext (iou, 'time_counter', 'long_name', long_name)
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
      ! Omega_A
      CALL defvar ('Omega_A', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'Omega_A', 'Aragonite saturation state', '1')
      CALL putatttext (iou, 'Omega_A', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('Omega_A', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), Omega_A(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'Omega_A.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'Omega_A.nc already exists'
   ENDIF

   ! If the output file does not exist, abort
   INQUIRE (file="Zsat_A.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file Zsat_A.nc not found...creating a new file..."
      CALL opennew ("Zsat_A.nc", iou)
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
      CALL putatttext (iou, 'time_counter', 'axis', axis)
      CALL putatttext (iou, 'time_counter', 'standard_name', standard_name)
      CALL putatttext (iou, 'time_counter', 'units', units)
      CALL putatttext (iou, 'time_counter', 'long_name', long_name)
      CALL putatttext (iou, 'time_counter', 'time_origin', time_origin)
      CALL putatttext (iou, 'time_counter', 'bounds', bounds)
      CALL defvar ('time_counter_bnds', iou, 2, (/id_tbnds, id_time/), 0., 0., ' ', 'D' &
             , '', '', '')
      ! Zsat_A
      CALL defvar ('Zsat_A', iou, 3, (/id_x, id_y, id_time/), 0., 0., ' ', 'F', &
                   'Zsat_A', 'Aragonite saturation horizon depth', 'm')
      CALL putatttext (iou, 'Zsat_A', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('Zsat_A', iou, imt*jmt*lm, (/1,1,1/), (/imt, jmt, lm/), zsat_a(:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'Zsat_A.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'Zsat_A.nc already exists'
   ENDIF

   ! If the output file does not exist, abort
   INQUIRE (file="Zsat_C.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file Zsat_C.nc not found...creating a new file..."
      CALL opennew ("Zsat_C.nc", iou)
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
      CALL putatttext (iou, 'time_counter', 'axis', axis)
      CALL putatttext (iou, 'time_counter', 'standard_name', standard_name)
      CALL putatttext (iou, 'time_counter', 'units', units)
      CALL putatttext (iou, 'time_counter', 'long_name', long_name)
      CALL putatttext (iou, 'time_counter', 'time_origin', time_origin)
      CALL putatttext (iou, 'time_counter', 'bounds', bounds)
      CALL defvar ('time_counter_bnds', iou, 2, (/id_tbnds, id_time/), 0., 0., ' ', 'D' &
             , '', '', '')
      ! Zsat_C
      CALL defvar ('Zsat_C', iou, 3, (/id_x, id_y, id_time/), 0., 0., ' ', 'F', &
                   'Zsat_C', 'Calcite saturation horizon depth', 'm')
      CALL putatttext (iou, 'Zsat_C', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('Zsat_C', iou, imt*jmt*lm, (/1,1,1/), (/imt, jmt, lm/), zsat_c(:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'Zsat_C.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'Zsat_C.nc already exists'
   ENDIF

   ! If the output file does not exist, abort
   INQUIRE (file="o2min.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file o2min.nc not found...creating a new file..."
      CALL opennew ("o2min.nc", iou)
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
      CALL putatttext (iou, 'time_counter', 'axis', axis)
      CALL putatttext (iou, 'time_counter', 'standard_name', standard_name)
      CALL putatttext (iou, 'time_counter', 'units', units)
      CALL putatttext (iou, 'time_counter', 'long_name', long_name)
      CALL putatttext (iou, 'time_counter', 'time_origin', time_origin)
      CALL putatttext (iou, 'time_counter', 'bounds', bounds)
      CALL defvar ('time_counter_bnds', iou, 2, (/id_tbnds, id_time/), 0., 0., ' ', 'D' &
             , '', '', '')
      ! o2min
      CALL defvar ('o2min', iou, 3, (/id_x, id_y, id_time/), 0., 0., ' ', 'F', &
                   'o2min', 'Minimum oxygen concentration', 'mol m^-3')
      CALL putatttext (iou, 'o2min', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('o2min', iou, imt*jmt*lm, (/1,1,1/), (/imt, jmt, lm/), o2min(:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'o2min.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'o2min.nc already exists'
   ENDIF

   ! If the output file does not exist, abort
   INQUIRE (file="zo2min.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file zo2min.nc not found...creating a new file..."
      CALL opennew ("zo2min.nc", iou)
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
      CALL putatttext (iou, 'time_counter', 'axis', axis)
      CALL putatttext (iou, 'time_counter', 'standard_name', standard_name)
      CALL putatttext (iou, 'time_counter', 'units', units)
      CALL putatttext (iou, 'time_counter', 'long_name', long_name)
      CALL putatttext (iou, 'time_counter', 'time_origin', time_origin)
      CALL putatttext (iou, 'time_counter', 'bounds', bounds)
      CALL defvar ('time_counter_bnds', iou, 2, (/id_tbnds, id_time/), 0., 0., ' ', 'D' &
             , '', '', '')
      ! zo2min
      CALL defvar ('zo2min', iou, 3, (/id_x, id_y, id_time/), 0., 0., ' ', 'F', &
                   'zo2min', 'Depth of minimum oxygen concentration', 'm')
      CALL putatttext (iou, 'zo2min', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('zo2min', iou, imt*jmt*lm, (/1,1,1/), (/imt, jmt, lm/), zo2min(:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'zo2min.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'zo2min.nc already exists'
   ENDIF



END PROGRAM nemo_diag_canoe

