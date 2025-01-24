PROGRAM nemo_diag_cmoc
   
   !!===============================================================
   !!                 ***    PROGRAM nemo_diag_cmoc    ***              
   !!                  CMIP6 nemo offline diagnostics
   !!===============================================================
   !! 2018-04 (D. Yang): Original code
   !! 2019-02 (J. Christian): biogeochemistry (CMOC) version
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
   !! 3D: [CO3--]sat, [CO3--], pH, [O2]sat, abiotic and natural pH and [CO3--], Omega_A, Omega_C
   !! 2D: calcite and aragonite saturation depth, minimum [O2], depth of minimum [O2]
   !!
   !! INPUT FILES
   !! grid_t
   !! ptrc_t
   !! orca_mesh_mask
   !!
   !! OUTPUT FILES
   !! CO3sata.nc, CO3satc.nc, o2sol.nc, CO3.nc, pH3D.nc, CO3abio.nc, pHabio.nc, CO3nat.nc, pHnat.nc, Omega_C.nc, Omega_C_nat.nc, Omega_C_abio.nc, Omega_A.nc, Omega_A_nat.nc, Omega_A_abio.nc, 
   !! Zsat_A.nc, Zsat_C.nc, o2min.nc, zo2min.nc
   !!---------------------------------------------------------------
   USE nemo_diag_glovars_cmoc     ! global variable declarations
   USE nemo_diag_cal_cmoc         ! diagnostics calculations
 
   IMPLICIT NONE

   CHARACTER(len=100) :: fname01, fname02, fname03, fname04, fname05
   CHARACTER(len=100) :: axis, standard_name, units, calendar, title
   CHARACTER(len=100) :: long_name, time_origin, bounds
   CHARACTER(len=100) :: dimnm
   INTEGER   :: iou, iou1, iou2, iou3, iou4, iou5
   INTEGER   :: ntrec, id_time, id_tbnds, id_l, id_s, id_x, id_y, id_z
   INTEGER   :: ntbnds, ndim, ntdim
   INTEGER   :: i, j, k, l
   LOGICAL   :: exists
   LOGICAL, PARAMETER :: process_abio = .false.
   !INTEGER   :: strlen
   INTEGER, DIMENSION(10)            :: ierr
   REAL, PARAMETER :: Ca=0.010280
   REAL, DIMENSION(:),       ALLOCATABLE :: time
   REAL, DIMENSION(:),       ALLOCATABLE :: ytime
   REAL, DIMENSION(:),       ALLOCATABLE :: x, y
   REAL, DIMENSION(:,:),     ALLOCATABLE :: nav_lon_t, nav_lat_t
   REAL, DIMENSION(:,:),     ALLOCATABLE :: time_bnds
   REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: alk_abio
   REAL                                  :: sss_glob_avg, r_sss_glob_avg
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
   ALLOCATE( TT(imt,jmt,km,lm), SS(imt,jmt,km,lm), CC(imt,jmt,km,lm), CAB(imt,jmt,km,lm), CNT(imt,jmt,km,lm), &
             AA(imt,jmt,km,lm), NO3(imt,jmt,km,lm), O2(imt,jmt,km,lm), asi3(imt,jmt,km,lm), STAT=ierr(4) )
   ALLOCATE( pHfull(imt,jmt,km,lm), CO3full(imt,jmt,km,lm), pHabio(imt,jmt,km,lm), CO3abio(imt,jmt,km,lm), &
             pHnat(imt,jmt,km,lm), CO3nat(imt,jmt,km,lm), alk_abio(imt,jmt,km,lm), STAT=ierr(5) )
   ALLOCATE( K_sp_cal(imt,jmt,km,lm), K_sp_arag(imt,jmt,km,lm), Omega_C(imt,jmt,km,lm), Omega_A(imt,jmt,km,lm), &
             Omega_C_abio(imt,jmt,km,lm), Omega_A_abio(imt,jmt,km,lm), Omega_C_nat(imt,jmt,km,lm), &
             Omega_A_nat(imt,jmt,km,lm), STAT=ierr(6) )
   ALLOCATE( zsat_c(imt,jmt,lm), zsat_a(imt,jmt,lm), o2min(imt,jmt,lm), zo2min(imt,jmt,lm), STAT=ierr(7) )
   ALLOCATE( prhop(imt,jmt,km,lm), co3_satc(imt,jmt,km,lm), co3_sata(imt,jmt,km,lm), o2sol(imt,jmt,km,lm), STAT=ierr(8) )
 
   IF (MAXVAL(ierr) /=0) THEN
      STOP 'Memory allocation error in cmip6_nemo_offl'
   ENDIF
      
   iou1=0
   iou2=0
   iou3=0
   iou4=0
   iou5=0

   !!--------------------
   !! Define NetCDF files   
   !!--------------------
   print*,'Reading data on NEMO grid...'
   fname01='orca_mesh_mask'
   fname02='grid_t'
   fname03='ptrc_t'
   fname04='si.nc'
   fname05='sss_glob_avg.nc'
   
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

   !!-------------------
   !! Get grid/mask data   
   !!-------------------
   CALL getvara ('thkcello', iou2, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/),e3t , 1., 0.)
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
   ! temperature
   CALL getvara ('thetao', iou2, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), TT, 1., 0.)
   ! salinity
   CALL getvara ('so', iou2, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), SS, 1., 0.)
   ! DIC
   CALL getvara ('DIC', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), CC, 1., 0.)
   ! abiotic DIC
   if (process_abio) CALL getvara ('DICabio', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), CAB, 1., 0.)
   ! natural DIC
   CALL getvara ('DICnat', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), CNT, 1., 0.)
   ! alkalinity
   CALL getvara ('Alkalini', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), AA, 1., 0.)
   ! Nitrate
   CALL getvara ('NO3', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), NO3, 1., 0.)
   ! Oxygen
   CALL getvara ('O2', iou3, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), O2, 1., 0.)
   ! Silicate
   CALL getvara ('Si', iou4, imt*jmt*km*lm, (/1,1,1,1/), (/imt,jmt,km,lm/), asi3, 1., 0.)
   ! Globally averaged Salinity 
   if (process_abio) CALL getvara ('sss_glob_avg', iou5, 1, (/1/), (/1/), sss_glob_avg, 1., 0.)
   print*, '-------------------'
   print*, 'Input data read OK!'
   print*, '-------------------'

   CALL closeall

   !!Calculate the abiotic alkalinity as specified in Orr et al. 2017
   if (process_abio) then
   r_sss_glob_avg = 1./sss_glob_avg
   DO i=1,imt ; DO j=1,jmt ; DO k=1,km ; DO l=1,lm
     alk_abio(i,j,k,l) = ( 2297. )*( SS(i,j,k,l)*r_sss_glob_avg ) ! Equation 27 of Orr et al. 2017 in micromol
   ENDDO ; ENDDO; ENDDO ; ENDDO
   endif
   !!---------------------------------------------------------
   !! Computations of solubility product, O2 solubility, pH, carbonate ion
   !!---------------------------------------------------------
   CALL density
   CALL cmip6_co3sat(K_sp_arag,K_sp_cal)
   CALL cmip6_o2sol
   CALL cmip6_zo2min
   CALL cmip6_cchem(CC,AA,K_sp_cal,K_sp_arag,CO3full,pHfull,Omega_C,Omega_A)
   if (process_abio) CALL cmip6_cchem(CAB,alk_abio,K_sp_cal,K_sp_arag,CO3abio,pHabio,Omega_C_abio,Omega_A_abio)
   CALL cmip6_cchem(CNT,AA,K_sp_cal,K_sp_arag,CO3nat,pHnat,Omega_C_nat,Omega_A_nat)
   CALL saturation_depth(Omega_C,Omega_A)

   DEALLOCATE( TT, SS, CC, AA, CAB, CNT, NO3, O2, asi3, alk_abio)

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
      CALL putvara ('CO3sata', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), co3_sata(:,:,:,:), 1., 0.)
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
      CALL putvara ('CO3satc', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), co3_satc(:,:,:,:), 1., 0.)
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
      CALL putvara ('CO3', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), CO3full(:,:,:,:), 1., 0.)
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
      CALL putvara ('pH3D', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), pHfull(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'pH3D.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'pH3D.nc already exists'
   ENDIF

   IF (process_abio) THEN
   ! If the output file does not exist, abort
   INQUIRE (file="CO3abio.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file CO3abio.nc not found...creating a new file..."
      CALL opennew ("CO3abio.nc", iou)
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
      ! CO3abio
      CALL defvar ('CO3abio', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'CO3abio', 'Carbonate ion concentration', 'mol m-3')
      CALL putatttext (iou, 'CO3abio', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_t(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_t(:,:), 1., 0.)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('CO3abio', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), CO3abio(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'CO3abio.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'CO3abio.nc already exists'
   ENDIF

   ! If the output file does not exist, abort
   INQUIRE (file="pHabio.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file pHabio.nc not found...creating a new file..."
      CALL opennew ("pHabio.nc", iou)
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
      ! pHabio
      CALL defvar ('pHabio', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'pHabio', 'pHabio', ' ')
      CALL putatttext (iou, 'pHabio', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_t(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_t(:,:), 1., 0.)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('pHabio', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), pHabio(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'pHabio.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'pHabio.nc already exists'
   ENDIF
   ENDIF

   ! If the output file does not exist, abort
   INQUIRE (file="CO3nat.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file CO3nat.nc not found...creating a new file..."
      CALL opennew ("CO3nat.nc", iou)
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
      ! CO3nat
      CALL defvar ('CO3nat', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'CO3nat', 'Carbonate ion concentration', 'mol m-3')
      CALL putatttext (iou, 'CO3nat', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_t(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_t(:,:), 1., 0.)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('CO3nat', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), CO3nat(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'CO3nat.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'CO3nat.nc already exists'
   ENDIF

   ! If the output file does not exist, abort
   INQUIRE (file="pHnat.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file pHnat.nc not found...creating a new file..."
      CALL opennew ("pHnat.nc", iou)
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
      ! pHnat
      CALL defvar ('pHnat', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'pHnat', 'pHnat', ' ')
      CALL putatttext (iou, 'pHnat', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      !CALL putvara ('nav_lon', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lon_t(:,:), 1., 0.)
      !CALL putvara ('nav_lat', iou, imt*jmt, (/1,1/), (/imt, jmt/), nav_lat_t(:,:), 1., 0.)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('pHnat', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), pHnat(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'pHnat.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'pHnat.nc already exists'
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
   INQUIRE (file="Omega_C_nat.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file Omega_C_nat.nc not found...creating a new file..."
      CALL opennew ("Omega_C_nat.nc", iou)
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
      ! Omega_C_nat
      CALL defvar ('Omega_C_nat', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'Omega_C_nat', 'Calcite saturation state for natural DIC', '1')
      CALL putatttext (iou, 'Omega_C_nat', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('Omega_C_nat', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), Omega_C_nat(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'Omega_C_nat.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'Omega_C_nat.nc already exists'
   ENDIF

   IF (process_abio) THEN
   ! If the output file does not exist, abort
   INQUIRE (file="Omega_C_abio.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file Omega_C_abio.nc not found...creating a new file..."
      CALL opennew ("Omega_C_abio.nc", iou)
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
      ! Omega_C_abio
      CALL defvar ('Omega_C_abio', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'Omega_C_abio', 'Calcite saturation state for abiotic DIC', '1')
      CALL putatttext (iou, 'Omega_C_abio', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('Omega_C_abio', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), Omega_C_abio(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'Omega_C_abio.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'Omega_C_abio.nc already exists'
   ENDIF
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
   INQUIRE (file="Omega_A_nat.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file Omega_A_nat.nc not found...creating a new file..."
      CALL opennew ("Omega_A_nat.nc", iou)
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
      ! Omega_A_nat
      CALL defvar ('Omega_A_nat', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'Omega_A_nat', 'Aragonite saturation state for natural DIC', '1')
      CALL putatttext (iou, 'Omega_A_nat', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('Omega_A_nat', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), Omega_A_nat(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'Omega_A_nat.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'Omega_A_nat.nc already exists'
   ENDIF

   IF (process_abio) THEN
   ! If the output file does not exist, abort
   INQUIRE (file="Omega_A_abio.nc", exist=exists)
   IF (.not. exists) THEN
      print*,"output file Omega_A_abio.nc not found...creating a new file..."
      CALL opennew ("Omega_A_abio.nc", iou)
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
      ! Omega_A_abio
      CALL defvar ('Omega_A_abio', iou, 4, (/id_x, id_y, id_z, id_time/), 0., 0., ' ', 'F', &
                   'Omega_A_abio', 'Aragonite saturation state for abiotic DIC', '1')
      CALL putatttext (iou, 'Omega_A_abio', 'coordinates', 'nav_lat nav_lon')
      CALL enddef (iou)
      CALL putvara ('deptht', iou, km, (/1/), (/km/), deptht(:), 1., 0.)
      CALL putvara ('time_counter', iou, lm, (/1/), (/lm/), time, 1., 0.)
      CALL putvara ('time_counter_bnds', iou, ntbnds*lm, (/1,1/), (/ntbnds,lm/), time_bnds, 1., 0.)
      CALL putvara ('Omega_A_abio', iou, imt*jmt*km*lm, (/1,1,1,1/), (/imt, jmt, km, lm/), Omega_A_abio(:,:,:,:), 1., 0.)
      print*, '---------------------'
      print*, 'Omega_A_abio.nc written OK!'
      print*, '---------------------'
      CALL closefile (iou)
   ELSE
      print*, 'Omega_A_abio.nc already exists'
   ENDIF
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

END PROGRAM nemo_diag_cmoc

