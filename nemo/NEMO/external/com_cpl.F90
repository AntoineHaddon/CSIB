!--- cpp defs (set via a combination of canesm.cfg, jobfile, and set_sizes)
#include "cppdef_config.h90"
#include "cppdef_sizes.h90"

module com_cpl
  !***************************************************************
  ! Common MPI initialization and communication routines and data
  !***************************************************************

#ifdef use_mpi
  use mpi
#endif

  implicit none

  private

  public :: impi
  public :: nlon_a, nlat_a, nlon_o, nlat_o, nlon_i, nlat_i
  public :: olap_a, olap_o, olap_i
  public :: nlon_canom, nlat_canom, olap_canom
  public :: find_cpl_vinfo
  public :: cpl_vinfo_t
  public :: copy_cpl_vinfo_t
  public :: print_cpl_vinfo_t
  public :: format_ibuf
  public :: cpl_var_list
  public :: add_cpl_var
  public :: define_cpl_var_by_name
  public :: event, nevents
  public :: overlap
  public :: uniq_names
  public :: print_var_stats, sprint_var_stats
  public :: cpl_var_new_tag
  public :: is_special_op
  public :: filter_group_ID
  public :: idx2d_from_1d
  public :: cpl_var_name_exists

#ifdef use_mpi
  public :: MPI_COMM_WORLD
  public :: MPI_REAL, MPI_REAL4, MPI_REAL8, MPI_DOUBLE_PRECISION
  public :: MPI_INTEGER, MPI_INTEGER4, MPI_INTEGER8
  public :: MPI_CHARACTER, MPI_UNDEFINED, MPI_STATUS_SIZE
  public :: define_group
  public :: com_mpi_init
  public :: cpl_initialize_events
  public :: recv_fld, send_fld
  public :: recv_data_rec, send_data_rec
  public :: bcastGroup
  public :: bcast_inter

  !--- Define generic interfaces for recv_data_rec and send_data_rec
  !--- to accomodate various input data array types and shapes
  interface recv_data_rec
    module procedure recv_data_rec_r1d
    module procedure recv_data_rec_r2d
    module procedure recv_data_rec_r3d
    module procedure recv_data_rec_i1d
    module procedure recv_data_rec_i1d_noibuf
  end interface

  interface send_data_rec
    module procedure send_data_rec_r1d
    module procedure send_data_rec_r1d_noibuf
    module procedure send_data_rec_r2d
    module procedure send_data_rec_r2d_noibuf
    module procedure send_data_rec_r3d
    module procedure send_data_rec_i1d
    module procedure send_data_rec_i1d_noibuf
  end interface

  !--- Define generic interfaces for recv_fld and send_fld
  !--- to accomodate various input data array types and shapes
  interface recv_fld
    module procedure recv_fld_r1d
    module procedure recv_fld_r2d
    module procedure recv_fld_i1d
  end interface

  interface send_fld
    module procedure send_fld_r1d
    module procedure send_fld_r2d
    module procedure send_fld_i1d
  end interface

  interface bcastGroup
    module procedure bcastGroup_Scalar_i4
    module procedure bcastGroup_1d_Array_i4
    module procedure bcastGroup_2d_Array_i4
    module procedure bcastGroup_Scalar_i8
    module procedure bcastGroup_1d_Array_i8
    module procedure bcastGroup_Scalar_r4
    module procedure bcastGroup_Scalar_r8
    module procedure bcastGroup_Scalar_logical
    module procedure bcastGroup_Scalar_char
    module procedure bcastGroup_1d_Array_char
  end interface

  interface bcast_inter
    module procedure bcast_inter_Scalar_char
    module procedure bcast_inter_1d_Array_char
    module procedure bcast_inter_2d_Array_char
    module procedure bcast_inter_3d_Array_char
    module procedure bcast_inter_Scalar_logical
    module procedure bcast_inter_1d_Array_logical
    module procedure bcast_inter_2d_Array_logical
    module procedure bcast_inter_3d_Array_logical
    module procedure bcast_inter_Scalar_i4
    module procedure bcast_inter_1d_Array_i4
    module procedure bcast_inter_2d_Array_i4
    module procedure bcast_inter_3d_Array_i4
    module procedure bcast_inter_Scalar_i8
    module procedure bcast_inter_1d_Array_i8
    module procedure bcast_inter_2d_Array_i8
    module procedure bcast_inter_3d_Array_i8
    module procedure bcast_inter_Scalar_r4
    module procedure bcast_inter_1d_Array_r4
    module procedure bcast_inter_2d_Array_r4
    module procedure bcast_inter_3d_Array_r4
    module procedure bcast_inter_Scalar_r8
    module procedure bcast_inter_1d_Array_r8
    module procedure bcast_inter_2d_Array_r8
    module procedure bcast_inter_3d_Array_r8
  end interface
#endif

  interface print_var_stats
    module procedure print_var_stats_data_nwrds
    module procedure print_var_stats_data
  end interface

  interface sprint_var_stats
    module procedure sprint_var_stats_data_nwrds
    module procedure sprint_var_stats_data
  end interface

#ifdef CrayFTN
#define QUOTEMAC(x) "x"
#else
#define QUOTEMAC(x) #x
#endif
#define QUOTEMACVAL(x) QUOTEMAC(x)
  !--- COUPLER_COMMIT_ID is the checksum associated with the version
  !--- of the coupler source code that is currently running
  !--- COUPLER_REPO_PATH is the path to this coupler source code
  !--- These macros are defined during the coupler build process
#ifdef COUPLER_COMMIT_ID
  character(64),  save, public :: cpl_build_commit_id = &
      QUOTEMACVAL(COUPLER_COMMIT_ID)
#else
  character(64),  save, public :: cpl_build_commit_id = " "
#endif
#ifdef COUPLER_REPO_PATH
  character(256), save, public :: cpl_build_repo_path = &
      QUOTEMACVAL(COUPLER_REPO_PATH)
#else
  character(256), save, public :: cpl_build_repo_path = " "
#endif

  !--- NEMO_COMMIT is the checksum associated with the version
  !--- of the nemo source code that is currently running
  !--- NEMO_REPO_PATH is the path to this nemo source code
  !--- These macros are defined during the nemo build process
#ifdef NEMO_COMMIT
  character(64),  save, public :: nemo_commit = &
      QUOTEMACVAL(NEMO_COMMIT)
#else
  character(64),  save, public :: nemo_commit = " "
#endif
#ifdef NEMO_REPO_PATH
  character(256), save, public :: nemo_repo_path = &
      QUOTEMACVAL(NEMO_REPO_PATH)
#else
  character(256), save, public :: nemo_repo_path = " "
#endif
#ifdef NEMO_CONFIG
  character(256), save, public :: nemo_config = &
      QUOTEMACVAL(NEMO_CONFIG)
#else
  character(256), save, public :: nemo_config = " "
#endif
  !--- The CanCPL commit and repo used to extract com_cpl.F90 during the NEMO build
#ifdef NEMO_CANCPL_COMMIT
  character(64),  save, public :: nemo_build_cancpl_commit = &
      QUOTEMACVAL(NEMO_CANCPL_COMMIT)
#else
  character(64),  save, public :: nemo_build_cancpl_commit = " "
#endif
#ifdef NEMO_CANCPL_REPO_PATH
  character(256), save, public :: nemo_build_cancpl_repo_path = &
      QUOTEMACVAL(NEMO_CANCPL_REPO_PATH)
#else
  character(256), save, public :: nemo_build_cancpl_repo_path = " "
#endif

  character(64),  save, public :: cpl_runtime_commit = " "
  character(256), save, public :: cpl_runtime_repo_path = " "
  character(64),  save, public :: nemo_runtime_commit = " "
  character(256), save, public :: nemo_runtime_repo_path = " "
  character(256), save, public :: nemo_runtime_config = " "

  !--- specified_bc_file is the name of file containing specified boundary conditions
  !--- In the case of multi-year boundary conditions this will be a prefix
  !--- If specified_bc_file is blank then read from the "AN" file

  !--- specified_bc is an integer flag to indicate the type of boundary data
  !--- to be read from the file named by the variable specified_bc_file
  !--- specified_bc = 1 ...Monthly climatology
  !--- specified_bc = 2 ...Daily climatology
  !--- specified_bc = 3 ...multi-year monthly data
  !--- specified_bc = 4 ...multi-year daily data

  !--- atm_forcing_from_file is a logical variable used in field_ops

  !--- Define default values for specified_bc, specified_bc_file and
  !--- atm_forcing_from_file based on the value of certain cpp tokens
#ifdef coupler_specified_bc
  integer, save, public :: specified_bc = coupler_specified_bc
# if coupler_specified_bc > 0
  logical, save, public :: atm_forcing_from_file=.true.
#  ifdef coupler_specified_bc_file
  !--- If the cpp macro coupler_specified_bc_file has a value then
  !--- assume it is the name of a file containing climatolgies
  !--- for GT, SIC and SICN.
  character(len=256), save, public :: specified_bc_file = QUOTEMACVAL(coupler_specified_bc_file)
#  else
  !--- Read climatologies for GT, SIC and SICN from the "AN" file
  character(len=256), save, public :: specified_bc_file = "AN"
#  endif
# else
  character(len=256), save, public :: specified_bc_file = ' '
  logical, save, public :: atm_forcing_from_file=.false.
# endif
#else
  integer, save, public :: specified_bc = 0
  character(len=256), save, public :: specified_bc_file = ' '
  logical, save, public :: atm_forcing_from_file=.false.
#endif

  !--- _IMPI_ is the word size used for integer parameters
  !--- passed in calls to mpi routines
#ifndef _IMPI_
#  define _IMPI_ 4
#endif

  !--- Word size for integer parameters used in mpi calls
  integer, parameter :: impi = _IMPI_

  !--- Define the horizontal resolution (nx,ny) for all component models
  !--- These dimension sizes DO NOT include any overlap or cyclic longitude
#ifndef _PAR_NLON_A
#  define _PAR_NLON_A 128
#endif

#ifndef _PAR_NLAT_A
#  define _PAR_NLAT_A  64
#endif

#ifndef _PAR_NLON_O
    !--- Using ORCA1 dimension sizes nlon_o = 360
#  define _PAR_NLON_O 360
#endif

#ifndef _PAR_NLAT_O
  !--- Using ORCA1 dimension sizes nlat_o = 292
#  define _PAR_NLAT_O 292
#endif

#ifndef _PAR_NLON_I
#  define _PAR_NLON_I 360
#endif
#ifndef _PAR_NLAT_I
#  define _PAR_NLAT_I 292
#endif

  integer, parameter :: nlon_a = _PAR_NLON_A
  integer, parameter :: nlat_a = _PAR_NLAT_A
  integer, parameter :: olap_a = 1
  integer, parameter :: nlon_o = _PAR_NLON_O
  integer, parameter :: nlat_o = _PAR_NLAT_O
  integer, parameter :: olap_o = 2
  integer, parameter :: nlon_i = _PAR_NLON_I
  integer, parameter :: nlat_i = _PAR_NLAT_I
  integer, parameter :: olap_i = 0

  !--- Temporary definitions
  integer, parameter :: nlon_canom = 256
  integer, parameter :: nlat_canom = 192
  integer, parameter :: olap_canom = 1
  integer, parameter :: nwds_canom = (nlon_canom+olap_canom) * nlat_canom

  !--- nwrds_a is the number of elements on the agcm grid
  !--- that will be passed in a single mpi send or receive
  !--- This includes any overlap longitude
  integer, parameter :: nwds_a = (nlon_a+olap_a) * nlat_a

  !--- nwrds_o is the number of elements on the ocean grid
  !--- that will be passed in a single mpi send or receive
  !--- This includes any overlap longitudes
  integer, parameter :: nwds_o = (nlon_o+olap_o) * nlat_o

  !--- nwrds_i is the number of elements on the ice grid
  !--- that will be passed in a single mpi send or receive
  !--- This includes any overlap longitudes
  integer, parameter :: nwds_i = (nlon_i+olap_i) * nlat_i

  !--- nwrds_d is the number of words that will be passed in a
  !--- single mpi send or receive when using the MSG array
  integer, parameter :: nwds_d = 8

  !--- Default coupling frequencies (in seconds) for each component model
  integer(8), save, public :: cpl_atm_freq = -1 !---  3hrs
  integer(8), save, public :: cpl_ocn_freq = -1 !---  3hrs
  integer(8), save, public :: cpl_ice_freq = -1 !--- 24hrs

  !--- Group IDs for each component model and other misc groups
  integer, parameter, public :: cpl_gid = 1
  integer, parameter, public :: ocn_gid = 2
  integer, parameter, public :: atm_gid = 3
  integer, parameter, public :: ice_gid = 4
  integer, parameter, public :: msg_gid = 5
  integer, parameter, public :: canom_gid = 6
  integer, parameter, public :: other_gid = 7

  !--- CTEM dimension sizes
  integer, parameter :: ictem=9 ! number of CTEM PFTs
  integer, parameter ::  ignd=3 ! number of SOIL LAYERS
  integer, parameter ::  ican=4 ! number of CLASS PFTs

  !--- Avoid extra calls to mpi_comm_rank
  integer, public, save :: curr_task = 0

  !--- A derived type used to hold details of a data transfer between
  !--- two component modules, sent through the coupler
  type event_t
    character(32) ::        get_model = " "
    character(32) ::        put_model = " "
    character(32) ::          get_var = " "
    character(32) ::          put_var = " "
    character(32) ::    regrid_method = " "
    character(32) ::  regrid_norm_opt = " "
    character(64) ::           action = " "
    integer       ::        freq = -1
    integer       ::         tag = -1
    integer       ::       ncopy = -1
    integer       ::  get_master = -1
    integer       ::  put_master = -1
    integer       ::   regrid_id = -1
    integer       :: alarm_index = -1
    logical       :: is_msg = .false.
  end type event_t

  !--- A list of data transfer events between component models
  integer, save :: nevents = 0
  type(event_t), pointer, save :: event(:)

  !--- A list of all field names that take part in inter-model transfers
  type uniq_names_t
    character(32) :: name
    character(128) :: grid
  end type uniq_names_t
  type(uniq_names_t), pointer, save :: uniq_names(:)

  !--- A derived type used to hold information about the type of remapping
  !--- that should be done for a particular variable
  type regrid_info_t
    !--- The type of interpolation to be done
    !--- Possible values are:
    !---     cv, conservative
    !---     bl, bilinear
    !---     bc, bicubic
    !---     nn, distwgt
    character(64) :: method = 'conservative'

    !--- The type of normalization to be done for a conservative remap
    !--- Possible values are (case insensitive):
    !---     fracArea
    !---     destArea
    character(64) :: norm_opt = 'fracarea'

    !--- Logical flag to indicate whether or not land/ocean masks supplied
    !--- in the grid description file will be used during the creation
    !--- of interpolation weights and associated map info
    !--- masked = T means mask found in the grid desc file will be used
    !--- masked = F means mask will not be used (ie all grid points will
    !---            participate in the weights calculation)
    logical :: masked = .false.

  end type regrid_info_t

  type(regrid_info_t), save :: &
    default_regrid = regrid_info_t('conservative', 'fracarea', .false.)

  !--- A derived type used to hold information associated with a variable
  !--- that is passed via MPI at run time
  type cpl_vinfo_t
    !--- A unique name used to identify this record
    character(len=32) :: name = " "

    !--- A string to identify the horizontal grid associated with this variable
    character(len=128) :: grid = " "

    !--- Integer tag used in mpi calls (tags must be unique for all records)
    integer(kind=impi) :: tag = -1

    !--- The number of words passed during each MPI send or receive
    !--- The size is defined in terms of nlon, nlat and olap as follows
    !---    size == (nlon+olap) * nlat
    integer :: size = -1

    !--- The number of longitudes in the associated grid
    !--- This will NOT include any overlap (cyclic) longitude
    integer :: nlon = -1

    !--- The number of latitudes in the associated grid
    integer :: nlat = -1

    !--- The number of overlapping longitudes
    integer :: olap = 0

    !--- The number of instances of this field to be passed sequentially
    integer :: ncopy = -1

    !--- The coupling frequency in seconds
    integer :: freq = -1

    !--- The 4 character CCCma variable name (ie ibuf3 as char)
    character(len=4) :: ccc_name = " "

    !--- Infomation about regridding this variable
    type(regrid_info_t) :: regrid = regrid_info_t('conservative', 'fracarea', .false.)

  end type cpl_vinfo_t

  type(regrid_info_t), pointer, save :: null_regrid => null()

  !--- An array of type cpl_vinfo_t to contain info about all known variables.
  !--- This information was originally based on info found in the taginfo.h file
  !--- but extra elements have been added to improve functionality.
  !--- NOTE: Every "name" (ie the first entry of each cpl_vinfo_t type assignment)
  !---       as well as every "tag" found in cpl_var_list MUST be unique.
  type(cpl_vinfo_t), pointer, save :: cpl_var_list(:)

  !--- The list of variables that will be sent from the agcm to the ocean
  !--- through the coupler
  integer :: n_atm_ocn_event = 0
  type(event_t), save :: atm_ocn_event(100)

  !--- The list of variables that will be sent from the ocean to the agcm
  !--- through the coupler
  integer :: n_ocn_atm_event = 0
  type(event_t), save :: ocn_atm_event(100)

  !--- The list of variables that will be sent from the ice model to the agcm
  !--- through the coupler
  integer :: n_ice_atm_event = 0
  type(event_t), save :: ice_atm_event(100)

  !--- The list of variables that will be sent from the agcm to the ice model
  !--- through the coupler
  integer :: n_atm_ice_event = 0
  type(event_t), save :: atm_ice_event(100)

  !--- The list of variables that will be sent from the ice model to the ocean
  !--- through the coupler
  integer :: n_ice_ocn_event = 0
  type(event_t), save :: ice_ocn_event(100)

  !--- The list of variables that will be sent from the ocean to the ice model
  !--- through the coupler
  integer :: n_ocn_ice_event = 0
  type(event_t), save :: ocn_ice_event(100)

  integer(kind=impi) :: ierr
  integer(kind=impi) :: curr_real_type = -1

  !--- nbuf is the size of the integer array that is prepended
  !--- to every MPI data transmission done below
  !--- In the past this has always been the 8 word ibuf array
  !--- but this could change in the future
  integer(kind=impi), parameter :: nbuf=8

  !--- The master task associated with each component model
  !--- The *_master variables indicate the rank in MPI_COMM_WORLD not the group rank
  integer(kind=impi), public, save :: cpl_master=-1
  integer(kind=impi), public, save :: atm_master=-1
  integer(kind=impi), public, save :: ocn_master=-1
  integer(kind=impi), public, save :: ice_master=-1

  type task_info_t
    integer(kind=4) :: master
    character(32)   :: group_name
    integer(kind=4) :: group_min_rank
    integer(kind=4) :: group_max_rank
  end type task_info_t

  type(task_info_t), pointer, save :: task_info(:)

  !--- model_group_comm is the group intra-communicator returned by mpi_comm_split
  integer(kind=impi), public, save :: model_group_comm = -1

  !--- Inter-communicators used for group to group exchanges
  integer(kind=impi), public, save :: intercomm_cpl_atm = -1
  integer(kind=impi), public, save :: intercomm_cpl_ocn = -1
  integer(kind=impi), public, save :: intercomm_atm_ocn = -1

  logical, save :: com_mpi_initialized = .false.

  real(kind=8), parameter, public :: fillvalue = 1e38
  real(kind=8), parameter, public :: pi_r8 = 3.141592653589793_8
  !--- The current coupler time as a string
  character(20), public, save :: cpl_time_string = " "

  !--- The total coupler elasped time in seconds
  integer(8), public, save :: cpl_elapsed_time_secs = 0_8

  !--- The coupler start time in seconds
  integer(8), public, save :: cpl_start_time_secs = 0_8

  !--- The coupler time step in seconds
  integer(8), public, save :: cpl_time_step_secs = 0_8

  !--- The current iteration of the coupler time step
  integer(8), public, save :: cpl_time_step_iter = 0_8

  !--- The current coupler year
  integer(4), public, save :: cpl_year = 1

  !--- The current coupler day of the year
  integer(4), public, save :: cpl_day_of_year = 1

  !--- nemo_nn_ice identifies the ice model used in NEMO
  !--- nemo_nn_ice = 0 ...no ice boundary condition
  !---             = 1 ...use observed ice-cover
  !---             = 2 ...lim2 ice-model used
  !---             = 3 ...lim3 ice-model used
  !---             = 4 ...cice ice-model used
  !--- nemo_nn_ice = nn_ice (a namelist parameter read in nemo)
  !--- It is defined in the nemo subroutine cpl_cancpl_init and broadcast to all
  !--- of MPI_COMM_WORLD when the subroutine define_group (defined below) is
  !--- called by component models
  !--- This variable must be typed with a specific size (it cannot be default integer)
  !--- since default integer could be of a different size on different mpi processes
  !--- and this causes problems when communicating between tasks (e.g. broadcasting)
  integer(kind=4), public :: nemo_nn_ice = -1

  !--- nemo_ncat is the number of ice catagories used by CICE
  integer(kind=4), public :: nemo_ncat = -1

  !--- nemo_rtd is the NEMO dynamics time step in seconds
  integer(kind=8), public :: nemo_rn_rdt = -1

  !--- nemo_nn_fsbc is the number of NEMO time steps between calls to the surface boundary
  !--- condition routines, including the sea ice model and coupling with the atmosphere
  integer(kind=8), public :: nemo_nn_fsbc = -1

  !--- nemo_nn_it000 is the index of the initial time step used by nemo
  integer(kind=8), public :: nemo_nn_it000 = -1

  !--- nemo_nn_itend is the index of the final time step used by nemo
  integer(kind=8), public :: nemo_nn_itend = -1

  !--- nemo_nn_date0 is the initial calendar date used by nemo
  integer(kind=8), public :: nemo_nn_date0 = -1

  !--- (nemo_jpiglo, nemo_jpjglo) is the shape of the NEMO global domain
  integer(kind=8), public :: nemo_jpiglo = -1
  integer(kind=8), public :: nemo_jpjglo = -1

  !--- nemo_namsbc_cpl_cldes will contain the description (aka cldes) associated
  !--- with variables that are sent/received by NEMO. These descriptions are
  !--- read from the namelist namsbc_cpl in the subroutine sbc_cpl_init
  !--- Currently the elements of nemo_namsbc_cpl_cldes are as follows
  !---   1) sn_snd_temp % cldes
  !---   2) sn_snd_alb % cldes
  !---   3) sn_snd_thick % cldes
  !---   4) sn_snd_crt % cldes
  !---   5) sn_snd_co2 % cldes
  !---   6) sn_rcv_w10m % cldes
  !---   7) sn_rcv_taumod % cldes
  !---   8) sn_rcv_tau % cldes
  !---   9) sn_rcv_dqnsdt % cldes
  !---  10) sn_rcv_qsr % cldes
  !---  11) sn_rcv_qns % cldes
  !---  12) sn_rcv_emp % cldes
  !---  13) sn_rcv_rnf % cldes
  !---  14) sn_rcv_cal % cldes
  !---  15) sn_rcv_iceflx % cldes
  !---  16) sn_rcv_co2 % cldes
  character(32), public ::  nemo_namsbc_cpl_cldes(16)

  !--- nemo_glamt, nemo_glamu, nemo_glamv, nemo_glamf are cell center longitudes 
  !--- (in degrees) on the global NEMO T,U,V,F grids respectively
  real(kind=8), pointer, public :: nemo_glamt(:,:)
  real(kind=8), pointer, public :: nemo_glamu(:,:)
  real(kind=8), pointer, public :: nemo_glamv(:,:)
  real(kind=8), pointer, public :: nemo_glamf(:,:)

  !--- nemo_gphit, nemo_gphiu, nemo_gphiv, nemo_gphif are cell center latitudes 
  !--- (in degrees) on the global NEMO T,U,V,F grids respectively
  real(kind=8), pointer, public :: nemo_gphit(:,:)
  real(kind=8), pointer, public :: nemo_gphiu(:,:)
  real(kind=8), pointer, public :: nemo_gphiv(:,:)
  real(kind=8), pointer, public :: nemo_gphif(:,:)

  !--- nemo_e1t, nemo_e1u, nemo_e1v, nemo_e1f
  !--- nemo_e2t, nemo_e2u, nemo_e2v, nemo_e2f
  !--- are horizontal scale factors on the global NEMO T,U,V,F grids respectively
  real(kind=8), pointer, public :: nemo_e1t(:,:)
  real(kind=8), pointer, public :: nemo_e1u(:,:)
  real(kind=8), pointer, public :: nemo_e1v(:,:)
  real(kind=8), pointer, public :: nemo_e1f(:,:)
  real(kind=8), pointer, public :: nemo_e2t(:,:)
  real(kind=8), pointer, public :: nemo_e2u(:,:)
  real(kind=8), pointer, public :: nemo_e2v(:,:)
  real(kind=8), pointer, public :: nemo_e2f(:,:)

  !--- nemo_tmask, nemo_umask, nemo_vmask, nemo_fmask
  !--- are the global ocean/land mask on the NEMO T,U,V,F grids respectively
  real(kind=8), pointer, public :: nemo_tmask(:,:)
  real(kind=8), pointer, public :: nemo_umask(:,:)
  real(kind=8), pointer, public :: nemo_vmask(:,:)
  real(kind=8), pointer, public :: nemo_fmask(:,:)

  !--- A list of variable names sent from NEMO to the coupler
  !--- This list is defined at run time in nemo subroutine cpl_cancpl_define
  integer(kind=4), save, public :: nemo_n_send_var=0
  character(32),   pointer, public :: nemo_send_var(:)

  !--- A list of variable names sent from the coupler to NEMO
  !--- This list is defined at run time in nemo subroutine cpl_cancpl_define
  integer(kind=4), save, public :: nemo_n_recv_var=0
  character(32),   pointer, public :: nemo_recv_var(:)

  !--- A list of variable names sent from the AGCM to the coupler
  !--- This list is defined at run time in the subroutine define_atm_send_var
  integer(kind=4), save, public :: atm_n_send_var=0
  character(32),   pointer, public :: atm_send_var(:)

  !--- A list of variable names sent from the coupler to AGCM
  !--- This list is defined at run time in the subroutine define_atm_recv_var
  integer(kind=4), save, public :: atm_n_recv_var=0
  character(32),   pointer, public :: atm_recv_var(:)

  !--- Ocean grid cell area
  real(kind=8), pointer, public :: ocn_grid_cell_area(:,:)

  !--- Atmosphere grid cell area
  real(kind=8), pointer, public :: atm_grid_cell_area(:,:)

  !--- The radius of the earth (meters) according to the agcm
  real(kind=8), public :: atm_earth_radius = -1.0_8

  !--- Number of atm longitudes as calculated in the AGCM at runtime
  !--- atm_nlon must be identical to nlon_a defined above
  integer(kind=8), public :: atm_nlon = -1

  !--- Number of atm latitudes as calculated in the AGCM at runtime
  !--- atm_nlat must be identical to nlat_a defined above
  integer(kind=8), public :: atm_nlat = -1

  !--- atm longitudes (degrees east)
  real(kind=8), pointer, public :: atm_lon(:)

  !--- atm latitudes (degrees north)
  real(kind=8), pointer, public :: atm_lat(:)

  !--- atm Gaussian weights
  real(kind=8), pointer, public :: atm_gwt(:)

  !--- atm_kount is the time step index used in the AGCM
  integer(kind=8), public :: atm_kount = -1

  !--- atm_kstartc is the variable KSTARTC used in the AGCM and sent to the coupler
  !--- at the end of each iteration of the ICOUNT loop in the AGCM driver via MSG
  integer(kind=8), public :: atm_kstartc = -1

  !--- atm_kount2 is the time step index used in the AGCM and sent to the coupler
  !--- at the end of each iteration of the ICOUNT loop in the AGCM driver via MSG
  integer(kind=8), public :: atm_kount2 = -1

  !--- atm_kstart is the initial time step index used in the AGCM
  integer(kind=8), public :: atm_kstart = -1

  !--- atm_ksteps is the number of AGCM time steps in each coupling interval
  integer(kind=8), public :: atm_ksteps = -1

  !--- atm_kfinal is the last time step index used in the AGCM
  integer(kind=8), public :: atm_kfinal = -1

  !--- atm_delt is the time step in seconds that is used in the agcm
  integer(kind=8), public :: atm_delt = -1

  !--- atm_lmask is a binary land mask on the AGCM grid (1=some land, 0=no land)
  real(kind=8), pointer, public :: atm_lmask(:)

  !--- atm_fland is the fractional land mask that is used in the agcm (0..1)
  real(kind=8), pointer, public :: atm_fland(:)

  !--- atm_cv_lmask is a "conservative" land mask the is 0 over cells that are
  !--- all or partial water and 2 over cells that are land only, determined
  !--- using the fractional land mask found in atm_fland
  real(kind=8), pointer, public :: atm_cv_lmask(:)

  !--- Define the coupler_par namelist

  !--- The coupler restart file names
  character(128), public :: cpl_rs_file_name_in   = "cpl_restart.nc"
  character(128), public :: cpl_rs_file_name_out  = "cpl_restart_out.nc"

  !--- The coupler history file name
  character(128), public :: cpl_hist_file_name  = "cpl_history.nc"

  !--- The name of an existing file containing coupler history data
  !--- This data may be used to overwrite fields that are sent to a component model
  !--- (ie ocean or agcm) ...currently only used for debugging
  character(128), public :: cpl_hist_file_name_to_read  = " "

  !--- Abort, or not, if the input coupler resart file is missing a field
  !--- If cpl_rs_abort_if_missing_field is false then a field of zeros will be supplied
  !--- when a missing field is encountered
  logical, public :: cpl_rs_abort_if_missing_field = .true.

  !--- When cpl_use_initial_conditions is true create an initial coupler restart file
  !--- if the initial restart file does not exist
  !--- Only do this at the start of a run
  logical, public :: cpl_use_initial_conditions = .false.

  !--- ..._grid_desc variables define names of files that contain the grid
  !--- description for individual component model
  character(128), public :: atm_grid_desc = " "
  character(128), public :: ocn_grid_desc = " "
  character(128), public :: ice_grid_desc = " "
  character(128), public :: canom4_grid_desc = " "

  !--- The name of a file containing NEMO "mesh_mask" data
  character(128), public :: nemo_mesh_mask_file = " "

  !--- The name of a file containing AGCM grid cell area data
  character(128), public :: agcm_grid_area_file = " "

  !--- ..._grid_idx variables, if > 0, are the index in grid_desc_list (defined
  !--- in module grid_desc) that contains grid information for each component
  integer(kind=4), public :: atm_grid_idx   = -1
  integer(kind=4), public :: ocn_grid_idx   = -1
  integer(kind=4), public :: ice_grid_idx   = -1
  integer(kind=4), public :: canom4_grid_idx = -1

  !--- The name of a CanESM4 coupler restart file that may be used to bootstrap
  !--- creation of a CanCPL restart file (used for development only)
  character(128), public :: canesm4_cpl_restart = " "

  !--- If the following variables are defined in the environment at run time
  !--- then they written to the namelist prior to execution
  character(128), public :: env_model = " "
  character(128), public :: env_start = " "
  character(128), public :: env_runid = " "

  !--- env_months will contain the value of the parmsub parameter "months" if it is
  !--- defined in the environment at run time
  integer(kind=4), public :: env_months = 1

  integer(kind=4), public :: env_run_start_year  = -1
  integer(kind=4), public :: env_run_start_month = -1
  integer(kind=4), public :: env_run_stop_year   = -1
  integer(kind=4), public :: env_run_stop_month  = -1

  !--- Is a bulk formulation used in the coupler to define NEMO fields
  logical, public :: bulk_in_cpl = .false.

  !--- Should O_Runoff be zeroed out before it is sent to NEMO
  logical, public :: zero_runoff_sent_to_ocean = .false.

  !--- couple_serial = T means run in serial mode, otherwise run in parallel mode
  !--- When coupled to an ocean and an atmosphere we run in serial mode where the
  !--- order of communication is
  !---     send to atm, recv from atm, send to ocn, recv from ocn
  !--- or in parallel mode, in which case the order of communication is determined
  !--- by the parameter coupler_parallel_type
  logical, public :: couple_serial = .false.

  !--- When coupled in parallel mode coupler_parallel_type determines the order
  !--- of communication between the coupler and component models
  !---   coupler_parallel_type = 1
  !---      ... send to atm, send to ocn, recv from atm, recv from ocn
  !---   coupler_parallel_type = 2
  !---      ... send to atm, send to ocn, recv from ocn, recv from atm
  integer(kind=4), public :: coupler_parallel_type = 2

  !--- A value to add to the year when reading AGCM specified boundary conditions
  !--- (GT, SIC, SICN) from a file using the subroutine read_spec_bc.
  !--- This is only used when reading multi-year boundary conditions and is required to
  !--- read files with data for years different from coupler time.
  !--- specified_bc_year_offset is added to the current year, as determined by the coupler,
  !--- so a negative value will be subtracted
  integer(kind=8), public :: specified_bc_year_offset = 0

  !--- The number of years to loop multi-year specified boundary conditions that are read
  !--- from a file in the subroutine read_spec_bc.
  !--- FIXME: NOT YET SUPPORTED
  integer(kind=8), public :: specified_bc_years_in_loop = 0
  
  ! Options used to reproduce the CMIP6 p2 or p1 variants of CanESM5
  logical, public       :: landfrac_bug = .false.
  character(32) , public :: windstress_remap = 'bilinear'

  public :: coupler_par
  namelist /coupler_par/ atm_grid_desc, ocn_grid_desc, ice_grid_desc, &
        canom4_grid_desc, canesm4_cpl_restart, nemo_mesh_mask_file, agcm_grid_area_file, &
        cpl_use_initial_conditions, env_model, env_start, env_months, env_runid, &
        env_run_start_year, env_run_start_month, env_run_stop_year, env_run_stop_month, &
        cpl_rs_abort_if_missing_field, bulk_in_cpl, zero_runoff_sent_to_ocean, atm_forcing_from_file, &
        specified_bc, specified_bc_file, couple_serial, coupler_parallel_type, &
        specified_bc_year_offset, specified_bc_years_in_loop, cpl_hist_file_name_to_read, &
        landfrac_bug, windstress_remap

  !-----------------------------------------------------------------------------
  !--- Everything below this line is defined for backward compatibility with
  !--- previous atm/ocn models that did not have a separate coupler.
  !--- These variables are depreciated and should not be used in any new code
  !-----------------------------------------------------------------------------

  !--- member is defined here for compatability with recv_fld and send_fld
  character(len=8) :: member(0:1000)

  !--- These arrays define integer tags associated with certain records
  !--- that are passed via MPI
  !--- These arrays are defined here for backward compatibility with older
  !--- versions of the AGCM and OGCM but are depreciated.
  !--- They are replaced by cpl_var_list defined above.

  !--- MSG record
  integer(kind=8), parameter, dimension(3), public :: &
        d_msg  = (/ 900, 900, nwds_d /)

  !--- Records on the atm grid
  integer(kind=8), parameter, dimension(3), public :: &
        a_mask   = (/ 170, 170, nwds_a /), &
        a_flnd   = (/ 171, 171, nwds_a /)

  !--- Records on the CanOM ocean grid
  integer(kind=8), parameter, dimension(3), public :: &
        o_data = (/ 200, 200, nwds_o /), &
        o_gtx  = (/ 201, 201, nwds_o /), &
        o_oufs = (/ 202, 202, nwds_o /), &
        o_ovfs = (/ 203, 203, nwds_o /), &
        o_obea = (/ 204, 204, nwds_o /), &
        o_obwa = (/ 205, 205, nwds_o /), &
        o_ogt  = (/ 206, 206, nwds_o /), &
        o_slt  = (/ 207, 207, nwds_o /), &
        o_uoc  = (/ 210, 210, nwds_o /), &
        o_voc  = (/ 211, 211, nwds_o /), &
        o_cnv  = (/ 212, 212, nwds_o /), &
        o_psib = (/ 213, 213, nwds_o /), &
        o_uiso = (/ 214, 214, nwds_o /), &
        o_viso = (/ 215, 215, nwds_o /), &
        o_obei = (/ 216, 216, nwds_o /), &
        o_odpx = (/ 217, 217, nwds_o /), &
        o_odpy = (/ 218, 218, nwds_o /), &
        o_tr01 = (/ 219, 219, nwds_o /), &
        o_sst2 = (/ 220, 220, nwds_o /), &
        o_hblm = (/ 221, 221, nwds_o /)

  contains

    !***************************************************************************
    !--- Normal exit or abort depending on the input value of N
    !---   -100 <= N <= -1  ...abort with non-zero exit status
    !---   otherwise normal exit
    !***************************************************************************
    subroutine err_exit(name,n)
      implicit none

      character(*) :: name
      integer, intent(in) :: n

      integer :: iou, lename
      integer(4) :: myrank, ierr, n4

      character(80), parameter :: dash = &
      '--------------------------------------------------------------------------------'
      character(80), parameter :: star = &
      '********************************************************************************'

#ifdef use_mpi
      call mpi_comm_rank ( MPI_COMM_WORLD, myrank, ierr )
#else
      myrank = 0
#endif
      if (myrank.eq.0 .or. (n.lt.0 .and. n.ge.-100) ) then
        lename=len_trim(name)
        if (lename.lt.1) lename=1
  
        if(n.ge.0) write(6,'(a,"  END  ",a,1x,a,i8)') &
          dash(1:8),name(1:lename),dash(1:81-lename),n

        if(n.lt.0) write(6,'(a,"  END  ",a,1x,a,i8)') &
          star(1:8),name(1:lename),star(1:81-lename),n

      endif
      if (myrank.eq.0) then
       print *, '0closing data i/o units'
       do iou=1,100
        if (iou.ne.5.and.iou.ne.6) close (iou)
       enddo
      endif
      call flush(6)

      if ( n.ge.0 .or. n.lt.-100 ) then

#ifdef use_mpi
        call system('sleep 3')
        print *, '0got to before mpi_finalize'
        call flush(6)
        call mpi_finalize(ierr)
        print *, '0got to after mpi_finalize'
        call flush(6)
#endif
        stop 
      else
#ifdef use_mpi
        n4=n
        call system('sleep 3')
        call mpi_abort(model_group_comm,n4,ierr)
#endif
        call abort
      endif
    end subroutine err_exit  

    !***************************************************************************
    !--- Convert all upper case A-Z characters in the input string to lower case
    !--- This function is also in the strings module but is included here
    !--- to minimize dependencies for this module
    !***************************************************************************
    function lowerc(strng) result(lc)
      character(*), intent(in) :: strng
      character( len=1024 ) :: lc
      integer :: idx, ich

      if ( len_trim(strng) > 1024 ) then
        write(6,*)'lowerc: Input string to long.'
        call err_exit("lowerc",-1)
      endif

      lc = " "
      lc = trim(strng)

      do idx=1,len_trim(strng)
        ich = iachar( lc(idx:idx) )
        if ( ich > 64 .and. ich < 91 ) then
          !--- This is an upper case character
          !--- Convert it to lower case
          lc(idx:idx) = achar( ich+32 )
        endif
      enddo

    end function lowerc

    !***************************************************************************
    !--- Convert all lower case a-z characters in the input string to upper case
    !--- This function is also in the strings module but is included here
    !--- to minimize dependencies for this module
    !***************************************************************************
    function upperc(strng) result(uc)
      character(*), intent(in) :: strng
      character( len=1024 ) :: uc
      integer :: idx, ich

      if ( len_trim(strng) > 1024 ) then
        write(6,*)'lowerc: Input string to long.'
        call err_exit("lowerc",-1)
      endif

      uc = " "
      uc = trim(strng)

      do idx=1,len_trim(strng)
        ich = iachar( uc(idx:idx) )
        if ( ich > 96 .and. ich < 123 ) then
          !--- This is a lower case character
          !--- Convert it to upper case
          uc(idx:idx) = achar( ich-32 )
        endif
      enddo

    end function upperc

    !***************************************************************************
    !--- Determine the number of overlap longitudes given the x and y grid
    !--- dimension (nx,ny) sizes as well as the total size of the transfer (sz)
    !***************************************************************************
    function overlap(nx,ny,sz) result(nolap)
      integer :: nx, ny, sz
      integer :: nolap

      if ( nx*ny == sz ) then
        !--- No overlap longitude
        nolap = 0
      else if ( (nx+1)*ny == sz ) then
        !--- One overlap longitude
        nolap = 1
      else if ( (nx+2)*ny == sz ) then
        !--- Two overlap longitudes
        nolap = 2
      else
        !--- Invalid
        nolap = -1
      endif

    end function overlap

    !***********************************************************************
    !--- Calculate the x and y indicies of a 2D array in column-major order
    !--- using
    !---   nx   ...the size of the first dimension of the 2D array
    !---   addr ...the index in a 1D array corresponding with the sequential
    !---           layout in memory of the 2D array
    !***********************************************************************
    subroutine idx2d_from_1d(ix, iy, addr, nx)
      implicit none
      integer, intent(out) :: ix, iy
      integer, intent(in)  :: addr, nx

      ix = 1 + mod(addr-1,nx)
      iy = (addr-1)/nx + 1
    end subroutine idx2d_from_1d

    !***************************************************************************
    !--- Print min,max,avg (and possibly other stats) for a given data array
    !***************************************************************************
    function sprint_var_stats_data_nwrds(var_data, nwrds, name, pfx) result(str)
      real(kind=8), intent(in) :: var_data(:)
      integer, intent(in) :: nwrds
      character(*), intent(in), optional :: name
      character(*), intent(in), optional :: pfx
      character(512) :: str

      if ( present(name) .and. present(pfx) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, name=trim(name), pfx=trim(pfx), str_out=str)
      else if ( present(name) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, name=trim(name), str_out=str)
      else if ( present(pfx) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, pfx=trim(pfx), str_out=str)
      else
        call print_var_stats_data_nwrds(var_data, nwrds, str_out=str)
      endif

    end function sprint_var_stats_data_nwrds

    function sprint_var_stats_data(var_data, name, pfx) result(str)
      real(kind=8), intent(in) :: var_data(:)
      character(*), intent(in), optional :: name
      character(*), intent(in), optional :: pfx
      character(512) :: str

      !--- Local
      integer :: nwrds

      !--- When nwrds is not supplied by the user, assume
      !--- the entire input array is to be used
      nwrds = size(var_data)

      if ( present(name) .and. present(pfx) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, name=trim(name), pfx=trim(pfx), str_out=str)
      else if ( present(name) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, name=trim(name), str_out=str)
      else if ( present(pfx) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, pfx=trim(pfx), str_out=str)
      else
        call print_var_stats_data_nwrds(var_data, nwrds, str_out=str)
      endif

    end function sprint_var_stats_data

    subroutine print_var_stats_data(var_data, name, pfx, gavg)
      real(kind=8), intent(in) :: var_data(:)
      character(*), intent(in), optional :: name
      character(*), intent(in), optional :: pfx
      real(kind=8), intent(in), optional :: gavg

      !--- Local
      integer :: nwrds

      !--- When nwrds is not supplied by the user, assume
      !--- the entire input array is to be used
      nwrds = size(var_data)

      if ( present(name) .and. present(pfx) .and. present(gavg) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, name=trim(name), pfx=trim(pfx), gavg=gavg)
      else if ( present(name) .and. present(pfx) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, name=trim(name), pfx=trim(pfx))
      else if ( present(name) .and. present(gavg) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, name=trim(name), gavg=gavg)
      else if ( present(pfx) .and. present(gavg) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, pfx=trim(pfx), gavg=gavg)
      else if ( present(name) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, name=trim(name))
      else if ( present(pfx) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, pfx=trim(pfx))
      else if ( present(gavg) ) then
        call print_var_stats_data_nwrds(var_data, nwrds, gavg=gavg)
      else
        call print_var_stats_data_nwrds(var_data, nwrds)
      endif

    end subroutine print_var_stats_data

    subroutine print_var_stats_data_nwrds(var_data, nwrds, name, pfx, gavg, str_out)
      real(kind=8), intent(in) :: var_data(:)
      integer, intent(in) :: nwrds
      character(*), intent(in), optional :: name
      character(*), intent(in), optional :: pfx
      real(kind=8), intent(in), optional :: gavg
      character(*), intent(inout), optional :: str_out

      !--- Local
      type(cpl_vinfo_t) :: curr_vinfo
      real(kind=8), pointer :: tmp_data(:)
      real(kind=8) :: vmin, vmax, vavg
      character(64) :: strng
      character(512) :: strng1_out, strng2_out
      character(20) :: zero_str, nan_str, spval_str, vmin_idx_str, vmax_idx_str
      logical :: avg_only_non_zero = .false.
      logical :: add_min_max_index
      integer :: spval_count, nan_count, zero_count, navg, ic1
      integer :: vmin_idx(1), vmax_idx(1), nx, ix, iy, status
      integer :: verbose=1

      !--- add_min_max_index = T means add min and max array index values to the output
      add_min_max_index = .true.

      if ( nwrds < 1 ) then
        if ( verbose > 0 ) then
          write(6,*)"print_var_stats: **WW** nwrds = ",nwrds
        endif
        return
      endif

      if ( nwrds > size(var_data) ) then
        write(6,*)"print_var_stats: Requesting ",nwrds, &
            " words from an array containing only ",size(var_data)," words."
        call err_exit("print_var_stats",-1)
      endif

      !--- Ignore elements containing fillvalue or NaN
      !--- NOTE: By definition NaN is not equal to anything, even itself
      nan_count = count( var_data(1:nwrds) /= var_data(1:nwrds) )
      spval_count = count( var_data(1:nwrds) == fillvalue )

      !--- Also keep track of the number of zeros in the data
      zero_count = count( var_data(1:nwrds) == 0.0_8 )

      if ( nan_count > 0 .or. spval_count > 0 ) then
        !--- Ignoring fillvalue and/or NaNs when calculating stats
        allocate( tmp_data(nwrds) )
        tmp_data(1:nwrds) = var_data(1:nwrds)

        !--- Set all bad values to the largest machine representable number
        !--- then find the minimum value
        if ( nan_count > 0 ) then
          !--- Set NANs to HUGE in the tmp array
          where ( tmp_data /= tmp_data ) tmp_data = HUGE(tmp_data)
        endif

        if ( spval_count > 0 ) then
          !--- Also set fillvalues to HUGE in the tmp array
          where ( tmp_data == fillvalue ) tmp_data = HUGE(tmp_data)
        endif
        vmin     = minval( tmp_data )
        vmin_idx = minloc( tmp_data )

        !--- Set all bad values to the smallest machine representable number
        !--- then find the maximum value
        where ( tmp_data == HUGE(tmp_data) ) tmp_data = -HUGE(tmp_data)
        vmax     = maxval( tmp_data )
        vmax_idx = maxloc( tmp_data )

        !--- Set all bad values to 0 then find the average value
        where ( tmp_data == -HUGE(tmp_data) ) tmp_data = 0.0
        if ( avg_only_non_zero ) then
          navg = nwrds - count( tmp_data(1:nwrds) == 0.0 )
        else
          navg = nwrds
        endif
        if ( navg == 0 ) then
          !--- All values are zero
          vavg = 0.0
        else
          !--- Average over non-fillvalue/non-NaN values
          vavg = sum( tmp_data )/real(navg,8)
        endif

      else
        !--- Calculate stats using all data values
        vmin     = minval( var_data(1:nwrds) )
        vmin_idx = minloc( var_data(1:nwrds) )
        vmax     = maxval( var_data(1:nwrds) )
        vmax_idx = maxloc( var_data(1:nwrds) )
        if ( avg_only_non_zero ) then
          navg = nwrds - zero_count
        else
          navg = nwrds
        endif
        if ( navg == 0 ) then
          !--- All values are zero
          vavg = 0.0
        else
          !--- Average over non-zero values
          vavg = sum( var_data(1:nwrds) )/real(navg,8)
        endif
      endif

      if ( present(gavg) ) then
        !--- Use the global average supplied by the user
        vavg = gavg
      endif

      !--- Always initialize these to blank
      strng1_out = " "
      strng2_out = " "

      if ( present(pfx) ) then
        ic1 = len_trim(pfx)
        if ( ic1 > 0 ) then
          strng = trim(pfx)
          ic1 = max(ic1,45)
          if ( present(str_out) ) then
            !--- Write the output string to a variable
            write(strng1_out,'(a)') adjustr(strng(1:ic1))
          else
            !--- Write the output string to stdout
            write(6,'(2a,$)') adjustr(strng(1:ic1))," "
          endif
        endif
      endif

      nan_str = " "
      if ( nan_count > 0 ) write(nan_str,'(a,i8)') '  NaNs=',nan_count

      zero_str = " "
      if ( zero_count > 0 ) write(zero_str,'(a,i8)') '  zeros=',zero_count

      spval_str = " "
      if ( spval_count > 0 ) write(spval_str,'(a,i8)') '  spvals=',spval_count

      vmin_idx_str = " "
      vmax_idx_str = " "
      if ( add_min_max_index ) then
        if ( present(name) ) then
          !--- Use the name supplied to determine the grid shape and then the
          !--- (ix,iy) indicies for the location of the min and max values
          curr_vinfo = find_cpl_vinfo( name=trim(adjustl(name)), status=status )
          if ( status == 0 ) then
            !--- The record was found
            nx = curr_vinfo%nlon + curr_vinfo%olap
            if ( nx > 0 ) then
              !--- The grid shape is known
              call idx2d_from_1d(ix, iy, vmin_idx(1), nx)
              write(vmin_idx_str,'(" (",i4,",",i4,")")') ix,iy
              call idx2d_from_1d(ix, iy, vmax_idx(1), nx)
              write(vmax_idx_str,'(" (",i4,",",i4,")")') ix,iy
            else
              write(vmin_idx_str,'(a,i9,a)') ' (',vmin_idx(1),')'
              write(vmax_idx_str,'(a,i9,a)') ' (',vmax_idx(1),')'
            endif
          else
            write(vmin_idx_str,'(a,i9,a)') ' (',vmin_idx(1),')'
            write(vmax_idx_str,'(a,i9,a)') ' (',vmax_idx(1),')'
          endif
        else
          write(vmin_idx_str,'(a,i9,a)') ' (',vmin_idx(1),')'
          write(vmax_idx_str,'(a,i9,a)') ' (',vmax_idx(1),')'
        endif
      endif

      if ( present(str_out) ) then
        !--- Write the output string to str_out
        write(strng2_out,'(3a,g16.8,3a,g16.8,a,g22.14,a,i8,3a)') &
              ' min',trim(adjustl(vmin_idx_str)),'=',vmin, &
              ' max',trim(adjustl(vmax_idx_str)),'=',vmax,' avg=',vavg, &
              '  size=',nwrds,trim(zero_str),trim(spval_str),trim(nan_str)
        strng1_out = trim(strng1_out)//" "//trim(strng2_out)
        if ( len(str_out) < len_trim(strng1_out) ) then
          write(6,*)"print_var_stats_data_nwrds: **WW** Input string is too short."
        endif
        ic1 = min( len_trim(strng1_out), len(str_out) )
        str_out = " "
        str_out(1:ic1) = strng1_out(1:ic1)
      else
        !--- Write the output string to stdout
        write(6,'(3a,g16.8,3a,g16.8,a,g22.14,a,i8,3a)') &
              ' min',trim(adjustl(vmin_idx_str)),'=',vmin, &
              ' max',trim(adjustl(vmax_idx_str)),'=',vmax,' avg=',vavg, &
              '  size=',nwrds,trim(zero_str),trim(spval_str),trim(nan_str)
        call flush(6)
      endif

    end subroutine print_var_stats_data_nwrds

    !***************************************************************************
    !--- Identify non-transfer events by their name
    !***************************************************************************
    function is_special_op(name) result(found)
      character(*), intent(in) :: name
      logical :: found

      found = .false.
      select case ( trim(adjustl(name)) )
        case ('start_cpl2atm')
          found = .true.
        case ('stop_cpl2atm')
          found = .true.
        case ('start_cpl2ocn')
          found = .true.
        case ('stop_cpl2ocn')
          found = .true.
        case ('start_atm2cpl')
          found = .true.
        case ('stop_atm2cpl')
          found = .true.
        case ('start_ocn2cpl')
          found = .true.
        case ('stop_ocn2cpl')
          found = .true.
      end select

    end function is_special_op

    !***************************************************************************
    !--- Reset the group ID for certain variable names that must always belong
    !--- to a particular group
    !***************************************************************************
    function filter_group_ID(name, gid_in) result(gid)
      character(*), intent(in) :: name
      integer, intent(in) :: gid_in
      integer :: gid

      gid = gid_in

      select case( trim(adjustl(name)) )

        case("GT_atm")
          if ( atm_forcing_from_file ) then
            !--- Force atm group when reading AGCM bcs from a file
            gid = atm_gid
          endif

        case("SIC_atm")
          if ( atm_forcing_from_file ) then
            !--- Force atm group when reading AGCM bcs from a file
            gid = atm_gid
          endif

        case("SICN_atm")
          if ( atm_forcing_from_file ) then
            !--- Force atm group when reading AGCM bcs from a file
            gid = atm_gid
          endif

        case("SNO_atm")
          if ( atm_forcing_from_file ) then
            !--- Force atm group when reading AGCM bcs from a file
            gid = atm_gid
          endif

        case("RES_atm")
          !--- RES_atm will always be in the agcm group
          gid = atm_gid
      end select

    end function filter_group_ID

    !***************************************************************************
    !--- Copy data between 2 event_t data types
    !***************************************************************************
    subroutine copy_event_t(from, to)
      type(event_t) :: from, to
      to%get_model       = from%get_model
      to%put_model       = from%put_model
      to%get_var         = from%get_var
      to%put_var         = from%put_var
      to%regrid_method   = from%regrid_method
      to%regrid_norm_opt = from%regrid_norm_opt
      to%action          = from%action
      to%freq            = from%freq
      to%tag             = from%tag
      to%ncopy           = from%ncopy
      to%get_master      = from%get_master
      to%put_master      = from%put_master
      to%regrid_id       = from%regrid_id
      to%alarm_index     = from%alarm_index
      to%is_msg          = from%is_msg
    end subroutine copy_event_t

    !***************************************************************************
    !--- Print data from an event_t data type
    !***************************************************************************
    subroutine print_event_t(event)
      type(event_t) :: event
      write(6,'(2a)')  '       get_model = ',trim(event%get_model)
      write(6,'(2a)')  '       put_model = ',trim(event%put_model)
      write(6,'(2a)')  '         get_var = ',trim(event%get_var)
      write(6,'(2a)')  '         get_var = ',trim(event%put_var)
      write(6,'(2a)')  ' regrid_method   = ',trim(event%regrid_method)
      write(6,'(2a)')  ' regrid_norm_opt = ',trim(event%regrid_norm_opt)
      write(6,'(2a)')  '          action = ',trim(event%action)
      write(6,'(a,i8)')'        freq = ',event%freq
      write(6,'(a,i8)')'         tag = ',event%tag
      write(6,'(a,i8)')'       ncopy = ',event%ncopy
      write(6,'(a,i8)')'  get_master = ',event%get_master
      write(6,'(a,i8)')'  put_master = ',event%put_master
      write(6,'(a,i8)')'   regrid_id = ',event%regrid_id
      write(6,'(a,i8)')' alarm_index = ',event%alarm_index
      write(6,'(a,l1)')'      is_msg = ',event%is_msg
      call flush(6)
    end subroutine print_event_t

    !***************************************************************************
    !--- Copy data between 2 cpl_vinfo_t data types
    !***************************************************************************
    subroutine copy_cpl_vinfo_t(from, to)
      type(cpl_vinfo_t) :: from, to
      to%name     = from%name
      to%grid     = from%grid
      to%tag      = from%tag
      to%size     = from%size
      to%nlon     = from%nlon
      to%nlat     = from%nlat
      to%olap     = from%olap
      to%ncopy    = from%ncopy
      to%freq     = from%freq
      to%ccc_name = from%ccc_name
      to%regrid%method   = from%regrid%method
      to%regrid%norm_opt = from%regrid%norm_opt
      to%regrid%masked   = from%regrid%masked
    end subroutine copy_cpl_vinfo_t

    !***************************************************************************
    !--- Print data from a cpl_vinfo_t data type
    !***************************************************************************
    subroutine print_cpl_vinfo_t(vinfo)
      type(cpl_vinfo_t) :: vinfo
      write(6,'(2a)')  '    name = ',trim(vinfo%name)
      write(6,'(2a)')  '    grid = ',trim(vinfo%grid)
      write(6,'(a,i8)')'     tag = ',vinfo%tag
      write(6,'(a,i8)')'    size = ',vinfo%size
      write(6,'(a,i8)')'    nlon = ',vinfo%nlon
      write(6,'(a,i8)')'    nlat = ',vinfo%nlat
      write(6,'(a,i8)')'    olap = ',vinfo%olap
      write(6,'(a,i8)')'   ncopy = ',vinfo%ncopy
      write(6,'(a,i8)')'    freq = ',vinfo%freq
      write(6,'(2a)')  'ccc_name = ',trim(vinfo%ccc_name)
      write(6,'(2a)')  'regrid%   method = ',trim(vinfo%regrid%method)
      write(6,'(2a)')  'regrid% norm_opt = ',trim(vinfo%regrid%norm_opt)
      write(6,'(a,l1)')'regrid%   masked = ',vinfo%regrid%masked
      call flush(6)
    end subroutine print_cpl_vinfo_t

    !***************************************************************************
    !--- Define the list of all fields that may be coupled
    !***************************************************************************
    subroutine define_cpl_var_list()
      !--- cpl_var_list is an array containing info about all known variables.
      !--- NOTE: Every "name" (ie the first entry of each cpl_vinfo_t type assignment)
      !---       as well as every "tag" found in cpl_var_list MUST be unique.

      !--- Temporary space for ATM variables
      type(cpl_vinfo_t), pointer :: atm_var_list(:)
      integer :: n_atm_vars

      !--- Temporary space for CCCma's old NCOM based ocean variables
      type(cpl_vinfo_t), pointer :: oocn_var_list(:)
      integer :: n_oocn_vars

      !--- Temporary space for NEMO variables
      !--- Names associated with these variables are OASIS names
      type(cpl_vinfo_t), pointer :: nemo_var_list(:)
      integer :: n_nemo_vars

      !--- Temporary space for ICE variables
      type(cpl_vinfo_t), pointer :: ice_var_list(:)
      integer :: n_ice_vars

      !--- Temporary space for MISC variables
      type(cpl_vinfo_t), pointer :: misc_var_list(:)
      integer :: n_misc_vars

      !--- Temporary space for MSG variables
      type(cpl_vinfo_t), pointer :: msg_var_list(:)
      integer :: n_msg_vars

      integer :: n_list, idx, n_cpl_var_list
      integer :: verbose=0

      if ( verbose > 1 ) then
        write(6,*)"define_cpl_var_list: IN"
        call flush(6)
      endif

      n_atm_vars = 54
      nullify(atm_var_list)
      allocate( atm_var_list(n_atm_vars) )
      atm_var_list = &
      (/ cpl_vinfo_t('DATA_atm',   'atm', 100, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "DATA", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('GC_atm',     'atm', 102, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "  GC", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('GT_atm',     'atm', 103, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "  GT", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SIC_atm',    'atm', 104, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, " SIC", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OUFS_atm',   'atm', 105, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "OUFS", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OVFS_atm',   'atm', 106, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "OVFS", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SNO_atm',    'atm', 108, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, " SNO", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('UOC_atm',    'atm', 110, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, " UOC", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('VOC_atm',    'atm', 111, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, " VOC", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('UICE_atm',   'atm', 112, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "UICE", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('VICE_atm',   'atm', 113, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "VICE", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('GTX_atm',    'atm', 114, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, " GTX", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SLTX_atm',   'atm', 115, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "SLTX", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('HADP_atm',   'atm', 116, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "HADP", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('MADP_atm',   'atm', 117, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "MADP", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('YGWO_atm',   'atm', 118, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "YGWO", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('YGWS_atm',   'atm', 119, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "YGWS", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('YSWO_atm',   'atm', 120, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "YSWO", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('YSWS_atm',   'atm', 121, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "YSWS", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('YSWI_atm',   'atm', 122, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "YSWI", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('YDEP_atm',   'atm', 123, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "YDEP", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('ROFO_atm',   'atm', 124, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "ROFO", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SICN_atm',   'atm', 125, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "SICN", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OBEI_atm',   'atm', 126, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "OBEI", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OFSG_atm',   'atm', 127, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "OFSG", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('PMSL_atm',   'atm', 128, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "PMSL", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SWMX_atm',   'atm', 129, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "SWMX", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('XSRF_atm',   'atm', 130, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "XSRF", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('XSFX_atm',   'atm', 131, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "XSFX", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('NAOD_atm',   'atm', 159, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "XXXX", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('BEGO_atm',   'atm', 162, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "BEGO", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('RAIN_atm',   'atm', 163, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "RAIN", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SNOW_atm',   'atm', 164, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "SNOW", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('HFLI_atm',   'atm', 165, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "HFLI", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('BWGO_atm',   'atm', 166, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "BWGO", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('HSEA_atm',   'atm', 167, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "HSEA", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('BEGI_atm',   'atm', 168, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "BEGI", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SLIM_atm',   'atm', 169, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "SLIM", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SLIMlw_atm',   'atm', 6101, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "SLIM", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SLIMsh_atm',   'atm', 6102, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "SLIM", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SLIMlh_atm',   'atm', 6103, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "SLIM", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('GTICE_atm',  'atm', 6104, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, " GTI", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('MASK_atm',   'atm', 170, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "MASK", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('FLND_atm',   'atm', 171, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "FLND", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('UFSO_atm',   'atm', 172, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "UFSO", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('VFSO_atm',   'atm', 173, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "VFSO", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('UFSI_atm',   'atm', 174, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "UFSI", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('VFSI_atm',   'atm', 175, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "VFSI", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('FSGO_atm',   'atm', 176, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "FSGO", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('FSGI_atm',   'atm', 177, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "FSGI", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('CO2_atm',    'atm', 178, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "XSRF", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('CO2flx_atm', 'atm', 179, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "XSFX", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('RIVO_atm',   'atm', 180, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "RIVO", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('BWGI_atm',   'atm', 181, nwds_a, nlon_a, nlat_a, olap_a, 1, -1, "BWGI", &
                     regrid_info_t('bilinear', 'fracarea', .false.) )  /)

      if ( verbose > 0 ) then
        write(6,*)"define_cpl_var_list: size(atm_var_list)=",size(atm_var_list)
        call flush(6)
      endif

      n_oocn_vars = 20
      nullify(oocn_var_list)
      allocate( oocn_var_list(n_oocn_vars) )
      oocn_var_list = &
      (/ cpl_vinfo_t('DATA_ocn',   'ocn', 200, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "DATA", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('GTX_ocn',    'ocn', 201, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, " GTX", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OUFS_ocn',   'ocn', 202, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "OUFS", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OVFS_ocn',   'ocn', 203, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "OVFS", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OBEA_ocn',   'ocn', 204, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "OBEA", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OBWA_ocn',   'ocn', 205, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "OBWA", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('GT_ocn',     'ocn', 206, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "  GT", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SLT_ocn',    'ocn', 207, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, " SLT", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('UOC_ocn',    'ocn', 210, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, " UOC", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('VOC_ocn',    'ocn', 211, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, " VOC", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('CNV_ocn',    'ocn', 212, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, " CNV", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('PSIB_ocn',   'ocn', 213, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "PSIB", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('UISO_ocn',   'ocn', 214, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "UISO", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('VISO_ocn',   'ocn', 215, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "VISO", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OBEI_ocn',   'ocn', 216, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "OBEI", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('ODPX_ocn',   'ocn', 217, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "ODPX", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('ODPY_ocn',   'ocn', 218, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "ODPY", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('TR01_ocn',   'ocn', 219, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "TR01", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SST2_ocn',   'ocn', 220, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "SST2", &
                     regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('HBLM_ocn',   'ocn', 221, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "HBLM", &
                     regrid_info_t('bilinear', 'fracarea', .false.) )  /)

      if ( verbose > 0 ) then
        write(6,*)"define_cpl_var_list: size(oocn_var_list)=",size(oocn_var_list)
        call flush(6)
      endif

      n_ice_vars = 1
      nullify(ice_var_list)
      allocate( ice_var_list(n_ice_vars) )
      ice_var_list = &
      (/ cpl_vinfo_t('DATA_ice',   'ice', 300, nwds_i, nlon_i, nlat_i, olap_i, 1, -1, "DATA", &
                     regrid_info_t('bilinear', 'fracarea', .false.) )  /)

      if ( verbose > 0 ) then
        write(6,*)"define_cpl_var_list: size(ice_var_list)=",size(ice_var_list)
        call flush(6)
      endif

      n_nemo_vars = 53
      nullify(nemo_var_list)
      allocate( nemo_var_list(n_nemo_vars) )
      nemo_var_list = &
      (/ cpl_vinfo_t('DATA_nemo', 'ocn', 400, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OIceFrc',   'ocn', 401, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_SSTSST',  'ocn', 402, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_TepIce',  'ocn', 403, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_TepMix',  'ocn', 404, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_AlbIce',  'ocn', 405, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_AlbMix',  'ocn', 406, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OIceTck',   'ocn', 407, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OSnwTck',   'ocn', 408, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_IVelx1',  'ocn', 409, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_IVely1',  'ocn', 410, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_IVelz1',  'ocn', 411, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_IVelx1',  'ocn', 412, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_IVely1',  'ocn', 413, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_IVelz1',  'ocn', 414, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_CO2FLX',  'ocn', 415, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_OTaux1',  'ocn', 416, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_OTauy1',  'ocn', 417, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_OTauz1',  'ocn', 418, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_OTaux2',  'ocn', 419, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_OTauy2',  'ocn', 420, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_OTauz2',  'ocn', 421, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_ITaux1',  'ocn', 422, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_ITauy1',  'ocn', 423, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_ITauz1',  'ocn', 424, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_ITaux2',  'ocn', 425, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_ITauy2',  'ocn', 426, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_ITauz2',  'ocn', 427, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_QsrOce',  'ocn', 428, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_QsrIce',  'ocn', 429, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_QsrMix',  'ocn', 430, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_QnsOce',  'ocn', 431, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_QnsIce',  'ocn', 432, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_QnsMix',  'ocn', 433, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OTotRain',  'ocn', 434, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OTotSnow',  'ocn', 435, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OTotEvap',  'ocn', 436, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OIceEvap',  'ocn', 437, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OSubMPre',  'ocn', 438, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OISubMSn',  'ocn', 439, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OOEvaMPr',  'ocn', 440, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_Wind10',  'ocn', 441, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_dQnsdT',  'ocn', 442, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_dQlwdT',  'ocn',6104, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_dQshdT',  'ocn',6104, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_dQlhdT',  'ocn',6106, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_Runoff',  'ocn', 443, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OCalving',  'ocn', 444, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_TauMod',  'ocn', 445, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_AtmCO2',  'ocn', 446, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OTopMlt',   'ocn', 447, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('OBotMlt',   'ocn', 448, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) ), &
         cpl_vinfo_t('O_MSLP',   'ocn', 449, nwds_o, nlon_o, nlat_o, olap_o, 1, -1, "NEMO", &
                     regrid_info_t('conservative', 'fracarea', .false.) )  /)

      if ( verbose > 0 ) then
        write(6,*)"define_cpl_var_list: size(nemo_var_list)=",size(nemo_var_list)
        call flush(6)
      endif

      n_msg_vars = 5
      nullify(msg_var_list)
      allocate( msg_var_list(n_msg_vars) )
      msg_var_list = &
      (/ cpl_vinfo_t('MSG',      'msg', 900, nwds_d, nwds_d, 1, 0, 1, -1, " MSG", default_regrid ), &
         cpl_vinfo_t('MSG_ocn',  'msg', 901, nwds_d, nwds_d, 1, 0, 1, -1, "MSGO", default_regrid ), &
         cpl_vinfo_t('MSG_atm',  'msg', 902, nwds_d, nwds_d, 1, 0, 1, -1, "MSGA", default_regrid ), &
         cpl_vinfo_t('MSG_ice',  'msg', 903, nwds_d, nwds_d, 1, 0, 1, -1, "MSGI", default_regrid ), &
         cpl_vinfo_t('MSG_nemo', 'msg', 904, nwds_d, nwds_d, 1, 0, 1, -1, "MSGN", default_regrid )  /)

      if ( verbose > 0 ) then
        write(6,*)"define_cpl_var_list: size(msg_var_list)=",size(msg_var_list)
        call flush(6)
      endif

      n_misc_vars = 34
      nullify(misc_var_list)
      allocate( misc_var_list(n_misc_vars) )
      misc_var_list = &
      (/ cpl_vinfo_t('Sp_Op1',    'msg', 801, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('Sp_Op2',    'msg', 802, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('Sp_Op3',    'msg', 803, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('Sp_Op4',    'msg', 804, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('Sp_Op5',    'msg', 805, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('Sp_Op6',    'msg', 806, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('Sp_Op7',    'msg', 807, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('Sp_Op8',    'msg', 808, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('Sp_Op9',    'msg', 809, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('Sp_Op10',   'msg', 810, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('SST_canom',  'canom', 811, nwds_canom, nlon_canom, nlat_canom, olap_canom, &
                     1, -1, " SST", regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SLT_canom',  'canom', 812, nwds_canom, nlon_canom, nlat_canom, olap_canom, &
                     1, -1, " SLT", regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OUFS_canom', 'canom', 813, nwds_canom, nlon_canom, nlat_canom, olap_canom, &
                     1, -1, "OUFS", regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OVFS_canom', 'canom', 814, nwds_canom, nlon_canom, nlat_canom, olap_canom, &
                     1, -1, "OVFS", regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OBEA_canom', 'canom', 815, nwds_canom, nlon_canom, nlat_canom, olap_canom, &
                     1, -1, "OBEA", regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('OBWA_canom', 'canom', 816, nwds_canom, nlon_canom, nlat_canom, olap_canom, &
                     1, -1, "OBWA", regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('QSOL_canom', 'canom', 817, nwds_canom, nlon_canom, nlat_canom, olap_canom, &
                     1, -1, "QSOL", regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('SICN_canom', 'canom', 818, nwds_canom, nlon_canom, nlat_canom, olap_canom, &
                     1, -1, "SICN", regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('FCOK_canom', 'canom', 819, nwds_canom, nlon_canom, nlat_canom, olap_canom, &
                     1, -1, "FCOK", regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('XCO2_canom', 'canom', 820, nwds_canom, nlon_canom, nlat_canom, olap_canom, &
                     1, -1, "XCO2", regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('P_canom',    'canom', 821, nwds_canom, nlon_canom, nlat_canom, olap_canom, &
                     1, -1, "   P", regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('XKW_canom',  'canom', 822, nwds_canom, nlon_canom, nlat_canom, olap_canom, &
                     1, -1, " XKW", regrid_info_t('bilinear', 'fracarea', .false.) ), &
         cpl_vinfo_t('start_cpl2atm', 'msg', 823, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('stop_cpl2atm',  'msg', 824, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('start_atm2cpl', 'msg', 825, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('stop_atm2cpl',  'msg', 826, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('start_cpl2ocn', 'msg', 827, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('stop_cpl2ocn',  'msg', 828, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('start_ocn2cpl', 'msg', 829, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('stop_ocn2cpl',  'msg', 830, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('start_cpl2ice', 'msg', 831, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('stop_cpl2ice',  'msg', 832, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('start_ice2cpl', 'msg', 833, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ), &
         cpl_vinfo_t('stop_ice2cpl',  'msg', 834, 1, 1, 1, 0, 1, -1, "SpOp", default_regrid ) /)

      if ( verbose > 0 ) then
        write(6,*)"define_cpl_var_list: size(misc_var_list)=",size(misc_var_list)
        call flush(6)
      endif

      n_cpl_var_list = n_atm_vars + n_oocn_vars + n_ice_vars + n_nemo_vars &
                     + n_msg_vars + n_misc_vars
      nullify(cpl_var_list)
      allocate( cpl_var_list(n_cpl_var_list) )

      n_list = 0
      do idx = 1,n_atm_vars
        n_list = n_list + 1
        call copy_cpl_vinfo_t( atm_var_list(idx), cpl_var_list(n_list) )
      enddo
      do idx = 1,n_oocn_vars
        n_list = n_list + 1
        call copy_cpl_vinfo_t( oocn_var_list(idx), cpl_var_list(n_list) )
      enddo
      do idx = 1,n_ice_vars
        n_list = n_list + 1
        call copy_cpl_vinfo_t( ice_var_list(idx), cpl_var_list(n_list) )
      enddo
      do idx = 1,n_nemo_vars
        n_list = n_list + 1
        call copy_cpl_vinfo_t( nemo_var_list(idx), cpl_var_list(n_list) )
      enddo
      do idx = 1,n_misc_vars
        n_list = n_list + 1
        call copy_cpl_vinfo_t( misc_var_list(idx), cpl_var_list(n_list) )
      enddo
      do idx = 1,n_msg_vars
        n_list = n_list + 1
        call copy_cpl_vinfo_t( msg_var_list(idx), cpl_var_list(n_list) )
      enddo

      if ( verbose > 0 ) then
        write(6,*)"define_cpl_var_list: size(cpl_var_list)=",size(cpl_var_list)
        call flush(6)
      endif

      if ( n_list /= n_cpl_var_list ) then
        write(6,*)'define_cpl_var_list: Size mismatch.'
        call err_exit("define_cpl_var_list",-1)
      endif

      deallocate( atm_var_list, oocn_var_list, nemo_var_list, ice_var_list, &
                  msg_var_list, misc_var_list )

    end subroutine define_cpl_var_list

    !***************************************************************************
    !--- Add an new entry to cpl_var_list
    !***************************************************************************
    subroutine add_cpl_var(new_var, exists)

      !--- new_var must be completely defined by the user
      type(cpl_vinfo_t), pointer, intent(in) :: new_var

      !--- If the input variable name is the same as an existing variable then
      !--- exists provides the user control over whether the existing
      !--- variable info is used (exists=0), overwritten with info supplied
      !--- by the user (exists=1) or the program aborts when there
      !--- is a name conflict (exists=2)
      integer, intent(in), optional :: exists

      !--- Local
      type(cpl_vinfo_t), pointer :: tmp_list(:)
      integer :: idx, n_curr
      logical :: abort_when_exists, overwrite_when_exists, found
      integer :: verbose=1

      if ( .not. associated(new_var) ) then
        write(6,*)'add_cpl_var: Disassociated new_var pointer in parameter list.'
        call err_exit("add_cpl_var",-1)
      endif

      if ( present(exists) ) then
        if ( exists == 0 ) then
          !--- Use the existing variable info when a name conflict occurs
          abort_when_exists = .false.
          overwrite_when_exists = .false.
        else if ( exists == 1 ) then
          !--- Overwrite the existing variable info when a name conflict occurs
          abort_when_exists = .false.
          overwrite_when_exists = .true.
        else if ( exists == 2 ) then
          !--- Abort when a name conflict occurs
          abort_when_exists = .true.
          overwrite_when_exists = .false.
        else
          !--- Abort when exists is out of range
          write(6,*)'add_cpl_var: exists=',exists,' is out of range.'
          call err_exit("add_cpl_var",-2)
        endif
      else
        !--- The default is to not abort when a name conflict occurs
        !--- and the existing variable info is used in this case
        abort_when_exists = .false.
        overwrite_when_exists = .false.
      endif

      if ( .not. associated(cpl_var_list) ) then
        !--- Allocate cpl_var_list with size 1 and copy the new variable info
        allocate( cpl_var_list(1) )
        call copy_cpl_vinfo_t( new_var, cpl_var_list(1) )
        if ( verbose > 0 ) then
          write(6,*)"add_cpl_var: New variable definition."
          call print_cpl_vinfo_t(new_var)
          call flush(6)
        endif
        return
      endif

      !--- Get the length of the existing variable list
      n_curr = size(cpl_var_list)

      !--- Find out if this variable name or tag already exists in the list
      found = .false.
      do idx = 1,n_curr
        !--- Ensure that the new name or tag does not conflict with an existing definition
        if ( trim(adjustl(new_var%name)) .eq. &
             trim(adjustl(cpl_var_list(idx)%name)) ) then
          !--- This name exists in cpl_var_list but all names must be unique
          write(6,*)'add_cpl_var: The name --> ',trim(new_var%name),' <--- already exists.'
          if ( abort_when_exists ) then
            !--- Abort when an existing name is found
            call flush(6)
            call err_exit("add_cpl_var",-2)
          else
            !--- Issue a warning but continue
            if ( overwrite_when_exists ) then
              write(6,*)'add_cpl_var: Overwritting variable definition for ',trim(new_var%name)
              call copy_cpl_vinfo_t( new_var, cpl_var_list(idx) )
              if ( verbose > 0 ) then
                !xxx write(6,*)"add_cpl_var: New variable definition."
                call print_cpl_vinfo_t(new_var)
                call flush(6)
              endif
            else
              write(6,*)'add_cpl_var: Using existing variable definition for ',trim(new_var%name)
            endif
            call flush(6)
            found = .true.
            !--- Do not exit the loop here so that all tags may be checked
          endif

        else if ( new_var%tag == cpl_var_list(idx)%tag ) then
          !--- This tag exists in cpl_var_list but all tags must be unique
          write(6,*)'add_cpl_var: The tag --> ',new_var%tag,' <--- already exists.'
          !--- Always abort when an existing tag is found
          call flush(6)
          call err_exit("add_cpl_var",-3)
        endif
      enddo

      !--- If the variable already exists in the list then do nothing more
      if ( found ) return

      allocate( tmp_list(n_curr+1) )

      !--- Populate tmp_list with the current contents of cpl_var_list
      do idx = 1,n_curr
        call copy_cpl_vinfo_t( cpl_var_list(idx), tmp_list(idx) )
      enddo

      !--- Add the new variable info as the last entry
      call copy_cpl_vinfo_t( new_var, tmp_list(n_curr+1) )

      !--- Destroy then reassign cpl_var_list with the additional info
      nullify( cpl_var_list )
      allocate( cpl_var_list(n_curr+1) )
      do idx = 1,n_curr+1
        call copy_cpl_vinfo_t( tmp_list(idx), cpl_var_list(idx) )
      enddo

      !--- Clean up
      nullify( tmp_list )

      if ( verbose > 0 ) then
        write(6,*)"add_cpl_var: New variable definition added to cpl_var_list."
        call print_cpl_vinfo_t(new_var)
        call flush(6)
      endif

    end subroutine add_cpl_var

    !***************************************************************************
    !--- Determine if a given variable name exists in cpl_var_list
    !***************************************************************************
    function cpl_var_name_exists(name) result(exists)
      !--- The name to find
      character(*), intent(in) :: name

      !--- The boolean return value
      logical :: exists

      !--- Local
      integer :: idx

      !--- Find out if this name exists in cpl_var_list
      exists = .false.
      if ( size(cpl_var_list) > 0 ) then
        do idx = 1,size(cpl_var_list)
          if ( trim(adjustl(name)) .eq. &
               trim(adjustl(cpl_var_list(idx)%name)) ) then
            !--- This name exists in cpl_var_list
            exists = .true.
            exit
          endif
        enddo
      endif

    end function cpl_var_name_exists

    !***************************************************************************
    !--- Find a tag value that is not currently used in cpl_var_list
    !***************************************************************************
    function cpl_var_new_tag(start_tag) result(new_tag)
      !--- start_tag is the minimum tag value that will be returned
      integer, intent(in) :: start_tag

      !--- The return value is a tag that does not yet exist in cpl_var_list
      integer :: new_tag

      !--- Local
      integer :: indx
      logical :: found

      new_tag = start_tag

      if ( associated(cpl_var_list) ) then
        !--- When cpl_var_list is not defined simply return the value start_tag
        found = .false.
        do while (.not. found)
          found = .true.
          do indx=1,size(cpl_var_list)
            if ( cpl_var_list(indx)%tag == new_tag ) then
              !--- This tag exists in cpl_var_list
              found = .false.
              exit
            endif
          enddo
          if (.not. found) then
            new_tag = new_tag + 1
            if ( new_tag > 10000 ) then
              write(6,*)"cpl_var_new_tag: tag=",new_tag," is out of range."
              call flush(6)
              call err_exit("cpl_var_new_tag",-1)
            endif
          endif
        enddo
      endif

    end function cpl_var_new_tag

    !***************************************************************************
    !--- Define a cpl_vinfo_t variable on a particular grid given the name of
    !--- the field and the name of the grid
    !***************************************************************************
    subroutine define_cpl_var_by_name(name, grid)
      character(*), intent(in) :: name, grid

      !--- A temporary variable used to create new entries in cpl_var_list
      type(cpl_vinfo_t), pointer :: tmp_vinfo
      integer :: n

      allocate( tmp_vinfo )
      n = min(len(name),4)
      select case ( trim(adjustl(grid)) )
        case ("atm")
          tmp_vinfo%name     = trim(adjustl(name))
          tmp_vinfo%grid     = trim(adjustl(grid))
          tmp_vinfo%tag      = cpl_var_new_tag(5000)
          tmp_vinfo%size     = (nlon_a+olap_a)*nlat_a
          tmp_vinfo%nlon     = nlon_a
          tmp_vinfo%nlat     = nlat_a
          tmp_vinfo%olap     = olap_a
          tmp_vinfo%ncopy    = 1
          tmp_vinfo%freq     = 0
          tmp_vinfo%ccc_name = name(1:n)
          tmp_vinfo%regrid%method   = 'conservative'
          tmp_vinfo%regrid%norm_opt = 'fracarea'
          tmp_vinfo%regrid%masked   = .false.
          !--- Add this variable defintion
          call add_cpl_var(tmp_vinfo, exists=0)
        case ("ocn")
          tmp_vinfo%name     = trim(adjustl(name))
          tmp_vinfo%grid     = trim(adjustl(grid))
          tmp_vinfo%tag      = cpl_var_new_tag(5000)
          tmp_vinfo%size     = (nlon_o+olap_o)*nlat_o
          tmp_vinfo%nlon     = nlon_o
          tmp_vinfo%nlat     = nlat_o
          tmp_vinfo%olap     = olap_o
          tmp_vinfo%ncopy    = 1
          tmp_vinfo%freq     = 0
          tmp_vinfo%ccc_name = name(1:n)
          tmp_vinfo%regrid%method   = 'conservative'
          tmp_vinfo%regrid%norm_opt = 'fracarea'
          tmp_vinfo%regrid%masked   = .false.
          !--- Add this variable defintion
          call add_cpl_var(tmp_vinfo, exists=0)
        case default
          write(6,*)"define_cpl_var_by_name: Unknown grid ",trim(grid)
          call err_exit("define_cpl_var_by_name",-1)
      end select
      deallocate( tmp_vinfo )
    end subroutine define_cpl_var_by_name

    !***************************************************************************
    !--- Add fields to the event list that are sent coupler to NEMO
    !***************************************************************************
    subroutine add_events_cpl_to_nemo(nrec, tmp_event)
      !--- nrec is the current number of variables defined in tmp_event
      integer :: nrec

      !--- The event list to be appended with cpl_to_nemo variables
      type(event_t), pointer :: tmp_event(:)

      !--- Local
      type(event_t) :: this_event
      integer :: indx

      if ( nemo_n_recv_var > 0 ) then
        !--- nemo_recv_var will contain the list of variables transfered cpl to NEMO

        !--- Set defaults for certain elements (some may be reset below)
        this_event%get_model = "cpl"
        this_event%put_model = "ocn"
        this_event%regrid_method = "conservative"
        this_event%regrid_norm_opt = "fracarea"
        this_event%is_msg = .false.
        this_event%action = " "
        this_event%freq = cpl_ocn_freq  !--- cpl_ocn_freq is defined above

        !--- Insert an event to indicate the start of transfers from the coupler to NEMO
        nrec = nrec + 1
        this_event%get_var = "start_cpl2ocn"
        this_event%put_var = "start_cpl2ocn"
        call copy_event_t( this_event, tmp_event(nrec) )
        tmp_event(nrec)%is_msg = .true.
        tmp_event(nrec)%action = "Begin Transfer cpl to NEMO"

        !--- Insert an event for each variable that is to be transfered cpl to NEMO
        do indx=1,nemo_n_recv_var
          nrec = nrec + 1
          this_event%get_var = trim(nemo_recv_var(indx))
          this_event%put_var = trim(nemo_recv_var(indx))
          call copy_event_t( this_event, tmp_event(nrec) )
        enddo

        !--- Insert an event to indicate the end of transfers from the coupler to NEMO
        nrec = nrec + 1
        this_event%get_var = "stop_cpl2ocn"
        this_event%put_var = "stop_cpl2ocn"
        call copy_event_t( this_event, tmp_event(nrec) )
        tmp_event(nrec)%is_msg = .true.
        tmp_event(nrec)%action = "End Transfer cpl to NEMO"

#ifdef use_NEMO
      else

        write(6,'(a)')'add_events_cpl_to_nemo: nemo_n_recv_var == 0'
        call flush(6)
        call err_exit("add_events_cpl_to_nemo",-1)
#endif

      endif

    end subroutine add_events_cpl_to_nemo

    !***************************************************************************
    !--- Add fields to the event list that are sent coupler to AGCM
    !***************************************************************************
    subroutine add_events_cpl_to_agcm(nrec, tmp_event)
      !--- nrec is the current number of variables defined in tmp_event
      integer :: nrec

      !--- The event list to be appended with cpl_to_agcm variables
      type(event_t), pointer :: tmp_event(:)

      !--- Local
      type(event_t) :: this_event
      integer :: indx

      if ( atm_n_recv_var > 0 ) then
        !--- atm_recv_var will contain the list of variables transfered cpl to agcm

        !--- Set defaults for certain elements (some may be reset below)
        this_event%get_model = "cpl"
        this_event%put_model = "atm"
        this_event%regrid_method   = "conservative"
        this_event%regrid_norm_opt = "fracarea"
        this_event%is_msg = .false.
        this_event%action = " "
        this_event%freq = cpl_atm_freq  !--- cpl_atm_freq is defined above

        !--- This message event must occur before any other cpl to AGCM event since
        !--- the AGCM receives the message at the top of the do count=1,9999999 meaning this
        !--- occurs before any other MPI transfer done in the AGCM at each coupling interval.
        !--- This includes the "start_cpl2atm" event which will broadcast coupler time to the
        !--- agcm and therefore block the MSG transfer if the time broadcast goes first.
        nrec = nrec + 1
        this_event%get_var = "MSG"
        this_event%put_var = "MSG"
        call copy_event_t( this_event, tmp_event(nrec) )
        tmp_event(nrec)%is_msg = .true.

        !--- Insert an event to indicate the start of transfers from the coupler to the agcm
        nrec = nrec + 1
        this_event%get_var = "start_cpl2atm"
        this_event%put_var = "start_cpl2atm"
        call copy_event_t( this_event, tmp_event(nrec) )
        tmp_event(nrec)%is_msg = .true.
        tmp_event(nrec)%action = "Begin Transfer cpl to atm"

        !--- Insert an event for each variable that is to be transfered cpl to agcm
        do indx=1,atm_n_recv_var
          nrec = nrec + 1
          this_event%get_var = trim(atm_recv_var(indx))
          this_event%put_var = trim(atm_recv_var(indx))
          call copy_event_t( this_event, tmp_event(nrec) )
        enddo

        !--- Insert an event to indicate the end of transfers from the coupler to the agcm
        nrec = nrec + 1
        this_event%get_var = "stop_cpl2atm"
        this_event%put_var = "stop_cpl2atm"
        call copy_event_t( this_event, tmp_event(nrec) )
        tmp_event(nrec)%is_msg = .true.
        tmp_event(nrec)%action = "End Transfer cpl to atm"

      else

        write(6,'(a)')'add_events_cpl_to_agcm: atm_n_recv_var == 0'
        call flush(6)
        call err_exit("add_events_cpl_to_agcm",-1)

      endif

    end subroutine add_events_cpl_to_agcm

    !***************************************************************************
    !--- Add fields to the event list that are sent NEMO to coupler
    !***************************************************************************
    subroutine add_events_nemo_to_cpl(nrec, tmp_event)
      !--- nrec is the current number of variables defined in tmp_event
      integer :: nrec

      !--- The event list to be appended with nemo to cpl variables
      type(event_t), pointer :: tmp_event(:)

      !--- Local
      type(event_t) :: this_event
      integer :: indx

      if ( nemo_n_send_var > 0 ) then
        !--- nemo_send_var will contain the list of variables transfered NEMO to cpl

        !--- Set defaults for certain elements (some may be reset below)
        this_event%get_model = "ocn"
        this_event%put_model = "cpl"
        this_event%regrid_method = "conservative"
        this_event%regrid_norm_opt = "fracarea"
        this_event%is_msg = .false.
        this_event%action = " "
        this_event%freq = cpl_ocn_freq  !--- cpl_ocn_freq is defined above

        !--- Insert an event to indicate the start of transfers from NEMO to coupler
        nrec = nrec + 1
        this_event%get_var = "start_ocn2cpl"
        this_event%put_var = "start_ocn2cpl"
        call copy_event_t( this_event, tmp_event(nrec) )
        tmp_event(nrec)%is_msg = .true.
        tmp_event(nrec)%action = "Begin Transfer NEMO to cpl"

        !--- Insert an event for each variable that is to be transfered NEMO to coupler
        do indx=1,nemo_n_send_var
          nrec = nrec + 1
          this_event%get_var = trim(nemo_send_var(indx))
          this_event%put_var = trim(nemo_send_var(indx))
          call copy_event_t( this_event, tmp_event(nrec) )
        enddo

        !--- Insert an event to indicate the end of transfers from NEMO to coupler
        nrec = nrec + 1
        this_event%get_var = "stop_ocn2cpl"
        this_event%put_var = "stop_ocn2cpl"
        call copy_event_t( this_event, tmp_event(nrec) )
        tmp_event(nrec)%is_msg = .true.
        tmp_event(nrec)%action = "End Transfer NEMO to cpl"

#ifdef use_NEMO
      else

        write(6,'(a)')'add_events_nemo_to_cpl: nemo_n_send_var == 0'
        call flush(6)
        call err_exit("add_events_nemo_to_cpl",-1)
#endif

      endif

    end subroutine add_events_nemo_to_cpl

    !***************************************************************************
    !--- Add fields to the event list that are sent AGCM to coupler
    !***************************************************************************
    subroutine add_events_agcm_to_cpl(nrec, tmp_event)
      !--- nrec is the current number of variables defined in tmp_event
      integer :: nrec

      !--- The event list to be appended with agcm to cpl variables
      type(event_t), pointer :: tmp_event(:)

      !--- Local
      type(event_t) :: this_event
      integer :: indx
      character(4) :: cccname
      logical :: create_missing_name

      !--- A temporary variable used to create new entries in cpl_var_list
      type(cpl_vinfo_t), pointer :: tmp_vinfo

      if ( atm_n_send_var > 0 ) then
        !--- atm_send_var will contain the list of variables transfered agcm to cpl

        !--- Set defaults for certain elements (some may be reset below)
        this_event%get_model = "atm"
        this_event%put_model = "cpl"
        this_event%regrid_method   = "conservative"
        this_event%regrid_norm_opt = "fracarea"
        this_event%is_msg = .false.
        this_event%action = " "
        this_event%freq = cpl_atm_freq  !--- cpl_atm_freq is defined above

        !--- Insert an event to indicate the start of transfers from the agcm to the coupler
        nrec = nrec + 1
        this_event%get_var = "start_atm2cpl"
        this_event%put_var = "start_atm2cpl"
        call copy_event_t( this_event, tmp_event(nrec) )
        tmp_event(nrec)%is_msg = .true.
        tmp_event(nrec)%action = "Begin Transfer atm to cpl"

        !--- Insert an event for each variable that is to be transfered agcm to cpl
        create_missing_name = .true.
        do indx=1,atm_n_send_var
          if ( create_missing_name ) then
            if ( .not. cpl_var_name_exists(atm_send_var(indx)) ) then
              allocate( tmp_vinfo )
              tmp_vinfo%name     = trim(atm_send_var(indx))
              tmp_vinfo%grid     = "atm"
              tmp_vinfo%tag      = cpl_var_new_tag(6000)
              tmp_vinfo%size     = nwds_a
              tmp_vinfo%nlon     = nlon_a
              tmp_vinfo%nlat     = nlat_a
              tmp_vinfo%olap     = olap_a
              tmp_vinfo%ncopy    = 1
              tmp_vinfo%freq     = cpl_atm_freq
              cccname = atm_send_var(indx)(1:4)
              tmp_vinfo%ccc_name = adjustr(cccname)
              tmp_vinfo%regrid%method   = this_event%regrid_method
              tmp_vinfo%regrid%norm_opt = this_event%regrid_norm_opt
              tmp_vinfo%regrid%masked   = .false.
              !--- Add this variable defintion
              call add_cpl_var(tmp_vinfo, exists=0)
              nullify( tmp_vinfo )
            endif
          endif
          nrec = nrec + 1
          this_event%get_var = trim(atm_send_var(indx))
          this_event%put_var = trim(atm_send_var(indx))
          call copy_event_t( this_event, tmp_event(nrec) )
        enddo

        !--- This message event must occur after all variables have been transfered
        nrec = nrec + 1
        this_event%get_var = "MSG"
        this_event%put_var = "MSG"
        call copy_event_t( this_event, tmp_event(nrec) )
        tmp_event(nrec)%is_msg = .true.

        !--- Insert an event to indicate the end of transfers from the agcm to the coupler
        nrec = nrec + 1
        this_event%get_var = "stop_atm2cpl"
        this_event%put_var = "stop_atm2cpl"
        call copy_event_t( this_event, tmp_event(nrec) )
        tmp_event(nrec)%is_msg = .true.
        tmp_event(nrec)%action = "End Transfer atm to cpl"

      else

        write(6,'(a)')'add_events_agcm_to_cpl: atm_n_send_var == 0'
        call flush(6)
        call err_exit("add_events_agcm_to_cpl",-1)

      endif

    end subroutine add_events_agcm_to_cpl

    !***************************************************************************
    !--- Initialize the default event list
    !***************************************************************************
    subroutine init_events_part1()
      integer :: indx, nrec, iev, iev_prev, tag
      integer(kind=impi) :: rank, ierr
      type(event_t), pointer :: tmp_event(:)
      type(cpl_vinfo_t) :: get_vinfo, put_vinfo
      integer :: verbose=0

      nrec = 0
      allocate( tmp_event(100) )

      write(6,*) "com_cpl.F90, init_events_part1(): couple_serial = ", couple_serial
      write(6,*) "com_cpl.F90, init_events_part1(): coupler_parallel_type = ", coupler_parallel_type


      if ( couple_serial ) then
        !***************************************************************************
        !--- The coupler communicates with component models serially
        !---   1) send data to agcm
        !---   2) receive data from agcm
        !---   3) send data to ocean
        !---   4) receive data from ocean
        !***************************************************************************

        !--- The coupler sends data to the AGCM
        !--- Append the event list with variables sent from the coupler to the AGCM
        call add_events_cpl_to_agcm(nrec, tmp_event)

        !--- The coupler receives data from the AGCM
        !--- Append the event list with variables sent from the AGCM to the coupler
        call add_events_agcm_to_cpl(nrec, tmp_event)

        !--- The coupler sends data to NEMO
        !--- Append the event list with variables sent from the coupler to NEMO
        call add_events_cpl_to_nemo(nrec, tmp_event)

        !--- The coupler receives data from NEMO
        !--- Append the event list with variables sent from NEMO to the coupler
        call add_events_nemo_to_cpl(nrec, tmp_event)

      else
        !***************************************************************************
        !--- The coupler communicates with component models in parallel
        !***************************************************************************
        write(6,*) "com_cpl.F90, init_events_part1(): selecting case for coupler_parallel_type = ", coupler_parallel_type
        select case (coupler_parallel_type)
          case (1)
            !***************************************************************************
            !--- The order of operation is as follows
            !---   1) send data to agcm
            !---   2) send data to ocean
            !---   3) receive data from agcm
            !---   4) receive data from ocean
            !***************************************************************************

            !--- The coupler sends data to the AGCM
            !--- Append the event list with variables sent from the coupler to the AGCM
            call add_events_cpl_to_agcm(nrec, tmp_event)

            !--- The coupler sends data to NEMO
            !--- Append the event list with variables sent from the coupler to NEMO
            call add_events_cpl_to_nemo(nrec, tmp_event)

            !--- The coupler receives data from the AGCM
            !--- Append the event list with variables sent from the AGCM to the coupler
            call add_events_agcm_to_cpl(nrec, tmp_event)

            !--- The coupler receives data from NEMO
            !--- Append the event list with variables sent from NEMO to the coupler
            call add_events_nemo_to_cpl(nrec, tmp_event)

          case (2)
            !***************************************************************************
            !--- The order of operation is as follows
            !---   1) send data to agcm
            !---   2) send data to ocean
            !---   3) receive data from ocean
            !---   4) receive data from agcm
            !***************************************************************************

            !--- The coupler sends data to the AGCM
            !--- Append the event list with variables sent from the coupler to the AGCM
            call add_events_cpl_to_agcm(nrec, tmp_event)

            !--- The coupler sends data to NEMO
            !--- Append the event list with variables sent from the coupler to NEMO
            call add_events_cpl_to_nemo(nrec, tmp_event)

            !--- The coupler receives data from NEMO
            !--- Append the event list with variables sent from NEMO to the coupler
            call add_events_nemo_to_cpl(nrec, tmp_event)

            !--- The coupler receives data from the AGCM
            !--- Append the event list with variables sent from the AGCM to the coupler
            call add_events_agcm_to_cpl(nrec, tmp_event)

          case default
            write(6,*)"init_events_part1: Invalid coupler_parallel_type = ", &
                      coupler_parallel_type
            call err_exit("init_events_part1",-1)

        end select
      endif

      !***************************************************************************
      !--- Assign other elements that are available here
      !--- Certain elements (e.g. get_master, put_master, alarm_index)
      !--- cannot be assigned until after commworlds, regrid operations and
      !--- alarms have been defined
      !***************************************************************************
      do iev=1,nrec
        !--- Define an MPI tag to be associated with this event
        !--- Some events will not use this tag but define one for all events so that
        !--- exceptions do not need to be determined here, mitigating error
        !--- None of these tags can be used elsewhere
        tag = 0
        if ( iev > 1 ) then
          !--- Determine if the current name has been used in a previous event
          !--- If so then set the tag to that of the previous event
          !--- This is required when there is more than one MSG event
          do iev_prev=1,iev-1
            !--- There is always a single tag associated with a given event and
            !--- this tag value will have been set on a previous loop iteration
            if ( trim(adjustl(tmp_event(iev     )%get_var)) .eq. &
                 trim(adjustl(tmp_event(iev_prev)%get_var)) ) then
              !--- Use the tag from the first of any duplicates
              tag = tmp_event(iev_prev)%tag
              exit
            endif
            if ( trim(adjustl(tmp_event(iev     )%put_var)) .eq. &
                 trim(adjustl(tmp_event(iev_prev)%put_var)) ) then
              !--- Use the tag from the first of any duplicates
              tag = tmp_event(iev_prev)%tag
              exit
            endif
          enddo
        endif
        if ( tag == 0 ) tag = 2300 + iev
        tmp_event(iev)%tag = tag

        !--- Ensure that both get_var and put_var appear in cpl_var_list and contain
        !--- the appropriate information for this run
        get_vinfo = find_cpl_vinfo(name=trim(adjustl(tmp_event(iev)%get_var)), list_index=indx)

        !--- Set the frequency for get_var in the master list (ie cpl_var_list)
        cpl_var_list(indx)%freq =  tmp_event(iev)%freq

        !--- Set the MPI tag used with get_var (the same tag must be used with put_var)
        cpl_var_list(indx)%tag =  tag

        put_vinfo = find_cpl_vinfo(name=trim(adjustl(tmp_event(iev)%put_var)), list_index=indx)

        !--- Set the frequency for put_var in the master list (ie cpl_var_list)
        cpl_var_list(indx)%freq =  tmp_event(iev)%freq

        !--- Set the MPI tag used with put_var (the same tag must be used with get_var)
        cpl_var_list(indx)%tag =  tag

        !--- Set the regrid method and normalization type for put_var in the master list
        cpl_var_list(indx)%regrid%method   = tmp_event(iev)%regrid_method
        cpl_var_list(indx)%regrid%norm_opt = tmp_event(iev)%regrid_norm_opt

        !--- Assign ncopy from get_var
        tmp_event(iev)%ncopy  = get_vinfo % ncopy

        if ( len_trim(tmp_event(iev)%action) <= 0 ) then
          !--- Assign an action here unless it was assigned above
          tmp_event(iev)%action = " "
          tmp_event(iev)%action = " Transfer "// &
              trim(tmp_event(iev)%get_model)//":"//trim(tmp_event(iev)%get_var)//" to "// &
              trim(tmp_event(iev)%put_model)//":"//trim(tmp_event(iev)%put_var)
        endif
      enddo

      !***************************************************************************
      !--- Set the total number of events
      !***************************************************************************
      nevents = nrec

      !***************************************************************************
      !--- Allocate and assign the event array
      !***************************************************************************
      if ( associated(event) ) nullify(event)
      allocate( event(nevents) )
      do iev=1,nevents
        call copy_event_t( tmp_event(iev), event(iev) )
        if ( verbose > 10 ) then
#ifdef use_mpi
          call mpi_comm_rank ( MPI_COMM_WORLD, rank, ierr )
          if ( rank == cpl_master .or. &
               rank == atm_master .or. &
               rank == ocn_master ) then
            !--- write(6,*)"init_events_part1: event number ",iev
            !--- call print_event_t( event(iev) )
            write(6,*)"init_events_part1: event=",iev,"  rank=",rank, &
                      "  ",trim(event(iev)%action), &
                      "  freq=",event(iev)%freq,"  tag=",event(iev)%tag
            call flush(6)
          endif
#else
          write(6,*)"init_events_part1: event number ",iev
          call print_event_t( event(iev) )
#endif
        endif
      enddo

      deallocate( tmp_event )

    end subroutine init_events_part1

    !***************************************************************************
    !--- Determine if a given name is to be coupled (ie is on the event stack)
    !***************************************************************************
    logical function is_coupled(name)
      !--- The name of the variable (this name must be known to the coupler)
      character(*), intent(in) :: name

      !--- Local
      integer :: idx1

      !--- Coupler events must be defined prior to invoking this routine
      if ( nevents <= 0 ) then
        write(6,*)"is_coupled: The events array is not defined."
        call err_exit("is_coupled",-1)
      endif

      is_coupled = .false.
      do idx1=1,nevents
        if ( trim(adjustl(event(idx1)%get_var)) .eq. trim(adjustl(name)) ) then
          !--- This variable is sent by at least one component model
          !--- If event(idx1)%get_master is defined (> -1) then we can determine if
          !--- the current task is sending it and set coupled accordingly
          !--- For the time being just indicate if the field is coupled or not
          is_coupled = .true.
          exit
        endif
        if ( trim(adjustl(event(idx1)%put_var)) .eq. trim(adjustl(name)) ) then
          !--- This variable is received by at least one component model
          !--- If event(idx1)%put_master is defined (> -1) then we can determine if
          !--- the current task is receiving it and set coupled accordingly
          !--- For the time being just indicate if the field is coupled or not
          is_coupled = .true.
          exit
        endif
      enddo

    end function is_coupled

    !***************************************************************************
    !--- Return a cpl_vinfo_t record from cpl_var_list that matches either
    !--- the "name" or the "tag" element
    !***************************************************************************
    function find_cpl_vinfo( name, tag, list_index, status ) result(vinfo)
      character(len=*), intent(in),  optional :: name
      integer,          intent(in),  optional :: tag
      !--- Optionally return the index in cpl_var_list of the matching record
      integer,          intent(out), optional :: list_index
      !--- Optionally return a status flag
      integer,          intent(out), optional :: status
      type(cpl_vinfo_t) :: vinfo

      !--- local
      integer :: idx, indx

      vinfo%name     = " "
      vinfo%grid     = " "
      vinfo%tag      = -1
      vinfo%size     = -1
      vinfo%nlon     = -1
      vinfo%nlat     = -1
      vinfo%ncopy    = -1
      vinfo%freq     = -1
      vinfo%ccc_name = " "
      vinfo%regrid%method   = " "
      vinfo%regrid%norm_opt = " "
      vinfo%regrid%masked   = .false.

      indx = -1
      if ( present(list_index) ) list_index=-1

      if ( present(status) ) status = 0

      if ( present(name) .and. present(tag) ) then
        write(6,'(3a,i6)')'find_cpl_vinfo: Both name and tag are present.  name=', &
            trim(name),'  tag=',tag
        call flush(6)
        call err_exit("FIND_CPL_VINFO",-1)
      endif

      if ( present(name) ) then
        do idx=1,size(cpl_var_list)
          if ( trim(adjustl(cpl_var_list(idx)%name)) .eq. trim(adjustl(name)) ) then
            indx = idx
            exit
          endif
        enddo

      else if ( present(tag) ) then
        do idx=1,size(cpl_var_list)
          if ( cpl_var_list(idx)%tag .eq. tag ) then
            indx = idx
            exit
          endif
        enddo

      else
        write(6,'(a)')'find_cpl_vinfo: Neither name nor tag are present.'
        call flush(6)
        call err_exit("FIND_CPL_VINFO",-2)
      endif

      if (indx == -1) then
        if ( present(name) ) then
          write(6,'(3a)')'find_cpl_vinfo: Unknown variable name --> ', &
                           trim(name),' <-- was requested.'
        else
          write(6,'(a,i6,a)')'find_cpl_vinfo: Unknown tag --> ', &
                           tag,' <-- was requested.'
        endif
        call flush(6)
        if ( present(status) ) then
          !--- Do not abort if the optional status flag is requested by the user
          status = 1
          write(6,'(3a)')'find_cpl_vinfo: Ignoring failed request.'
          call flush(6)
        else
          !--- Abort when the record cannot be found
          call err_exit("FIND_CPL_VINFO",-3)
        endif

      else
        call copy_cpl_vinfo_t(cpl_var_list(indx), vinfo)
        if ( present(list_index) ) list_index=indx

      endif

    end function find_cpl_vinfo

    !***************************************************************************
    !--- Return a character string containing a formatted version of ibuf(1:8)
    !***************************************************************************
    function format_ibuf(ibuf) result(fmt_ibuf)
      integer(kind=8), intent(in) :: ibuf(8)
      character(len=82) :: fmt_ibuf
      character(len=4) :: cckind, ccname
      integer :: idx

      cckind = " "
      if ( ibuf(1).gt.9999 ) then
        !--- Assume ibuf1 represents a 4 character string
        cckind = transfer(ibuf(1),"XXXX")
        !--- Replace any nulls found in cckind with blanks
        do idx=1,len(cckind)
          if ( cckind(idx:idx) .eq. achar(0) ) cckind(idx:idx) = " "
        enddo
      else if ( ibuf(1).ge.0 ) then
        !--- Write the integer value of ibuf1 (up to 4 digits) to cckind
        write(cckind,'(i4)')ibuf(1)
      else
        !--- This could be an error, but carry on anyway
        cckind = " "
      endif

      ccname = " "
      if ( ibuf(3).gt.9999 ) then
        !--- Assume ibuf3 represents a 4 character string
        ccname = transfer(ibuf(3),"XXXX")
        !--- Replace any nulls found in ccname with blanks
        do idx=1,len(ccname)
          if ( ccname(idx:idx) .eq. achar(0) ) ccname(idx:idx) = " "
        enddo
      else if ( ibuf(3).ge.0 ) then
        !--- Write the integer value of ibuf3 (up to 4 digits) to ccname
        write(ccname,'(i4)')ibuf(3)
      else
        !--- This could be an error, but carry on anyway
        ccname = " "
      endif

      fmt_ibuf = " "
      if ( len_trim(cckind) <= 0 .and. len_trim(ccname) <= 0 ) then
        write(fmt_ibuf,'(i10,1x,i10,1x,i10,5i10)')ibuf(1:8)
      else if ( len_trim(cckind) <= 0 ) then
        write(fmt_ibuf,'(i10,1x,i10,1x,6x,a4,5i10)')ibuf(1),ibuf(2),ccname,ibuf(4:8)
      else if ( len_trim(ccname) <= 0 ) then
        write(fmt_ibuf,'(6x,a4,1x,i10,1x,i10,5i10)')cckind,ibuf(2),ibuf(3),ibuf(4:8)
      else
        write(fmt_ibuf,'(6x,a4,1x,i10,1x,6x,a4,5i10)')cckind,ibuf(2),ccname,ibuf(4:8)
      endif
    end function format_ibuf

    !-----------------------------------------------------------------------------
    !--- String hash function
    !--- Developed by Daniel J. Bernstein
    !--- https://fortrandev.wordpress.com/2013/07/06/fortran-hashing-algorithm/
    !-----------------------------------------------------------------------------
    function djb_hash(str) result(hash)
      character(len=*), intent(in) :: str
      integer(kind=8) :: hash
      integer :: i = 0

      hash = 5381
      do i=1,len(str)
        hash = (ishft(hash,5) + hash) + ichar(str(i:i))
      end do

    end function djb_hash

    !-----------------------------------------------------------------------------
    !--- Another string hash function
    !--- Origin unknown
    !-----------------------------------------------------------------------------
    function modish_hash( key ) result(index)
      character(len=*), intent(in) :: key
      integer :: index

      !--- local
      integer :: idx
      integer, parameter :: hash_size  = 4993
      integer, parameter :: multiplier = 31

      index = 0

      do idx = 1,len(key)
        index = multiplier * index + ichar(key(idx:idx))
      enddo

      index = 1 + mod( index-1, hash_size )
    end function modish_hash

    !-----------------------------------------------------------------------------
    !--- Bob Jenkin's "One at a time" hash
    !--- https://en.wikipedia.org/wiki/Jenkins_hash_function#one-at-a-time
    !-----------------------------------------------------------------------------
!???? type elemental has problems on xlf
!????    elemental function bj_hash(str)
    function bj_hash(str)
      !! Bob Jenkins' one-at-a-time hash function
      !! Map a given character string onto an integer. For best performance,
      !! this should be:
      !!     1. fast
      !!     2. deterministic
      !!     3. uniformly distributed over the integers
      !!     4. chaotic
      !! See <http://burtleburtle.net/bob/hash/doobs.html>, and
      !! <http://www.stanford.edu/class/ee380/Abstracts/121017-slides.pdf>
      !! (slides from talk of google cityhash developers to Stanford)
      character(len=*), intent(in)   :: str
      !! string to hash
      integer(selected_int_kind(8))  :: bj_hash
      !! Numeric bj_hash
      character(len=len(str))        :: str2
      integer(selected_int_kind(8))  :: i, stlen
      str2 = adjustl(str)
      stlen = len(trim(str2))
      bj_hash  = 0
      do i = 1,stlen
         bj_hash = bj_hash + ichar(str2(i:i))
         bj_hash = bj_hash + ishft(bj_hash,10)
         bj_hash = ieor(bj_hash,ishft(bj_hash,-6))
      end do
      bj_hash = bj_hash + ishft(bj_hash,3)
      bj_hash = ieor(bj_hash,ishft(bj_hash,-11))
      bj_hash = abs(bj_hash + ishft(bj_hash,15))
    end function bj_hash

#ifdef use_mpi
    !--- Everthing below this line requires mpi

    !***************************************************************************
    !--- Initialize certain variables
    !--- This subroutine must be called prior to using any variables, procedures
    !--- etc from com_cpl
    !***************************************************************************
    subroutine com_mpi_init()
      integer(kind=impi) :: myrank, ierr
      integer :: verbose=0
      logical(kind=impi) :: mpi_init_was_called

      !--- Do nothing if already initialized
      if ( com_mpi_initialized ) return

      !--- Call MPI_init if MPI is not already initialized
      !--- The atm and the coupler will have already called mpi_init
      !--- but NEMO will not call mpi_init when it is coupled
      call mpi_initialized ( mpi_init_was_called, ierr )
      if ( .not. mpi_init_was_called ) then
        call mpi_init( ierr )
      else
        if ( verbose > 0 ) then
          call mpi_comm_rank ( MPI_COMM_WORLD, myrank, ierr )
          write(6,'(a,i4)')'com_mpi_init: ** II ** MPI is already initialized on task ',myrank
          call flush(6)
        endif
      endif

      call mpi_comm_rank ( MPI_COMM_WORLD, myrank, ierr )

      !--- Assign a global variable containing the current task number
      curr_task = myrank

      !--- Allocate and populate cpl_var_list
      call define_cpl_var_list()

      !--- Define the mpi real type according to word size on the execution machine
      if ( kind(0.0) == selected_real_kind(13) ) then
        !--- default real is MPI double precisison
        curr_real_type = MPI_DOUBLE_PRECISION
      else if ( kind(0.0) == selected_real_kind(6) ) then
        !--- default real is 32 bits
        curr_real_type = MPI_REAL
      else
        write(6,*)'com_mpi_init: Unable to determine real type to use with MPI.'
        call flush(6)
        call err_exit("com_mpi_init",-2)
      end if

      !--- Verify that the horizontal grid resolution parameters are reasonable
      select case (nlon_a)
        case (64, 96, 128, 256)
          !--- These are expected values
        case default
          write(6,'(a,i8,a,i4)') &
            'com_mpi_init: Invalid number of longitudes nlon_a = ',nlon_a, &
            '  on rank ',myrank
          call flush(6)
          call err_exit("com_mpi_init",-3)
      end select

      select case (nlat_a)
        case (32, 48, 64, 128)
          !--- These are expected values
        case default
          write(6,'(a,i8,a,i4)') &
            'com_mpi_init: Invalid number of latitudes nlat_a = ',nlat_a, &
            '  on rank ',myrank
          call flush(6)
          call err_exit("com_mpi_init",-4)
      end select

      select case (nlon_o)
        case (192, 256, 360)
          !--- These are expected values
        case default
          write(6,'(a,i8,a,i4)') &
            'com_mpi_init: Invalid number of longitudes nlon_o = ',nlon_o, &
            '  on rank ',myrank
          call flush(6)
          call err_exit("com_mpi_init",-5)
      end select

      select case (nlat_o)
        case (96, 192, 292)
          !--- These are expected values
        case default
          write(6,'(a,i8,a,i4)') &
            'com_mpi_init: Invalid number of latitudes nlat_o = ',nlat_o, &
            '  on rank ',myrank
          call flush(6)
          call err_exit("com_mpi_init",-6)
      end select

      !--- WARNING ABOUT THE XLF COMPILER AND STRING INITIALIZATION ---
      !--- It seems that the xlf compiler (and perhaps others) cannot initialize character
      !--- variables correctly when those variables are part of an array.
      !--- If the first element in the list is initialized with a string of a given length
      !--- then all other elements in the list will be truncated at that length.
      !--- e.g. if the first element is initialized as "GT" and the second element is
      !--- initialized as "SIC" then the second element will be truncated
      !--- to 2 characters "SI"
      !--- WARNING ABOUT THE XLF COMPILER AND STRING INITIALIZATION ---

      if ( verbose > 0 ) then
        write(6,'(a,i4)')"com_mpi_init: ** II ** Initialization complete on task ",myrank
        call flush(6)
      endif

      !--- Set a logical flag to indicate that com_cpl has been initialized
      com_mpi_initialized = .true.

    end subroutine com_mpi_init

    !***************************************************************************
    !--- Initialize the event array
    !--- A complete definition of each event will be done later
    !***************************************************************************
    subroutine cpl_initialize_events()
      integer(kind=impi) :: myrank, ierr
      integer :: idx1, verbose=1

      !--- Broadcast values of certain atm and ocn variables to all mpi tasks
      !--- cpl_define_global_parameters must be called after MPI has been initialized
      !--- and after atm_master, ocn_master and cpl_master have been assigned
      !--- (ie after the call to define_group)
      !--- nemo_nn_ice and nemo_ncat are two of these values and are required
      !--- to build the event array
      call cpl_define_global_parameters()

      !--- Define the nominal coupling frequency for all component models
      !--- atm_ksteps and atm_delt are defined in the agcm and broadcast
      !--- to all component models in cpl_define_global_parameters
      !--- Currently these 2 parameters define the coupling frequency for both
      !--- the atm and the ocean.
      !--- It is possible to set a different coupling frequency for the ocean but
      !--- as it stands now, NEMO will simply store the values after each transfer and
      !--- use the same data over again if the NEMO sbc time step (dependent on nn_fsbc)
      !--- is shorter than the coupling frequency
      cpl_atm_freq = atm_ksteps * atm_delt
      cpl_ocn_freq = cpl_atm_freq

      !--- Initialize the list of all transfer events
      !--- This partially defines the event list, the assignment is completed later
      call init_events_part1()

      if ( nevents <= 0 ) then
        write(6,*)"cpl_initialize_events: No events."
        call err_exit("cpl_initialize_events",-1)
      endif

      if (verbose > 0) then
        call mpi_comm_rank ( MPI_COMM_WORLD, myrank, ierr )
        do idx1=1,nevents
          write(6,'(a,i3,8a,a,i8,a,i4,a)')' com_mpi_init: event ',idx1, &
              '  transfer ',trim(event(idx1)%get_model),":",trim(event(idx1)%get_var), &
              ' to ',trim(event(idx1)%put_model),":",trim(event(idx1)%put_var), &
              " freq=",event(idx1)%freq," (on task ",myrank,")"
        enddo
        call flush(6)
      endif

    end subroutine cpl_initialize_events

    !***************************************************************************
    !--- Define groups for coupler, atm, ocn, ice
    !--- For each group assign/determine the master task and group communicator
    !--- This will define: cpl_master, atm_master, ocn_master, ice_master
    !***************************************************************************
    subroutine define_group(group_name, group_comm, group_member)

      !--- A name used to identify this group
      !--- Only the first 3 characters are used
      !--- Valid values are : cpl, atm, ocn, ice
      character(*), intent(in) :: group_name

      !--- The mpi communicator that is assigned to this group
      integer(kind=impi), intent(inout) :: group_comm

      character(len=8), intent(out), optional :: group_member(0:1000)

      !--- local
      integer :: cpl_rank_min, cpl_rank_max
      integer :: atm_rank_min, atm_rank_max
      integer :: ocn_rank_min, ocn_rank_max
      integer :: ice_rank_min, ice_rank_max
      integer :: ivar
      integer :: verbose=0
      character(len=3) :: name
      integer(kind=impi) :: rank, nproc, proc, color, key

      if ( .not. com_mpi_initialized ) call com_mpi_init()
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      if ( group_name(1:3).ne.'cpl' .and. &
           group_name(1:3).ne.'atm' .and. &
           group_name(1:3).ne.'ocn' .and. &
           group_name(1:3).ne.'ice' ) then
        write(6,'(2a)') 'define_group: Invalid group name ',group_name(1:3)
        write(6,'(a)')  'define_group: Valid names are: atm cpl ocn ice'
        call flush(6)
        call err_exit('define_group',-1)
      endif

      !---Determine the rank of the calling process in MPI_COMM_WORLD
      call mpi_comm_rank ( MPI_COMM_WORLD, rank, ierr )

      !---Determine the number of processes in MPI_COMM_WORLD
      call mpi_comm_size ( MPI_COMM_WORLD, nproc, ierr )

      !--- Allocate the task info array for the current process
      if ( associated(task_info) ) nullify(task_info)
      allocate( task_info(0:nproc-1) )

      if (verbose > 1) then
        write(6,'(4a,3i4)') &
          'define_group: group_name=',group_name(1:3), &
          '  MPI_COMM_WORLD,rank,size ',MPI_COMM_WORLD,rank,nproc
        call flush(6)
      endif

      atm_rank_min = nproc
      atm_rank_max = 0
      ocn_rank_min = nproc
      ocn_rank_max = 0
      cpl_rank_min = nproc
      cpl_rank_max = 0
      ice_rank_min = nproc
      ice_rank_max = 0
      key = 1

      !--- Broadcast the name assigned to each process to MPI_COMM_WORLD and
      !--- use the responses to determine process ranges for each sub group
      do proc=0,nproc-1
        !--- Default integer ivar is required for the min/max functions below
        ivar = proc

        if (proc.eq.rank) then
          !--- Broadcast the current group name (passed in by caller)
          name = group_name(1:3)
        else
          name = 'XXX'
        endif

        call mpi_bcast(name, 3_impi, MPI_CHARACTER, proc, MPI_COMM_WORLD, ierr)

        select case( trim(adjustl(name)) )
          case ('atm')
            atm_rank_min = min(atm_rank_min, ivar)
            atm_rank_max = max(atm_rank_max, ivar)
          case ('ocn')
            ocn_rank_min = min(ocn_rank_min, ivar)
            ocn_rank_max = max(ocn_rank_max, ivar)
          case ('cpl')
            cpl_rank_min = min(cpl_rank_min, ivar)
            cpl_rank_max = max(cpl_rank_max, ivar)
          case ('ice')
            ice_rank_min = min(ice_rank_min, ivar)
            ice_rank_max = max(ice_rank_max, ivar)
          case default
            !--- This should never happen, but ...
            write(6,'(4a)')'define_group: group name ',group_name(1:3),'  Invalid name ',name
            write(6,'(a)') 'define_group: Valid names are: atm cpl ocn ice'
            call err_exit('define_group',-2)
        end select

      enddo

      if (verbose > 1) then
        ! call mpi_barrier(MPI_COMM_WORLD, ierr)
        write(6,'(2a,4(a,i3,a,i3))') &
          'define_group: group_name=',group_name(1:3), &
          '  cpl_rank_min=',cpl_rank_min,' cpl_rank_max=',cpl_rank_max, &
          '  atm_rank_min=',atm_rank_min,' atm_rank_max=',atm_rank_max, &
          '  ocn_rank_min=',ocn_rank_min,' ocn_rank_max=',ocn_rank_max, &
          '  ice_rank_min=',ice_rank_min,' ice_rank_max=',ice_rank_max
        call flush(6)
      endif

      !--- Split MPI_COMM_WORLD into sub groups, one group for each of cpl, atm, ocn, ice
      color = MPI_UNDEFINED
      if ((rank.ge.cpl_rank_min).and.(rank.le.cpl_rank_max)) then
        !--- coupler
        color = 1
      endif
      if ((rank.ge.atm_rank_min).and.(rank.le.atm_rank_max)) then
        !--- atmosphere
        color = 2
      endif
      if ((rank.ge.ocn_rank_min).and.(rank.le.ocn_rank_max)) then
        !--- ocean
        color = 3
      endif
      if ((rank.ge.ice_rank_min).and.(rank.le.ice_rank_max)) then
        !--- sea ice
        color = 4
      endif

      !--- Build intra-communicator (group_comm) for local sub-groups
      group_comm = MPI_UNDEFINED
      call mpi_comm_split(MPI_COMM_WORLD,color,key,group_comm,ierr)

      !--- The group communicator returned by mpi_comm_split will have the same value
      !--- for every process involved (ie every process in MPI_COMM_WORLD)
      !--- Use this group communicator for collective operations that are to be
      !--- isolated to a particular group
      !--- Save the value of the group communicator for use elswhere
      model_group_comm = group_comm
      if ( verbose > 2 ) then
        write(6,*)"define_group: Task ",rank,"  group_name=",trim(group_name),"  group_comm=",group_comm
        call flush(6)
      endif

      !--- Assign the master task associated with each group
      !--- The master task is always the task with the lowest rank
      !--- Note that if a group is not used (e.g. no ocean) then the master task
      !--- associated with that group will be nproc, which is outside the range of
      !--- mpi tasks in MPI_COMM_WORLD and therefore cannot be used without error
      cpl_master = cpl_rank_min
      atm_master = atm_rank_min
      ocn_master = ocn_rank_min
      ice_master = ice_rank_min

      !--- Assign a complete copy of the task info array for each process
      do proc=0,nproc-1
       if ((proc.le.cpl_rank_max).and.(proc.ge.cpl_rank_min)) then
         task_info(proc)%group_name     = "cpl"
         task_info(proc)%group_min_rank = cpl_rank_min
         task_info(proc)%group_max_rank = cpl_rank_max
         task_info(proc)%master         = task_info(proc)%group_min_rank
       endif
       if ((proc.le.atm_rank_max).and.(proc.ge.atm_rank_min)) then
         task_info(proc)%group_name     = "atm"
         task_info(proc)%group_min_rank = atm_rank_min
         task_info(proc)%group_max_rank = atm_rank_max
         task_info(proc)%master         = task_info(proc)%group_min_rank
       endif
       if ((proc.le.ocn_rank_max).and.(proc.ge.ocn_rank_min)) then
         task_info(proc)%group_name     = "ocn"
         task_info(proc)%group_min_rank = ocn_rank_min
         task_info(proc)%group_max_rank = ocn_rank_max
         task_info(proc)%master         = task_info(proc)%group_min_rank
       endif
       if ((proc.le.ice_rank_max).and.(proc.ge.ice_rank_min)) then
         task_info(proc)%group_name     = "ice"
         task_info(proc)%group_min_rank = ice_rank_min
         task_info(proc)%group_max_rank = ice_rank_max
         task_info(proc)%master         = task_info(proc)%group_min_rank
       endif
      enddo

      !--- Assign the member array
      !--- Depreciated. This is only used with send_fld and recv_fld
      do proc=0,nproc-1
       if ((proc.le.cpl_rank_max).and.(proc.ge.cpl_rank_min)) member(proc)='coupler '
       if ((proc.le.atm_rank_max).and.(proc.ge.atm_rank_min)) member(proc)='atmos   '
       if ((proc.le.ocn_rank_max).and.(proc.ge.ocn_rank_min)) member(proc)='ocean   '
      enddo

      if ( present(group_member) ) then
        group_member(:) = member(:)
      endif

      !--- Wait for all mpi communications to complete
      call mpi_barrier(MPI_COMM_WORLD, ierr)

#undef dev_inter
#ifdef dev_inter
!---TODO--- The ocean is not used when reading atm boundary conditions from a file and this
!---TODO--- needs to be accounted for in the definition of group inter-communicators
!---TODO--- Do not compile the following code until this is sorted out
      if (verbose.gt.2) then
        write(6,*)"define_group: Task ",rank," about to assign inter-communicators."
        call flush(6)
      endif
      !--- Define group inter-communicators
!---TODO--- This is not quite there yet
!---TODO--- (use only 2 inter-communicators which will be assigned appropriate values on each process)
      !--- With 3 groups (cpl, atm, ocn), each group will require 2 inter-communicators
      !--- and there are 3 unique inter-communicators (cpl_atm, cpl_ocn, atm_ocn)
      !---
      !---             cpl_atm
      !---            intercomm
      !---                |
      !---           CPL-----ATM
      !---             \     /
      !---   cpl_ocn -- \   /-- atm_ocn
      !--- intercomm     OCN    intercomm
      !---
      !--- local_leader is the master in MPI_COMM_WORLD
      local_leader = 0
      if ( color == 1 ) then
        !--- coupler
        !--- cpl communicates with atm
        remote_leader = atm_master
        tag = 11
        call mpi_intercomm_create(group_comm, local_leader, MPI_COMM_WORLD, &
                                  remote_leader, tag, intercomm_cpl_atm, ierr)
        !--- cpl communicates with ocn
        remote_leader = ocn_master
        tag = 12
        call mpi_intercomm_create(group_comm, local_leader, MPI_COMM_WORLD, &
                                  remote_leader, tag, intercomm_cpl_ocn, ierr)
      endif
      if ( color == 2 ) then
        !--- atmosphere
        !--- cpl communicates with atm
        remote_leader = cpl_master
        tag = 11
        call mpi_intercomm_create(group_comm, local_leader, MPI_COMM_WORLD, &
                                  remote_leader, tag, intercomm_cpl_atm, ierr)
        !--- atm communicates with ocn
        remote_leader = ocn_master
        tag = 13
        call mpi_intercomm_create(group_comm, local_leader, MPI_COMM_WORLD, &
                                  remote_leader, tag, intercomm_atm_ocn, ierr)
      endif
      if ( color == 3 ) then
        !--- ocean
        !--- cpl communicates with ocn
        remote_leader = cpl_master
        tag = 12
        call mpi_intercomm_create(group_comm, local_leader, MPI_COMM_WORLD, &
                                  remote_leader, tag, intercomm_cpl_ocn, ierr)
        !--- atm communicates with ocn
        remote_leader = atm_master
        tag = 13
        call mpi_intercomm_create(group_comm, local_leader, MPI_COMM_WORLD, &
                                  remote_leader, tag, intercomm_atm_ocn, ierr)
      endif

      if (verbose.gt.1) then
        call mpi_comm_rank(group_comm, group_rank, ierr)
        call mpi_comm_size(group_comm, group_size, ierr)
        !call mpi_barrier(MPI_COMM_WORLD, ierr)

        write(6,'(a,9i8)')"color,rank,group_comm,group_rank,group_size: ", &
                           color,rank,group_comm,group_rank,group_size
        !call mpi_barrier(MPI_COMM_WORLD, ierr)
        call flush(6)
      endif

      if ( intercomm_cpl_atm > 0 .and. verbose > 1 ) then
        call MPI_COMM_TEST_INTER(intercomm_cpl_atm, is_inter, ierr)
        write(6,*)"define_group: is_inter = ",is_inter, &
                  "  intercomm_cpl_atm = ",intercomm_cpl_atm
        call flush(6)
      endif

      if ( intercomm_cpl_ocn > 0 .and. verbose > 1 ) then
        call MPI_COMM_TEST_INTER(intercomm_cpl_ocn, is_inter, ierr)
        write(6,*)"define_group: is_inter = ",is_inter, &
                  "  intercomm_cpl_ocn = ",intercomm_cpl_ocn
        call flush(6)
      endif

      if ( intercomm_atm_ocn > 0 .and. verbose > 1 ) then
        call MPI_COMM_TEST_INTER(intercomm_atm_ocn, is_inter, ierr)
        write(6,*)"define_group: is_inter = ",is_inter, &
                  "  intercomm_atm_ocn = ",intercomm_atm_ocn
        call flush(6)
      endif

      if ( verbose > 2 ) then
        write(6,*)"define_group: task=",rank,"  intercomm_cpl_atm=",intercomm_cpl_atm
        write(6,*)"define_group: task=",rank,"  intercomm_cpl_ocn=",intercomm_cpl_ocn
        write(6,*)"define_group: task=",rank,"  intercomm_atm_ocn=",intercomm_atm_ocn
        call flush(6)
      endif

      if (verbose.gt.2) then
        write(6,*)"define_group: Task ",rank," done assigning inter-communicators."
        call flush(6)
      endif
#endif

      !--- Ensure all processes have arrived before proceeding
      if (verbose.gt.2) then
        write(6,*)"define_group: Task ",rank," waiting for all."
        call flush(6)
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      if (verbose.gt.2) then
        write(6,*)"define_group: Task ",rank," moving on now."
        call flush(6)
      endif

      if ( verbose > 1 ) then
        write(6,*)"define_group: ",trim(group_name)," group_comm,rank: ",group_comm,rank
        call flush(6)
      endif

      if ( verbose > 2 ) then
        do proc=0,nproc-1
          write(6,*)"define_group: rank=",rank,"  name=",trim(task_info(proc)%group_name), &
          "  min_rank, max_rank, master: ", task_info(proc)%group_min_rank, &
          task_info(proc)%group_max_rank, task_info(proc)%master
        enddo
        call flush(6)
      endif

    end subroutine define_group

    !***************************************************************************
    !--- Broadcast values of certain variables to all mpi tasks
    !***************************************************************************
    subroutine cpl_define_global_parameters()
      integer(kind=impi) :: rank, ierr
      integer :: verbose = 2
      integer(kind=8) :: ibuf(8)

      !--- Local
      integer :: idx1

      call mpi_comm_rank ( MPI_COMM_WORLD, rank, ierr )

      if (  rank == cpl_master ) then
        if ( verbose > 1 ) then
          write(6,'(4(a,i6))') "cpl_define_global_parameters: rank=",rank, &
          "  cpl_master=",cpl_master,"  atm_master=",atm_master,"  ocn_master=",ocn_master
          call flush(6)
        endif
      endif

      if (  rank == cpl_master ) then
        !--- Query the environment for the runtime value of certain variables
        call get_environment_variable("cancpl_ver",cpl_runtime_commit)
        call get_environment_variable("cancpl_repo",cpl_runtime_repo_path)
        call get_environment_variable("nemo_ver",nemo_runtime_commit)
        call get_environment_variable("nemo_repo",nemo_runtime_repo_path)
        call get_environment_variable("nemo_config",nemo_runtime_config)
      endif

      if ( verbose > 10 ) then
        write(6,*) "cpl_define_global_parameters: rank=",rank, &
        "  After reading env variables"
        call flush(6)
      endif

      !--- Variables initially defined on the coupler master task
      call bcastGroup(cpl_build_commit_id,   cpl_master, MPI_COMM_WORLD)
      call bcastGroup(cpl_build_repo_path,   cpl_master, MPI_COMM_WORLD)
      call bcastGroup(cpl_runtime_commit,    cpl_master, MPI_COMM_WORLD)
      call bcastGroup(cpl_runtime_repo_path, cpl_master, MPI_COMM_WORLD)
      call bcastGroup(atm_forcing_from_file, cpl_master, MPI_COMM_WORLD)
      call bcastGroup(couple_serial,         cpl_master, MPI_COMM_WORLD)
      call bcastGroup(coupler_parallel_type, cpl_master, MPI_COMM_WORLD)

      if ( verbose > 10 ) then
        write(6,*) "cpl_define_global_parameters: rank=",rank, &
        "  After broadcasting env variables"
        call flush(6)
      endif

#ifndef use_NEMO
      if ( atm_forcing_from_file ) then
        !--- Hard code nemo_recv_var when reading AGCM boundary conditions from a file
        !--- so that these variables will be created in the coupler and are saved to
        !--- the coupler history file
        if ( trim(adjustl(lowerc(task_info(rank)%group_name))) .eq. "cpl" ) then
          nemo_n_recv_var = 15
          if ( associated(nemo_recv_var) ) deallocate(nemo_recv_var)
          allocate( nemo_recv_var(nemo_n_recv_var) )
          nemo_recv_var = (/ "O_OTaux1                        ", &
                             "O_OTauy1                        ", &
                             "O_ITaux1                        ", &
                             "O_ITauy1                        ", &
                             "O_QsrOce                        ", &
                             "O_QsrIce                        ", &
                             "O_QnsOce                        ", &
                             "O_QnsIce                        ", &
                             "OIceEvap                        ", &
                             "OSubMPre                        ", &
                             "OISubMSn                        ", &
                             "OOEvaMPr                        ", &
                             "O_Wind10                        ", &
                             "O_dQnsdT                        ", &
                             "O_Runoff                        " /)
        endif
        if ( verbose > 2 ) then
          write(6,*)"cpl_define_global_parameters: atm_forcing_from_file=",atm_forcing_from_file
          write(6,*)"cpl_define_global_parameters: nemo_n_recv_var=",nemo_n_recv_var
          if ( associated(nemo_recv_var) ) then
             write(6,*)"cpl_define_global_parameters: nemo_recv_var:"
             write(6,'(10(a,2x))')(trim(nemo_recv_var(idx1)),idx1=1,nemo_n_recv_var)
          endif
          call flush(6)
        endif
      endif
#endif

#ifdef use_NEMO
      !--- Variables initially defined on the ocean master task
      call bcastGroup(nemo_nn_ice,           ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_ncat,             ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_rn_rdt,           ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_nn_fsbc,          ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_nn_it000,         ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_nn_itend,         ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_nn_date0,         ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_namsbc_cpl_cldes, ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_commit,           ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_repo_path,        ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_config,           ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_build_cancpl_commit,    ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_build_cancpl_repo_path, ocn_master, MPI_COMM_WORLD)

      if ( verbose > 2 ) then
        write(6,*) "cpl_define_global_parameters: rank=",rank, &
            "  nemo_namsbc_cpl_cldes"
        write(6,'(6x,a)')nemo_namsbc_cpl_cldes
        call flush(6)
      endif

      !--- Fields send from NEMO to the coupler
      !--- These field names are listed in the order they will be transferred
      call bcastGroup(nemo_n_send_var, ocn_master, MPI_COMM_WORLD)
      if ( nemo_n_send_var > 0 ) then
        if ( trim(adjustl(lowerc(task_info(rank)%group_name))) .ne. "ocn" ) then
          if ( associated(nemo_send_var) ) deallocate(nemo_send_var)
          allocate( nemo_send_var(nemo_n_send_var) )
        endif
        call bcastGroup(nemo_send_var, ocn_master, MPI_COMM_WORLD)
      endif

      !--- Fields received by NEMO from the coupler
      !--- These field names are listed in the order they will be transferred
      call bcastGroup(nemo_n_recv_var, ocn_master, MPI_COMM_WORLD)
      if ( nemo_n_recv_var > 0 ) then
        if ( trim(adjustl(lowerc(task_info(rank)%group_name))) .ne. "ocn" ) then
          if ( associated(nemo_recv_var) ) deallocate(nemo_recv_var)
          allocate( nemo_recv_var(nemo_n_recv_var) )
        endif
        call bcastGroup(nemo_recv_var, ocn_master, MPI_COMM_WORLD)
      endif

      call bcastGroup(nemo_jpiglo,     ocn_master, MPI_COMM_WORLD)
      call bcastGroup(nemo_jpjglo,     ocn_master, MPI_COMM_WORLD)
      if ( verbose > 2 ) then
        write(6,*)"cpl_define_global_parameters: nemo_jpiglo,nemo_jpjglo ", &
                  nemo_jpiglo,nemo_jpjglo
        call flush(6)
      endif
      if ( nemo_jpiglo > 0 .and. nemo_jpjglo > 0 ) then
        if ( trim(adjustl(lowerc(task_info(rank)%group_name))) .eq. "cpl" ) then
          !--- This is a coupler node

          !--- NEMO coordinates
          if ( associated(nemo_glamt) ) deallocate(nemo_glamt)
          allocate( nemo_glamt(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_glamu) ) deallocate(nemo_glamu)
          allocate( nemo_glamu(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_glamv) ) deallocate(nemo_glamv)
          allocate( nemo_glamv(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_glamf) ) deallocate(nemo_glamf)
          allocate( nemo_glamf(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_gphit) ) deallocate(nemo_gphit)
          allocate( nemo_gphit(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_gphiu) ) deallocate(nemo_gphiu)
          allocate( nemo_gphiu(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_gphiv) ) deallocate(nemo_gphiv)
          allocate( nemo_gphiv(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_gphif) ) deallocate(nemo_gphif)
          allocate( nemo_gphif(nemo_jpiglo,nemo_jpjglo) )

          !--- NEMO grid scale lengths
          if ( associated(nemo_e1t) ) deallocate(nemo_e1t)
          allocate( nemo_e1t(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_e1u) ) deallocate(nemo_e1u)
          allocate( nemo_e1u(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_e1v) ) deallocate(nemo_e1v)
          allocate( nemo_e1v(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_e1f) ) deallocate(nemo_e1f)
          allocate( nemo_e1f(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_e2t) ) deallocate(nemo_e2t)
          allocate( nemo_e2t(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_e2u) ) deallocate(nemo_e2u)
          allocate( nemo_e2u(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_e2v) ) deallocate(nemo_e2v)
          allocate( nemo_e2v(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_e2f) ) deallocate(nemo_e2f)
          allocate( nemo_e2f(nemo_jpiglo,nemo_jpjglo) )

          !--- NEMO masks
          if ( associated(nemo_tmask) ) deallocate(nemo_tmask)
          allocate( nemo_tmask(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_umask) ) deallocate(nemo_umask)
          allocate( nemo_umask(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_vmask) ) deallocate(nemo_vmask)
          allocate( nemo_vmask(nemo_jpiglo,nemo_jpjglo) )
          if ( associated(nemo_fmask) ) deallocate(nemo_fmask)
          allocate( nemo_fmask(nemo_jpiglo,nemo_jpjglo) )
        endif
        if ( rank == ocn_master .or. &
             trim(adjustl(lowerc(task_info(rank)%group_name))) .eq. "cpl" ) then
          !--- This is a coupler node or the ocean master

          !--- Coupler receives NEMO coordinates from ocean
          call bcast_inter(nemo_glamt, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_glamt shape = ",shape(nemo_glamt)
            write(6,*)"cpl_define_global_parameters: nemo_glamt   min = ",minval(nemo_glamt)
            write(6,*)"cpl_define_global_parameters: nemo_glamt   max = ",maxval(nemo_glamt)
            write(6,*)"cpl_define_global_parameters: nemo_glamt  NaNs = ",count(nemo_glamt /= nemo_glamt)
          endif
          call bcast_inter(nemo_glamu, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_glamu shape = ",shape(nemo_glamu)
            write(6,*)"cpl_define_global_parameters: nemo_glamu   min = ",minval(nemo_glamu)
            write(6,*)"cpl_define_global_parameters: nemo_glamu   max = ",maxval(nemo_glamu)
            write(6,*)"cpl_define_global_parameters: nemo_glamu  NaNs = ",count(nemo_glamu /= nemo_glamu)
          endif
          call bcast_inter(nemo_glamv, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_glamv shape = ",shape(nemo_glamv)
            write(6,*)"cpl_define_global_parameters: nemo_glamv   min = ",minval(nemo_glamv)
            write(6,*)"cpl_define_global_parameters: nemo_glamv   max = ",maxval(nemo_glamv)
            write(6,*)"cpl_define_global_parameters: nemo_glamv  NaNs = ",count(nemo_glamv /= nemo_glamv)
          endif
          call bcast_inter(nemo_glamf, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_glamf shape = ",shape(nemo_glamf)
            write(6,*)"cpl_define_global_parameters: nemo_glamf   min = ",minval(nemo_glamf)
            write(6,*)"cpl_define_global_parameters: nemo_glamf   max = ",maxval(nemo_glamf)
            write(6,*)"cpl_define_global_parameters: nemo_glamf  NaNs = ",count(nemo_glamf /= nemo_glamf)
          endif
          call bcast_inter(nemo_gphit, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_gphit shape = ",shape(nemo_gphit)
            write(6,*)"cpl_define_global_parameters: nemo_gphit   min = ",minval(nemo_gphit)
            write(6,*)"cpl_define_global_parameters: nemo_gphit   max = ",maxval(nemo_gphit)
            write(6,*)"cpl_define_global_parameters: nemo_gphit  NaNs = ",count(nemo_gphit /= nemo_gphit)
          endif
          call bcast_inter(nemo_gphiu, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_gphiu shape = ",shape(nemo_gphiu)
            write(6,*)"cpl_define_global_parameters: nemo_gphiu   min = ",minval(nemo_gphiu)
            write(6,*)"cpl_define_global_parameters: nemo_gphiu   max = ",maxval(nemo_gphiu)
            write(6,*)"cpl_define_global_parameters: nemo_gphiu  NaNs = ",count(nemo_gphiu /= nemo_gphiu)
          endif
          call bcast_inter(nemo_gphiv, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_gphiv shape = ",shape(nemo_gphiv)
            write(6,*)"cpl_define_global_parameters: nemo_gphiv   min = ",minval(nemo_gphiv)
            write(6,*)"cpl_define_global_parameters: nemo_gphiv   max = ",maxval(nemo_gphiv)
            write(6,*)"cpl_define_global_parameters: nemo_gphiv  NaNs = ",count(nemo_gphiv /= nemo_gphiv)
          endif
          call bcast_inter(nemo_gphif, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_gphif shape = ",shape(nemo_gphif)
            write(6,*)"cpl_define_global_parameters: nemo_gphif   min = ",minval(nemo_gphif)
            write(6,*)"cpl_define_global_parameters: nemo_gphif   max = ",maxval(nemo_gphif)
            write(6,*)"cpl_define_global_parameters: nemo_gphif  NaNs = ",count(nemo_gphif /= nemo_gphif)
          endif

          !--- Coupler receives NEMO grid cell scale lengths from ocean
          call bcast_inter(nemo_e1t, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_e1t shape = ",shape(nemo_e1t)
            write(6,*)"cpl_define_global_parameters: nemo_e1t   min = ",minval(nemo_e1t)
            write(6,*)"cpl_define_global_parameters: nemo_e1t   max = ",maxval(nemo_e1t)
            write(6,*)"cpl_define_global_parameters: nemo_e1t  NaNs = ",count(nemo_e1t /= nemo_e1t)
          endif
          call bcast_inter(nemo_e1u, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_e1u shape = ",shape(nemo_e1u)
            write(6,*)"cpl_define_global_parameters: nemo_e1u   min = ",minval(nemo_e1u)
            write(6,*)"cpl_define_global_parameters: nemo_e1u   max = ",maxval(nemo_e1u)
            write(6,*)"cpl_define_global_parameters: nemo_e1u  NaNs = ",count(nemo_e1u /= nemo_e1u)
          endif
          call bcast_inter(nemo_e1v, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_e1v shape = ",shape(nemo_e1v)
            write(6,*)"cpl_define_global_parameters: nemo_e1v   min = ",minval(nemo_e1v)
            write(6,*)"cpl_define_global_parameters: nemo_e1v   max = ",maxval(nemo_e1v)
            write(6,*)"cpl_define_global_parameters: nemo_e1v  NaNs = ",count(nemo_e1v /= nemo_e1v)
          endif
          call bcast_inter(nemo_e1f, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_e1f shape = ",shape(nemo_e1f)
            write(6,*)"cpl_define_global_parameters: nemo_e1f   min = ",minval(nemo_e1f)
            write(6,*)"cpl_define_global_parameters: nemo_e1f   max = ",maxval(nemo_e1f)
            write(6,*)"cpl_define_global_parameters: nemo_e1f  NaNs = ",count(nemo_e1f /= nemo_e1f)
          endif
          call bcast_inter(nemo_e2t, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_e2t shape = ",shape(nemo_e2t)
            write(6,*)"cpl_define_global_parameters: nemo_e2t   min = ",minval(nemo_e2t)
            write(6,*)"cpl_define_global_parameters: nemo_e2t   max = ",maxval(nemo_e2t)
            write(6,*)"cpl_define_global_parameters: nemo_e2t  NaNs = ",count(nemo_e2t /= nemo_e2t)
          endif
          call bcast_inter(nemo_e2u, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_e2u shape = ",shape(nemo_e2u)
            write(6,*)"cpl_define_global_parameters: nemo_e2u   min = ",minval(nemo_e2u)
            write(6,*)"cpl_define_global_parameters: nemo_e2u   max = ",maxval(nemo_e2u)
            write(6,*)"cpl_define_global_parameters: nemo_e2u  NaNs = ",count(nemo_e2u /= nemo_e2u)
          endif
          call bcast_inter(nemo_e2v, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_e2v shape = ",shape(nemo_e2v)
            write(6,*)"cpl_define_global_parameters: nemo_e2v   min = ",minval(nemo_e2v)
            write(6,*)"cpl_define_global_parameters: nemo_e2v   max = ",maxval(nemo_e2v)
            write(6,*)"cpl_define_global_parameters: nemo_e2v  NaNs = ",count(nemo_e2v /= nemo_e2v)
          endif
          call bcast_inter(nemo_e2f, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_e2f shape = ",shape(nemo_e2f)
            write(6,*)"cpl_define_global_parameters: nemo_e2f   min = ",minval(nemo_e2f)
            write(6,*)"cpl_define_global_parameters: nemo_e2f   max = ",maxval(nemo_e2f)
            write(6,*)"cpl_define_global_parameters: nemo_e2f  NaNs = ",count(nemo_e2f /= nemo_e2f)
          endif

          !--- Coupler receives NEMO masks from ocean
          call bcast_inter(nemo_tmask, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_tmask shape = ",shape(nemo_tmask)
            write(6,*)"cpl_define_global_parameters: nemo_tmask zeros = ",count(nemo_tmask == 0.0_8)
            write(6,*)"cpl_define_global_parameters: nemo_tmask  ones = ",count(nemo_tmask == 1.0_8)
            write(6,*)"cpl_define_global_parameters: nemo_tmask  NaNs = ",count(nemo_tmask /= nemo_tmask)
          endif
          call bcast_inter(nemo_umask, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_umask shape = ",shape(nemo_umask)
            write(6,*)"cpl_define_global_parameters: nemo_umask zeros = ",count(nemo_umask == 0.0_8)
            write(6,*)"cpl_define_global_parameters: nemo_umask  ones = ",count(nemo_umask == 1.0_8)
            write(6,*)"cpl_define_global_parameters: nemo_umask  NaNs = ",count(nemo_umask /= nemo_umask)
          endif
          call bcast_inter(nemo_vmask, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            write(6,*)"cpl_define_global_parameters: nemo_vmask shape = ",shape(nemo_vmask)
            write(6,*)"cpl_define_global_parameters: nemo_vmask zeros = ",count(nemo_vmask == 0.0_8)
            write(6,*)"cpl_define_global_parameters: nemo_vmask  ones = ",count(nemo_vmask == 1.0_8)
            write(6,*)"cpl_define_global_parameters: nemo_vmask  NaNs = ",count(nemo_vmask /= nemo_vmask)
          endif
          call bcast_inter(nemo_fmask, ocn_master, "cpl")
          if ( verbose > 1 .and. rank == cpl_master ) then
            !--- Note: fmask may have values that are not 0/1 (e.g. 0.5,2,3)
            write(6,*)"cpl_define_global_parameters: nemo_fmask shape = ",shape(nemo_fmask)
            write(6,*)"cpl_define_global_parameters: nemo_fmask zeros = ",count(nemo_fmask == 0.0_8)
            write(6,*)"cpl_define_global_parameters: nemo_fmask  ones = ",count(nemo_fmask == 1.0_8)
            write(6,*)"cpl_define_global_parameters: nemo_fmask  NaNs = ",count(nemo_fmask /= nemo_fmask)
          endif
        endif
      endif
#endif

      !--- Variables initially defined on the atm master task
      call bcastGroup(atm_delt,    atm_master, MPI_COMM_WORLD)
      call bcastGroup(atm_kstart,  atm_master, MPI_COMM_WORLD)
      call bcastGroup(atm_ksteps,  atm_master, MPI_COMM_WORLD)
      call bcastGroup(atm_kfinal,  atm_master, MPI_COMM_WORLD)
      call bcastGroup(atm_kount,   atm_master, MPI_COMM_WORLD)

      if ( verbose > 5 ) then
        write(6,*) "cpl_define_global_parameters: rank=",rank, &
        "  After broadcasting agcm variables delt, kstart, ksteps, kfinal, kount"
        call flush(6)
      endif

      if ( rank == atm_master .or. &
           trim(adjustl(lowerc(task_info(rank)%group_name))) .eq. "cpl" ) then
        !--- This is a coupler node or the agcm master

        !--- Coupler receives the AGCM earth radius (meters)
        call bcast_inter(atm_earth_radius, atm_master, "cpl")
        if ( verbose > 1 .and. rank == cpl_master ) then
          write(6,*)"cpl_define_global_parameters: atm_earth_radius = ",atm_earth_radius
        endif

        !--- Coupler receives grid shape from agcm
        call bcast_inter(atm_nlon, atm_master, "cpl")
        call bcast_inter(atm_nlat, atm_master, "cpl")
        if (  rank == cpl_master ) then
          !--- Verify AGCM dimension sizes at run time
          if ( atm_nlon /= nlon_a ) then
            write(6,*)"cpl_define_global_parameters: *** ERROR *** atm_nlon does not equal nlon_a"
            write(6,*)"  atm_nlon=",atm_nlon,"  nlon_a=",nlon_a
            call flush(6)
            call err_exit("cpl_define_global_parameters",-1)
          endif
          if ( atm_nlat /= nlat_a ) then
            write(6,*)"cpl_define_global_parameters: *** ERROR *** atm_nlat does not equal nlat_a"
            write(6,*)"  atm_nlat=",atm_nlat,"  nlat_a=",nlat_a
            call flush(6)
            call err_exit("cpl_define_global_parameters",-2)
          endif
        endif

        !--- Coupler receives grid cell areas from agcm
        if ( trim(adjustl(lowerc(task_info(rank)%group_name))) .eq. "cpl" ) then
          !--- This is any coupler task
          if ( associated(atm_grid_cell_area) ) deallocate(atm_grid_cell_area)
          allocate( atm_grid_cell_area(atm_nlon,atm_nlat) )
        endif
        call bcast_inter(atm_grid_cell_area, atm_master, "cpl")
        if ( verbose > 1 .and. rank == cpl_master ) then
          write(6,*)"cpl_define_global_parameters: atm_grid_cell_area shape = ",shape(atm_grid_cell_area)
          write(6,*)"cpl_define_global_parameters: atm_grid_cell_area   min = ",minval(atm_grid_cell_area)
          write(6,*)"cpl_define_global_parameters: atm_grid_cell_area   max = ",maxval(atm_grid_cell_area)
          write(6,*)"cpl_define_global_parameters: atm_grid_cell_area  NaNs = ", &
                    count(atm_grid_cell_area /= atm_grid_cell_area)
          write(6,*)"cpl_define_global_parameters: atm_grid_cell_area column 1 = "
          write(6,'(5g24.16)')(atm_grid_cell_area(1,idx1),idx1=1,atm_nlat)
          write(6,*)"cpl_define_global_parameters: atm_grid_cell_area row 1 = "
          write(6,'(5g24.16)')(atm_grid_cell_area(idx1,1),idx1=1,atm_nlon)
        endif

        if ( trim(adjustl(lowerc(task_info(rank)%group_name))) .eq. "cpl" ) then
          !--- This is any coupler task
          if ( associated(atm_lon) ) deallocate(atm_lon)
          allocate( atm_lon(atm_nlon) )
          if ( associated(atm_lat) ) deallocate(atm_lat)
          allocate( atm_lat(atm_nlat) )
          if ( associated(atm_gwt) ) deallocate(atm_gwt)
          allocate( atm_gwt(atm_nlat) )
        endif
        call bcast_inter(atm_lon, atm_master, "cpl")
        if ( verbose > 1 .and. rank == cpl_master ) then
          write(6,*)"cpl_define_global_parameters: atm_lon shape = ",shape(atm_lon)
          write(6,*)"cpl_define_global_parameters: atm_lon   min = ",minval(atm_lon)
          write(6,*)"cpl_define_global_parameters: atm_lon   max = ",maxval(atm_lon)
          write(6,*)"cpl_define_global_parameters: atm_lon  NaNs = ",count(atm_lon /= atm_lon)
          write(6,*)"cpl_define_global_parameters: atm_lon = "
          write(6,'(5g24.16)')atm_lon
        endif
        call bcast_inter(atm_lat, atm_master, "cpl")
        if ( verbose > 1 .and. rank == cpl_master ) then
          write(6,*)"cpl_define_global_parameters: atm_lat shape = ",shape(atm_lat)
          write(6,*)"cpl_define_global_parameters: atm_lat   min = ",minval(atm_lat)
          write(6,*)"cpl_define_global_parameters: atm_lat   max = ",maxval(atm_lat)
          write(6,*)"cpl_define_global_parameters: atm_lat  NaNs = ",count(atm_lat /= atm_lat)
          write(6,*)"cpl_define_global_parameters: atm_lat = "
          write(6,'(5g24.16)')atm_lat
        endif
        call bcast_inter(atm_gwt, atm_master, "cpl")
        if ( verbose > 1 .and. rank == cpl_master ) then
          write(6,*)"cpl_define_global_parameters: atm_gwt shape = ",shape(atm_gwt)
          write(6,*)"cpl_define_global_parameters: atm_gwt   min = ",minval(atm_gwt)
          write(6,*)"cpl_define_global_parameters: atm_gwt   max = ",maxval(atm_gwt)
          write(6,*)"cpl_define_global_parameters: atm_gwt  NaNs = ",count(atm_gwt /= atm_gwt)
          write(6,*)"cpl_define_global_parameters: atm_gwt = "
          write(6,'(5g24.16)')atm_gwt
        endif

      endif

      !--- Fields send from the AGCM to the coupler
      call bcastGroup(atm_n_send_var, atm_master, MPI_COMM_WORLD)
      if ( atm_n_send_var > 0 ) then
        if ( trim(adjustl(lowerc(task_info(rank)%group_name))) .ne. "atm" ) then
          !--- These field names are defined in the agcm so only allocate
          !--- when we are not in the agcm
          if ( associated(atm_send_var) ) deallocate(atm_send_var)
          allocate( atm_send_var(atm_n_send_var) )
        endif
        call bcastGroup(atm_send_var, atm_master, MPI_COMM_WORLD)
      endif

      !--- Fields send from the coupler to the AGCM
      call bcastGroup(atm_n_recv_var, atm_master, MPI_COMM_WORLD)
      if ( atm_n_recv_var > 0 ) then
        if ( trim(adjustl(lowerc(task_info(rank)%group_name))) .ne. "atm" ) then
          !--- These field names are defined in the agcm so only allocate
          !--- when we are not in the agcm
          if ( associated(atm_recv_var) ) deallocate(atm_recv_var)
          allocate( atm_recv_var(atm_n_recv_var) )
        endif
        call bcastGroup(atm_recv_var, atm_master, MPI_COMM_WORLD)
      endif

      if ( verbose > 1 ) then
        write(6,*) "cpl_define_global_parameters: rank=",rank, &
            "  atm_n_send_var=",atm_n_send_var,"  atm_n_recv_var=",atm_n_recv_var
        if ( atm_n_send_var > 0 ) then
          write(6,*) "cpl_define_global_parameters: atm_send_var"
          write(6,'(5(2x,a))')(atm_send_var(idx1),idx1=1,atm_n_send_var)
          call flush(6)
        endif
        if ( atm_n_recv_var > 0 ) then
          write(6,*) "cpl_define_global_parameters: atm_recv_var"
          write(6,'(5(2x,a))')(atm_recv_var(idx1),idx1=1,atm_n_recv_var)
          call flush(6)
        endif
      endif

      !--- Land mask from the AGCM
      if (  rank == cpl_master ) then
        !--- Only the coupler gets the atm land mask
        !--- NOTE: atm_lmask (and all related fields) contain 1 overlap longitude

        if ( associated(atm_lmask) ) deallocate(atm_lmask)
        allocate( atm_lmask((nlon_a+olap_a)*nlat_a) )

        !--- Receive the fractional land mask from the atm
        if ( associated(atm_fland) ) deallocate(atm_fland)
        allocate( atm_fland(size(atm_lmask)) )
        call recv_data_rec(atm_fland, ibuf, atm_master, "FLND_atm", dbg=1)
        if ( verbose > 1 ) then
          write(6,*)"cpl_define_global_parameters: atm_fland shape   = ",shape(atm_fland)
          write(6,*)"cpl_define_global_parameters: min/max atm_fland = ", &
              minval(atm_fland),maxval(atm_fland)
          call flush(6)
        endif

        !--- atm_lmask is zero in grid cells where there is no land
        !--- and one in grid cells where there is some land
        atm_lmask = 1.0_8
        where (atm_fland <= 0.0_8) atm_lmask = 0.0_8

        !--- Assign atm_cv_lmask based on the fractional land mask recieved from the atm
        if ( associated(atm_cv_lmask) ) deallocate(atm_cv_lmask)
        allocate( atm_cv_lmask(size(atm_lmask)) )
        atm_cv_lmask = 0
        where ( atm_fland >= 1.0_8 ) atm_cv_lmask = 2

      endif

      if ( verbose > 2 ) then
        write(6,*) "cpl_define_global_parameters: rank=",rank, &
          " nn_ice, ncat, rn_rdt, nn_it000, nn_itend, nn_date0: ", &
          nemo_nn_ice, nemo_ncat, nemo_rn_rdt, nemo_nn_it000, nemo_nn_itend, nemo_nn_date0

        write(6,*)"cpl_define_global_parameters: rank=",rank, &
          "  atm_kstart, atm_ksteps, atm_kfinal, atm_delt: ", &
          atm_kstart, atm_ksteps, atm_kfinal, atm_delt

        if ( nemo_n_send_var > 0 ) then
          write(6,*)"cpl_define_global_parameters: rank=",rank,"  nemo_n_send_var=",nemo_n_send_var
          write(6,'(5(2x,a))')nemo_send_var(1:nemo_n_send_var)
        endif
        if ( nemo_n_recv_var > 0 ) then
          write(6,*)"cpl_define_global_parameters: rank=",rank,"  nemo_n_recv_var=",nemo_n_recv_var
          write(6,'(5(2x,a))')nemo_recv_var(1:nemo_n_recv_var)
        endif
        call flush(6)
      endif

      if ( verbose > 2 ) then
        write(6,*) "cpl_define_global_parameters: rank=",rank, &
        "  About to return to caller"
        call flush(6)
      endif

    end subroutine cpl_define_global_parameters

    !***************************************************************************
    !--- Determine data record size and tag, given the record name
    !--- Run some sanity checks and define task ID strings
    !***************************************************************************
    subroutine verify_name(rec_name, rem_proc, curr_size, curr_tag, my_taskid, rem_taskid)

      implicit none

      !--- The name associated with the record to be transferred
      !--- This name must exist in cpl_var_list (defined in this module)
      character(*), intent(in) :: rec_name

      !--- Rank of the remote MPI process from or to which data will be transferred
      integer(kind=impi), intent(in)  :: rem_proc

      !--- Length of the data record to be transferred (this does not include ibuf)
      integer, intent(out) :: curr_size

      !--- The MPI message tag associated with this record
      integer, intent(out) :: curr_tag

      !--- A string to contain a task ID for the current process
      character(*), intent(out) :: my_taskid

      !--- A string to contain a task ID for the remote process
      character(*), intent(out) :: rem_taskid

      !--- Local
      integer(kind=impi) :: myrank
      type(cpl_vinfo_t)  :: cpl_vinfo

      !-------------------------------------------------------------------------

      if ( .not. com_mpi_initialized ) call com_mpi_init()

      !--- Determine the tag and data array size to be used with this field
      cpl_vinfo = find_cpl_vinfo( name=trim(adjustl(rec_name)) )
      curr_size = cpl_vinfo%size
      curr_tag  = cpl_vinfo%tag

      !--- Sanity checks on input/output parameters

      if ( rem_proc .lt. 0 ) then
        write(6,'(a,i8)')'verify_name: Invalid rem_proc task id = ',rem_proc
        call flush(6)
        call err_exit("VERIFY_NAME",-1)
      endif

      if ( curr_size .lt. 0 ) then
        write(6,'(a,i8)')'verify_name: Invalid data size = ',curr_size
        call flush(6)
        call err_exit("VERIFY_NAME",-2)
      endif

      if ( curr_tag .lt. 0 ) then
        write(6,'(a,i8)')'verify_name: Invalid tag = ',curr_tag
        call flush(6)
        call err_exit("VERIFY_NAME",-3)
      endif
 
      !--- Verify that we are not receiving from self
      call mpi_comm_rank( MPI_COMM_WORLD, myrank, ierr )
      if (myrank == rem_proc) call err_exit('VERIFY_NAME',-4)

      !--- Communication with another executable is only done from the master task
      !--- Verify that this is one of the master tasks
      if ( myrank.ne.cpl_master .and. &
           myrank.ne.atm_master .and. &
           myrank.ne.ocn_master .and. &
           myrank.ne.ice_master ) then
        write(6,'(a,i4,a)')'verify_name: The current mpi task --> ', &
          myrank,' <-- is not a master task.'
        call err_exit("VERIFY_NAME",-5)
      endif

      !--- Verify that rem_proc is one of the master tasks
      if ( rem_proc.ne.cpl_master .and. &
           rem_proc.ne.atm_master .and. &
           rem_proc.ne.ocn_master .and. &
           rem_proc.ne.ice_master ) then
        write(6,'(a,i4,a)')'verify_name: The remote mpi task --> ', &
          rem_proc,' <-- is not a master task.'
        call err_exit("VERIFY_NAME",-6)
      endif

      !--- Define a string containing the current task ID
      my_taskid = " "
      if ( myrank == atm_master ) then
        my_taskid = "ATM"
      else if ( myrank == ocn_master ) then
        my_taskid = "OCN"
      else if ( myrank == cpl_master ) then
        my_taskid = "CPL"
      else if ( myrank == ice_master ) then
        my_taskid = "ICE"
      else
        my_taskid = "UNK"
      endif

      !--- Define a string containing the remote task ID
      rem_taskid = " "
      if ( rem_proc == atm_master ) then
        rem_taskid = "ATM"
      else if ( rem_proc == ocn_master ) then
        rem_taskid = "OCN"
      else if ( rem_proc == cpl_master ) then
        rem_taskid = "CPL"
      else if ( rem_proc == ice_master ) then
        rem_taskid = "ICE"
      else
        rem_taskid = "UNK"
      endif

    end subroutine verify_name

    !***************************************************************************
    !--- Define interface routines for generic recv_data_rec and send_data_rec
    !--- Each interface will do what is required to get the user supplied data
    !--- array as real and call recv_data_rec_r1d or send_data_rec_r1d
    !***************************************************************************
    subroutine recv_data_rec_i1d (i_data, ibuf, origin, rec_name, dbg)
      !-------------------------------------------------------------------------
      !--- Receive 1D integer data array
      !-------------------------------------------------------------------------
      integer(kind=8),     intent(out) :: i_data(:)
      integer(kind=8),     intent(out) :: ibuf(:)
      integer (kind=impi), intent(in)  :: origin
      character(*),        intent(in)  :: rec_name
      integer,             intent(in), optional :: dbg
      !--- Local
      type(cpl_vinfo_t)         :: cpl_vinfo
      integer                   :: curr_size
      real(kind=8), allocatable :: wrk(:)
      !-------------------------------------------------------------------------
      cpl_vinfo = find_cpl_vinfo( name=trim(adjustl(rec_name)) )
      curr_size = cpl_vinfo%size
      if ( size(i_data) .lt. curr_size ) then
        write(6,'(a)')'recv_data_rec: Input data array is to small.'
        write(6,'(a,i10)')'           input size=',size(i_data)
        write(6,'(a,i10)')'  current record size=',curr_size
        call flush(6)
        call err_exit("RECV_DATA_REC_I1D",-1)
      endif
      allocate( wrk(curr_size) )
      if ( present(dbg) ) then
        call recv_data_rec_r1d(wrk, ibuf, origin, rec_name, dbg=dbg)
      else
        call recv_data_rec_r1d(wrk, ibuf, origin, rec_name)
      endif
      i_data(1:curr_size) = wrk(1:curr_size)
      deallocate( wrk )
    end subroutine recv_data_rec_i1d

    subroutine recv_data_rec_i1d_noibuf (i_data, origin, rec_name, dbg)
      !-------------------------------------------------------------------------
      !--- Receive 1D integer data array
      !-------------------------------------------------------------------------
      integer(kind=8),     intent(out) :: i_data(:)
      integer (kind=impi), intent(in)  :: origin
      character(*),        intent(in)  :: rec_name
      integer,             intent(in), optional :: dbg
      !--- Local
      integer(kind=8)           :: ibuf(8)
      type(cpl_vinfo_t)         :: cpl_vinfo
      integer                   :: curr_size, nc
      character(4)              :: recv_name
      real(kind=8), allocatable :: wrk(:)
      !-------------------------------------------------------------------------
      cpl_vinfo = find_cpl_vinfo( name=trim(adjustl(rec_name)) )
      curr_size = cpl_vinfo%size
      if ( size(i_data) .lt. curr_size ) then
        write(6,'(a)')'recv_data_rec: Input data array is to small.'
        write(6,'(a,i10)')'           input size=',size(i_data)
        write(6,'(a,i10)')'  current record size=',curr_size
        call flush(6)
        call err_exit("RECV_DATA_REC_I1D",-1)
      endif

      !--- Use the first 4 chars (at most) of the field name for ibuf3 and
      !--- ensure the ibuf3 name is right justified if less than 4 chars
      nc = min(4,len_trim(rec_name))
      if ( nc <= 0 ) then
        write(6,*)"recv_data_rec_i1d_noibuf: Missing field name."
        call err_exit("RECV_DATA_REC_I1D_NOIBUF",-1)
      endif
      recv_name = " "
      recv_name(1:nc) = rec_name(1:nc)
      recv_name = adjustr(recv_name)

      !--- Define ibuf on the fly for this transfer
      ibuf(1) = transfer("GRID",1_8)
      ibuf(2) = 0
      ibuf(3) = transfer(recv_name,1_8)
      ibuf(4) = 1
      ibuf(5) = cpl_vinfo%nlon + cpl_vinfo%olap
      ibuf(6) = cpl_vinfo%nlat
      ibuf(7) = 0
      ibuf(8) = 1

      allocate( wrk(curr_size) )
      if ( present(dbg) ) then
        call recv_data_rec_r1d(wrk, ibuf, origin, rec_name, dbg=dbg)
      else
        call recv_data_rec_r1d(wrk, ibuf, origin, rec_name)
      endif
      i_data(1:curr_size) = wrk(1:curr_size)
      deallocate( wrk )

    end subroutine recv_data_rec_i1d_noibuf

    subroutine recv_data_rec_r2d (r_data, ibuf, origin, rec_name, dbg)
      !-------------------------------------------------------------------------
      !--- Receive 2D real data array
      !-------------------------------------------------------------------------
      real(kind=8),        intent(out) :: r_data(:,:)
      integer(kind=8),     intent(out) :: ibuf(:)
      integer (kind=impi), intent(in)  :: origin
      character(*),        intent(in)  :: rec_name
      integer,             intent(in), optional :: dbg
      !--- Local
      type(cpl_vinfo_t) :: cpl_vinfo
      integer           :: curr_size, dims(2)
      real(kind=8), allocatable :: wrk(:)
      !-------------------------------------------------------------------------
      cpl_vinfo = find_cpl_vinfo( name=trim(adjustl(rec_name)) )
      curr_size = cpl_vinfo%size
      dims(1:2) = shape( r_data )
      if ( curr_size .ne. dims(1)*dims(2) ) then
        write(6,'(2a)') &
          'recv_data_rec_r2d: Input 2D array size is inconsistent with ',trim(rec_name)
        write(6,'(a,i8,a,2i8)') &
          '  record size: ',curr_size,'  array shape: ',dims(1),dims(2)
        call err_exit("RECV_DATA_REC_R2D",-1)
      endif
      allocate( wrk(curr_size) )
      if ( present(dbg) ) then
        call recv_data_rec_r1d(wrk, ibuf, origin, rec_name, dbg=dbg)
      else
        call recv_data_rec_r1d(wrk, ibuf, origin, rec_name)
      endif
      r_data(1:dims(1),1:dims(2)) = &
        reshape( wrk(1:dims(1)*dims(2)), (/dims(1),dims(2)/) )
      deallocate( wrk )
    end subroutine recv_data_rec_r2d

    subroutine recv_data_rec_r3d (r_data, ibuf, origin, rec_name, dbg)
!!! recv_data_rec_r3d is Untested !!!
      !-------------------------------------------------------------------------
      !--- Receive 3D real data array
      !-------------------------------------------------------------------------
      real(kind=8),        intent(out) :: r_data(:,:,:)
      integer(kind=8),     intent(out) :: ibuf(:)
      integer (kind=impi), intent(in)  :: origin
      character(*),        intent(in)  :: rec_name
      integer,             intent(in), optional :: dbg
      !--- Local
      type(cpl_vinfo_t) :: cpl_vinfo
      integer           :: curr_size, dims(3)
      real(kind=8), allocatable :: wrk(:)
      !-------------------------------------------------------------------------
      cpl_vinfo = find_cpl_vinfo( name=trim(adjustl(rec_name)) )
      curr_size = cpl_vinfo%size
      dims(1:3) = shape( r_data )
      if ( curr_size .ne. dims(1)*dims(2)*dims(3) ) then
        write(6,'(2a)') &
          'recv_data_rec_r3d: Input 3D array size is inconsistent with ',trim(rec_name)
        write(6,'(a,i8,a,2i8)') &
          '  record size: ',curr_size,'  array shape: ',dims(1),dims(2),dims(3)
        call err_exit("RECV_DATA_REC_R3D",-1)
      endif
      allocate( wrk(curr_size) )
      if ( present(dbg) ) then
        call recv_data_rec_r1d(wrk, ibuf, origin, rec_name, dbg=dbg)
      else
        call recv_data_rec_r1d(wrk, ibuf, origin, rec_name)
      endif
      r_data(1:dims(1),1:dims(2),1:dims(3)) = &
        reshape( wrk(1:dims(1)*dims(2)*dims(3)), (/dims(1),dims(2),dims(3)/) )
      deallocate( wrk )
    end subroutine recv_data_rec_r3d

    subroutine send_data_rec_i1d (i_data, ibuf, dest, rec_name, dbg)
      !-------------------------------------------------------------------------
      !--- Send 1D integer data array
      !-------------------------------------------------------------------------
      integer(kind=8),     intent(in) :: i_data(:)
      integer(kind=8),     intent(in) :: ibuf(:)
      integer (kind=impi), intent(in) :: dest
      character(*),        intent(in) :: rec_name
      integer,             intent(in), optional :: dbg
      !--- Local
      type(cpl_vinfo_t) :: cpl_vinfo
      integer           :: curr_size
      real(kind=8), allocatable :: wrk(:)
      !-------------------------------------------------------------------------
      cpl_vinfo = find_cpl_vinfo( name=trim(adjustl(rec_name)) )
      curr_size = cpl_vinfo%size
      if ( size(i_data) .lt. curr_size ) then
        write(6,'(a)')'send_data_rec: Input data array is to small.'
        write(6,'(a,i10)')'           input size=',size(i_data)
        write(6,'(a,i10)')'  current record size=',curr_size
        call flush(6)
        call err_exit("SEND_DATA_REC_I1D",-1)
      endif
      allocate( wrk(curr_size) )
      wrk(1:curr_size) = i_data(1:curr_size)
      if ( present(dbg) ) then
        call send_data_rec_r1d(wrk, ibuf, dest, rec_name, dbg=dbg)
      else
        call send_data_rec_r1d(wrk, ibuf, dest, rec_name)
      endif
      deallocate( wrk )
    end subroutine send_data_rec_i1d

    subroutine send_data_rec_i1d_noibuf (i_data, dest, rec_name, dbg)
      !-------------------------------------------------------------------------
      !--- Send 1D integer data array --- No user supplied ibuf
      !-------------------------------------------------------------------------
      integer(kind=8),     intent(in) :: i_data(:)
      integer (kind=impi), intent(in) :: dest
      character(*),        intent(in) :: rec_name
      integer,             intent(in), optional :: dbg
      !--- Local
      integer(kind=8)   :: ibuf(8)
      type(cpl_vinfo_t) :: curr_vinfo
      integer           :: curr_size, nc
      character(len=4)  :: send_name
      real(kind=8), allocatable :: wrk(:)
      !-------------------------------------------------------------------------
      curr_vinfo = find_cpl_vinfo( name=trim(adjustl(rec_name)) )
      curr_size = curr_vinfo%size
      if ( size(i_data) .lt. curr_size ) then
        write(6,'(a)')'send_data_rec: Input data array is to small.'
        write(6,'(a,i10)')'           input size=',size(i_data)
        write(6,'(a,i10)')'  current record size=',curr_size
        call flush(6)
        call err_exit("SEND_DATA_REC_I1D",-1)
      endif
      allocate( wrk(curr_size) )
      wrk(1:curr_size) = i_data(1:curr_size)

      !--- Use the first 4 chars (at most) of the field name for ibuf3 and
      !--- ensure the ibuf3 name is right justified if less than 4 chars
      nc = min(4,len_trim(rec_name))
      if ( nc <= 0 ) then
        write(6,*)"send_data_rec_i1d_noibuf: Missing field name."
        call err_exit("SEND_DATA_REC_I1D_NOIBUF",-1)
      endif
      send_name = " "
      send_name(1:nc) = rec_name(1:nc)
      send_name = adjustr(send_name)

      !--- Define ibuf on the fly for this transfer
      ibuf(1) = transfer("GRID",1_8)
      ibuf(2) = 0
      ibuf(3) = transfer(send_name,1_8)
      ibuf(4) = 1
      ibuf(5) = curr_vinfo%nlon + curr_vinfo%olap
      ibuf(6) = curr_vinfo%nlat
      ibuf(7) = 0
      ibuf(8) = 1

      if ( present(dbg) ) then
        call send_data_rec_r1d(wrk, ibuf, dest, rec_name, dbg=dbg)
      else
        call send_data_rec_r1d(wrk, ibuf, dest, rec_name)
      endif
      deallocate( wrk )
    end subroutine send_data_rec_i1d_noibuf

    subroutine send_data_rec_r2d (r_data, ibuf, dest, rec_name, dbg)
      !-------------------------------------------------------------------------
      !--- Send 2D real data array
      !-------------------------------------------------------------------------
      real(kind=8),        intent(in) :: r_data(:,:)
      integer(kind=8),     intent(in) :: ibuf(:)
      integer (kind=impi), intent(in) :: dest
      character(*),        intent(in) :: rec_name
      integer,             intent(in), optional :: dbg
      !--- Local
      type(cpl_vinfo_t) :: cpl_vinfo
      integer           :: curr_size, dims(2)
      real(kind=8), allocatable :: wrk(:)
      !-------------------------------------------------------------------------
      cpl_vinfo = find_cpl_vinfo( name=trim(adjustl(rec_name)) )
      curr_size = cpl_vinfo%size
      dims(1:2) = shape( r_data )
      if ( curr_size .ne. dims(1)*dims(2) ) then
        write(6,'(2a)') &
          'send_data_rec_r2d: Input 2D array size is inconsistent with ',trim(rec_name)
        write(6,'(a,i8,a,2i8)') &
          '  record size: ',curr_size,'  array shape: ',dims(1),dims(2)
        call err_exit("SEND_DATA_REC_R2D",-1)
      endif
      allocate( wrk(curr_size) )
      wrk(1:dims(1)*dims(2)) = &
        reshape( r_data(1:dims(1),1:dims(2)), (/dims(1)*dims(2)/) )
      if ( present(dbg) ) then
        call send_data_rec_r1d(wrk, ibuf, dest, rec_name, dbg=dbg)
      else
        call send_data_rec_r1d(wrk, ibuf, dest, rec_name)
      endif
      deallocate( wrk )
    end subroutine send_data_rec_r2d

    subroutine send_data_rec_r2d_noibuf (r_data, dest, rec_name, dbg)
      !-------------------------------------------------------------------------
      !--- Send 2D real data array --- No user supplied ibuf
      !-------------------------------------------------------------------------
      real(kind=8),        intent(in) :: r_data(:,:)
      integer (kind=impi), intent(in) :: dest
      character(*),        intent(in) :: rec_name
      integer,             intent(in), optional :: dbg
      !--- Local
      integer(kind=8)   :: ibuf(8)
      type(cpl_vinfo_t) :: curr_vinfo
      integer           :: curr_size, dims(2), nc
      character(len=4)  :: send_name
      real(kind=8), allocatable :: wrk(:)
      !-------------------------------------------------------------------------
      curr_vinfo = find_cpl_vinfo( name=trim(adjustl(rec_name)) )
      curr_size = curr_vinfo%size
      dims(1:2) = shape( r_data )
      if ( curr_size .ne. dims(1)*dims(2) ) then
        write(6,'(2a)') &
          'send_data_rec_r2d_noibuf: Input 2D array size is inconsistent with ',trim(rec_name)
        write(6,'(a,i8,a,2i8)') &
          '  record size: ',curr_size,'  array shape: ',dims(1),dims(2)
        call err_exit("SEND_DATA_REC_R2D_noibuf",-1)
      endif
      allocate( wrk(curr_size) )
      wrk(1:dims(1)*dims(2)) = &
        reshape( r_data(1:dims(1),1:dims(2)), (/dims(1)*dims(2)/) )

      !--- Use the first 4 chars (at most) of the field name for ibuf3 and
      !--- ensure the ibuf3 name is right justified if less than 4 chars
      nc = min(4,len_trim(rec_name))
      if ( nc <= 0 ) then
        write(6,*)"send_data_rec_i2d_noibuf: Missing field name."
        call err_exit("SEND_DATA_REC_I2D_NOIBUF",-1)
      endif
      send_name = " "
      send_name(1:nc) = rec_name(1:nc)
      send_name = adjustr(send_name)

      !--- Define ibuf on the fly for this transfer
      ibuf(1) = transfer("GRID",1_8)
      ibuf(2) = 0
      ibuf(3) = transfer(send_name,1_8)
      ibuf(4) = 1
      ibuf(5) = curr_vinfo%nlon + curr_vinfo%olap
      ibuf(6) = curr_vinfo%nlat
      ibuf(7) = 0
      ibuf(8) = 1

      if ( present(dbg) ) then
        call send_data_rec_r1d(wrk, ibuf, dest, rec_name, dbg=dbg)
      else
        call send_data_rec_r1d(wrk, ibuf, dest, rec_name)
      endif
      deallocate( wrk )
    end subroutine send_data_rec_r2d_noibuf

    subroutine send_data_rec_r3d (r_data, ibuf, dest, rec_name, dbg)
!!! send_data_rec_r3d is Untested !!!
      !-------------------------------------------------------------------------
      !--- Send 3D real data array
      !-------------------------------------------------------------------------
      real(kind=8),        intent(in) :: r_data(:,:,:)
      integer(kind=8),     intent(in) :: ibuf(:)
      integer (kind=impi), intent(in) :: dest
      character(*),        intent(in) :: rec_name
      integer,             intent(in), optional :: dbg
      !--- Local
      type(cpl_vinfo_t) :: curr_vinfo
      integer           :: curr_size, dims(3)
      real(kind=8), allocatable :: wrk(:)
      !-------------------------------------------------------------------------
      curr_vinfo = find_cpl_vinfo( name=trim(adjustl(rec_name)) )
      curr_size = curr_vinfo%size
      dims(1:3) = shape( r_data )
      if ( curr_size .ne. dims(1)*dims(2)*dims(3) ) then
        write(6,'(2a)') &
          'send_data_rec_r3d: Input 3D array size is inconsistent with ',trim(rec_name)
        write(6,'(a,i8,a,2i8)') &
          '  record size: ',curr_size,'  array shape: ',dims(1),dims(2),dims(3)
        call err_exit("SEND_DATA_REC_R3D",-1)
      endif
      allocate( wrk(curr_size) )
      wrk(1:dims(1)*dims(2)*dims(3)) = &
        reshape( r_data(1:dims(1),1:dims(2),1:dims(3)), (/dims(1)*dims(2)*dims(3)/) )
      if ( present(dbg) ) then
        call send_data_rec_r1d(wrk, ibuf, dest, rec_name, dbg=dbg)
      else
        call send_data_rec_r1d(wrk, ibuf, dest, rec_name)
      endif
      deallocate( wrk )
    end subroutine send_data_rec_r3d

    !***************************************************************************
    !--- Receive a single logical record (ibuf + data) via MPI
    !***************************************************************************
    subroutine recv_data_rec_r1d (rdata, ibuf, origin, rec_name, dbg)

      implicit none

      !--- Real data array to be communicated (received) via MPI
      real(kind=8), intent(out) :: rdata(:)

      !--- Integer array containing info about the field to be received
      integer(kind=8), intent(out) :: ibuf(:)

      !--- Rank of the MPI process from which the data was sent
      integer (kind=impi), intent(in) :: origin

      !--- The name associated with the record to be received
      character(*), intent(in) :: rec_name

      !--- An optional integer flag used to determine the amount
      !--- of diagnostic output written to stdout
      integer, intent(in), optional :: dbg

      !--- Local
      integer :: rdata_size
      integer :: tag
      integer (kind=impi) :: tag_impi
      !--- buffer to contain data that will be received
      real, dimension(:), allocatable :: recv_buf
      integer :: verbose
      character(len=512) :: strnga, strngb
      integer (kind=impi) :: recv_size
      integer (kind=impi) :: status(MPI_status_size)
      character(len=16) :: recvid, sendid

      !-----------------------------------------------------------

      if ( present(dbg) ) then
        verbose = dbg
      else
        verbose = 0
      endif

      if ( .not. com_mpi_initialized ) call com_mpi_init()

!xxx      !--- Determine if this field is on the event stack
!xxx      if ( .not. is_coupled(rec_name) ) then
!xxx        !--- Do not receive the data
!xxx        write(6,'(3a)')"recv_data_rec: The variable ",trim(rec_name), &
!xxx                       " is not coupled. No data was received."
!xxx        call flush(6)
!xxx        return
!xxx      endif

      !--- Determine data record size and tag, given the record name
      !--- rdata_size is the length of the data to be received (without ibuf)
      !--- tag is the MPI tag associated with this record
      call verify_name(rec_name, origin, rdata_size, tag, recvid, sendid)

      if ( size(rdata) .lt. rdata_size ) then
        write(6,'(a,i8)')'recv_data_rec: rdata is to small.  size=',size(rdata)
        call flush(6)
        call err_exit("RECV_DATA_REC",-1)
      endif

      if ( size(ibuf) .lt. nbuf ) then
        write(6,'(a,i8)')'recv_data_rec: ibuf is too small.  size=',size(ibuf)
        call flush(6)
        call err_exit("RECV_DATA_REC",-2)
      endif

      !--- The record to be received will contain ibuf and the data array
      recv_size = rdata_size + nbuf
      allocate( recv_buf(recv_size) )

      if ( verbose > 10 ) then
        write(6,'(3a,i8,a,i6)')"recv_data_rec_r1d: Receiving ",trim(rec_name), &
                               "  tag=",tag,"  from ",origin
        call flush(6)
      endif

      !--- Receive data
      tag_impi = tag
      curr_real_type=MPI_REAL4
      call mpi_recv( recv_buf, recv_size, curr_real_type, origin, tag_impi, &
                     MPI_COMM_WORLD, status, ierr )

      !-- Extract ibuf from the array that was received
      ibuf(1:nbuf) = recv_buf(1:nbuf)

      !-- Extract data from the array that was received
      rdata(1:rdata_size) = recv_buf(nbuf+1:nbuf+rdata_size)

      !--- Clean up
      deallocate(recv_buf)

      if ( verbose > -10 ) then
        !--- Write diagnostic info to stdout
        strnga = " "
        strngb = " "
        strngb = trim(adjustl(rec_name))
        write(strnga,'(9a,i4)') &
          '=== ',cpl_time_string,' ',trim(recvid),' recvs from ', &
          trim(sendid),' --> ',strngb(1:10),' <-- tag=',tag
!xxx        strngb = trim( format_ibuf(ibuf) )
        strngb = sprint_var_stats(rdata, rdata_size, name=rec_name)
        write(6,'(a,2x,a)') trim(strnga),trim(strngb)

        if ( verbose.gt.1 ) then
          write(6,'(3a,1x,5g14.5)') &
            '=== ',trim(adjustl(rec_name)),' data(1:5)=',rdata(1:5)
        endif
        call flush(6)
      endif

    end subroutine recv_data_rec_r1d

    !***************************************************************************
    !--- Send a single logical record via MPI --- No user supplied ibuf
    !--- Define ibuf here then call send_data_rec_r1d
    !***************************************************************************
    subroutine send_data_rec_r1d_noibuf (rdata, dest, rec_name, dbg)

      implicit none

      !--- Real data array to be communicated (received) via MPI
      real(kind=8), intent(in) :: rdata(:)

      !--- Rank of the MPI process from which the data was sent
      integer (kind=impi), intent(in)  :: dest

      !--- The name associated with the record to be received
      character(*), intent(in) :: rec_name

      !--- An optional integer flag used to determine the amount
      !--- of diagnostic output written to stdout
      integer, intent(in), optional :: dbg

      !--- Local
      integer(kind=8)   :: ibuf(8)
      type(cpl_vinfo_t) :: curr_vinfo
      integer           :: nc
      character(len=4)  :: send_name
      !-------------------------------------------------------------------------

      !--- Access information about this field
      curr_vinfo = find_cpl_vinfo( name=trim(adjustl(rec_name)) )

      !--- Use the first 4 chars (at most) of the field name for ibuf3 and
      !--- ensure the ibuf3 name is right justified if less than 4 chars
      nc = min(4,len_trim(rec_name))
      if ( nc <= 0 ) then
        write(6,*)"send_data_rec_r1d_noibuf: Missing field name."
        call err_exit("SEND_DATA_REC_R1D_NOIBUF",-1)
      endif
      send_name = " "
      send_name(1:nc) = rec_name(1:nc)
      send_name = adjustr(send_name)

      !--- Define ibuf on the fly for this transfer
      ibuf(1) = transfer("GRID",1_8)
      ibuf(2) = 0
      ibuf(3) = transfer(send_name,1_8)
      ibuf(4) = 1
      ibuf(5) = curr_vinfo%nlon + curr_vinfo%olap
      ibuf(6) = curr_vinfo%nlat
      ibuf(7) = 0
      ibuf(8) = 1

      if ( present(dbg) ) then
        call send_data_rec_r1d(rdata, ibuf, dest, rec_name, dbg=dbg)
      else
        call send_data_rec_r1d(rdata, ibuf, dest, rec_name)
      endif

    end subroutine send_data_rec_r1d_noibuf

    !***************************************************************************
    !--- Send a single logical record (ibuf + data) via MPI
    !***************************************************************************
    subroutine send_data_rec_r1d (rdata, ibuf, destination, rec_name, dbg)

      implicit none

      !--- Real data array to be communicated (received) via MPI
      real(kind=8), intent(in) :: rdata(:)

      !--- Integer array containing info about the field to be received
      integer(kind=8), intent(in) :: ibuf(:)

      !--- Rank of the MPI process from which the data was sent
      integer (kind=impi), intent(in)  :: destination

      !--- The name associated with the record to be received
      character(*), intent(in) :: rec_name

      !--- An optional integer flag used to determine the amount
      !--- of diagnostic output written to stdout
      integer, intent(in), optional :: dbg

      !--- Local
      integer :: rdata_size
      integer :: tag
      integer (kind=impi) :: tag_impi
      !--- buffer to contain data that will be sent
      real, dimension(:), allocatable :: send_buf
      integer :: verbose
      character(len=512) :: strnga, strngb
      integer (kind=impi) :: send_size
      character(len=16) :: recvid, sendid

      !-----------------------------------------------------------

      if ( present(dbg) ) then
        verbose = dbg
      else
        verbose = 0
      endif

      if ( .not. com_mpi_initialized ) call com_mpi_init()

!xxx      !--- Determine if this field is on the event stack
!xxx      if ( .not. is_coupled(rec_name) ) then
!xxx        !--- Do not send the data
!xxx        write(6,'(3a)')"send_data_rec: The variable ",trim(rec_name), &
!xxx                       " is not coupled. No data was sent."
!xxx        call flush(6)
!xxx        return
!xxx      endif

      !--- Determine data record size and tag, given the record name
      !--- rdata_size is the length of the data to be sent (without ibuf)
      !--- tag is the MPI tag associated with this record
      call verify_name(rec_name, destination, rdata_size, tag, sendid, recvid)

      if ( size(rdata) .lt. rdata_size ) then
        write(6,'(a,i8)')'send_data_rec: rdata is to small.  size=',size(rdata)
        call flush(6)
        call err_exit("SEND_DATA_REC",-1)
      endif

      if ( size(ibuf) .lt. nbuf ) then
        write(6,'(a,i8)')'send_data_rec: ibuf is too small.  size=',size(ibuf)
        call flush(6)
        call err_exit("SEND_DATA_REC",-2)
      endif

      !--- The record to be sent will contain ibuf and the data array
      send_size = rdata_size + nbuf
      allocate( send_buf(send_size) )

      !-- Insert ibuf into the array that will be sent
      send_buf(1:nbuf) = ibuf(1:nbuf)

      !-- Insert the real data into the array that will be sent
      send_buf(nbuf+1:nbuf+rdata_size) = rdata(1:rdata_size)

      if ( verbose > 10 ) then
        write(6,'(3a,i8,a,i6)')"send_data_rec_r1d: Sending ",trim(rec_name), &
            "  tag=",tag,"  to ",destination
        call flush(6)
      endif

      !--- Send data
      tag_impi = tag
      curr_real_type=MPI_REAL4
      call mpi_ssend( send_buf, send_size, curr_real_type, destination, &
                      tag_impi, MPI_COMM_WORLD, ierr )

      !--- Clean up
      deallocate(send_buf) 

      if ( verbose > -10 ) then
        !--- Write diagnostic info to stdout
        strnga = " "
        strngb = " "
        strngb = trim(adjustl(rec_name))
        write(strnga,'(9a,i4)') &
          '=== ',cpl_time_string,' ',trim(sendid),' sends  to  ', &
          trim(recvid),' --> ',strngb(1:10),' <-- tag=',tag
!xxx        strngb = trim( format_ibuf(ibuf) )
        strngb = sprint_var_stats(rdata, rdata_size, name=rec_name)
        write(6,'(a,2x,a)') trim(strnga),trim(strngb)

        if ( verbose.gt.1 ) then
          write(6,'(3a,1x,5g14.5)') &
            '=== ',trim(adjustl(rec_name)),' data(1:5)=',rdata(1:5)
        endif
        call flush(6)
      endif

    end subroutine send_data_rec_r1d

    !***************************************************************************
    !--- recv_fld with a real 1D data array as output
    !***************************************************************************
    subroutine recv_fld_r1d (r_data, ibuf, origin, info_arr)

      !--- Data to be communicated (received) via MPI
      real(kind=8), intent(out) :: r_data(:)

      !--- 8 word integer array containing info about the field to be received
      integer(kind=8), intent(out) :: ibuf(:)

      !--- Rank of the MPI process from which the data was sent
      integer(kind=impi), intent(in)  :: origin

      !--- info_arr(1) = tagid
      !--- info_arr(2) = NOT USED -- currently the same as info_arr(1)
      !--- info_arr(3) = max size of the array
      integer(kind=8), intent(in)  :: info_arr(3)

      !--- Local
      type(cpl_vinfo_t) :: cpl_vinfo
      integer :: tag
      character(len=32) :: rec_name
      !-------------------------------------------------------------------------
      tag = info_arr(1)
      cpl_vinfo = find_cpl_vinfo( tag=tag )
      rec_name = trim(cpl_vinfo%name)

      call recv_data_rec_r1d(r_data, ibuf, origin, rec_name)

    end subroutine recv_fld_r1d

    !***************************************************************************
    !--- recv_fld with a real 2D data array as output
    !***************************************************************************
    subroutine recv_fld_r2d (r_data2d, ibuf, origin, info_arr)

      !--- Data to be communicated (received) via MPI
      real(kind=8), intent(out) :: r_data2d(:,:)

      !--- 8 word integer array containing info about the field to be received
      integer(kind=8), intent(out) :: ibuf(:)

      !--- Rank of the MPI process from which the data was sent
      integer(kind=impi), intent(in)  :: origin

      !--- info_arr(1) = tagid
      !--- info_arr(2) = NOT USED -- currently the same as info_arr(1)
      !--- info_arr(3) = max size of the array
      integer(kind=8), intent(in)  :: info_arr(3)

      !--- Local
      type(cpl_vinfo_t) :: cpl_vinfo
      integer :: tag
      character(len=32) :: rec_name
      !-------------------------------------------------------------------------
      tag = info_arr(1)
      cpl_vinfo = find_cpl_vinfo( tag=tag )
      rec_name = trim(cpl_vinfo%name)

      call recv_data_rec_r2d(r_data2d, ibuf, origin, rec_name)

    end subroutine recv_fld_r2d

    !***************************************************************************
    !--- recv_fld with a 1D integer data array as output
    !***************************************************************************
    subroutine recv_fld_i1d (i_data, ibuf, origin, info_arr)

      implicit none

      !--- Data to be communicated (received) via MPI
      integer(kind=8), intent(out) :: i_data(:)

      !--- 8 word integer array containing info about the field to be received
      integer(kind=8), intent(out) :: ibuf(:)

      !--- Rank of the MPI process from which the data was sent
      integer(kind=impi), intent(in)  :: origin

      !--- info_arr(1) = tagid
      !--- info_arr(2) = NOT USED -- currently the same as info_arr(1)
      !--- info_arr(3) = max size of the array
      integer(kind=8), intent(in)  :: info_arr(3)

      !--- Local
      type(cpl_vinfo_t) :: cpl_vinfo
      integer :: tag
      character(len=32) :: rec_name
      !-------------------------------------------------------------------------
      tag = info_arr(1)
      cpl_vinfo = find_cpl_vinfo( tag=tag )
      rec_name = trim(cpl_vinfo%name)

      call recv_data_rec_i1d(i_data, ibuf, origin, rec_name)

    end subroutine recv_fld_i1d

    !***************************************************************************
    !--- send_fld with a real 1D data array as input
    !***************************************************************************
    subroutine send_fld_r1d (r_data, ibuf, destid, info_arr)

      !--- Data to be communicated (sent) via MPI
      real(kind=8), intent(in) :: r_data(:)

      !--- 8 word integer array containing info about the field to be sent
      integer(kind=8), intent(in) :: ibuf(:)

      !--- Rank of the MPI process to which the data will be sent
      integer (kind=impi), intent(in) :: destid

      !--- info_arr(1) = tagid
      !--- info_arr(2) = NOT USED -- currently the same as info_arr(1)
      !--- info_arr(3) = max size of the array
      integer(kind=8), intent(in) :: info_arr(3)

      !--- Local
      type(cpl_vinfo_t) :: cpl_vinfo
      integer :: tag
      character(len=32) :: rec_name
      !-------------------------------------------------------------------------
      tag = info_arr(1)
      cpl_vinfo = find_cpl_vinfo( tag=tag )
      rec_name = trim(cpl_vinfo%name)

      call send_data_rec_r1d(r_data, ibuf, destid, rec_name)

    end subroutine send_fld_r1d

    !***************************************************************************
    !--- send_fld with a real 2D data array as input
    !***************************************************************************
    subroutine send_fld_r2d (r_data2d, ibuf, destid, info_arr)

      !--- Data to be communicated (sent) via MPI
      real(kind=8), intent(in) :: r_data2d(:,:)

      !--- 8 word integer array containing info about the field to be sent
      integer(kind=8), intent(in) :: ibuf(:)

      !--- Rank of the MPI process to which the data will be sent
      integer (kind=impi), intent(in) :: destid

      !--- info_arr(1) = tagid
      !--- info_arr(2) = NOT USED -- currently the same as info_arr(1)
      !--- info_arr(3) = max size of the array
      integer(kind=8), intent(in) :: info_arr(3)

      !--- Local
      type(cpl_vinfo_t) :: cpl_vinfo
      integer :: tag
      character(len=32) :: rec_name
      !-------------------------------------------------------------------------
      tag = info_arr(1)
      cpl_vinfo = find_cpl_vinfo( tag=tag )
      rec_name = trim(cpl_vinfo%name)

      call send_data_rec_r2d(r_data2d, ibuf, destid, rec_name)

    end subroutine send_fld_r2d

    !***************************************************************************
    !--- send_fld with a 1D integer data array as input
    !***************************************************************************
    subroutine send_fld_i1d (i_data, ibuf, destid, info_arr)

      !--- Data to be communicated (sent) via MPI
      integer(kind=8), intent(in) :: i_data(:)

      !--- 8 word integer array containing info about the field to be sent
      integer(kind=8), intent(in) :: ibuf(:)

      !--- Rank of the MPI process to which the data will be sent
      integer (kind=impi), intent(in) :: destid

      !--- info_arr(1) = tagid
      !--- info_arr(2) = NOT USED -- currently the same as info_arr(1)
      !--- info_arr(3) = max size of the array
      integer(kind=8), intent(in) :: info_arr(3)

      !--- Local
      type(cpl_vinfo_t) :: cpl_vinfo
      integer :: tag
      character(len=32) :: rec_name
      !-------------------------------------------------------------------------
      tag = info_arr(1)
      cpl_vinfo = find_cpl_vinfo( tag=tag )
      rec_name = trim(cpl_vinfo%name)

      !--- Send the desired record (as a real array) to destid
      call send_data_rec_i1d(i_data, ibuf, destid, rec_name)

    end subroutine send_fld_i1d

    !-----------------------------------------------------------------------------
    !------------------- Values to be broadcast are integer(4) -------------------
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    !--- Broadcast a single integer(4) value to all tasks in a given comm world
    !--- from a specified task
    !--- If comm_world is an intra-communicator then it must be MPI_COMM_WORLD
    !--- If 
    !-----------------------------------------------------------------------------
    subroutine bcastGroup_Scalar_i4(val, root, comm_world)
 
      !--- The integer to be broadcast
      integer(kind=4), intent(inout) :: val

      !--- The process to broadcast from
      integer(impi), intent(in) :: root

      !--- The comm world to broadcast to
      integer(impi), intent(in) :: comm_world
 
      !--- Local
      integer(impi), parameter :: one=1
      integer(impi) :: ierr
 
      call mpi_bcast(val, one, MPI_INTEGER4, root, comm_world, ierr)
      call mpi_barrier(comm_world, ierr)

    end subroutine bcastGroup_Scalar_i4

    !-----------------------------------------------------------------------------
    !--- Broadcast an array of integer(4) value to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcastGroup_1d_Array_i4(val, root, comm_world, start, count)
 
      !--- The integer array to be broadcast
      integer(kind=4), intent(inout) :: val(:)

      !--- The process to broadcast from
      integer(impi), intent(in) :: root

      !--- The comm world to broadcast to
      integer(impi), intent(in) :: comm_world

      !--- Starting array element in val (default is the lower bound of val)
      integer(kind=8), intent(in), optional :: start

      !--- The number of elements to send (default is the entire array)
      integer(kind=8), intent(in), optional :: count
 
      !--- Local
      integer :: lb1, ub1, idx1, idx2
      integer(impi) :: send_size, ierr

      lb1 = lbound(val, dim=1)
      ub1 = ubound(val, dim=1)

      if ( present(start) ) then
        idx1 = start
        if ( idx1 < lb1 .or. idx1 > ub1 ) then
          write(6,*)'bcastGroup_1d_Array_i4: Start index ',idx1,' is out of range.'
          write(6,*)'   Valid bounds: ',lb1,ub1
          call err_exit("bcastGroup_1d_Array_i4",-1)
        endif
      else
        idx1 = lb1
      endif

      if ( present(count) ) then
        send_size = count
      else
        send_size = ub1 - idx1 + 1
      endif

      idx2 = idx1 + send_size - 1
      if ( idx2 < lb1 .or. idx2 > ub1 ) then
        write(6,*)'bcastGroup_1d_Array_i4: End index ',idx2,' is out of range.'
        write(6,*)'   Valid bounds: ',lb1,ub1,'  Send size=',send_size
        call err_exit("bcastGroup_1d_Array_i4",-2)
      endif

      call mpi_bcast(val(idx1:idx2), send_size, MPI_INTEGER4, root, comm_world, ierr)
      call mpi_barrier(comm_world, ierr)

    end subroutine bcastGroup_1d_Array_i4

    !-----------------------------------------------------------------------------
    !--- Broadcast an array of integer(4) value to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcastGroup_2d_Array_i4(val, root, comm_world, start, count)
 
      !--- The integer array to be broadcast
      integer(kind=4), intent(inout) :: val(:,:)

      !--- The process to broadcast from
      integer(impi), intent(in) :: root

      !--- The comm world to broadcast to
      integer(impi), intent(in) :: comm_world

      !--- Starting array element in val (default is the lower bound of val)
      integer(kind=8), intent(in), optional :: start(:)

      !--- The number of elements to send (default is the entire array)
      integer(kind=8), intent(in), optional :: count(:)
 
      !--- Local
      integer :: lb1, ub1, lb2, ub2
      integer :: ix(2), iy(2)
      integer(impi) :: send_size, ierr

      lb1 = lbound(val, dim=1)
      ub1 = ubound(val, dim=1)
      lb2 = lbound(val, dim=2)
      ub2 = ubound(val, dim=2)

      if ( present(start) ) then
        ix(1) = start(1)
        if ( ix(1) < lb1 .or. ix(1) > ub1 ) then
          write(6,*)'bcastGroup_2d_Array_i4: Start index 1 --> ',ix(1),' is out of range.'
          write(6,*)'   Valid bounds: ',lb1,ub1
          call err_exit("bcastGroup_2d_Array_i4",-1)
        endif
        iy(1) = start(2)
        if ( ix(2) < lb2 .or. ix(2) > ub2 ) then
          write(6,*)'bcastGroup_2d_Array_i4: Start index 2 --> ',iy(1),' is out of range.'
          write(6,*)'   Valid bounds: ',lb2,ub2
          call err_exit("bcastGroup_2d_Array_i4",-2)
        endif
      else
        ix(1) = lb1
        iy(1) = lb2
      endif

      if ( present(count) ) then
        ix(2) = ix(1) + count(1) - 1
        iy(2) = iy(1) + count(2) - 1
        send_size = count(1)*count(2)
      else
        ix(2) = ub1
        iy(2) = ub2
        send_size = (ub1 - ix(1) + 1) * (ub2 - iy(1) + 1)
      endif

      if ( ix(2) < lb1 .or. ix(2) > ub1 ) then
        write(6,*)'bcastGroup_2d_Array_i4: End index ',ix(2),' is out of range.'
        write(6,*)'   Valid bounds: ',lb1,ub1,'  Send size=',send_size
        call err_exit("bcastGroup_2d_Array_i4",-3)
      endif
      if ( iy(2) < lb2 .or. iy(2) > ub2 ) then
        write(6,*)'bcastGroup_2d_Array_i4: End index ',iy(2),' is out of range.'
        write(6,*)'   Valid bounds: ',lb2,ub2,'  Send size=',send_size
        call err_exit("bcastGroup_2d_Array_i4",-4)
      endif
      call mpi_bcast(val(ix(1):ix(2),iy(1):iy(2)), send_size, MPI_INTEGER4, &
                     root, comm_world, ierr)
      call mpi_barrier(comm_world, ierr)

    end subroutine bcastGroup_2d_Array_i4

    !-----------------------------------------------------------------------------
    !------------------- Values to be broadcast are integer(8) -------------------
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    !--- Broadcast a single integer(8) value to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcastGroup_Scalar_i8(val, root, comm_world)
 
      !--- The integer to be broadcast
      integer(kind=8), intent(inout) :: val

      !--- The process to broadcast from
      integer(impi), intent(in) :: root

      !--- The comm world to broadcast to
      integer(impi), intent(in) :: comm_world
 
      !--- Local
      integer(impi), parameter :: one=1
      integer(impi) :: ierr
 
      call mpi_bcast(val, one, MPI_INTEGER8, root, comm_world, ierr)
      call mpi_barrier(comm_world, ierr)

    end subroutine bcastGroup_Scalar_i8

    !-----------------------------------------------------------------------------
    !--- Broadcast an array of integer(8) values to all tasks in a given
    !--- comm world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcastGroup_1d_Array_i8(val, root, comm_world, start, count)
 
      !--- The integer to be broadcast
      integer(kind=8), intent(inout) :: val(:)

      !--- The process to broadcast from
      integer(impi), intent(in) :: root

      !--- The comm world to broadcast to
      integer(impi), intent(in) :: comm_world

      !--- Starting array element in val (default is the lower bound of val)
      integer(kind=8), intent(in), optional :: start

      !--- The number of elements to send (default is the entire array)
      integer(kind=8), intent(in), optional :: count
 
      !--- Local
      integer :: lb1, ub1, idx1, idx2
      integer(impi) :: send_size, ierr

      lb1 = lbound(val, dim=1)
      ub1 = ubound(val, dim=1)

      if ( present(start) ) then
        idx1 = start
        if ( idx1 < lb1 .or. idx1 > ub1 ) then
          write(6,*)'bcastGroup_1d_Array_i8: Start index ',idx1,' is out of range.'
          write(6,*)'   Valid bounds: ',lb1,ub1
          call err_exit("bcastGroup_1d_Array_i8",-1)
        endif
      else
        idx1 = lb1
      endif

      if ( present(count) ) then
        send_size = count
      else
        send_size = ub1 - idx1 + 1
      endif

      idx2 = idx1 + send_size - 1
      if ( idx2 < lb1 .or. idx2 > ub1 ) then
        write(6,*)'bcastGroup_1d_Array_i8: End index ',idx2,' is out of range.'
        write(6,*)'   Valid bounds: ',lb1,ub1,'  Send size=',send_size
        call err_exit("bcastGroup_1d_Array_i8",-2)
      endif
 
      call mpi_bcast(val(idx1:idx2), send_size, MPI_INTEGER8, root, comm_world, ierr)
      call mpi_barrier(comm_world, ierr)

    end subroutine bcastGroup_1d_Array_i8

    !-----------------------------------------------------------------------------
    !------------------- Values to be broadcast are real(4) ----------------------
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    !--- Broadcast a single real(4) value to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcastGroup_Scalar_r4(val, root, comm_world)
 
      !--- The real variable to be broadcast
      real(kind=4), intent(inout) :: val

      !--- The process to broadcast from
      integer(impi), intent(in) :: root

      !--- The comm world to broadcast to
      integer(impi), intent(in) :: comm_world
 
      !--- Local
      integer(impi), parameter :: one=1
      integer(impi) :: ierr
 
      call MPI_BCAST(val, one, MPI_REAL4, root, comm_world, ierr)
      call MPI_BARRIER(comm_world, ierr)

    end subroutine bcastGroup_Scalar_r4

    !-----------------------------------------------------------------------------
    !------------------- Values to be broadcast are real(8) ----------------------
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    !--- Broadcast a single real(8) value to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcastGroup_Scalar_r8(val, root, comm_world)
 
      !--- The real variable to be broadcast
      real(kind=8), intent(inout) :: val

      !--- The process to broadcast from
      integer(impi), intent(in) :: root

      !--- The comm world to broadcast to
      integer(impi), intent(in) :: comm_world
 
      !--- Local
      integer(impi), parameter :: one=1
      integer(impi) :: ierr
 
      call MPI_BCAST(val, one, MPI_REAL8, root, comm_world, ierr)
      call MPI_BARRIER(comm_world, ierr)

    end subroutine bcastGroup_Scalar_r8

    !-----------------------------------------------------------------------------
    !------------------- Values to be broadcast are character --------------------
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    !--- Broadcast a single character variable to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcastGroup_Scalar_char(val, root, comm_world)
 
      !--- The character variable to be broadcast
      character(*), intent(inout) :: val

      !--- The process to broadcast from
      integer(impi), intent(in) :: root

      !--- The comm world to broadcast to
      integer(impi), intent(in) :: comm_world
 
      !--- Local
      integer(impi), parameter :: one=1
      integer(impi) :: csize, ierr
 
      csize = len(val)
      call mpi_bcast(val, csize, MPI_CHARACTER, root, comm_world, ierr)

      call mpi_barrier(comm_world, ierr)

    end subroutine bcastGroup_Scalar_char

    !-----------------------------------------------------------------------------
    !--- Broadcast a 1d array of character variables to all tasks in a given
    !--- comm world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcastGroup_1d_Array_char(val, root, comm_world)
 
      !--- The character array to be broadcast
      character(*), intent(inout) :: val(:)

      !--- The process to broadcast from
      integer(impi), intent(in) :: root

      !--- The comm world to broadcast to
      integer(impi), intent(in) :: comm_world
 
      !--- Local
      integer(impi), parameter :: one=1
      integer(impi) :: csize, ierr, idx

      !--- There is surely a better way to do this, but efficiently is not
      !--- paramount here
      do idx=1,size(val)
        csize = len(val(idx))
        call mpi_bcast(val(idx), csize, MPI_CHARACTER, root, comm_world, ierr)
        call mpi_barrier(comm_world, ierr)
      enddo

    end subroutine bcastGroup_1d_Array_char

    !-----------------------------------------------------------------------------
    !------------------- Values to be broadcast are logical ----------------------
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    !--- Broadcast a single logical value to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcastGroup_Scalar_logical(val, root, comm_world)
 
      !--- The logical variable to be broadcast
      logical, intent(inout) :: val

      !--- The process to broadcast from
      integer(impi), intent(in) :: root

      !--- The comm world to broadcast to
      integer(impi), intent(in) :: comm_world
 
      !--- Local
      integer(impi), parameter :: one=1
      integer(impi) :: ierr
      integer(kind=4) :: i4val
 
      if (val) then
        i4val = 1
      else
        i4val = 0
      endif

      call MPI_BCAST(i4val, one, MPI_INTEGER4, root, comm_world, ierr)
      call MPI_BARRIER(comm_world, ierr)

      if ( i4val == 1 ) then
        val = .true.
      else
        val = .false.
      endif

    end subroutine bcastGroup_Scalar_logical

!==============================================================================
! Begin definition of routines associated with the generic inteface bcast_inter
!==============================================================================

    !-----------------------------------------------------------------------------
    !--- Given the name of a group return the rank of the group leader (master)
    !-----------------------------------------------------------------------------
    function group_leader(group_name) result(master)

      !--- The name of the group
      character(*), intent(in) :: group_name

      !--- The master task rank in MPI_COMM_WORLD that is associated with this group
      integer :: master

      !--- Local
      integer :: idx

      if ( .not. associated(task_info) ) then
        write(6,*)"group_leader: task_info is not defined."
        call err_exit("group_leader",-1)
      endif

      master = -1
      do idx=1,size(task_info)
        if ( trim(adjustl(lowerc(task_info(idx)%group_name))) .eq. &
             trim(adjustl(lowerc(group_name))) ) then
          !--- All group members will have the same leader so we use the first one
          master = task_info(idx)%master
          exit
        endif
      enddo

      if (master == -1 ) then
        write(6,*)"group_leader: Invalid group name ",trim(group_name)
        call err_exit("group_leader",-2)
      endif

    end function group_leader

    !-----------------------------------------------------------------------------
    !--- Broadcast a single character variable to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_Scalar_char(val, src_rank, remote_group)
 
      !--- The character variable to be broadcast
      character(*), intent(inout) :: val

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast scalar character value
      call bcast_inter_Any(src_rank, remote_group, val_schar=val)

    end subroutine bcast_inter_Scalar_char

    !-----------------------------------------------------------------------------
    !--- Broadcast a 1D character variable array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_1D_Array_char(val, src_rank, remote_group)
 
      !--- The character variable to be broadcast
      character(*), intent(inout) :: val(:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 1D character array
      call bcast_inter_Any(src_rank, remote_group, val_a1char=val)

    end subroutine bcast_inter_1D_array_char

    !-----------------------------------------------------------------------------
    !--- Broadcast a 3D character variable array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_3D_Array_char(val, src_rank, remote_group)
 
      !--- The character variable to be broadcast
      character(*), intent(inout) :: val(:,:,:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 3D character array
      call bcast_inter_Any(src_rank, remote_group, val_a3char=val)

    end subroutine bcast_inter_3D_array_char

    !-----------------------------------------------------------------------------
    !--- Broadcast a 2D character variable array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_2D_Array_char(val, src_rank, remote_group)
 
      !--- The character variable to be broadcast
      character(*), intent(inout) :: val(:,:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 2D character array
      call bcast_inter_Any(src_rank, remote_group, val_a2char=val)

    end subroutine bcast_inter_2D_array_char

    !-----------------------------------------------------------------------------
    !--- Broadcast a single logical variable to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_Scalar_logical(val, src_rank, remote_group)
 
      !--- The logical variable to be broadcast
      logical, intent(inout) :: val

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast scalar character value
      call bcast_inter_Any(src_rank, remote_group, val_sbool=val)

    end subroutine bcast_inter_Scalar_logical

    !-----------------------------------------------------------------------------
    !--- Broadcast a 1D logical array to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_1d_Array_logical(val, src_rank, remote_group)
 
      !--- The logical variable to be broadcast
      logical, intent(inout) :: val(:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 1d logical array
      call bcast_inter_Any(src_rank, remote_group, val_a1bool=val)

    end subroutine bcast_inter_1d_Array_logical

    !-----------------------------------------------------------------------------
    !--- Broadcast a 2D logical array to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_2d_Array_logical(val, src_rank, remote_group)
 
      !--- The logical variable to be broadcast
      logical, intent(inout) :: val(:,:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 2d logical array
      call bcast_inter_Any(src_rank, remote_group, val_a2bool=val)

    end subroutine bcast_inter_2d_Array_logical

    !-----------------------------------------------------------------------------
    !--- Broadcast a 3D logical array to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_3d_Array_logical(val, src_rank, remote_group)
 
      !--- The logical variable to be broadcast
      logical, intent(inout) :: val(:,:,:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 3d logical array
      call bcast_inter_Any(src_rank, remote_group, val_a3bool=val)

    end subroutine bcast_inter_3d_Array_logical

    !-----------------------------------------------------------------------------
    !--- Broadcast a single integer*4 variable to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_Scalar_i4(val, src_rank, remote_group)
 
      !--- The integer variable to be broadcast
      integer(4), intent(inout) :: val

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast scalar integer kind=4
      call bcast_inter_Any(src_rank, remote_group, val_si4=val)

    end subroutine bcast_inter_Scalar_i4

    !-----------------------------------------------------------------------------
    !--- Broadcast a 1D integer*4 array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_1D_Array_i4(val, src_rank, remote_group)
 
      !--- The integer variable to be broadcast
      integer(4), intent(inout) :: val(:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 1D integer kind=4 array
      call bcast_inter_Any(src_rank, remote_group, val_a1i4=val)

    end subroutine bcast_inter_1D_array_i4

    !-----------------------------------------------------------------------------
    !--- Broadcast a 2D integer*4 array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_2D_Array_i4(val, src_rank, remote_group)
 
      !--- The integer variable to be broadcast
      integer(4), intent(inout) :: val(:,:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 2D integer kind=4 array
      call bcast_inter_Any(src_rank, remote_group, val_a2i4=val)

    end subroutine bcast_inter_2D_array_i4


    !-----------------------------------------------------------------------------
    !--- Broadcast a 3D integer*4 array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_3D_Array_i4(val, src_rank, remote_group)
 
      !--- The integer variable to be broadcast
      integer(4), intent(inout) :: val(:,:,:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 3D integer kind=4 array
      call bcast_inter_Any(src_rank, remote_group, val_a3i4=val)

    end subroutine bcast_inter_3D_array_i4

    !-----------------------------------------------------------------------------
    !--- Broadcast a single integer*8 variable to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_Scalar_i8(val, src_rank, remote_group)
 
      !--- The integer variable to be broadcast
      integer(8), intent(inout) :: val

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast scalar integer kind=8
      call bcast_inter_Any(src_rank, remote_group, val_si8=val)

    end subroutine bcast_inter_Scalar_i8

    !-----------------------------------------------------------------------------
    !--- Broadcast a 1D integer*8 array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_1D_Array_i8(val, src_rank, remote_group)
 
      !--- The integer variable to be broadcast
      integer(8), intent(inout) :: val(:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 1D integer kind=8 array
      call bcast_inter_Any(src_rank, remote_group, val_a1i8=val)

    end subroutine bcast_inter_1D_array_i8

    !-----------------------------------------------------------------------------
    !--- Broadcast a 2D integer*8 array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_2D_Array_i8(val, src_rank, remote_group)
 
      !--- The integer variable to be broadcast
      integer(8), intent(inout) :: val(:,:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 2D integer kind=8 array
      call bcast_inter_Any(src_rank, remote_group, val_a2i8=val)

    end subroutine bcast_inter_2D_array_i8

    !-----------------------------------------------------------------------------
    !--- Broadcast a 3D integer*8 array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_3D_Array_i8(val, src_rank, remote_group)
 
      !--- The integer variable to be broadcast
      integer(8), intent(inout) :: val(:,:,:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 3D integer kind=8 array
      call bcast_inter_Any(src_rank, remote_group, val_a3i8=val)

    end subroutine bcast_inter_3D_array_i8

    !-----------------------------------------------------------------------------
    !--- Broadcast a single real*4 variable to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_Scalar_r4(val, src_rank, remote_group)
 
      !--- The real variable to be broadcast
      real(4), intent(inout) :: val

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast scalar real kind=4
      call bcast_inter_Any(src_rank, remote_group, val_sr4=val)

    end subroutine bcast_inter_Scalar_r4

    !-----------------------------------------------------------------------------
    !--- Broadcast a 1D real*4 array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_1D_Array_r4(val, src_rank, remote_group)
 
      !--- The real variable to be broadcast
      real(4), intent(inout) :: val(:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 1D real kind=4 array
      call bcast_inter_Any(src_rank, remote_group, val_a1r4=val)

    end subroutine bcast_inter_1D_array_r4

    !-----------------------------------------------------------------------------
    !--- Broadcast a 2D real*4 array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_2D_Array_r4(val, src_rank, remote_group)
 
      !--- The real variable to be broadcast
      real(4), intent(inout) :: val(:,:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 2D real kind=4 array
      call bcast_inter_Any(src_rank, remote_group, val_a2r4=val)

    end subroutine bcast_inter_2D_array_r4

    !-----------------------------------------------------------------------------
    !--- Broadcast a 3D real*4 array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_3D_Array_r4(val, src_rank, remote_group)
 
      !--- The real variable to be broadcast
      real(4), intent(inout) :: val(:,:,:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 3D real kind=4 array
      call bcast_inter_Any(src_rank, remote_group, val_a3r4=val)

    end subroutine bcast_inter_3D_array_r4

    !-----------------------------------------------------------------------------
    !--- Broadcast a single real*8 variable to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_Scalar_r8(val, src_rank, remote_group)
 
      !--- The real variable to be broadcast
      real(8), intent(inout) :: val

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast scalar real kind=4
      call bcast_inter_Any(src_rank, remote_group, val_sr8=val)

    end subroutine bcast_inter_Scalar_r8

    !-----------------------------------------------------------------------------
    !--- Broadcast a 1D real*8 array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_1D_Array_r8(val, src_rank, remote_group)
 
      !--- The real variable to be broadcast
      real(8), intent(inout) :: val(:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 1D real kind=4 array
      call bcast_inter_Any(src_rank, remote_group, val_a1r8=val)

    end subroutine bcast_inter_1D_array_r8

    !-----------------------------------------------------------------------------
    !--- Broadcast a 2D real*8 array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_2D_Array_r8(val, src_rank, remote_group)
 
      !--- The real variable to be broadcast
      real(8), intent(inout) :: val(:,:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 2D real kind=4 array
      call bcast_inter_Any(src_rank, remote_group, val_a2r8=val)

    end subroutine bcast_inter_2D_array_r8

    !-----------------------------------------------------------------------------
    !--- Broadcast a 3D real*8 array to all tasks in a given comm
    !--- world from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_3D_Array_r8(val, src_rank, remote_group)
 
      !--- The real variable to be broadcast
      real(8), intent(inout) :: val(:,:,:)

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group

      !--- Broadcast 3D real kind=4 array
      call bcast_inter_Any(src_rank, remote_group, val_a3r8=val)

    end subroutine bcast_inter_3D_array_r8

    !-----------------------------------------------------------------------------
    !--- Broadcast a variable to all tasks in a given comm world
    !--- from a specified task
    !-----------------------------------------------------------------------------
    subroutine bcast_inter_Any(src_rank, remote_group, &
                 val_schar,  val_sbool,  val_si4,  val_si8,  val_sr4,  val_sr8,  &
                 val_a1char, val_a1bool, val_a1i4, val_a1i8, val_a1r4, val_a1r8, &
                 val_a2char, val_a2bool, val_a2i4, val_a2i8, val_a2r4, val_a2r8, &
                 val_a3char, val_a3bool, val_a3i4, val_a3i8, val_a3r4, val_a3r8  )

      !--- The task that will send the data (this cannot be in remote_group)
      integer(impi), intent(in) :: src_rank

      !--- The name of the group to receive the broadcast
      character(*), intent(in) :: remote_group
 
      !--- The variable to be broadcast (only one type allowed per call)
      character(*), intent(inout), optional :: val_schar
      logical,      intent(inout), optional :: val_sbool
      integer(4),   intent(inout), optional :: val_si4
      integer(8),   intent(inout), optional :: val_si8
      real(4),      intent(inout), optional :: val_sr4
      real(8),      intent(inout), optional :: val_sr8
      character(*), intent(inout), optional :: val_a1char(:)
      logical,      intent(inout), optional :: val_a1bool(:)
      integer(4),   intent(inout), optional :: val_a1i4(:)
      integer(8),   intent(inout), optional :: val_a1i8(:)
      real(4),      intent(inout), optional :: val_a1r4(:)
      real(8),      intent(inout), optional :: val_a1r8(:)
      character(*), intent(inout), optional :: val_a2char(:,:)
      logical,      intent(inout), optional :: val_a2bool(:,:)
      integer(4),   intent(inout), optional :: val_a2i4(:,:)
      integer(8),   intent(inout), optional :: val_a2i8(:,:)
      real(4),      intent(inout), optional :: val_a2r4(:,:)
      real(8),      intent(inout), optional :: val_a2r8(:,:)
      character(*), intent(inout), optional :: val_a3char(:,:,:)
      logical,      intent(inout), optional :: val_a3bool(:,:,:)
      integer(4),   intent(inout), optional :: val_a3i4(:,:,:)
      integer(8),   intent(inout), optional :: val_a3i8(:,:,:)
      real(4),      intent(inout), optional :: val_a3r4(:,:,:)
      real(8),      intent(inout), optional :: val_a3r8(:,:,:)
 
      !--- Local
      integer(impi) :: root, send_size, rank, ierr, i4val
      integer(impi) :: remote_leader, remote_comm_world, tag
      integer(impi) :: status(MPI_status_size)
      integer(impi), pointer :: i4val_a1(:)
      integer(impi), pointer :: i4val_a2(:,:)
      integer(impi), pointer :: i4val_a3(:,:,:)
      integer :: idx, vid
      logical :: in_group
      logical :: use_var(24)

      !--- Identify the incomming variable type
      use_var(:) = .false.
      if ( present(val_schar)  ) use_var( 1) = .true.
      if ( present(val_sbool)  ) use_var( 2) = .true.
      if ( present(val_si4)    ) use_var( 3) = .true.
      if ( present(val_si8)    ) use_var( 4) = .true.
      if ( present(val_sr4)    ) use_var( 5) = .true.
      if ( present(val_sr8)    ) use_var( 6) = .true.
      if ( present(val_a1char) ) use_var( 7) = .true.
      if ( present(val_a1bool) ) use_var( 8) = .true.
      if ( present(val_a1i4)   ) use_var( 9) = .true.
      if ( present(val_a1i8)   ) use_var(10) = .true.
      if ( present(val_a1r4)   ) use_var(11) = .true.
      if ( present(val_a1r8)   ) use_var(12) = .true.
      if ( present(val_a2char) ) use_var(13) = .true.
      if ( present(val_a2bool) ) use_var(14) = .true.
      if ( present(val_a2i4)   ) use_var(15) = .true.
      if ( present(val_a2i8)   ) use_var(16) = .true.
      if ( present(val_a2r4)   ) use_var(17) = .true.
      if ( present(val_a2r8)   ) use_var(18) = .true.
      if ( present(val_a3char) ) use_var(19) = .true.
      if ( present(val_a3bool) ) use_var(20) = .true.
      if ( present(val_a3i4)   ) use_var(21) = .true.
      if ( present(val_a3i8)   ) use_var(22) = .true.
      if ( present(val_a3r4)   ) use_var(23) = .true.
      if ( present(val_a3r8)   ) use_var(24) = .true.

      !--- Ensure that there is exactly one incomming variable
      idx = count(use_var)
      if ( idx == 0 ) then
        write(6,*)"bcast_inter_Any: User supplied data to be broadcast is missing."
        call err_exit("bcast_inter_Any",-1)
      endif
      if ( idx > 1 ) then
        write(6,*)"bcast_inter_Any: Only one variable type is allowed."
        call err_exit("bcast_inter_Any",-2)
      endif

      !--- Determine the variable ID (aka the index in use_var)
      vid = 0
      do idx=1,size(use_var)
        if ( use_var(idx) ) then
          vid = idx
          !--- There is only 1 variable allowed so we exit the loop
          exit
        endif
      enddo
      if (vid == 0) then
        write(6,*)"bcast_inter_Any: Unable to determine variable ID."
        call err_exit("bcast_inter_Any",-3)
      endif

      !--- Assign the tag associated with the point to point transfer
      !--- With 18 different variable type this will use tags 1001-1096 inclusive
      !--- No other procedures should use tags in this range
!--TODO-- Need a better way to assign these with a unique tag value
      select case( trim(adjustl(lowerc(remote_group))) )
        case ("cpl")
          tag = 1000 + vid
        case ("atm")
          tag = 1000 + size(use_var) + vid
        case ("ocn")
          tag = 1000 + 2*size(use_var) + vid
        case ("ice")
          tag = 1000 + 3*size(use_var) + vid
        case default
          write(6,*)"bcast_inter_Any: Invalid group name ",trim(remote_group)
          call err_exit("bcast_inter_Any",-4)
      end select

      !---Determine the rank of the calling process in MPI_COMM_WORLD
      call mpi_comm_rank ( MPI_COMM_WORLD, rank, ierr )

      !--- Is this process in the destination group
      in_group = trim(adjustl(lowerc(task_info(rank)%group_name))) .eq. &
                 trim(adjustl(lowerc(remote_group)))

      if ( rank == src_rank ) then
        !--- Send the data from src_rank
        if ( in_group ) then
          write(6,*)"bcast_inter_Any: src_rank=",src_rank, &
                    " is in the remote group ",trim(remote_group)
          call err_exit("bcast_inter_Any",-5)
        endif

        !--- Send the data to the master task in the destination group
        remote_leader = group_leader(remote_group)

        !--- vid is the index of the variable type to be sent
        select case(vid)
          case (1)  !--- Scalar character
            send_size = len(val_schar)
            call mpi_ssend(val_schar, send_size, MPI_CHARACTER, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (2)  !--- Scalar logical
            send_size = 1
            i4val = 0
            if (val_sbool) i4val = 1
            call mpi_ssend(i4val, send_size, MPI_INTEGER4, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (3)  !--- Scalar integer kind=4
            send_size = 1
            call mpi_ssend(val_si4, send_size, MPI_INTEGER4, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (4)  !--- Scalar integer kind=8
            send_size = 1
            call mpi_ssend(val_si8, send_size, MPI_INTEGER8, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (5)  !--- Scalar real kind=4
            send_size = 1
            call mpi_ssend(val_sr4, send_size, MPI_REAL4, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (6)  !--- Scalar real kind=8
            send_size = 1
            call mpi_ssend(val_sr8, send_size, MPI_REAL8, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (7)  !--- 1D character array
            send_size = len(val_a1char(1)) * size(val_a1char)
            call mpi_ssend(val_a1char, send_size, MPI_CHARACTER, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (8)  !--- 1D logical array
            send_size = size(val_a1bool)
            if ( associated(i4val_a1) ) deallocate(i4val_a1)
            allocate( i4val_a1(send_size) )
            i4val_a1(:) = 0
            where ( val_a1bool ) i4val_a1 = 1
            call mpi_ssend(i4val_a1, send_size, MPI_INTEGER4, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
            deallocate(i4val_a1)
          case (9)  !--- 1D integer kind=4 array
            send_size = size(val_a1i4)
            call mpi_ssend(val_a1i4, send_size, MPI_INTEGER4, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (10)  !--- 1D integer kind=8 array
            send_size = size(val_a1i8)
            call mpi_ssend(val_a1i8, send_size, MPI_INTEGER8, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (11)  !--- 1D real kind=4 array
            send_size = size(val_a1r4)
            call mpi_ssend(val_a1r4, send_size, MPI_REAL4, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (12)  !--- 1D real kind=8 array
            send_size = size(val_a1r8)
            call mpi_ssend(val_a1r8, send_size, MPI_REAL8, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (13)  !--- 2D character array
            send_size = len(val_a2char(1,1)) * size(val_a2char)
            call mpi_ssend(val_a2char, send_size, MPI_CHARACTER, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (14)  !--- 2D logical array
            send_size = size(val_a2bool)
            if ( associated(i4val_a2) ) deallocate(i4val_a2)
            allocate( i4val_a2(size(val_a2bool,1),size(val_a2bool,2)) )
            i4val_a2(:,:) = 0
            where ( val_a2bool ) i4val_a2 = 1
            call mpi_ssend(i4val_a2, send_size, MPI_INTEGER4, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
            deallocate(i4val_a2)
          case (15)  !--- 2D integer kind=4 array
            send_size = size(val_a2i4)
            call mpi_ssend(val_a2i4, send_size, MPI_INTEGER4, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (16)  !--- 2D integer kind=8 array
            send_size = size(val_a2i8)
            call mpi_ssend(val_a2i8, send_size, MPI_INTEGER8, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (17)  !--- 2D real kind=4 array
            send_size = size(val_a2r4)
            call mpi_ssend(val_a2r4, send_size, MPI_REAL4, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (18)  !--- 2D real kind=8 array
            send_size = size(val_a2r8)
            call mpi_ssend(val_a2r8, send_size, MPI_REAL8, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (19)  !--- 3D character array
            send_size = len(val_a3char(1,1,1)) * size(val_a3char)
            call mpi_ssend(val_a3char, send_size, MPI_CHARACTER, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (20)  !--- 3D logical array
            send_size = size(val_a3bool)
            if ( associated(i4val_a3) ) deallocate(i4val_a3)
            allocate( i4val_a3(size(val_a3bool,1),size(val_a3bool,2),size(val_a3bool,3)) )
            i4val_a3(:,:,:) = 0
            where ( val_a3bool ) i4val_a3 = 1
            call mpi_ssend(i4val_a3, send_size, MPI_INTEGER4, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
            deallocate(i4val_a3)
          case (21)  !--- 3D integer kind=4 array
            send_size = size(val_a3i4)
            call mpi_ssend(val_a3i4, send_size, MPI_INTEGER4, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (22)  !--- 3D integer kind=8 array
            send_size = size(val_a3i8)
            call mpi_ssend(val_a3i8, send_size, MPI_INTEGER8, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (23)  !--- 3D real kind=4 array
            send_size = size(val_a3r4)
            call mpi_ssend(val_a3r4, send_size, MPI_REAL4, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case (24)  !--- 3D real kind=8 array
            send_size = size(val_a3r8)
            call mpi_ssend(val_a3r8, send_size, MPI_REAL8, &
                           remote_leader, tag, MPI_COMM_WORLD, ierr )
          case default
            write(6,*)": Invalid variable ID ",vid
            call err_exit("bcast_inter_Any",-6)
        end select
      endif

      if ( in_group ) then
        !--- Do nothing more unless this task is in the destination group

        remote_leader = task_info(rank)%master
        remote_comm_world = model_group_comm

        !--- The master task in the destination group will receive the data sent by src_rank
        if ( rank == remote_leader ) then
          !--- vid is the index of the variable to be sent
          select case(vid)
            case (1)  !--- Scalar character
              send_size = len(val_schar)
              call mpi_recv(val_schar, send_size, MPI_CHARACTER, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (2)  !--- Scalar logical
              send_size = 1
              call mpi_recv(i4val, send_size, MPI_INTEGER4, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
              val_sbool = .false.
              if (i4val == 1) val_sbool = .true.
            case (3)  !--- Scalar integer kind=4
              send_size = 1
              call mpi_recv(val_si4, send_size, MPI_INTEGER4, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (4)  !--- Scalar integer kind=8
              send_size = 1
              call mpi_recv(val_si8, send_size, MPI_INTEGER8, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (5)  !--- Scalar real kind=4
              send_size = 1
              call mpi_recv(val_sr4, send_size, MPI_REAL4, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (6)  !--- Scalar real kind=8
              send_size = 1
              call mpi_recv(val_sr8, send_size, MPI_REAL8, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (7)  !--- 1D character array
              send_size = len(val_a1char(1)) * size(val_a1char)
              call mpi_recv(val_a1char, send_size, MPI_CHARACTER, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (8)  !--- 1D logical array
              send_size = size(val_a1bool)
              if ( associated(i4val_a1) ) deallocate(i4val_a1)
              allocate( i4val_a1(send_size) )
              call mpi_recv(i4val_a1, send_size, MPI_INTEGER4, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
              val_a1bool(:) = .false.
              where ( i4val_a1 == 1 ) val_a1bool = .true.
              deallocate(i4val_a1)
            case (9)  !--- 1D integer kind=4 array
              send_size = size(val_a1i4)
              call mpi_recv(val_a1i4, send_size, MPI_INTEGER4, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (10)  !--- 1D integer kind=8 array
              send_size = size(val_a1i8)
              call mpi_recv(val_a1i8, send_size, MPI_INTEGER8, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (11)  !--- 1D real kind=4 array
              send_size = size(val_a1r4)
              call mpi_recv(val_a1r4, send_size, MPI_REAL4, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (12)  !--- 1D real kind=8 array
              send_size = size(val_a1r8)
              call mpi_recv(val_a1r8, send_size, MPI_REAL8, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (13)  !--- 2D character array
              send_size = len(val_a2char(1,1)) * size(val_a2char)
              call mpi_recv(val_a2char, send_size, MPI_CHARACTER, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (14)  !--- 2D logical array
              send_size = size(val_a2bool)
              if ( associated(i4val_a2) ) deallocate(i4val_a2)
              allocate( i4val_a2(size(val_a2bool,1),size(val_a2bool,2)) )
              call mpi_recv(i4val_a2, send_size, MPI_INTEGER4, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
              val_a2bool(:,:) = .false.
              where ( i4val_a2 == 1 )  val_a2bool = .true.
              deallocate(i4val_a2)
            case (15)  !--- 2D integer kind=4 array
              send_size = size(val_a2i4)
              call mpi_recv(val_a2i4, send_size, MPI_INTEGER4, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (16)  !--- 2D integer kind=8 array
              send_size = size(val_a2i8)
              call mpi_recv(val_a2i8, send_size, MPI_INTEGER8, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (17)  !--- 2D real kind=4 array
              send_size = size(val_a2r4)
              call mpi_recv(val_a2r4, send_size, MPI_REAL4, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (18)  !--- 2D real kind=8 array
              send_size = size(val_a2r8)
              call mpi_recv(val_a2r8, send_size, MPI_REAL8, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (19)  !--- 3D character array
              send_size = len(val_a3char(1,1,1)) * size(val_a3char)
              call mpi_recv(val_a3char, send_size, MPI_CHARACTER, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (20)  !--- 3D logical array
              send_size = size(val_a3bool)
              if ( associated(i4val_a3) ) deallocate(i4val_a3)
              allocate( i4val_a3(size(val_a3bool,1),size(val_a3bool,2),size(val_a3bool,3)) )
              call mpi_recv(i4val_a3, send_size, MPI_INTEGER4, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
              val_a3bool(:,:,:) = .false.
              where ( i4val_a3 == 1 )  val_a3bool = .true.
              deallocate(i4val_a3)
            case (21)  !--- 3D integer kind=4 array
              send_size = size(val_a3i4)
              call mpi_recv(val_a3i4, send_size, MPI_INTEGER4, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (22)  !--- 3D integer kind=8 array
              send_size = size(val_a3i8)
              call mpi_recv(val_a3i8, send_size, MPI_INTEGER8, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (23)  !--- 3D real kind=4 array
              send_size = size(val_a3r4)
              call mpi_recv(val_a3r4, send_size, MPI_REAL4, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case (24)  !--- 3D real kind=8 array
              send_size = size(val_a3r8)
              call mpi_recv(val_a3r8, send_size, MPI_REAL8, src_rank, &
                            tag, MPI_COMM_WORLD, status, ierr )
            case default
              write(6,*)": Invalid variable ID ",vid
              call err_exit("bcast_inter_Any",-7)
          end select
        endif
 
        !--- Broadcast the data within the remote group
        root = 0   !--- root is the master task in the destination group comm world

        !--- vid is the index of the variable to be sent
        select case(vid)
          case (1)  !--- Scalar character
            send_size = len(val_schar)
            call mpi_bcast(val_schar, send_size, MPI_CHARACTER, root, &
                           remote_comm_world, ierr)
          case (2)  !--- Scalar logical
            send_size = 1
            i4val = 0
            if (val_sbool) i4val = 1
            call mpi_bcast(i4val, send_size, MPI_INTEGER4, root, &
                           remote_comm_world, ierr)
            val_sbool = .false.
            if (i4val == 1) val_sbool = .true.
          case (3)  !--- Scalar integer kind=4
            send_size = 1
            call mpi_bcast(val_si4, send_size, MPI_INTEGER4, root, &
                           remote_comm_world, ierr)
          case (4)  !--- Scalar integer kind=8
            send_size = 1
            call mpi_bcast(val_si8, send_size, MPI_INTEGER8, root, &
                           remote_comm_world, ierr)
          case (5)  !--- Scalar real kind=4
            send_size = 1
            call mpi_bcast(val_sr4, send_size, MPI_REAL4, root, &
                           remote_comm_world, ierr)
          case (6)  !--- Scalar real kind=8
            send_size = 1
            call mpi_bcast(val_sr8, send_size, MPI_REAL8, root, &
                           remote_comm_world, ierr)
          case (7)  !--- 1D character array
            send_size = len(val_a1char(1)) * size(val_a1char)
            call mpi_bcast(val_a1char, send_size, MPI_CHARACTER, root, &
                           remote_comm_world, ierr)
          case (8)  !--- 1D logical array
            send_size = size(val_a1bool)
            if ( associated(i4val_a1) ) deallocate(i4val_a1)
            allocate( i4val_a1(send_size) )
            i4val_a1(:) = 0
            where ( val_a1bool ) i4val_a1 = 1
            call mpi_bcast(i4val_a1, send_size, MPI_INTEGER4, root, &
                           remote_comm_world, ierr)
            val_a1bool(:) = .false.
            where ( i4val_a1 == 1 )  val_a1bool = .true.
            deallocate(i4val_a1)
          case (9)  !--- 1D integer kind=4 array
            send_size = size(val_a1i4)
            call mpi_bcast(val_a1i4, send_size, MPI_INTEGER4, root, &
                           remote_comm_world, ierr)
          case (10)  !--- 1D integer kind=8 array
            send_size = size(val_a1i8)
            call mpi_bcast(val_a1i8, send_size, MPI_INTEGER8, root, &
                           remote_comm_world, ierr)
          case (11)  !--- 1D real kind=4 array
            send_size = size(val_a1r4)
            call mpi_bcast(val_a1r4, send_size, MPI_REAL4, root, &
                           remote_comm_world, ierr)
          case (12)  !--- 1D real kind=8 array
            send_size = size(val_a1r8)
            call mpi_bcast(val_a1r8, send_size, MPI_REAL8, root, &
                           remote_comm_world, ierr)
          case (13)  !--- 2D character array
            send_size = len(val_a2char(1,1)) * size(val_a2char)
            call mpi_bcast(val_a2char, send_size, MPI_CHARACTER, root, &
                           remote_comm_world, ierr)
          case (14)  !--- 2D logical array
            send_size = size(val_a2bool)
            if ( associated(i4val_a2) ) deallocate(i4val_a2)
            allocate( i4val_a2(size(val_a2bool,1),size(val_a2bool,2)) )
            i4val_a2(:,:) = 0
            where ( val_a2bool ) i4val_a2 = 1
            call mpi_bcast(i4val_a2, send_size, MPI_INTEGER4, root, &
                           remote_comm_world, ierr)
            val_a2bool(:,:) = .false.
            where ( i4val_a2 == 1 )  val_a2bool = .true.
            deallocate(i4val_a2)
          case (15)  !--- 2D integer kind=4 array
            send_size = size(val_a2i4)
            call mpi_bcast(val_a2i4, send_size, MPI_INTEGER4, root, &
                           remote_comm_world, ierr)
          case (16)  !--- 2D integer kind=8 array
            send_size = size(val_a2i8)
            call mpi_bcast(val_a2i8, send_size, MPI_INTEGER8, root, &
                           remote_comm_world, ierr)
          case (17)  !--- 2D real kind=4 array
            send_size = size(val_a2r4)
            call mpi_bcast(val_a2r4, send_size, MPI_REAL4, root, &
                           remote_comm_world, ierr)
          case (18)  !--- 2D real kind=8 array
            send_size = size(val_a2r8)
            call mpi_bcast(val_a2r8, send_size, MPI_REAL8, root, &
                           remote_comm_world, ierr)
          case (19)  !--- 3D character array
            send_size = len(val_a3char(1,1,1)) * size(val_a3char)
            call mpi_bcast(val_a3char, send_size, MPI_CHARACTER, root, &
                           remote_comm_world, ierr)
          case (20)  !--- 3D logical array
            send_size = size(val_a3bool)
            if ( associated(i4val_a3) ) deallocate(i4val_a3)
            allocate( i4val_a3(size(val_a3bool,1),size(val_a3bool,2),size(val_a3bool,3)) )
            i4val_a3(:,:,:) = 0
            where ( val_a3bool ) i4val_a3 = 1
            call mpi_bcast(i4val_a3, send_size, MPI_INTEGER4, root, &
                           remote_comm_world, ierr)
            val_a3bool(:,:,:) = .false.
            where ( i4val_a3 == 1 )  val_a3bool = .true.
            deallocate(i4val_a3)
          case (21)  !--- 3D integer kind=4 array
            send_size = size(val_a3i4)
            call mpi_bcast(val_a3i4, send_size, MPI_INTEGER4, root, &
                           remote_comm_world, ierr)
          case (22)  !--- 3D integer kind=8 array
            send_size = size(val_a3i8)
            call mpi_bcast(val_a3i8, send_size, MPI_INTEGER8, root, &
                           remote_comm_world, ierr)
          case (23)  !--- 3D real kind=4 array
            send_size = size(val_a3r4)
            call mpi_bcast(val_a3r4, send_size, MPI_REAL4, root, &
                           remote_comm_world, ierr)
          case (24)  !--- 3D real kind=8 array
            send_size = size(val_a3r8)
            call mpi_bcast(val_a3r8, send_size, MPI_REAL8, root, &
                           remote_comm_world, ierr)
          case default
            write(6,*)": Invalid variable ID ",vid
            call err_exit("bcast_inter_Any",-8)
        end select

        !--- Everybody waits here until all tasks have caught up
        call mpi_barrier(remote_comm_world, ierr)
      endif

    end subroutine bcast_inter_Any

!==============================================================================
! End definition of routines associated with the generic inteface bcast_inter
!==============================================================================

#endif

end module com_cpl
