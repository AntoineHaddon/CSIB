MODULE cpl_cancpl
  !!=======================================================================
  !!                    ***  MODULE cpl_cancpl  ***
  !!
  !! Coupled ocean-atmosphere using the coupler CanCPL
  !!
  !! This module is intended to serve the same function as cpl_oasis3 does
  !! in nemo v3.4.1, except (of course) it will communicate with CanCPL
  !!
  !! Larry Solheim Aug,2015
  !!=======================================================================
#if defined key_cancpl
  !!-----------------------------------------------------------------------
  !!   'key_cancpl'                    coupled Ocean/Atmosphere via CanCPL
  !!-----------------------------------------------------------------------
  !!   cpl_cancpl_init     : initialize coupled mode communication
  !!   cpl_cancpl_define   : define all fields to be sent of received
  !!   cpl_cancpl_snd      : send fields in coupled mode
  !!   cpl_cancpl_rcv      : receive fields in coupled mode
  !!   cpl_cancpl_freq     : return the coupling frequency of a given field
  !!   cpl_cancpl_finalize : finalize the coupled mode communication
  !!-----------------------------------------------------------------------

  !--- CanCPL variables/routines common to all component models
  use com_cpl

  use par_oce                      ! ocean parameters
  use dom_oce                      ! ocean space and time domain
  use in_out_manager               ! I/O manager
  use lbclnk                       ! ocean lateral boundary conditions (or mpp link)
  use sbc_oce, only: nn_ice
#if defined key_cice
  use ice_domain_size, only: ncat
#endif
  use lib_mpp, only: mpi_comm_opa

  implicit none
  private

  public :: cpl_cancpl_init
  public :: cpl_cancpl_define
  public :: cpl_cancpl_snd
  public :: cpl_cancpl_rcv
  public :: cpl_cancpl_freq
  public :: cpl_cancpl_finalize

  logical, public, parameter ::   lk_cpl = .true.   !: coupled flag
  integer, public, save      ::   oasis_idle = 0    !: return code if no send or recv
  integer, public, save      ::   oasis_rcv  = 1    !: return code if field received
  integer, public, save      ::   oasis_snd  = 2    !: return code if field sent

  !--- These are defined in com_cpl
  public :: cpl_vinfo_t, find_cpl_vinfo
  public :: bcast_inter, cpl_time_string, cpl_elapsed_time_secs, cpl_master

  integer            ::   nerror        ! return error code
  integer, parameter ::   nmaxfld=40    ! maximum number of coupling fields

  !--- Derived type used to contain coupling field information
  type, public :: fld_cpl
     logical               ::   laction   ! To be coupled or not
     character(len = 8)    ::   clname    ! Name of the coupling field   
     character(len = 1)    ::   clgrid    ! Grid type  
     real(wp)              ::   nsgn      ! Control of the sign change
     integer, dimension(9) ::   nid = 0  ! Id of the field (no more than 9 categories)
     integer               ::   nct       ! Number of categories in field
  end type fld_cpl

  !--- Lists of all fields that are to be coupled
  type(fld_cpl), save, dimension(nmaxfld), public ::   srcv, ssnd

  TYPE :: FLD_C
     CHARACTER(len = 32) ::   cldes                  ! desciption of the coupling strategy
     CHARACTER(len = 32) ::   clcat                  ! multiple ice categories strategy
     CHARACTER(len = 32) ::   clvref                 ! reference of vector ('spherical' or 'cartesian')
     CHARACTER(len = 32) ::   clvor                  ! orientation of vector fields ('eastward-northward' or 'local grid')
     CHARACTER(len = 32) ::   clvgrd                 ! grids on which is located the vector fields
  END TYPE FLD_C
  ! Send to the atmosphere                           !
  TYPE(FLD_C) ::   sn_snd_temp, sn_snd_alb, sn_snd_thick, sn_snd_crt, sn_snd_co2                        
  ! Received from the atmosphere                     !
  TYPE(FLD_C) ::   sn_rcv_w10m, sn_rcv_taumod, sn_rcv_tau, sn_rcv_dqnsdt, sn_rcv_qsr, sn_rcv_qns, sn_rcv_emp, sn_rcv_rnf
  TYPE(FLD_C) ::   sn_rcv_cal, sn_rcv_iceflx, sn_rcv_co2

  NAMELIST /namsbc_cpl/ sn_snd_temp, sn_snd_alb   , sn_snd_thick, sn_snd_crt   , sn_snd_co2,   &
                        sn_rcv_w10m, sn_rcv_taumod, sn_rcv_tau  , sn_rcv_dqnsdt, sn_rcv_qsr,   &
                        sn_rcv_qns , sn_rcv_emp   , sn_rcv_rnf  , sn_rcv_cal   , sn_rcv_iceflx  , sn_rcv_co2

  integer :: nn_fsbc, nn_ice, nn_fwb
  logical :: ln_ana, ln_flx, ln_blk_clio, ln_blk_core, ln_cpl, ln_blk_mfs, ln_apr_dyn, ln_dm2dc, ln_rnf, ln_ssr, ln_cdgw
  NAMELIST/namsbc/ nn_fsbc, ln_ana, ln_flx, ln_blk_clio, ln_blk_core, ln_cpl,   &
                   ln_blk_mfs, ln_apr_dyn, nn_ice, ln_dm2dc, ln_rnf, ln_ssr, nn_fwb, ln_cdgw

  !--- tmp space for use with MPI gather/scatter operations
  real(wp), allocatable, save, dimension(:,:,:), private :: png

  !--- tmp char space
  character(512), save :: strng

  !--- Work space used to temporarily hold fields passed via MPI
  integer, parameter    :: maxx=8392704  !---4098x2048 = (2+2^12)x(2^11)
  real(kind=8), save    :: wrk(maxx)
  integer(kind=8), save :: ibuf(8)

  !--- A derived type holding coupler related information
  !--- cpl_vinfo_t is defined in the com_cpl module
  type(cpl_vinfo_t), save :: cpl_vinfo

  !--- NOTE: nproc is not equal to the MPI task in MPI_COMM_WORLD since it will always
  !--- be one of 0,1,2,...(jpnij-1) and the AGCM gets the first set of MPI tasks.
  !--- However nproc == 0 should still correspond with the ocn_master task
  !--- nproc is use associated through the module dom_oce

contains

  subroutine cpl_cancpl_init( kl_comm )
     !!-------------------------------------------------------------------
     !!             ***  ROUTINE cpl_cancpl_init  ***
     !!
     !! ** Purpose :   Initialize coupled mode communication for ocean
     !!    exchange between AGCM, OGCM and COUPLER.
     !!--------------------------------------------------------------------------
     integer, intent(out) ::   kl_comm   ! MPI group ID for the ocean comm_world
     !!--------------------------------------------------------------------------

     integer(kind=impi) :: rank, ierr
     integer(kind=impi) :: local_ocn_comm

     !!============================================
     !! WARNING: No write in numout in this routine
     !!============================================

     !--- Initialize groups for cpl, atm, ocn, ice, ...
     !--- This will, among other things, define cpl_master, atm_master, ocn_master
     !--- and return an ocean intra-communicator as local_ocn_comm
     !--- NOTE: ocn_master is the rank in MPI_COMM_WORLD not the rank in local_ocn_comm
     local_ocn_comm = -1
     call define_group('ocn', local_ocn_comm)

     kl_comm  = local_ocn_comm

     !---Determine the rank of the calling process in MPI_COMM_WORLD
     call mpi_comm_rank ( MPI_COMM_WORLD, rank, ierr )

     !--- ocn_master is defined during the call to define_group
     if ( rank == ocn_master ) then
       !--- Write to stdout (unit 6) since numout is not attached yet
       write(6,*) 'cpl_cancpl_init: ocean comm_world = ',kl_comm
       call flush(6)
     endif

  end subroutine cpl_cancpl_init


  subroutine cpl_cancpl_define( krcv, ksnd )
     !!-------------------------------------------------------------------
     !!             ***  ROUTINE cpl_cancpl_define  ***
     !!
     !! ** Purpose :   Define grid and field information for ocean
     !!    exchange between AGCM, OGCM and COUPLER.
     !!--------------------------------------------------------------------
     integer, intent(in) :: krcv   ! Number of all possible fields received
     integer, intent(in) :: ksnd   ! Number of all possible fields sent

     !--- Local
     integer :: ji,jc,jx
     character(len=8) :: zclname
     integer :: ldbg=1
     integer(kind=impi) :: rank, ierr
     integer :: verbose=1
     integer :: min_rank, min_index

     !--- var_list_info will be assigned enough info about each list of variables
     !--- sent or received to allow a simple reordering of these lists
     type var_list_info_t
       character(32) :: name
       integer       :: index
       integer       :: rank
       logical       :: used
     end type var_list_info_t
     type(var_list_info_t) :: var_list_info(50)

     !--- The order that fields are sent from NEMO to the coupler is determined
     !--- in the subroutine sbc_cpl_snd. In terms of the index in ssnd, this order is
     !--- jps_toce=2, jps_tice=3, jps_tmix=4, jps_albice=5, jps_albmix=6, jps_fice=1,
     !--- jps_hice=7, jps_hsnw=8, jps_co2=15, jps_ocx1=9, jps_ocy1=10, jps_ocz1=11,
     !--- jps_ivx1=12, jps_ivy1=13, jps_ivz1=14
     integer, dimension(15) :: send_order = &
         (/ 2, 3, 4, 5, 6, 1, 7, 8, 15, 9, 10, 11, 12, 13, 14 /)
     !!--------------------------------------------------------------------

     !--- Determine the rank of the calling process in MPI_COMM_WORLD
     call mpi_comm_rank ( MPI_COMM_WORLD, rank, ierr )

     if ( rank == ocn_master ) then
       write(numout,*)
       write(numout,*) 'cpl_cancpl_define : initialize coupled ocean/atmosphere'
       write(numout,*) '~~~~~~~~~~~~~~~~~'
       write(numout,*)
       call flush(numout)
     endif

     !--- Allocate temporary space used with MPI gather/scatter ops below
     allocate( png(jpi,jpj,jpnij), stat=nerror )
     if( nerror > 0 ) then
       call ctl_stop("STOP", " cpl_cancpl_define", "Problem allocating png")
     endif

     ! -----------------------------------------------------------------
     ! ... Assign variable name lists for ssnd and srcv variables
     ! -----------------------------------------------------------------
     !
     !--- Send variables
     !
     !--- nemo_n_send_var is the number of fields in nemo_send_var
     !--- These variables are available from the com_cpl module
     nemo_n_send_var=0
     var_list_info(:)%name  = " "
     var_list_info(:)%index = 0
     var_list_info(:)%rank  = 0
     var_list_info(:)%used  = .false.
     do ji = 1, ksnd
        if ( ssnd(ji)%laction ) then 
           do jc = 1, ssnd(ji)%nct
              if ( ssnd(ji)%nct .gt. 1 ) then
                 write(zclname,'( a7, i1)') ssnd(ji)%clname,jc
              else
                 zclname=ssnd(ji)%clname
              endif

              nemo_n_send_var = nemo_n_send_var + 1
              if ( nemo_n_send_var > 40 ) then
                write(numout,*) "cpl_cancpl_define: Too many send variables"
                write(6,*) "cpl_cancpl_define: Too many send variables"
                call flush(6)
                call ctl_stop("STOP", " cpl_cancpl_define", "Too many send variables")
              endif

              var_list_info(nemo_n_send_var)%name  = trim(zclname)
              var_list_info(nemo_n_send_var)%index = ji
              do jx=1,size(send_order)
                if ( var_list_info(nemo_n_send_var)%index == send_order(jx) ) then
                  var_list_info(nemo_n_send_var)%rank = jx
                  exit
                endif
              enddo
              if ( var_list_info(nemo_n_send_var)%rank == 0 ) then
                write(6,*)"cpl_cancpl_define: Unable to determine rank for ", &
                          trim(var_list_info(nemo_n_send_var)%name)
                call flush(6)
                call ctl_stop("STOP", " cpl_cancpl_define", "Unable to determine send rank")
              endif

           end do
        endif
     end do
     if ( nemo_n_send_var > 0 ) then
       if ( associated(nemo_send_var) ) deallocate(nemo_send_var)
       allocate( nemo_send_var(nemo_n_send_var) )
       nemo_send_var(1:nemo_n_send_var) = var_list_info(1:nemo_n_send_var)%name
     endif

     if ( nemo_n_send_var > 1 ) then
       !--- nemo_send_var needs to be reordered
       !--- The fields in this list must be in the same order as the the data
       !--- that is sent to the coupler.
       !--- The subroutine sbc_cpl_snd determines the send order
       !--- There is no clear way to determine this order on the fly so
       !--- it must be hard coded here (this is bad).
       !--- If there are any changes in sbc_cpl_snd that alter this order
       !--- then there must also be changes here.
       do ji=1,nemo_n_send_var
         min_rank  = 1000
         min_index = -1
         !--- Find the index of the variable with the lowest rank
         do jx=1,nemo_n_send_var
           !--- Ignore names already in the ordered list
           if ( var_list_info(jx)%used ) cycle
           if ( var_list_info(jx)%rank < min_rank ) then
             min_rank = var_list_info(jx)%rank
             min_index = jx
           endif
         enddo
         if ( min_index < 1 ) then
           write(6,*)"cpl_cancpl_define: Unable to find min rank at send list element ",ji
           call flush(6)
           call ctl_stop("STOP", " cpl_cancpl_define", "Unable to find send list min rank")
         endif
         var_list_info(min_index)%used = .true.
         !--- Overwrite nemo_send_var with the properly ordered names
         nemo_send_var(ji) = var_list_info(min_index)%name
       enddo
     endif

     !--- Receive variables
     !
     !--- nemo_n_recv_var is the number of fields in nemo_recv_var
     !--- These variables are available from the com_cpl module
     nemo_n_recv_var=0
     var_list_info(:)%name  = " "
     var_list_info(:)%index = 0
     var_list_info(:)%rank  = 0
     var_list_info(:)%used  = .false.
     do ji = 1, krcv
        if ( srcv(ji)%laction ) then 
           do jc = 1, srcv(ji)%nct
              if ( srcv(ji)%nct .gt. 1 ) then
                 write(zclname,'( a7, i1)') srcv(ji)%clname,jc
              else
                 zclname=srcv(ji)%clname
              endif

              nemo_n_recv_var = nemo_n_recv_var + 1
              if ( nemo_n_recv_var > 40 ) then
                write(numout,*) "cpl_cancpl_define: Too many receive variables"
                write(6,*) "cpl_cancpl_define: Too many receive variables"
                call flush(6)
                call ctl_stop("STOP", " cpl_cancpl_define", "Too many receive variables")
              endif

              var_list_info(nemo_n_recv_var)%name  = trim(zclname)
              var_list_info(nemo_n_recv_var)%index = ji
              !--- The rank in srcv also indicates the order data is received
              var_list_info(nemo_n_recv_var)%rank = ji

           end do
        endif
     end do
     if ( nemo_n_recv_var > 0 ) then
       if ( associated(nemo_recv_var) ) deallocate(nemo_recv_var)
       allocate( nemo_recv_var(nemo_n_recv_var) )
       !--- The fields in this list must be in the same order as the the data
       !--- that is received from the coupler
       !--- nemo_recv_var should already be in the correct order which is
       !--- the order the fields appear in srcv
       !--- If this ever changes then this list will need to be reordered
       nemo_recv_var(1:nemo_n_recv_var) = var_list_info(1:nemo_n_recv_var)%name
     endif

     if ( rank == ocn_master .and. verbose > 1 ) then
       write(6,*)"cpl_cancpl_define: nemo_n_send_var=",nemo_n_send_var
       write(6,'(5(2x,a))')nemo_send_var(1:nemo_n_send_var)
       write(6,*)"cpl_cancpl_define: nemo_n_recv_var=",nemo_n_recv_var
       write(6,'(5(2x,a))')nemo_recv_var(1:nemo_n_recv_var)
       call flush(6)
     endif

     !--- The following values will be broadcast to all mpi tasks
     !--- in the call to cpl_initialize_events

     !--- Set a value for nemo_nn_ice, defined in com_cpl
     nemo_nn_ice = nn_ice

     !--- Set a value for nemo_rn_rdt, defined in com_cpl
     !--- rn_rdt is defined in the module dom_oce
     nemo_rn_rdt = nint(rn_rdt,8)

     REWIND( numnam )                    ! ... read namlist namsbc
     READ  ( numnam, namsbc )

     !--- Set a value for nemo_nn_fsbc, defined in com_cpl
     !--- nn_fsbc is defined in the namelist namsbc
     nemo_nn_fsbc = nn_fsbc

     !--- Set a value for nemo_nn_it000, defined in com_cpl
     !--- nn_it000 is defined in the module in_out_manager
     nemo_nn_it000 = nn_it000

     !--- Set a value for nemo_nn_it000, defined in com_cpl
     !--- nn_itend is defined in the module in_out_manager
     nemo_nn_itend = nn_itend

     !--- Set a value for nemo_nn_it000, defined in com_cpl
     !--- nn_date0 is defined in the module in_out_manager
     nemo_nn_date0 = nn_date0

     !--- Set a value for nemo_ncat, defined in com_cpl
     !--- ncat is defined in the CICE module ice_domain_size
#if defined key_cice
     nemo_ncat = ncat
#else
     nemo_ncat = 1
#endif

     !--- Set values for nemo_jpiglo and nemo_jpjglo, defined in com_cpl
     !--- jpiglo and jpjglo are defined in the module par_oce
     nemo_jpiglo = jpiglo
     nemo_jpjglo = jpjglo

     !--- Assign nemo_namsbc_cpl_cldes with namelist parameters read into namsbc_cpl
     !--- These values will be used by the coupler
     !--- Set defaults
     sn_snd_temp   = FLD_C( 'weighted oce and ice',    'no'    ,     ''      ,         ''           ,   ''   ) 
     sn_snd_alb    = FLD_C( 'weighted ice'        ,    'no'    ,     ''      ,         ''           ,   ''   ) 
     sn_snd_thick  = FLD_C( 'none'                ,    'no'    ,     ''      ,         ''           ,   ''   ) 
     sn_snd_crt    = FLD_C( 'none'                ,    'no'    , 'spherical' , 'eastward-northward' ,  'T'   )     
     sn_snd_co2    = FLD_C( 'none'                ,    'no'    ,     ''      ,         ''           ,   ''   )     
     sn_rcv_w10m   = FLD_C( 'none'                ,    'no'    ,     ''      ,         ''          ,   ''    )
     sn_rcv_taumod = FLD_C( 'coupled'             ,    'no'    ,     ''      ,         ''          ,   ''    )
     sn_rcv_tau    = FLD_C( 'oce only'            ,    'no'    , 'cartesian' , 'eastward-northward',  'U,V'  )  
     sn_rcv_dqnsdt = FLD_C( 'coupled'             ,    'no'    ,     ''      ,         ''          ,   ''    )
     sn_rcv_qsr    = FLD_C( 'oce and ice'         ,    'no'    ,     ''      ,         ''          ,   ''    )
     sn_rcv_qns    = FLD_C( 'oce and ice'         ,    'no'    ,     ''      ,         ''          ,   ''    )
     sn_rcv_emp    = FLD_C( 'conservative'        ,    'no'    ,     ''      ,         ''          ,   ''    )
     sn_rcv_rnf    = FLD_C( 'coupled'             ,    'no'    ,     ''      ,         ''          ,   ''    )
     sn_rcv_cal    = FLD_C( 'coupled'             ,    'no'    ,     ''      ,         ''          ,   ''    )
     sn_rcv_iceflx = FLD_C( 'none'                ,    'no'    ,     ''      ,         ''          ,   ''    )
     sn_rcv_co2    = FLD_C( 'none'                ,    'no'    ,     ''      ,         ''          ,   ''    )

     REWIND( numnam )                    ! ... read namlist namsbc_cpl
     READ  ( numnam, namsbc_cpl )

     nemo_namsbc_cpl_cldes(:) = " "
     nemo_namsbc_cpl_cldes(1) = trim(sn_snd_temp%cldes)
     nemo_namsbc_cpl_cldes(2) = trim(sn_snd_alb%cldes)
     nemo_namsbc_cpl_cldes(3) = trim(sn_snd_thick%cldes)
     nemo_namsbc_cpl_cldes(4) = trim(sn_snd_crt%cldes)
     nemo_namsbc_cpl_cldes(5) = trim(sn_snd_co2%cldes)
     nemo_namsbc_cpl_cldes(6) = trim(sn_rcv_w10m%cldes)
     nemo_namsbc_cpl_cldes(7) = trim(sn_rcv_taumod%cldes)
     nemo_namsbc_cpl_cldes(8) = trim(sn_rcv_tau%cldes)
     nemo_namsbc_cpl_cldes(9) = trim(sn_rcv_dqnsdt%cldes)
     nemo_namsbc_cpl_cldes(10) = trim(sn_rcv_qsr%cldes)
     nemo_namsbc_cpl_cldes(11) = trim(sn_rcv_qns%cldes)
     nemo_namsbc_cpl_cldes(12) = trim(sn_rcv_emp%cldes)
     nemo_namsbc_cpl_cldes(13) = trim(sn_rcv_rnf%cldes)
     nemo_namsbc_cpl_cldes(14) = trim(sn_rcv_cal%cldes)
     nemo_namsbc_cpl_cldes(15) = trim(sn_rcv_iceflx%cldes)
     nemo_namsbc_cpl_cldes(16) = trim(sn_rcv_co2%cldes)

     !--- Gather tmask at the surface into the temporary global array png
     !--- tmask is found in module dom_oce
     call mppsync
     call mppgather (tmask(:,:,1),0,png)
     call mppsync

     if ( rank == ocn_master ) then
       !--- Allocate space for nemo_tmask, which is defined in com_cpl
       if ( associated(nemo_tmask) ) deallocate(nemo_tmask)
       allocate( nemo_tmask(nemo_jpiglo,nemo_jpjglo) )

       !--- Assign the global array containing the surface tmask to nemo_tmask
       !--- This data is then sent to the coupler in cpl_initialize_events
       !--- tmask is assigned in module dommsk using data from the array mbathy
       call copy_3d_to_2d_global(nemo_tmask, png)
     endif

     !--- Gather umask at the surface into the temporary global array png
     !--- umask is found in module dom_oce
     call mppsync
     call mppgather (umask(:,:,1),0,png)
     call mppsync

     if ( rank == ocn_master ) then
       !--- Allocate space for nemo_umask, which is defined in com_cpl
       if ( associated(nemo_umask) ) deallocate(nemo_umask)
       allocate( nemo_umask(nemo_jpiglo,nemo_jpjglo) )

       !--- Assign the global array containing the surface umask to nemo_umask
       !--- This data is then sent to the coupler in cpl_initialize_events
       !--- umask is defined in module dommsk in terms of tmask as follows
       !---    umask(i,j) = tmask(i,j) * tmask(i+1,j)
       call copy_3d_to_2d_global(nemo_umask, png)
     endif

     !--- Gather vmask at the surface into the temporary global array png
     !--- vmask is found in module dom_oce
     call mppsync
     call mppgather (vmask(:,:,1),0,png)
     call mppsync

     if ( rank == ocn_master ) then
       !--- Allocate space for nemo_vmask, which is defined in com_cpl
       if ( associated(nemo_vmask) ) deallocate(nemo_vmask)
       allocate( nemo_vmask(nemo_jpiglo,nemo_jpjglo) )

       !--- Assign the global array containing the surface vmask to nemo_vmask
       !--- This data is then sent to the coupler in cpl_initialize_events
       !--- vmask is defined in module dommsk in terms of tmask as follows
       !---    vmask(i,j) = tmask(i,j) * tmask(i,j+1)
       call copy_3d_to_2d_global(nemo_vmask, png)
     endif

     !--- Gather fmask at the surface into the temporary global array png
     !--- fmask is found in module dom_oce
     call mppsync
     call mppgather (fmask(:,:,1),0,png)
     call mppsync

     if ( rank == ocn_master ) then
       !--- Allocate space for nemo_fmask, which is defined in com_cpl
       if ( associated(nemo_fmask) ) deallocate(nemo_fmask)
       allocate( nemo_fmask(nemo_jpiglo,nemo_jpjglo) )

       !--- Assign the global array containing the surface fmask to nemo_fmask
       !--- This data is then sent to the coupler in cpl_initialize_events
       !--- fmask is defined in module dommsk in terms of tmask as follows
       !---    fmask(i,j) = tmask(i,j) * tmask(i+1,j) * tmask(i,j+1) * tmask(i+1,j+1)
       !--- then further modified for lateral boundary conditions on velocity
       !--- and to increase lateral friction near certain straights
       !--- Note: this means fmask will have values other than 0/1 (e.g. 0.5,2,3)
       call copy_3d_to_2d_global(nemo_fmask, png)
     endif

     !--- Initialize coupler events and broadcast global variables
     call cpl_initialize_events()

     !--- Broadcast the initial date and time from the coupler to all tasks
     !--- cpl_time_string is defined in com_cpl
     call bcastGroup(cpl_time_string, cpl_master, MPI_COMM_WORLD)

     ! -----------------------------------------------------------------
     ! ... Assign MPI tags to ssnd and srcv variables
     !--- This must be done after the call to cpl_initialize_events
     !--- because it will define these tags
     ! -----------------------------------------------------------------
     if ( nemo_n_send_var > 0 ) then
       do ji=1,nemo_n_send_var
         cpl_vinfo = find_cpl_vinfo( name=trim(nemo_send_var(ji)) )
         do jx=1,40
           if ( trim(adjustl(nemo_send_var(ji))) .eq. trim(adjustl(ssnd(jx)%clname)) ) then
             !--- jx is the index in ssnd for this name
             do jc=1,ssnd(jx)%nct
               ssnd(jx)%nid(jc) = cpl_vinfo%tag
             enddo
           endif
         enddo
         if ( rank == ocn_master ) then
           !--- Write to NEMO's ocean.output file
           write(numout,*) "cpl_cancpl_define: Send field ",ji, &
               "  name=",trim(nemo_send_var(ji))," tag=",cpl_vinfo%tag
           call flush(numout)

           !--- Also write to stdout (unit 6)
           write(6,*) "cpl_cancpl_define: Send field ",ji, &
               "  name=",trim(nemo_send_var(ji))," tag=",cpl_vinfo%tag
           call flush(6)
         endif
       enddo
     endif

     if ( nemo_n_recv_var > 0 ) then
       do ji=1,nemo_n_recv_var
         cpl_vinfo = find_cpl_vinfo( name=trim(nemo_recv_var(ji)) )
         do jx=1,40
           if ( trim(adjustl(nemo_recv_var(ji))) .eq. trim(adjustl(srcv(jx)%clname)) ) then
             !--- jx is the index in srcv for this name
             do jc=1,srcv(jx)%nct
               srcv(jx)%nid(jc) = cpl_vinfo%tag
             enddo
           endif
         enddo
         if ( rank == ocn_master ) then
           !--- Write to NEMO's ocean.output file
           write(numout,*) "cpl_cancpl_define: Recv field ",ji, &
               "  name=",trim(nemo_recv_var(ji))," tag=",cpl_vinfo%tag
           call flush(numout)

           !--- Also write to stdout (unit 6)
           write(6,*) "cpl_cancpl_define: Recv field ",ji, &
               "  name=",trim(nemo_recv_var(ji))," tag=",cpl_vinfo%tag
           call flush(6)
         endif
       enddo
     endif

  end subroutine cpl_cancpl_define

  subroutine copy_1d_to_3d_global(wrk, png)
    !------------------------------------------------------------------------
    !--- Copy values from a global 1D wrk array containing data received from
    !--- the agcm to global 3D array suitable for use with mppscatter
    !------------------------------------------------------------------------
    real(kind=8), intent(in) :: wrk(:)
    real(wp), intent(out) :: png(jpi,jpj,jpnij)

    !--- Local
    real(kind=8) :: glob_arr(jpiglo, jpjglo)
!xxx    integer :: ji, jj, jn, ji_glob, jj_glob

    glob_arr = 0.0_8
    glob_arr(1:jpiglo,1:jpjglo) = reshape( wrk(1:jpiglo*jpjglo), (/ jpiglo,jpjglo /) )

    call copy_2d_to_3d_global(glob_arr, png)

!xxx    do jn = 1,jpnij
!xxx      !--- jn loops over all subdomains
!xxx      png(:,:,jn) = 0.0_8
!xxx      do ji=nldit(jn),nleit(jn)
!xxx        do jj=nldjt(jn),nlejt(jn)
!xxx          !--- nimppt(jn),njmppt(jn) are the global indicies corresponding to the
!xxx          !--- (1,1) grid cell in the local index space of the current subdomain
!xxx          ji_glob = ji + nimppt(jn) - 1
!xxx          jj_glob = jj + njmppt(jn) - 1
!xxx          if ( ji_glob < 1      .or. jj_glob < 1 .or. &
!xxx               ji_glob > jpiglo .or. jj_glob > jpjglo ) then
!xxx            write(6,*)'copy_1d_to_3d_global: Global index is out of range.'
!xxx            write(6,*)'jn, ji, jj, ji_glob, jj_glob: ',jn, ji, jj, ji_glob, jj_glob
!xxx            call ctl_stop("STOP", "copy_1d_to_3d_global", "Global index is out of range")
!xxx          endif
!xxx          png(ji,jj,jn) = glob_arr(ji_glob,jj_glob)
!xxx        enddo
!xxx      enddo
!xxx    enddo

  end subroutine copy_1d_to_3d_global

  subroutine copy_2d_to_3d_global(glob_a2d, png)
    !------------------------------------------------------------------------
    !--- Copy values from a global 2D array dimensioned (jpiglo, jpjglo)
    !--- to a global 3D array suitable for use with mppscatter
    !------------------------------------------------------------------------
    real(kind=8) :: glob_a2d(jpiglo, jpjglo)
    real(wp), intent(out) :: png(jpi,jpj,jpnij)

    !--- Local
    integer :: ji, jj, jn, ji_glob, jj_glob

    do jn = 1,jpnij
      !--- jn loops over all subdomains
      png(:,:,jn) = 0.0_8
      do ji=nldit(jn),nleit(jn)
        do jj=nldjt(jn),nlejt(jn)
          !--- nimppt(jn),njmppt(jn) are the global indicies corresponding to the
          !--- (1,1) grid cell in the local index space of the current subdomain
          ji_glob = ji + nimppt(jn) - 1
          jj_glob = jj + njmppt(jn) - 1
          if ( ji_glob < 1      .or. jj_glob < 1 .or. &
               ji_glob > jpiglo .or. jj_glob > jpjglo ) then
            write(6,*)'copy_2d_to_3d_global: Global index is out of range.'
            write(6,*)'jn, ji, jj, ji_glob, jj_glob: ',jn, ji, jj, ji_glob, jj_glob
            call ctl_stop("STOP", "copy_2d_to_3d_global", "Global index is out of range")
          endif
          png(ji,jj,jn) = glob_a2d(ji_glob,jj_glob)
        enddo
      enddo
    enddo

  end subroutine copy_2d_to_3d_global

  subroutine copy_3d_to_1d_global(wrk, png)
    !------------------------------------------------------------------------
    !--- Copy values from a global 3D array containing data from a recent
    !--- call to mppgather to a global 1D wrk array to be sent to the agcm
    !------------------------------------------------------------------------
    real(kind=8), intent(out) :: wrk(:)
    real(wp), intent(in) :: png(jpi,jpj,jpnij)

    !--- Local
    real(kind=8) :: glob_arr(jpiglo, jpjglo)
!xxx    integer :: ji, jj, jn, ji_glob, jj_glob

!xxx    glob_arr = 0.0_8
!xxx    do jn = 1,jpnij
!xxx      !--- jn loops over all subdomains
!xxx      do ji=nldit(jn),nleit(jn)
!xxx        do jj=nldjt(jn),nlejt(jn)
!xxx          !--- nimppt(jn),njmppt(jn) are the global indicies corresponding to the
!xxx          !--- (1,1) grid cell in the local index space of the current subdomain
!xxx          ji_glob = ji + nimppt(jn) - 1
!xxx          jj_glob = jj + njmppt(jn) - 1
!xxx          if ( ji_glob < 1      .or. jj_glob < 1 .or. &
!xxx               ji_glob > jpiglo .or. jj_glob > jpjglo ) then
!xxx            write(6,*)'copy_3d_to_1d_global: Global index is out of range.'
!xxx            write(6,*)'jn, ji, jj, ji_glob, jj_glob: ',jn, ji, jj, ji_glob, jj_glob
!xxx            call ctl_stop("STOP", "copy_3d_to_1d_global", "Global index is out of range")
!xxx          endif
!xxx          glob_arr(ji_glob,jj_glob) = png(ji,jj,jn)
!xxx        enddo
!xxx      enddo
!xxx    enddo

    call  copy_3d_to_2d_global(glob_arr, png)

    wrk = 0.0_8
    wrk(1:jpiglo*jpjglo) = reshape( glob_arr(1:jpiglo,1:jpjglo), (/ jpiglo*jpjglo /) )

  end subroutine copy_3d_to_1d_global

  subroutine copy_3d_to_2d_global(glob_a2d, png)
    !------------------------------------------------------------------------
    !--- Copy values from a global 3D array containing data from a recent
    !--- call to mppgather to a global 2D array to be sent to the agcm
    !------------------------------------------------------------------------
    real(kind=8), intent(out) :: glob_a2d(jpiglo, jpjglo)
    real(wp), intent(in) :: png(jpi,jpj,jpnij)

    !--- Local
    integer :: ji, jj, jn, ji_glob, jj_glob

    glob_a2d = 0.0_8
    do jn = 1,jpnij
      !--- jn loops over all subdomains
      do ji=nldit(jn),nleit(jn)
        do jj=nldjt(jn),nlejt(jn)
          !--- nimppt(jn),njmppt(jn) are the global indicies corresponding to the
          !--- (1,1) grid cell in the local index space of the current subdomain
          ji_glob = ji + nimppt(jn) - 1
          jj_glob = jj + njmppt(jn) - 1
          if ( ji_glob < 1      .or. jj_glob < 1 .or. &
               ji_glob > jpiglo .or. jj_glob > jpjglo ) then
            write(6,*)'copy_3d_to_2d_global: Global index is out of range.'
            write(6,*)'jn, ji, jj, ji_glob, jj_glob: ',jn, ji, jj, ji_glob, jj_glob
            call ctl_stop("STOP", "copy_3d_to_2d_global", "Global index is out of range")
          endif
          glob_a2d(ji_glob,jj_glob) = png(ji,jj,jn)
        enddo
      enddo
    enddo

  end subroutine copy_3d_to_2d_global


  subroutine cpl_cancpl_snd( kid, kstep, pdata, kinfo )
     !!---------------------------------------------------------------------
     !!              ***  ROUTINE cpl_cancpl_snd  ***
     !!
     !! ** Purpose : - At each coupling time step,this routine sends the field
     !!                associated with the index kid to the coupler
     !!----------------------------------------------------------------------
     !--- Index in the array ssnd of the variable to be sent to the coupler
     integer,  intent(in)  :: kid

     !--- Ocean time in seconds
     integer,  intent(in)  :: kstep

     !--- pdata contains data on the local domain (local MPI task) to be sent
     !--- It will be dimensioned pdata(jpi, jpj, ssnd(kid)%nct)
     real(wp), intent(in)  :: pdata(:,:,:)

     !--- Integer flag to indicate if srcv(kid) was sent or not
     !--- kinfo = OASIS_idle means the field was not sent to the coupler
     !--- kinfo = OASIS_snd  means the field was sent to the coupler
     integer,  intent(out) :: kinfo

     !-- Local
     integer :: jc
     integer :: ldbg=1
     integer :: freq
     integer(kind=impi) :: rank, ierr
     integer :: idx, nwrds
     integer :: verbose=3
     !!--------------------------------------------------------------------

     !--- Determine the rank of the calling process in MPI_COMM_WORLD
     call mpi_comm_rank ( MPI_COMM_WORLD, rank, ierr )

     !--- Loop over all "catagories" (normally only 1) for this variable
     !--- sending separate data to the coupler for each catagory
     do jc = 1, ssnd(kid)%nct
       !--- The MPI tag associated with this transfer is ssnd(kid)%nid(jc)

       kinfo = OASIS_idle

       !--- Ensure that this variable was configured
       if ( ssnd(kid)%nid(jc) <= 0 ) then
         write(numout,*)'cpl_cancpl_snd: ',trim(ssnd(kid)%clname), &
           ' is not sent in this configuration.  kid = ',kid
         call ctl_stop("STOP", "cpl_cancpl_snd", &
                       "Invalid send variable "//trim(ssnd(kid)%clname))
       endif

       !--- Determine if this variable should be coupled now
       cpl_vinfo = find_cpl_vinfo( tag=ssnd(kid)%nid(jc) )
       freq = cpl_vinfo%freq

       !--- Ignore the rest of this loop if this is not a coupling time step
       if ( mod(kstep,freq) /= 0 ) then
         if ( rank == ocn_master .and. verbose > 2 ) then
           write(numout,'(a,i8,3a,i8)')'cpl_cancpl_snd: ',kstep, &
             ' is not a coupling time step for ',trim(ssnd(kid)%clname),'  freq=',freq
           call flush(numout)
         endif
         cycle
       endif

       if ( rank == ocn_master .and. verbose > 2 ) then
         write(numout,*)'cpl_cancpl_snd: NEMO sending ',trim(ssnd(kid)%clname), &
             ' from task ',rank,' to task ',cpl_master,'  kstep=',kstep,'  freq=',freq, &
             ' catagory=',jc
         call flush(numout)
       endif

       kinfo = OASIS_Snd

       if ( ln_ctl .and. verbose > 1 ) then
         !--- Write info for each sub-domain to the ocean output file
         write(numout,*) '****************'
         write(numout,*) 'cpl_cancpl_snd: Outgoing ', ssnd(kid)%clname
         write(numout,*) 'cpl_cancpl_snd:      tag ', ssnd(kid)%nid(jc)
         write(numout,*) 'cpl_cancpl_snd:    kstep ', kstep
         write(numout,*) 'cpl_cancpl_snd: mpi task ', rank
         write(numout,*) '      - minimum value is ', minval(pdata(:,:,jc))
         write(numout,*) '      - maximum value is ', maxval(pdata(:,:,jc))
         write(numout,*) '      -     sum value is ', sum(pdata(:,:,jc))
         write(numout,*) '****************'
       endif

       !--- Gather data into the global array png
       call mppsync
       call mppgather (pdata(:,:,jc),0,png)
       call mppsync

       !--- Skip the rest of this loop unless this is the master task
       if ( rank /= ocn_master ) cycle

       if ( verbose > 2 ) then
         !--- Count the number of NaNs in the global png array
         idx = count( png /= png )
         write(numout,*)'cpl_cancpl_snd: ',trim(ssnd(kid)%clname),'  Nans in png = ',idx
         write(numout,*)'cpl_cancpl_snd: ',trim(ssnd(kid)%clname),'  min,max,avg = ', &
             minval(png),maxval(png),sum(png)/real(size(png),kind=8)
         call flush(numout)
       endif

       !--- Map the the global 3D array png onto the 1D wrk array
       call copy_3d_to_1d_global(wrk, png)

       if ( verbose > 2 ) then
         !--- Count the number of NaNs in the wrk array
         nwrds = jpiglo*jpjglo
         idx = count( wrk(1:nwrds) /= wrk(1:nwrds) )
         write(numout,*)'cpl_cancpl_snd: ',trim(ssnd(kid)%clname),'  Nans in wrk = ',idx
         write(numout,*)'cpl_cancpl_snd: ',trim(ssnd(kid)%clname),'  min,max,avg = ', &
             minval(wrk(1:nwrds)),maxval(wrk(1:nwrds)),sum(wrk(1:nwrds))/real(nwrds,kind=8)
         call flush(numout)
       endif

       if ( verbose > 3 ) then
         write(numout,*) &
           'cpl_cancpl_snd: kstep=',kstep,'  Sending ', trim(ssnd(kid)%clname), &
           ' from task ',rank,' to task ',cpl_master
         call flush(numout)
       endif

       !--- Send the global array to the coupler
       call send_data_rec(wrk, cpl_master, trim(ssnd(kid)%clname), dbg=ldbg)

       if ( verbose > 0 ) then
         strng = sprint_var_stats(wrk, cpl_vinfo%size, name=ssnd(kid)%clname)
         write(numout,*)'cpl_cancpl_snd: kstep=',kstep,'  Sent ', &
             trim(ssnd(kid)%clname),'  ',trim(strng)
         call flush(numout)
       endif

     enddo

  end subroutine cpl_cancpl_snd


  subroutine cpl_cancpl_rcv( kid, kstep, pdata, kinfo )
     !!---------------------------------------------------------------------
     !!              ***  ROUTINE cpl_cancpl_rcv  ***
     !!
     !! ** Purpose : - At each coupling time-step,this routine receives fields
     !!      like stresses and fluxes from the coupler or remote application.
     !!----------------------------------------------------------------------

     !--- Index in the array srcv of the variable to be received from the coupler
     integer,  intent(in   ) :: kid

     !--- ocean time in seconds (ie number of seconds since time step nit000)
     integer,  intent(in   ) :: kstep

     !--- Space for data to be received
     !--- IN to keep the value if nothing is done
     !--- pdata contains data on the local domain (local MPI task) to be received
     !--- It will be dimensioned pdata(jpi, jpj, srcv(kid)%nct)
     real(wp), intent(inout) :: pdata(:,:,:)

     !--- Integer flag to indicate if srcv(kid) was recieved or not
     !--- kinfo = OASIS_idle means the field was not received from the coupler
     !--- kinfo = OASIS_rcv  means the field was received from the coupler
     integer,  intent(  out) :: kinfo

     !--- Local
     integer :: jc
     integer :: ldbg=1
     integer :: freq, idx, nwrds
     integer(kind=impi) :: rank, ierr, sz, tag
     integer :: verbose=3
     integer (kind=impi) :: status(MPI_status_size)
     !!--------------------------------------------------------------------

     !---Determine the rank of the calling process in MPI_COMM_WORLD
     call mpi_comm_rank ( MPI_COMM_WORLD, rank, ierr )

     !--- Loop over all "catagories" (normally only 1) for this variable
     !--- receiving separate data from the coupler for each catagory
     do jc = 1, srcv(kid)%nct
       !--- The MPI tag associated with this transfer is srcv(kid)%nid(jc)

       kinfo = OASIS_idle

       !--- This routine is called for every variable that could be coupled
       !--- so ignore variables that are not to be coupled
       !--- srcv(:)%nid(:) is initialized to zero then defined in
       !--- cpl_cancpl_define for variables that are to be coupled and so
       !--- it will only be non-zero for fields that are coupled
       if ( srcv(kid)%nid(jc) <= 0 ) cycle

       !--- Determine if this variable should be coupled now
       cpl_vinfo = find_cpl_vinfo( tag=srcv(kid)%nid(jc) )
       freq = cpl_vinfo%freq

       !--- Ignore the rest of this loop if this is not a coupling time step
       if ( mod(kstep,freq) /= 0 ) then
         if ( rank == ocn_master .and. verbose > 2 ) then
           write(numout,'(a,i8,3a,i8)')'cpl_cancpl_rvc: ',kstep, &
             ' is not a coupling time step for ',trim(srcv(kid)%clname),'  freq=',freq
           call flush(numout)
         endif
         cycle
       endif

       kinfo = OASIS_Rcv

       if ( rank == ocn_master ) then
         !--- This is the ocean master task
         if ( verbose > 2 ) then
           write(numout,'(3a,i2,a,i4,a,i8,a,i8)') 'cpl_cancpl_rcv: NEMO receiving ', &
             trim(srcv(kid)%clname),' from task ',cpl_master, &
             '  tag=',srcv(kid)%nid(jc),'  kstep=',kstep,'  freq=',freq
           call flush(numout)
         endif

         !--- Receive the global array from the coupler
         call recv_data_rec(wrk, ibuf, cpl_master, trim(srcv(kid)%clname), dbg=ldbg)

         !--- Map the 1D wrk array onto the global 3D array png
         call copy_1d_to_3d_global(wrk, png)
       endif

       !--- Scatter the global array onto each NEMO task
       call mppsync
       call mppscatter (png,0,pdata(:,:,jc)) 
       call mppsync

       if ( rank == ocn_master .and. verbose > 2 ) then
         !--- Count the number of NaNs in the global png array
         idx = count( png /= png )
         write(numout,*)'cpl_cancpl_rcv: ',trim(srcv(kid)%clname),'  Before lbc_lnk'
         write(numout,*)'cpl_cancpl_rcv: ',trim(srcv(kid)%clname),'  Nans in png = ',idx
         write(numout,*)'cpl_cancpl_rcv: ',trim(srcv(kid)%clname),'  min,max,avg = ', &
             minval(png),maxval(png),sum(png)/real(size(png),kind=8)
         call flush(numout)
       endif

       !--- Fill overlap areas and extra hallows and check periodicity
       call lbc_lnk( pdata(:,:,jc), srcv(kid)%clgrid, srcv(kid)%nsgn )

       if ( rank == ocn_master .and. verbose > 2 ) then
         !--- Count the number of NaNs in the global png array
         idx = count( png /= png )
           ! write(numout,*)'cpl_cancpl_rcv: ',trim(srcv(kid)%clname),'  After lbc_lnk'
           write(numout,*)'cpl_cancpl_rcv: ',trim(srcv(kid)%clname),'  Nans in png = ',idx
           ! write(numout,*)'cpl_cancpl_rcv: ',trim(srcv(kid)%clname),'  min,max,avg = ', &
           !     minval(png),maxval(png),sum(png)/real(size(png),kind=8)
           call flush(numout)
       endif

       if ( ln_ctl .and. verbose > 1 ) then
         !--- Write info for each sub-domain to the ocean output file
         write(numout,*) '****************'
         write(numout,*) 'cpl_cancpl_rcv: Incoming ', srcv(kid)%clname
         write(numout,*) 'cpl_cancpl_rcv: ivarid ',   srcv(kid)%nid(jc)
         write(numout,*) 'cpl_cancpl_rcv:   kstep', kstep
         write(numout,*) 'cpl_cancpl_rcv:   info ', kinfo
         write(numout,*) '     - minimum value is ', minval(pdata(:,:,jc))
         write(numout,*) '     - maximum value is ', maxval(pdata(:,:,jc))
         write(numout,*) '     -     sum value is ', sum(pdata(:,:,jc))
         write(numout,*) '****************'
       endif

       if ( rank == ocn_master .and. verbose > 0 ) then
         strng = sprint_var_stats(wrk, cpl_vinfo%size, name=srcv(kid)%clname)
         write(numout,*)'cpl_cancpl_rcv: kstep=',kstep,'  Received ', &
             trim(srcv(kid)%clname),'  ',trim(strng)
         call flush(numout)
       endif
         
     enddo

  end subroutine cpl_cancpl_rcv


  function cpl_cancpl_freq( vname ) result(freq)
    !!---------------------------------------------------------------------
    !!              ***  ROUTINE cpl_cancpl_freq  ***
    !!
    !! ** Purpose : - send back the coupling frequency for a particular field
    !!----------------------------------------------------------------------
    !--- vname is the name of the variable whose frequency is requested
    !--- This is the name known to the coupler e.g. srcv(:)%clname
    character(*), intent(in) :: vname
    integer :: freq
    !!----------------------------------------------------------------------

!xxx    !--- Find the coupling frequency for this variable
!xxx    cpl_vinfo = find_cpl_vinfo( name=trim(vname) )
!xxx    freq = cpl_vinfo%freq

    !--- cpl_ocn_freq is defined in the com_cpl module
    freq = cpl_ocn_freq

    if ( freq < 1 ) then
      write(6,*)"cpl_cancpl_freq: Invalid coupling frequency ",freq," for ",trim(vname)
      call flush(6)
      call ctl_stop("STOP", " cpl_cancpl_freq", "Invalid coupling frequency")
    endif

  end function cpl_cancpl_freq


  subroutine cpl_cancpl_finalize
    !!---------------------------------------------------------------------
    !!              ***  ROUTINE cpl_cancpl_finalize  ***
    !!
    !! ** Purpose : - Finalizes the coupling. If MPI_init has not been
    !!      called explicitly before cpl_cancpl_init it will also close
    !!      MPI communication.
    !!----------------------------------------------------------------------

    if ( allocated(png) ) DEALLOCATE( png )
    call mppstop
    !--- TODO --- Also tell coupler that the ocean has stopped

  end subroutine cpl_cancpl_finalize

#else
   !!----------------------------------------------------------------------
   !!  Serial mode (no MPI)       Dummy module      Forced Ocean/Atmosphere
   !!----------------------------------------------------------------------
   use in_out_manager               ! i/o manager
   logical, public, parameter :: lk_cpl = .false.   !: coupled flag
   public cpl_cancpl_init
   public cpl_cancpl_finalize
contains
   subroutine cpl_cancpl_init (kl_comm) 
      integer, intent(out)   :: kl_comm
      kl_comm = -1
      write(numout,*) 'cpl_cancpl_init: Called in serial mode, MPI not in use.'
   end subroutine cpl_cancpl_init
   subroutine cpl_cancpl_finalize
      write(numout,*) 'cpl_cancpl_finalize: Called in serial mode, MPI not in use.'
   end subroutine cpl_cancpl_finalize
#endif

   !!=====================================================================
end module cpl_cancpl
