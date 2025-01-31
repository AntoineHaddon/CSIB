MODULE cpl_cancpl
!DIR$ NOOPTIMIZE
  !!=======================================================================
  !!                    ***  MODULE cpl_cancpl  ***
  !!
  !! Coupled ocean-atmosphere using the coupler CanCPL
  !!
  !! This module is intended to serve the same function as cpl_oasis3 does
  !! in nemo v3.4.1, except (of course) it will communicate with CanCPL
  !!
  !! Larry Solheim Aug,2015
#if defined key_cancpl
  !!=======================================================================
  !!-----------------------------------------------------------------------
  !!   'lk_cancpl'                    coupled Ocean/Atmosphere via CanCPL
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
  use iom            ! I/O library
  use in_out_manager               ! I/O manager
  use lbclnk                       ! ocean lateral boundary conditions (or mpp link)
  use timing
  use lib_mpp, only : mpi_comm_oce, ctl_stop, mppgather, mppsync, mppscatter, mppstop
  use lib_mpp, only : reconstruct_global_2d
  use cpl_types, only : srcv, ssnd, FLD_C, FLD_CPL, nmaxfld
  USE iom            ! I/O library

  implicit none
  private

  public :: cpl_cancpl_init
  public :: cpl_cancpl_define
  public :: cpl_cancpl_snd
  public :: cpl_cancpl_rcv
  public :: cpl_cancpl_freq
  public :: cpl_cancpl_finalize
  public :: set_cancpl_params
  public :: check_value2d, check_value3d
  public :: query_start_cpl2ocn

  integer, public, save      ::   oasis_idle = 0    !: return code if no send or recv
  integer, public, save      ::   oasis_rcv  = 1    !: return code if field received
  integer, public, save      ::   oasis_snd  = 2    !: return code if field sent

  !--- These are defined in com_cpl
  public :: cpl_vinfo_t, find_cpl_vinfo
  public :: bcast_inter, cpl_time_string, cpl_elapsed_time_secs, cpl_master

  integer            ::   nerror        ! return error code

  integer :: nn_fsbc, nn_ice

  !--- tmp char space
  character(512), save :: strng

  !--- Work space used to temporarily hold fields passed via MPI
  integer, parameter    :: maxx=8392704  !---4098x2048 = (2+2^12)x(2^11)
  real(kind=8), save    :: wrk(maxx)
  integer(kind=8), save :: ibuf(8)

  !--- A derived type holding coupler related information
  !--- cpl_vinfo_t is defined in the com_cpl module
  type(cpl_vinfo_t), save :: cpl_vinfo

  !--- NOTE: nproc is not equal to the MPI task in model_communicator since it will always
  !--- be one of 0,1,2,...(jpnij-1) and the AGCM gets the first set of MPI tasks.
  !--- However nproc == 0 should still correspond with the ocn_master task
  !--- nproc is use associated through the module dom_oce
   
  !--- Coupler dimension limits
  integer :: cpl_jpiglo, cpl_jpjglo

  !real(kind=8), pointer, public :: work(:,:) 

contains


  !DBG
  subroutine check_value2d(name,var,kt,msg,mode)
    character(*) :: name
    real(wp) :: var(:,:)
    integer :: kt
    character(*), optional :: msg
    integer, optional :: mode
    real(wp) :: absmaxval
    integer :: lmode, iu
    logical :: exists
    character(256) :: strng
    if ( present(mode) ) then
      lmode = mode
    else
      lmode = 0
    endif
    absmaxval = max(abs(maxval(var)), abs(minval(var)))
    if ( absmaxval /= 0.0_wp ) then
      if ( present(msg) ) then
         write(numout,*)"** EE ** kt=",kt,"  absolute max value of ", &
              trim(name)," on tile ",narea," is non-zero.   ",absmaxval,"  ",trim(msg)
      else
         write(numout,*)"** EE ** kt=",kt,"  absolute max value of ", &
              trim(name)," on tile ",narea," is non-zero.   ",absmaxval
      endif
      if ( lmode == 0 ) then
        call ctl_stop("STOP", " check_value2d", trim(name)//" is out of range")
      endif
    endif
    if ( lmode == 2 ) then
      strng=" "
      write(strng,'("out.",a,"_",i4.4)')trim(name),narea-1
      inquire(file=trim(strng),exist=exists)
      if ( .not. exists ) then
        !--- Only write the first time this variable is passed
        iu = 827+narea
        open(iu,file=trim(strng),form="unformatted")
        write(iu) var
        close(iu)
      endif
    endif
  end subroutine check_value2d

  subroutine check_value3d(name,var,kt,msg,mode)
    character(*) :: name
    real(wp) :: var(:,:,:)
    integer :: kt
    character(*), optional :: msg
    integer, optional :: mode
    real(wp) :: absmaxval
    integer :: lmode, iu
    logical :: exists
    character(256) :: strng
    if ( present(mode) ) then
      lmode = mode
    else
      lmode = 0
    endif
    absmaxval = max(abs(maxval(var)), abs(minval(var)))
    absmaxval = max(abs(maxval(var)), abs(minval(var)))
    if ( absmaxval /= 0.0_wp ) then
      if ( present(msg) ) then
         write(numout,*)"** EE ** kt=",kt,"  absolute max value of ", &
              trim(name)," on tile ",narea," is non-zero.   ",absmaxval,"  ",trim(msg)
      else
         write(numout,*)"** EE ** kt=",kt,"  absolute max value of ", &
              trim(name)," on tile ",narea," is non-zero.   ",absmaxval
      endif
      if ( lmode == 0 ) then
        call ctl_stop("STOP", " check_value2d", trim(name)//" is out of range")
      endif
    endif
    if ( lmode == 2 ) then
      strng=" "
      write(strng,'("out.",a,"_",i4.4)')trim(name),narea-1
      inquire(file=trim(strng),exist=exists)
      if ( .not. exists ) then
        !--- Only write the first time this variable is passed
        iu = 927+narea
        open(iu,file=trim(strng),form="unformatted")
        write(iu) var
        close(iu)
      endif
    endif
  end subroutine check_value3d

  subroutine dump_array1d(name,var,pos)
    character(*) :: name
    real(wp) :: var(:)
    integer, optional :: pos
    integer :: lpos, iu
    logical :: exists
    character(256) :: strng
    if ( present(pos) ) then
      lpos = pos
    else
      lpos = 0
    endif
    strng=" "
    write(strng,'("out.",a,"_",i4.4)')trim(name),narea-1
    iu = 727+narea
    if ( lpos == 0 ) then
      !--- Only write the first time this name is used
      inquire(file=trim(strng),exist=exists)
      if ( .not. exists ) then
        open(iu,file=trim(strng),form="unformatted")
        write(iu) var
        close(iu)
      endif
    else if ( lpos == 1 ) then
      !--- Overwrite the file each time a write is requested
      open(iu,file=trim(strng),form="unformatted")
      rewind(iu)
      write(iu) var
      close(iu)
    else if ( lpos == 2 ) then
      !--- Append to the file each time a write is requested
      open(iu,file=trim(strng),form="unformatted",position="append")
      write(iu) var
      close(iu)
    endif
  end subroutine dump_array1d

  subroutine dump_array2d(name,var,pos)
    character(*) :: name
    real(wp) :: var(:,:)
    integer, optional :: pos
    integer :: lpos, iu
    logical :: exists
    character(256) :: strng
    if ( present(pos) ) then
      lpos = pos
    else
      lpos = 0
    endif
    strng=" "
    write(strng,'("out.",a,"_",i4.4)')trim(name),narea-1
    iu = 727+narea
    if ( lpos == 0 ) then
      !--- Only write the first time this name is used
      inquire(file=trim(strng),exist=exists)
      if ( .not. exists ) then
        open(iu,file=trim(strng),form="unformatted")
        write(iu) var
        close(iu)
      endif
    else if ( lpos == 1 ) then
      !--- Overwrite the file each time a write is requested
      open(iu,file=trim(strng),form="unformatted")
      rewind(iu)
      write(iu) var
      close(iu)
    else if ( lpos == 2 ) then
      !--- Append to the file each time a write is requested
      open(iu,file=trim(strng),form="unformatted",position="append")
      write(iu) var
      close(iu)
    endif
  end subroutine dump_array2d

  subroutine dump_array3d(name,var,pos)
    character(*) :: name
    real(wp) :: var(:,:,:)
    integer, optional :: pos
    integer :: lpos, iu
    logical :: exists
    character(256) :: strng
    if ( present(pos) ) then
      lpos = pos
    else
      lpos = 0
    endif
    strng=" "
    write(strng,'("out.",a,"_",i4.4)')trim(name),narea-1
    iu = 727+narea
    if ( lpos == 0 ) then
      !--- Only write the first time this name is used
      inquire(file=trim(strng),exist=exists)
      if ( .not. exists ) then
        open(iu,file=trim(strng),form="unformatted")
        write(iu) var
        close(iu)
      endif
    else if ( lpos == 1 ) then
      !--- Overwrite the file each time a write is requested
      open(iu,file=trim(strng),form="unformatted")
      rewind(iu)
      write(iu) var
      close(iu)
    else if ( lpos == 2 ) then
      !--- Append to the file each time a write is requested
      open(iu,file=trim(strng),form="unformatted",position="append")
      write(iu) var
      close(iu)
    endif
  end subroutine dump_array3d

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
     !    This will, among other things, define the communicator for all model components,
     !      as well as cpl_master, atm_master, ocn_master
     !      It also returns an local ocean only communicator as local_ocn_comm
     !    NOTE: ocn_master is the rank in model_communicator not the rank in local_ocn_comm
     local_ocn_comm = -1
     call define_group('ocn', local_ocn_comm)
     call init_common_coupler_parameters()

     kl_comm  = local_ocn_comm

     !---Determine the rank of the calling process in model_communicator
     call mpi_comm_rank ( model_communicator, rank, ierr )

     !--- ocn_master is defined during the call to define_group
     if ( rank == ocn_master ) then
       !--- Write to stdout (unit 6) since numout is not attached yet
       write(6,*) 'cpl_cancpl_init: ocean comm_world = ',kl_comm
       call flush(6)
     endif

  end subroutine cpl_cancpl_init


  subroutine cpl_cancpl_define( krcv, ksnd, kcplmodel )
     !!-------------------------------------------------------------------
     !!             ***  ROUTINE cpl_cancpl_define  ***
     !!
     !! ** Purpose :   Define grid and field information for ocean
     !!    exchange between AGCM, OGCM and COUPLER.
     !!--------------------------------------------------------------------
     integer, intent(in) :: krcv   ! Number of all possible fields received
     integer, intent(in) :: ksnd   ! Number of all possible fields sent
     integer, intent(in) :: kcplmodel ! Number of models to send too. Note this is a dummy argument for now
                                      ! so that the interface matches the oasis equivalent

     !--- Local
     integer :: ji,jc,jx
     character(len=8) :: zclname
     integer :: ldbg=1
     integer(kind=impi) :: rank, ierr
     integer :: verbose=2
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

     !--- Determine the rank of the calling process in model_communicator
     call mpi_comm_rank ( model_communicator, rank, ierr )

     if ( rank == ocn_master ) then
       write(numout,*)
       write(numout,*) 'cpl_cancpl_define: initialize coupled ocean/atmosphere'
       write(numout,*) '~~~~~~~~~~~~~~~~~'
       write(numout,*)
       call flush(numout)
     endif

     if ( rank == ocn_master .and. verbose > 1 ) then
       write(numout,*) 'cpl_cancpl_define: Assign sent variables'
       call flush(numout)
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
              if ( nemo_n_send_var > nmaxfld ) then
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
              elseif ( rank == ocn_master ) then
                write(6,*)"cpl_cancpl_define: rank for ", &
                          trim(var_list_info(nemo_n_send_var)%name),' is ',var_list_info(nemo_n_send_var)%rank
                call flush(6)
              endif

           end do
        endif
     end do
     if ( rank == ocn_master .and. verbose > 1 ) then
       write(numout,*) 'cpl_cancpl_define: nemo_n_send_var=',nemo_n_send_var
       call flush(numout)
     endif
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

     if ( rank == ocn_master .and. verbose > 1 ) then
       write(numout,*) 'cpl_cancpl_define : Assign received variables'
       call flush(numout)
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
              if ( nemo_n_recv_var > nmaxfld ) then
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
     if ( rank == ocn_master .and. verbose > 1 ) then
       write(numout,*) 'cpl_cancpl_define: nemo_n_recv_var=',nemo_n_recv_var
       call flush(numout)
     endif
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

     !--- Set a value for nemo_rn_rdt, defined in com_cpl
     !--- rn_Dt is defined in the module dom_oce
     nemo_rn_rdt = nint(rn_Dt,8)

     !--- Set a value for nemo_nn_ice, defined in com_cpl
     nemo_nn_ice = nn_ice

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
     !--- Ni0glo and Nj0glo are defined in the namelist ("interior" grid, no halos)
     !--- NOTE: the nemo_* variables are sent to the coupler
     cpl_jpiglo = Ni0glo
     cpl_jpjglo = Nj0glo 
     nemo_jpiglo = cpl_jpiglo
     nemo_jpjglo = cpl_jpjglo

     if ( rank == ocn_master .and. verbose > 1 ) then
       write(6,*)"cpl_cancpl_define: READ namsbc_cpl namelist"
       call flush(6)
     endif

     ! lat/lon and grid size read from domcfg to have the vavlues even over land eliminated processors
     if ( rank == ocn_master ) then
         write(6,*)"cpl_cancpl_define: READ domain variables from ",trim(cn_domcfg)
         call flush(6)
         IF (.NOT. ASSOCIATED(nemo_glamt)) ALLOCATE(nemo_glamt(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_glamu)) ALLOCATE(nemo_glamu(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_glamv)) ALLOCATE(nemo_glamv(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_glamf)) ALLOCATE(nemo_glamf(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_gphit)) ALLOCATE(nemo_gphit(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_gphiu)) ALLOCATE(nemo_gphiu(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_gphiv)) ALLOCATE(nemo_gphiv(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_gphif)) ALLOCATE(nemo_gphif(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_e1t)) ALLOCATE(nemo_e1t(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_e1u)) ALLOCATE(nemo_e1u(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_e1v)) ALLOCATE(nemo_e1v(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_e1f)) ALLOCATE(nemo_e1f(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_e2t)) ALLOCATE(nemo_e2t(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_e2u)) ALLOCATE(nemo_e2u(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_e2v)) ALLOCATE(nemo_e2v(nemo_jpiglo,nemo_jpjglo))
         IF (.NOT. ASSOCIATED(nemo_e2f)) ALLOCATE(nemo_e2f(nemo_jpiglo,nemo_jpjglo))
         CALL hgr_read_glo( nemo_glamt , nemo_glamu , nemo_glamv  , nemo_glamf  ,   &    ! gridpoints position (required)
             &              nemo_gphit , nemo_gphiu , nemo_gphiv  , nemo_gphif  ,   &     
             &              nemo_e1t  , nemo_e1u  , nemo_e1v   , nemo_e1f   ,   &    ! scale factors       (required)
             &              nemo_e2t  , nemo_e2u  , nemo_e2v   , nemo_e2f   )
     endif

     ! Mask fields not in the domcfg file, must relu on the calcuated done in dommsk
     !--- Gather tmask into nemo_tmask (found in com_cpl)
     !--- tmask_i is found in module dom_oce
     call cpl_gather("tmask_i", rank)

     !--- Gather umask (level 1) into nemo_umask (found in com_cpl)
     !--- umask is found in module dom_oce
     call cpl_gather("umask", rank)

     !--- Gather vmask (level 1) into nemo_vmask (found in com_cpl)
     !--- vmask is found in module dom_oce
     call cpl_gather("vmask", rank)

     !--- Gather fmask (level 1) into nemo_fmask (found in com_cpl)
     !--- fmask is found in module dom_oce
     call cpl_gather("fmask", rank)

     if ( rank == ocn_master .and. verbose > 1 ) then
       write(6,*)"cpl_cancpl_define: call cpl_initialize_events()"
       call flush(6)
     endif

     !--- Initialize coupler events and broadcast global variables
     call cpl_initialize_events()

     if ( rank == ocn_master .and. verbose > 1 ) then
       write(6,*)"cpl_cancpl_define: call bcastGroup(cpl_time_string, cpl_master, model_communicator)"
       write(6,*)"cpl_cancpl_define: cpl_master = ",cpl_master
       write(6,*)"cpl_cancpl_define: model_communicator = ",model_communicator
       call flush(6)
     endif

     !--- Broadcast the initial date and time from the coupler to all tasks
     !--- cpl_time_string is defined in com_cpl
     call bcastGroup(cpl_time_string, cpl_master, model_communicator)

     ! -----------------------------------------------------------------
     ! ... Assign MPI tags to ssnd and srcv variables
     !--- This must be done after the call to cpl_initialize_events
     !--- because it will define these tags
     ! -----------------------------------------------------------------
     if ( nemo_n_send_var > 0 ) then
       do ji=1,nemo_n_send_var
         cpl_vinfo = find_cpl_vinfo( name=trim(nemo_send_var(ji)) )
         do jx=1,nmaxfld
           if ( trim(adjustl(nemo_send_var(ji))) .eq. trim(adjustl(ssnd(jx)%clname)) ) then
             !--- jx is the index in ssnd for this name
             do jc=1,ssnd(jx)%nct
               ssnd(jx)%nid(jc,1) = cpl_vinfo%tag
             enddo
           endif
         enddo
         if ( rank == ocn_master .and. verbose > -1 ) then
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
         do jx=1,nmaxfld
           if ( trim(adjustl(nemo_recv_var(ji))) .eq. trim(adjustl(srcv(jx)%clname)) ) then
             !--- jx is the index in srcv for this name
             do jc=1,srcv(jx)%nct
               srcv(jx)%nid(jc,1) = cpl_vinfo%tag
             enddo
           endif
         enddo
         if ( rank == ocn_master .and. verbose > -1 ) then
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

   SUBROUTINE hgr_read_glo( plamt , plamu , plamv  , plamf  ,   &    ! gridpoints position (required)
      &                     pphit , pphiu , pphiv  , pphif  ,   &     
      &                     pe1t  , pe1u  , pe1v   , pe1f   ,   &    ! scale factors       (required)
      &                     pe2t  , pe2u  , pe2v   , pe2f   )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE hgr_read  ***
      !!
      !! ** Purpose :   Read a mesh_mask file in NetCDF format using IOM
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   plamt, plamu, plamv, plamf   ! longitude outputs 
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pphit, pphiu, pphiv, pphif   ! latitude outputs
      REAL(wp), DIMENSION(:,:), INTENT(out)  :: pe1v! i-scale factors
      REAL(wp), DIMENSION(:,:), INTENT(out)  :: pe1t, pe1u, pe1f! i-scale factors
      REAL(wp), DIMENSION(:,:), INTENT(out)  :: pe2u! j-scale factors
      REAL(wp), DIMENSION(:,:), INTENT(out)  :: pe2t, pe2v, pe2f! j-scale factors
      !
      INTEGER  ::   inum                  ! logical unit
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   hgr_read_glo : read the global horizontal coordinates in mesh_mask (cpl_cancpl_define)'
         WRITE(numout,*) '   ~~~~~~~~      jpiglo = ', jpiglo, ' jpjglo = ', jpjglo, ' jpk = ', jpk
      ENDIF
      !
      CALL iom_open( cn_domcfg, inum )
      !
      CALL iom_get( inum, jpdom_unknown, 'glamt', plamt(1:Ni0glo,1:Nj0glo), cd_type = 'T', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_unknown, 'glamu', plamu(1:Ni0glo,1:Nj0glo), cd_type = 'U', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_unknown, 'glamv', plamv(1:Ni0glo,1:Nj0glo), cd_type = 'V', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_unknown, 'glamf', plamf(1:Ni0glo,1:Nj0glo), cd_type = 'F', psgn = 1._wp, kfill = jpfillcopy )
      !
      CALL iom_get( inum, jpdom_unknown, 'gphit', pphit(1:Ni0glo,1:Nj0glo), cd_type = 'T', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_unknown, 'gphiu', pphiu(1:Ni0glo,1:Nj0glo), cd_type = 'U', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_unknown, 'gphiv', pphiv(1:Ni0glo,1:Nj0glo), cd_type = 'V', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_unknown, 'gphif', pphif(1:Ni0glo,1:Nj0glo), cd_type = 'F', psgn = 1._wp, kfill = jpfillcopy )
      !
      CALL iom_get( inum, jpdom_unknown, 'e1t'  , pe1t(1:Ni0glo,1:Nj0glo) , cd_type = 'T', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_unknown, 'e1u'  , pe1u(1:Ni0glo,1:Nj0glo) , cd_type = 'U', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_unknown, 'e1v'  , pe1v(1:Ni0glo,1:Nj0glo) , cd_type = 'V', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_unknown, 'e1f'  , pe1f(1:Ni0glo,1:Nj0glo) , cd_type = 'F', psgn = 1._wp, kfill = jpfillcopy )
      !
      CALL iom_get( inum, jpdom_unknown, 'e2t'  , pe2t(1:Ni0glo,1:Nj0glo) , cd_type = 'T', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_unknown, 'e2u'  , pe2u(1:Ni0glo,1:Nj0glo) , cd_type = 'U', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_unknown, 'e2v'  , pe2v(1:Ni0glo,1:Nj0glo) , cd_type = 'V', psgn = 1._wp, kfill = jpfillcopy )
      CALL iom_get( inum, jpdom_unknown, 'e2f'  , pe2f(1:Ni0glo,1:Nj0glo) , cd_type = 'F', psgn = 1._wp, kfill = jpfillcopy )
      !
      CALL iom_close( inum )
      !
   END SUBROUTINE hgr_read_glo

  subroutine cpl_gather(vname, rank)
    character(*), intent(in) :: vname
    integer(kind=impi), intent(in) :: rank
    real(wp), dimension(:, :), POINTER :: work
    integer inum

    IF (.NOT. ASSOCIATED(work)) ALLOCATE(work(Ni0glo,Nj0glo))
 
    ! initialize work array with 0. of the fields in cn_domcfg if present
    work=0.
    IF( ln_read_cfg ) then
        CALL iom_open( cn_domcfg, inum )
        if (iom_varid( inum, vname,ldstop=.false. ) > 0) CALL iom_get( inum, jpdom_unknown,vname, work)
        CALL iom_close(inum)
    ENDIF
    !--- Gather the variable vname into work, and then
    !       store it in global array named nemo_vname. work will contain the
    !       Nemo northfold, but we will avoid copying this into the nemo_* arrays
    !--- This data is then sent to the coupler in cpl_initialize_events
    select case (trim(adjustl(vname)))
      case ("glamt")
        call reconstruct_global_2d(glamt,0,work)
        IF (.NOT. ASSOCIATED(nemo_glamt)) ALLOCATE(nemo_glamt(nemo_jpiglo,nemo_jpjglo))
        nemo_glamt(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("glamu")
        call reconstruct_global_2d(glamu,0,work)
        IF (.NOT. ASSOCIATED(nemo_glamu)) ALLOCATE(nemo_glamu(nemo_jpiglo,nemo_jpjglo))
        nemo_glamu(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("glamv")
        call reconstruct_global_2d(glamv,0,work)
        IF (.NOT. ASSOCIATED(nemo_glamv)) ALLOCATE(nemo_glamv(nemo_jpiglo,nemo_jpjglo))
        nemo_glamv(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("glamf")
        call reconstruct_global_2d(glamf,0,work)
        IF (.NOT. ASSOCIATED(nemo_glamf)) ALLOCATE(nemo_glamf(nemo_jpiglo,nemo_jpjglo))
        nemo_glamf(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("gphit")
        call reconstruct_global_2d(gphit,0,work)
        IF (.NOT. ASSOCIATED(nemo_gphit)) ALLOCATE(nemo_gphit(nemo_jpiglo,nemo_jpjglo))
        nemo_gphit(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("gphiu")
        call reconstruct_global_2d(gphiu,0,work)
        IF (.NOT. ASSOCIATED(nemo_gphiu)) ALLOCATE(nemo_gphiu(nemo_jpiglo,nemo_jpjglo))
        nemo_gphiu(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("gphiv")
        call reconstruct_global_2d(gphiv,0,work)
        IF (.NOT. ASSOCIATED(nemo_gphiv)) ALLOCATE(nemo_gphiv(nemo_jpiglo,nemo_jpjglo))
        nemo_gphiv(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("gphif")
        call reconstruct_global_2d(gphif,0,work)
        IF (.NOT. ASSOCIATED(nemo_gphif)) ALLOCATE(nemo_gphif(nemo_jpiglo,nemo_jpjglo))
        nemo_gphif(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("e1t")
        call reconstruct_global_2d(e1t,0,work)
        IF (.NOT. ASSOCIATED(nemo_e1t)) ALLOCATE(nemo_e1t(nemo_jpiglo,nemo_jpjglo))
        nemo_e1t(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("e1u")
        call reconstruct_global_2d(e1u,0,work)
        IF (.NOT. ASSOCIATED(nemo_e1u)) ALLOCATE(nemo_e1u(nemo_jpiglo,nemo_jpjglo))
        nemo_e1u(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("e1v")
        call reconstruct_global_2d(e1v,0,work)
        IF (.NOT. ASSOCIATED(nemo_e1v)) ALLOCATE(nemo_e1v(nemo_jpiglo,nemo_jpjglo))
        nemo_e1v(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("e1f")
        call reconstruct_global_2d(e1f,0,work)
        IF (.NOT. ASSOCIATED(nemo_e1f)) ALLOCATE(nemo_e1f(nemo_jpiglo,nemo_jpjglo))
        nemo_e1f(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("e2t")
        call reconstruct_global_2d(e2t,0,work)
        IF (.NOT. ASSOCIATED(nemo_e2t)) ALLOCATE(nemo_e2t(nemo_jpiglo,nemo_jpjglo))
        nemo_e2t(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("e2u")
        call reconstruct_global_2d(e2u,0,work)
        IF (.NOT. ASSOCIATED(nemo_e2u)) ALLOCATE(nemo_e2u(nemo_jpiglo,nemo_jpjglo))
        nemo_e2u(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("e2v")
        call reconstruct_global_2d(e2v,0,work)
        IF (.NOT. ASSOCIATED(nemo_e2v)) ALLOCATE(nemo_e2v(nemo_jpiglo,nemo_jpjglo))
        nemo_e2v(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("e2f")
        call reconstruct_global_2d(e2f,0,work)
        IF (.NOT. ASSOCIATED(nemo_e2f)) ALLOCATE(nemo_e2f(nemo_jpiglo,nemo_jpjglo))
        nemo_e2f(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("tmask_i")
        call reconstruct_global_2d(tmask_i,0,work)
        IF (.NOT. ASSOCIATED(nemo_tmask)) ALLOCATE(nemo_tmask(nemo_jpiglo,nemo_jpjglo))
        nemo_tmask(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("tmask")
        call reconstruct_global_2d(tmask(:,:,1),0,work)
        IF (.NOT. ASSOCIATED(nemo_tmask)) ALLOCATE(nemo_tmask(nemo_jpiglo,nemo_jpjglo))
        nemo_tmask(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("umask")
        call reconstruct_global_2d(umask(:,:,1),0,work)
        IF (.NOT. ASSOCIATED(nemo_umask)) ALLOCATE(nemo_umask(nemo_jpiglo,nemo_jpjglo))
        nemo_umask(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("vmask")
        call reconstruct_global_2d(vmask(:,:,1),0,work)
        IF (.NOT. ASSOCIATED(nemo_vmask)) ALLOCATE(nemo_vmask(nemo_jpiglo,nemo_jpjglo))
        nemo_vmask(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case ("fmask")
        call reconstruct_global_2d(fmask(:,:,1),0,work)
        IF (.NOT. ASSOCIATED(nemo_fmask)) ALLOCATE(nemo_fmask(nemo_jpiglo,nemo_jpjglo))
        nemo_fmask(1:nemo_jpiglo,1:nemo_jpjglo) = work(1:nemo_jpiglo,1:nemo_jpjglo)

      case default
        write(6,*)"cpl_gather: Invalid variable name ",trim(vname)
        call flush(6)

    end select
    deallocate(work)
  end subroutine cpl_gather

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
     real(wp), intent(in   )  :: pdata(:,:,:)

     !--- Integer flag to indicate if srcv(kid) was sent or not
     !--- kinfo = OASIS_idle means the field was not sent to the coupler
     !--- kinfo = OASIS_snd  means the field was sent to the coupler
     integer,  intent(out) :: kinfo

     !-- Local
     integer :: jc,ij,jj
     integer :: ldbg=1
     integer :: freq
     integer(kind=impi) :: rank, ierr
     integer :: idx, nwrds
     integer :: verbose=1
     real(wp), dimension(Ni0glo,Nj0glo) :: global_array
     !!--------------------------------------------------------------------

     !--- Determine the rank of the calling process in model_communicator
     call mpi_comm_rank ( model_communicator, rank, ierr )

     !--- Loop over all "catagories" (normally only 1) for this variable
     !--- sending separate data to the coupler for each catagory
     do jc = 1, ssnd(kid)%nct
       !--- The MPI tag associated with this transfer is ssnd(kid)%nid(jc)

       kinfo = OASIS_idle

       !--- Ensure that this variable was configured
       if ( ssnd(kid)%nid(jc,1) <= 0 ) then
         write(numout,*)'cpl_cancpl_snd: ',trim(ssnd(kid)%clname), &
           ' is not sent in this configuration.  kid = ',kid
         call ctl_stop("STOP", "cpl_cancpl_snd", &
                       "Invalid send variable "//trim(ssnd(kid)%clname))
       endif

       !--- Determine if this variable should be coupled now
       cpl_vinfo = find_cpl_vinfo( tag=ssnd(kid)%nid(jc,1) )
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

       if ( sn_cfctl%l_prtctl .and. verbose > 1 ) then
         !--- Write info for each sub-domain to the ocean output file
         write(numout,*) '****************'
         write(numout,*) 'cpl_cancpl_snd: Outgoing ', ssnd(kid)%clname
         write(numout,*) 'cpl_cancpl_snd:      tag ', ssnd(kid)%nid(jc,1)
         write(numout,*) 'cpl_cancpl_snd:    kstep ', kstep
         write(numout,*) 'cpl_cancpl_snd: mpi task ', rank
         write(numout,*) '      - minimum value is ', minval(pdata(:,:,jc))
         write(numout,*) '      - maximum value is ', maxval(pdata(:,:,jc))
         write(numout,*) '      -     sum value is ', sum(pdata(:,:,jc))
         write(numout,*) '****************'
       endif

       IF( ln_timing )   call timing_start('cancpl_snd_gather')
       !--- Gather data into the global array global_array
       global_array=0. 
       call reconstruct_global_2d(pdata(:,:,jc),0,global_array)
       IF( ln_timing )   call timing_stop('cancpl_snd_gather')

       !--- Skip the rest of this loop unless this is the master task
       if ( rank /= ocn_master ) cycle

       if ( verbose > 2 ) then
         !--- Count the number of NaNs in the global global_array array
         idx = count( global_array /= global_array )
         write(numout,*)'cpl_cancpl_snd: ',trim(ssnd(kid)%clname),'  Nans in global_array = ',idx
         write(numout,*)'cpl_cancpl_snd: ',trim(ssnd(kid)%clname),'  min,max,avg = ', &
             minval(global_array),maxval(global_array),sum(global_array)/real(size(global_array),kind=8)
         call flush(numout)
       endif

       !--- Map the the global 3D array global_array onto the 1D wrk array
       wrk(1:Ni0glo*Nj0glo) = reshape(global_array,[Ni0glo*Nj0glo])

       if ( verbose > 2 ) then
         !--- Count the number of NaNs in the wrk array
         nwrds = Ni0glo*Nj0glo
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


  subroutine cpl_cancpl_rcv( kid, kstep, pdata, pmask, kinfo )
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

     REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pmask     ! coupling mask

     !--- Integer flag to indicate if srcv(kid) was recieved or not
     !--- kinfo = OASIS_idle means the field was not received from the coupler
     !--- kinfo = OASIS_rcv  means the field was received from the coupler
     integer,  intent(  out) :: kinfo

     !--- Local
     integer :: jc
     integer :: ldbg=1
     integer :: freq, idx, nwrds
     integer(kind=impi) :: rank, ierr, sz, tag
     integer :: verbose=1
     integer (kind=impi) :: status(MPI_status_size)
     type(FLD_CPL), pointer :: cpl_ptr
     real, dimension(Ni0glo,Nj0glo) :: global_array
     !!--------------------------------------------------------------------

     !---Determine the rank of the calling process in model_communicator
     call mpi_comm_rank ( model_communicator, rank, ierr )
     cpl_ptr => srcv(kid)
     !--- Loop over all "catagories" (normally only 1) for this variable
     !--- receiving separate data from the coupler for each catagory

     do jc = 1, srcv(kid)%nct
       !--- The MPI tag associated with this transfer is srcv(kid)%nid(jc)

       kinfo = OASIS_idle

       ! The following is no longer true. srcv(kid)%laction determines whether
       ! it should be received
      !  !--- This routine is called for every variable that could be coupled
      !  !--- so ignore variables that are not to be coupled
      !  !--- srcv(:)%nid(:) is initialized to zero then defined in
      !  !--- cpl_cancpl_define for variables that are to be coupled and so
      !  !--- it will only be non-zero for fields that are coupled
      !  if ( srcv(kid)%nid(jc,1) <= 0 ) cycle

       !--- Determine if this variable should be coupled now
       cpl_vinfo = find_cpl_vinfo( tag=srcv(kid)%nid(jc,1) )
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
             '  tag=',srcv(kid)%nid(jc,1),'  kstep=',kstep,'  freq=',freq
           call flush(numout)
         endif

         !--- Receive the global array from the coupler
         call recv_data_rec(wrk, ibuf, cpl_master, trim(srcv(kid)%clname), dbg=ldbg)

       endif

       !--- Scatter the global array onto each NEMO task
       IF( ln_timing )   call timing_start('cancpl_rcv_scatter')
       global_array = RESHAPE(wrk,[Ni0glo,Nj0glo])
       call mppsync
       call mppscatter(global_array, 0, pdata(:,:,jc))
       call mppsync
       IF( ln_timing )   call timing_stop('cancpl_rcv_scatter')

       if ( rank == ocn_master .and. verbose > 2 ) then
         !--- Count the number of NaNs in the global_array
         idx = count( global_array /= global_array )
         write(numout,*)'cpl_cancpl_rcv: ',trim(srcv(kid)%clname),'  Before lbc_lnk'
         write(numout,*)'cpl_cancpl_rcv: ',trim(srcv(kid)%clname),'  Nans in global_array = ',idx
         write(numout,*)'cpl_cancpl_rcv: ',trim(srcv(kid)%clname),'  min,max,avg = ', &
             minval(global_array),maxval(global_array),sum(global_array)/real(size(global_array),kind=8)
         call flush(numout)
       endif

       !--- Fill overlap areas and extra hallows and check periodicity
       call lbc_lnk( 'cpl_cancpl_rcv', pdata(:,:,jc), srcv(kid)%clgrid, srcv(kid)%nsgn )

       if ( rank == ocn_master .and. verbose > 2 ) then
         !--- Count the number of NaNs in the global global_array array
         idx = count( global_array /= global_array )
           write(numout,*)'cpl_cancpl_rcv: ',trim(srcv(kid)%clname),'  After lbc_lnk'
           write(numout,*)'cpl_cancpl_rcv: ',trim(srcv(kid)%clname),'  Nans in global_array = ',idx
           write(numout,*)'cpl_cancpl_rcv: ',trim(srcv(kid)%clname),'  min,max,avg = ', &
               minval(global_array),maxval(global_array),sum(global_array)/real(size(global_array),kind=8)
           call flush(numout)
       endif

       if ( sn_cfctl%l_prtctl .and. verbose > 1 ) then
         !--- Write info for each sub-domain to the ocean output file
         write(numout,*) '****************'
         write(numout,*) 'cpl_cancpl_rcv: Incoming ', srcv(kid)%clname
         write(numout,*) 'cpl_cancpl_rcv: ivarid ',   srcv(kid)%nid(jc,1)
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

    call mppstop
    !--- TODO --- Also tell coupler that the ocean has stopped

  end subroutine cpl_cancpl_finalize

  !> Propagate parameters from NEMO to CanCPL
  subroutine set_cancpl_params( coupled_fields, num_ice_steps, num_boundary_calls )
    type(FLD_C), dimension(:), intent(in) :: coupled_fields !< Contains all the field descriptors for coupled fields
    integer,                   intent(in) :: num_ice_steps  !< How many dynamics timesteps per ice timestep
    integer,                   intent(in) :: num_boundary_calls !< How many dynamics timesteps per surface boundary call

    integer :: i

    ! Set the coupling strategy
    !! Note that for now this retains the need for a static mapping between the incoming fields and what cancpl expects
    nemo_namsbc_cpl_cldes(:) = " "
    do i = 1,SIZE(coupled_fields,1)
      nemo_namsbc_cpl_cldes(i) = trim(coupled_fields(i)%cldes)
    enddo

    ! Set other
    nn_ice = num_ice_steps
    nn_fsbc = num_boundary_calls

  end subroutine set_cancpl_params

  !> Subroutine to detect when the ocean should get ready to receive fields from teh coupler
  subroutine query_start_cpl2ocn( isec )
    integer, intent(in) :: isec

    type(cpl_vinfo_t) :: cpl_vinfo

      cpl_vinfo = find_cpl_vinfo( name="start_cpl2ocn" )
      if ( mod(isec,cpl_vinfo%freq) == 0 ) then
        !--- Receive cpl_time_string from the coupler
        !--- cpl_time_string is found in the com_cpl module
        call bcast_inter(cpl_time_string, cpl_master, "ocn")
        !--- Send elapsed coupler time in seconds from the coupler to the ocean
        !--- cpl_elapsed_time_secs is found in the com_cpl module
        call bcast_inter(cpl_elapsed_time_secs, cpl_master, "ocn")
        write(numout,*)"isec=",isec,"  cpl_elapsed_time_secs=",cpl_elapsed_time_secs
        call flush(numout)
      endif
  end subroutine query_start_cpl2ocn

#else


  use par_kind, only : wp
  use cpl_types, only : FLD_C

  public :: cpl_cancpl_init
  public :: cpl_cancpl_define
  public :: cpl_cancpl_snd
  public :: cpl_cancpl_rcv
  public :: cpl_cancpl_freq
  public :: cpl_cancpl_finalize
  public :: set_cancpl_params
  public :: query_start_cpl2ocn

contains

  subroutine cpl_cancpl_init( kl_comm )
     !!-------------------------------------------------------------------
     !!             ***  ROUTINE cpl_cancpl_init  ***
     !!
     !!--------------------------------------------------------------------------
     integer, intent(out) ::   kl_comm   ! MPI group ID for the ocean comm_world
     kl_comm = 0
     WRITE(numout,*) 'cpl_cancpl_init: Error you sould not be there...'
  end subroutine cpl_cancpl_init
  subroutine cpl_cancpl_define( krcv, ksnd, kcplmodel )
     !!-------------------------------------------------------------------
     !!             ***  ROUTINE cpl_cancpl_define  ***
     !!
     !!--------------------------------------------------------------------
     integer, intent(in) :: krcv   ! Number of all possible fields received
     integer, intent(in) :: ksnd   ! Number of all possible fields sent
     integer, intent(in) :: kcplmodel ! Number of models to send too. Note this is a dummy argument for now
     WRITE(numout,*) 'cpl_cancpl_define: Error you sould not be there...'
  end subroutine cpl_cancpl_define
  subroutine cpl_cancpl_snd( kid, kstep, pdata, kinfo )
     !!---------------------------------------------------------------------
     !!              ***  ROUTINE cpl_cancpl_snd  ***
     !!
     !!----------------------------------------------------------------------
     !--- Index in the array ssnd of the variable to be sent to the coupler
     integer,  intent(in)  :: kid
     integer,  intent(in)  :: kstep
     real(wp), intent(in   )  :: pdata(:,:,:)
     integer,  intent(out) :: kinfo
     kinfo=0 
     WRITE(numout,*) 'cpl_cancpl_snd: Error you sould not be there...'
  end subroutine cpl_cancpl_snd
  subroutine cpl_cancpl_rcv( kid, kstep, pdata, pmask, kinfo )
     !!---------------------------------------------------------------------
     !!              ***  ROUTINE cpl_cancpl_rcv  ***
     !!
     !!----------------------------------------------------------------------
     integer,  intent(in   ) :: kid
     integer,  intent(in   ) :: kstep
     real(wp), intent(inout) :: pdata(:,:,:)
     REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pmask     ! coupling mask
     integer,  intent(  out) :: kinfo
     kinfo=0; pdtata=0. 
     WRITE(numout,*) 'cpl_cancpl_rcv: Error you sould not be there...'
  end subroutine cpl_cancpl_rcv
  function cpl_cancpl_freq( vname ) result(freq)
    !!---------------------------------------------------------------------
    !!              ***  ROUTINE cpl_cancpl_freq  ***
    !!
    !!----------------------------------------------------------------------
    character(*), intent(in) :: vname
    integer :: freq
    freq=0
     WRITE(numout,*) 'cpl_cancpl_freq: Error you sould not be there...'
  end function cpl_cancpl_freq
  subroutine cpl_cancpl_finalize
    !!---------------------------------------------------------------------
    !!              ***  ROUTINE cpl_cancpl_finalize  ***
    !!
    !!----------------------------------------------------------------------
     WRITE(numout,*) 'cpl_cancpl_finalize: Error you sould not be there...'
  end subroutine cpl_cancpl_finalize
  subroutine set_cancpl_params( coupled_fields, num_ice_steps, num_boundary_calls )
    type(FLD_C), dimension(:), intent(in) :: coupled_fields !< Contains all the field descriptors for coupled fields
    integer,                   intent(in) :: num_ice_steps  !< How many dynamics timesteps per ice timestep
    integer,                   intent(in) :: num_boundary_calls !< How many dynamics timesteps per surface boundary call
    WRITE(numout,*) 'set_cancpl_params: Error you sould not be there...'
  end subroutine set_cancpl_params
  subroutine query_start_cpl2ocn( isec )
    integer, intent(in) :: isec
    WRITE(numout,*) 'query_start_cpl2ocn: Error you sould not be there...'
  end subroutine query_start_cpl2ocn

#endif

end module cpl_cancpl
