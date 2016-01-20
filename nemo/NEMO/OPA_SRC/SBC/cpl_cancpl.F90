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

  !--- The mpi communicator assigned to the ocean
  !--- This is defined in cpl_cancpl_init and so use of ocn_comm
  !--- must occur after the call to cpl_cancpl_init
  integer(kind=impi), public :: ocn_comm

  logical, public, parameter ::   lk_cpl = .true.   !: coupled flag
  integer, public            ::   oasis_idle = 0    !: return code if no send or recv
  integer, public            ::   oasis_rcv  = 1    !: return code if field received
  integer, public            ::   oasis_snd  = 2    !: return code if field sent

  !--- These are defined in com_cpl
  public :: cpl_vinfo_t, find_cpl_vinfo
  public :: bcast_inter, cpl_time_string, cpl_time_secs, cpl_master

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

  !--- tmp space for use with MPI gather/scatter operations
  real(wp), allocatable, save, dimension(:,:,:), private :: png

  !--- tmp char space
  character(512) :: strng

  !--- Work space used to temporarily hold fields passed via MPI
  integer, parameter :: maxx=8392704  !---4098x2048 = (2+2^12)x(2^11)
  real(kind=8)        :: wrk(maxx)
  integer(kind=8)     :: ibuf(8)

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

     !!============================================
     !! WARNING: No write in numout in this routine
     !!============================================

     !--- Initialize groups for cpl, atm, ocn, ice, ...
     !--- This will, among other things, define cpl_master, atm_master, ocn_master
     !--- NOTE: ocn_master is the rank in MPI_COMM_WORLD not the rank in ocn_comm
     call define_group('ocn', ocn_comm)

     kl_comm = ocn_comm

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
     !
     integer :: ji,jc          ! local loop indices
     character(len=8) :: zclname
     integer :: ldbg=1
     integer(kind=impi) :: rank, ierr
     integer :: verbose=1
     character(32), allocatable :: var_list(:)
     type var_order_t
       character(32) :: name
       integer       :: index
     end type var_order_t
     type(var_order_t), pointer, dimension(:) :: var_order => null()
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
       !--- TODO --- Also tell coupler that there is a problem
       call ctl_stop("STOP", " cpl_cancpl_define", "Problem allocating png")
     endif

     ! -----------------------------------------------------------------
     ! ... Assign MPI tags to ssnd and srcv variables
     ! -----------------------------------------------------------------
     !
     !--- Send variables
     !
     !--- nemo_n_send_var is the number of fields in nemo_send_var
     !--- These are available from the com_cpl module
     nemo_n_send_var=0
     nemo_send_var(:) = " "
     do ji = 1, ksnd
        if ( ssnd(ji)%laction ) then 
           do jc = 1, ssnd(ji)%nct
              if ( ssnd(ji)%nct .gt. 1 ) then
                 write(zclname,'( a7, i1)') ssnd(ji)%clname,jc
              else
                 zclname=ssnd(ji)%clname
              endif

              !--- Assign the mpi tag associated with this variable to ssnd
              cpl_vinfo = find_cpl_vinfo( name=trim(zclname) )
              ssnd(ji)%nid(jc) = cpl_vinfo%tag

              nemo_n_send_var = nemo_n_send_var + 1
              if ( nemo_n_send_var > nemo_send_var_max ) then
                write(numout,*) "cpl_cancpl_define: Too many send variables"
                write(6,*) "cpl_cancpl_define: Too many send variables"
                call flush(6)
                call ctl_stop("STOP", " cpl_cancpl_define", "Too many send variables")
              endif
              nemo_send_var(nemo_n_send_var) = trim(zclname)

              if ( rank == ocn_master ) then
                !--- Write to NEMO's ocean.output file
                write(numout,*) "cpl_cancpl_define: Send field ",ji, &
                    "  name=",trim(zclname)," tag=",cpl_vinfo%tag
                call flush(numout)

                !--- Also write to stdout (unit 6)
                write(6,*) "cpl_cancpl_define: Send field ",ji, &
                    "  name=",trim(zclname)," tag=",cpl_vinfo%tag
                call flush(6)
              endif
           end do
        endif
     end do
     !
     !--- Receive variables
     !
     !--- nemo_n_recv_var is the number of fields in nemo_recv_var
     !--- These are available from the com_cpl module
     nemo_n_recv_var=0
     nemo_recv_var(:) = " "
     do ji = 1, krcv
        if ( srcv(ji)%laction ) then 
           do jc = 1, srcv(ji)%nct
              if ( srcv(ji)%nct .gt. 1 ) then
                 write(zclname,'( a7, i1)') srcv(ji)%clname,jc
              else
                 zclname=srcv(ji)%clname
              endif

              !--- Assign the mpi tag associated with this variable to srcv
              cpl_vinfo = find_cpl_vinfo( name=trim(zclname) )
              srcv(ji)%nid(jc) = cpl_vinfo%tag

              nemo_n_recv_var = nemo_n_recv_var + 1
              if ( nemo_n_recv_var > nemo_recv_var_max ) then
                write(numout,*) "cpl_cancpl_define: Too many receive variables"
                write(6,*) "cpl_cancpl_define: Too many receive variables"
                call flush(6)
                call ctl_stop("STOP", " cpl_cancpl_define", "Too many receive variables")
              endif
              nemo_recv_var(nemo_n_recv_var) = trim(zclname)

              if ( rank == ocn_master ) then
                !--- Write to NEMO's ocean.output file
                write(numout,*) "cpl_cancpl_define: Recv field ",ji, &
                    "  name=",trim(zclname)," tag=",cpl_vinfo%tag
                call flush(numout)

                !--- Also write to stdout (unit 6)
                write(6,*) "cpl_cancpl_define: Recv field ",ji, &
                    "  name=",trim(zclname)," tag=",cpl_vinfo%tag
                call flush(6)
              endif
           end do
        endif
     end do

!---TODO--- nemo_send_var and nemo_recv_var need to be reordered
!---TODO--- The fields in these lists must be in the same order as the the data
!---TODO--- is sent or received.
!---TODO--- See subroutine sbc_cpl_snd for the send order
!---TODO--- nemo_recv_var should already be in the correct order which is
!---TODO--- the order the fields appear in srcv
     if ( nemo_n_send_var > 0 ) then
       !--- nemo_send_var needs to be reordered
       if (allocated(var_list) ) deallocate(var_list)
       allocate( var_list(nemo_n_send_var) )
       !--- The subroutine sbc_cpl_snd determines the send order
       !--- Since there is no clear way to determine this order on the fly
       !--- it must be hard coded here (this is bad). Therefore if there
       !--- are any changes in sbc_cpl_snd that alter this order then there
       !--- must also be changes here
       do ji=1,nemo_n_send_var
       enddo
     endif

     if ( rank == ocn_master .and. verbose > 0 ) then
       write(6,*)"cpl_cancpl_define: nemo_n_send_var=",nemo_n_send_var
       write(6,'(5(2x,a))')nemo_send_var(1:nemo_n_send_var)
       write(6,*)"cpl_cancpl_define: nemo_n_recv_var=",nemo_n_recv_var
       write(6,'(5(2x,a))')nemo_recv_var(1:nemo_n_recv_var)
       call flush(6)
     endif

     !---DBG--- !--- Verify the relationship between nproc, narea and rank
     !---DBG--- call mppsync
     !---DBG--- write(numout,*) "cpl_cancpl_define: nproc,narea,rank: ",nproc,narea,rank
     !---DBG--- call flush(numout)
     !---DBG--- write(6,*) "cpl_cancpl_define: nproc,narea,rank: ",nproc,narea,rank
     !---DBG--- call flush(6)

     !--- These values will be broadcast to all mpi tasks
     !--- in the following call to cpl_initialize_events

     !--- Set a value for nemo_nn_ice, defined in com_cpl
     nemo_nn_ice = nn_ice

     !--- Set a value for nemo_rn_rdt, defined in com_cpl
     !--- rn_rdt is defined in the module dom_oce
     nemo_rn_rdt = nint(rn_rdt,8)

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

     !--- Initialize coupler events
     call cpl_initialize_events()

     !--- Broadcast the initial date and time from the coupler to all tasks
     !--- cpl_time_string is defined in com_cpl
     call bcastGroup(cpl_time_string, cpl_master, MPI_COMM_WORLD)

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
    integer :: ji, jj, jn, ji_glob, jj_glob

    glob_arr = 0.0_8
    glob_arr(1:jpiglo,1:jpjglo) = reshape( wrk(1:jpiglo*jpjglo), (/ jpiglo,jpjglo /) )

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
            write(6,*)'copy_1d_to_3d_global: Global index is out of range.'
            write(6,*)'jn, ji, jj, ji_glob, jj_glob: ',jn, ji, jj, ji_glob, jj_glob
            call ctl_stop("STOP", "copy_1d_to_3d_global", "Global index is out of range")
          endif
          png(ji,jj,jn) = glob_arr(ji_glob,jj_glob)
        enddo
      enddo
    enddo

  end subroutine copy_1d_to_3d_global


  subroutine copy_3d_to_1d_global(wrk, png)
    !------------------------------------------------------------------------
    !--- Copy values from a global 3D array containing data from a recent
    !--- call to mppgather to a global 1D wrk array to be sent to the agcm
    !------------------------------------------------------------------------
    real(kind=8), intent(out) :: wrk(:)
    real(wp), intent(in) :: png(jpi,jpj,jpnij)

    !--- Local
    real(kind=8) :: glob_arr(jpiglo, jpjglo)
    integer :: ji, jj, jn, ji_glob, jj_glob

    glob_arr = 0.0_8
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
            write(6,*)'copy_3d_to_1d_global: Global index is out of range.'
            write(6,*)'jn, ji, jj, ji_glob, jj_glob: ',jn, ji, jj, ji_glob, jj_glob
            call ctl_stop("STOP", "copy_3d_to_1d_global", "Global index is out of range")
          endif
          glob_arr(ji_glob,jj_glob) = png(ji,jj,jn)
        enddo
      enddo
    enddo

    wrk = 0.0_8
    wrk(1:jpiglo*jpjglo) = reshape( glob_arr(1:jpiglo,1:jpjglo), (/ jpiglo*jpjglo /) )

  end subroutine copy_3d_to_1d_global


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
     integer :: verbose=1
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

       if ( rank == ocn_master .and. verbose > -2 ) then
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
     integer :: verbose=1
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
