module ncdf
  !***************************************************************
  ! Common netcdf related routines
  !***************************************************************

  use netcdf

  implicit none

  public

  interface ncdf_write_var
    module procedure ncdf_write_var_1d_i4
    module procedure ncdf_write_var_1d_i8
    module procedure ncdf_write_var_1d_r4
    module procedure ncdf_write_var_1d_r8
    module procedure ncdf_write_var_2d_i4
    module procedure ncdf_write_var_2d_i8
    module procedure ncdf_write_var_2d_r4
    module procedure ncdf_write_var_2d_r8
    module procedure ncdf_write_var_3d_i4
    module procedure ncdf_write_var_3d_i8
    module procedure ncdf_write_var_3d_r4
    module procedure ncdf_write_var_3d_r8
    module procedure ncdf_write_var_4d_i4
    module procedure ncdf_write_var_4d_i8
    module procedure ncdf_write_var_4d_r4
    module procedure ncdf_write_var_4d_r8
    module procedure ncdf_write_var_1d_char
  end interface

  private :: ncdf_write_var_any
  private :: ncdf_write_var_1d_i4
  private :: ncdf_write_var_1d_i8
  private :: ncdf_write_var_1d_r4
  private :: ncdf_write_var_1d_r8
  private :: ncdf_write_var_2d_i4
  private :: ncdf_write_var_2d_i8
  private :: ncdf_write_var_2d_r4
  private :: ncdf_write_var_2d_r8
  private :: ncdf_write_var_3d_i4
  private :: ncdf_write_var_3d_i8
  private :: ncdf_write_var_3d_r4
  private :: ncdf_write_var_3d_r8
  private :: ncdf_write_var_4d_i4
  private :: ncdf_write_var_4d_i8
  private :: ncdf_write_var_4d_r4
  private :: ncdf_write_var_4d_r8
  private :: ncdf_write_var_1d_char

  interface ncdf_read_var
    module procedure ncdf_read_var_1d_i4
    module procedure ncdf_read_var_1d_i8
    module procedure ncdf_read_var_2d_i4
    module procedure ncdf_read_var_2d_i8
    module procedure ncdf_read_var_3d_i4
    module procedure ncdf_read_var_3d_i8
    module procedure ncdf_read_var_4d_i4
    module procedure ncdf_read_var_4d_i8

    module procedure ncdf_read_var_1d_r4
    module procedure ncdf_read_var_1d_r8
    module procedure ncdf_read_var_2d_r4
    module procedure ncdf_read_var_2d_r8
    module procedure ncdf_read_var_3d_r4
    module procedure ncdf_read_var_3d_r8
    module procedure ncdf_read_var_4d_r4
    module procedure ncdf_read_var_4d_r8

    module procedure ncdf_read_var_byname_1d_i4
    module procedure ncdf_read_var_byname_1d_i8
    module procedure ncdf_read_var_byname_2d_i4
    module procedure ncdf_read_var_byname_2d_i8
    module procedure ncdf_read_var_byname_3d_i4
    module procedure ncdf_read_var_byname_3d_i8
    module procedure ncdf_read_var_byname_4d_i4
    module procedure ncdf_read_var_byname_4d_i8

    module procedure ncdf_read_var_byname_1d_r4
    module procedure ncdf_read_var_byname_1d_r8
    module procedure ncdf_read_var_byname_2d_r4
    module procedure ncdf_read_var_byname_2d_r8
    module procedure ncdf_read_var_byname_3d_r4
    module procedure ncdf_read_var_byname_3d_r8
    module procedure ncdf_read_var_byname_4d_r4
    module procedure ncdf_read_var_byname_4d_r8
  end interface

  private :: ncdf_read_var_1d_i4
  private :: ncdf_read_var_1d_i8
  private :: ncdf_read_var_2d_i4
  private :: ncdf_read_var_2d_i8
  private :: ncdf_read_var_3d_i4
  private :: ncdf_read_var_3d_i8
  private :: ncdf_read_var_4d_i4
  private :: ncdf_read_var_4d_i8

  private :: ncdf_read_var_1d_r4
  private :: ncdf_read_var_1d_r8
  private :: ncdf_read_var_2d_r4
  private :: ncdf_read_var_2d_r8
  private :: ncdf_read_var_3d_r4
  private :: ncdf_read_var_3d_r8
  private :: ncdf_read_var_4d_r4
  private :: ncdf_read_var_4d_r8

  private :: ncdf_read_var_byname_1d_i4
  private :: ncdf_read_var_byname_1d_i8
  private :: ncdf_read_var_byname_2d_i4
  private :: ncdf_read_var_byname_2d_i8
  private :: ncdf_read_var_byname_3d_i4
  private :: ncdf_read_var_byname_3d_i8
  private :: ncdf_read_var_byname_4d_i4
  private :: ncdf_read_var_byname_4d_i8

  private :: ncdf_read_var_byname_1d_r4
  private :: ncdf_read_var_byname_1d_r8
  private :: ncdf_read_var_byname_2d_r4
  private :: ncdf_read_var_byname_2d_r8
  private :: ncdf_read_var_byname_3d_r4
  private :: ncdf_read_var_byname_3d_r8
  private :: ncdf_read_var_byname_4d_r4
  private :: ncdf_read_var_byname_4d_r8

  interface ncdf_add_coord
    module procedure ncdf_add_coord_none
    module procedure ncdf_add_coord_i4
    module procedure ncdf_add_coord_i8
    module procedure ncdf_add_coord_r4
    module procedure ncdf_add_coord_r8
    module procedure ncdf_add_coord_char
  end interface

  private :: ncdf_add_coord_any
  private :: ncdf_add_coord_none
  private :: ncdf_add_coord_i4
  private :: ncdf_add_coord_i8
  private :: ncdf_add_coord_r4
  private :: ncdf_add_coord_r8
  private :: ncdf_add_coord_char

  interface ncdf_add_var
    module procedure ncdf_add_var_no_vals
    module procedure ncdf_add_var_1d_i4
    module procedure ncdf_add_var_1d_i8
    module procedure ncdf_add_var_1d_r4
    module procedure ncdf_add_var_1d_r8
    module procedure ncdf_add_var_2d_i4
    module procedure ncdf_add_var_2d_i8
    module procedure ncdf_add_var_2d_r4
    module procedure ncdf_add_var_2d_r8
    module procedure ncdf_add_var_3d_i4
    module procedure ncdf_add_var_3d_i8
    module procedure ncdf_add_var_3d_r4
    module procedure ncdf_add_var_3d_r8
    module procedure ncdf_add_var_4d_i4
    module procedure ncdf_add_var_4d_i8
    module procedure ncdf_add_var_4d_r4
    module procedure ncdf_add_var_4d_r8
  end interface

  private :: ncdf_add_var_no_vals
  private :: ncdf_add_var_1d_i4
  private :: ncdf_add_var_1d_i8
  private :: ncdf_add_var_1d_r4
  private :: ncdf_add_var_1d_r8
  private :: ncdf_add_var_2d_i4
  private :: ncdf_add_var_2d_i8
  private :: ncdf_add_var_2d_r4
  private :: ncdf_add_var_2d_r8
  private :: ncdf_add_var_3d_i4
  private :: ncdf_add_var_3d_i8
  private :: ncdf_add_var_3d_r4
  private :: ncdf_add_var_3d_r8
  private :: ncdf_add_var_4d_i4
  private :: ncdf_add_var_4d_i8
  private :: ncdf_add_var_4d_r4
  private :: ncdf_add_var_4d_r8

  interface ncdf_get_att
    module procedure ncdf_get_att_char
    module procedure ncdf_get_att_si4
    module procedure ncdf_get_att_ai4
    module procedure ncdf_get_att_si8
    module procedure ncdf_get_att_ai8
    module procedure ncdf_get_att_sr4
    module procedure ncdf_get_att_ar4
    module procedure ncdf_get_att_sr8
    module procedure ncdf_get_att_ar8
  end interface

  private :: ncdf_get_att_char
  private :: ncdf_get_att_si4
  private :: ncdf_get_att_ai4
  private :: ncdf_get_att_si8
  private :: ncdf_get_att_ai8
  private :: ncdf_get_att_sr4
  private :: ncdf_get_att_ar4
  private :: ncdf_get_att_sr8
  private :: ncdf_get_att_ar8

  interface ncdf_add_att
    module procedure ncdf_add_att_char
    module procedure ncdf_add_att_si4
    module procedure ncdf_add_att_ai4
    module procedure ncdf_add_att_si8
    module procedure ncdf_add_att_ai8
    module procedure ncdf_add_att_sr4
    module procedure ncdf_add_att_ar4
    module procedure ncdf_add_att_sr8
    module procedure ncdf_add_att_ar8
  end interface

  private :: ncdf_add_att_char
  private :: ncdf_add_att_si4
  private :: ncdf_add_att_ai4
  private :: ncdf_add_att_si8
  private :: ncdf_add_att_ai8
  private :: ncdf_add_att_sr4
  private :: ncdf_add_att_ar4
  private :: ncdf_add_att_sr8
  private :: ncdf_add_att_ar8

  interface ncdf_inquire_vars
    module procedure ncdf_inquire_vars_by_file_name
    module procedure ncdf_inquire_vars_by_ID
  end interface

  private :: ncdf_inquire_vars_by_file_name
  private :: ncdf_inquire_vars_by_ID

  interface ncdf_get_time
    module procedure ncdf_get_time_by_file_name
    module procedure ncdf_get_time_by_ID
  end interface

  private :: ncdf_get_time_by_file_name
  private :: ncdf_get_time_by_ID

  type ncdf_file_name_t
    character(128) :: name = " "
    integer        :: fid  = -1
    integer        :: mode = -1
  end type ncdf_file_name_t

  !--- Space for name/fid/mode info about active netcdf files
  type(ncdf_file_name_t), save :: ncdf_file_name(101)

  contains

    !***************************************************************************
    !--- Handle error conditions returned from calls to netcdf routines.
    !***************************************************************************
    subroutine nc_error_handler(fid, status, eid, msg)

      implicit none

      !--- Input
      integer(kind=4) :: fid     !--- netcdf file id
      integer(kind=4) :: status  !--- return status from netcdf routine
      integer :: eid
      character(*), optional :: msg

      !--- Local
      integer(kind=4) :: lstatus

      !--- If there is no error then simply return
      if (status == NF90_NOERR) return

      !--- See what the netcdf libs have to say about this error
      write(6,'(a)') trim( nf90_strerror(status) )

      !--- Close the netcdf file
      lstatus = nf90_close(fid)

      !--- Print the user supplied message, if any
      if ( present(msg) ) write(6,'(a)') trim(msg)

      !--- Abort
      if (eid.gt.0) eid = -eid
      call xit("nc_error", eid)

    end subroutine nc_error_handler

    subroutine strng_rm_nulls(strng)
      !--- Replace any nulls found in the input string with blanks
      implicit none
      character(*), intent(inout) :: strng
      integer :: idx

      do idx=1,len(strng)
        if ( strng(idx:idx) .eq. achar(0) ) strng(idx:idx) = " "
      enddo

    end subroutine strng_rm_nulls

    !***************************************************************************
    !--- Create a new netcdf file or open an existing netcdf file
    !***************************************************************************
    subroutine ncdf_open(fid, fname, mode)
      !--- Return the netcdf file ID
      integer(kind=4), intent(out) :: fid

      !--- The name of the netcdf file to create or open
      character(*), intent(in) :: fname

      !--- File creation/open mode
      character(*), intent(in), optional :: mode

      !--- Local
      integer(kind=4) :: ncstat, idx, imode
      character(32) :: local_mode
      character(64) :: msg

      if ( len_trim(fname) < 1 ) then
        write(6,*)"ncdf_open: Missing file name."
        call xit("ncdf_open",-1)
      endif

      if ( present(mode) ) then
        local_mode = trim(mode)
      else
        !--- Use the safe default
        local_mode = "read"
      endif

      select case(trim(adjustl(local_mode)))
        case ("create", "c")
          !--- Create a netcdf4 file, overwriting any existing file by that name
          !--- The new file is in define mode and ready to write
          imode = 1
          ncstat = nf90_create(fname, ior(NF90_NETCDF4,NF90_CLOBBER), fid)
          msg = "ncdf_open: nf90_create("//trim(fname)//",..."
          call nc_error_handler(fid, ncstat, 1, msg=msg)

        case ("create_noclobber", "cnc")
          !--- Create a netcdf4 file, abort if any file by that name exists
          !--- The new file is in define mode and ready to write
          imode = 2
          ncstat = nf90_create(fname, ior(NF90_NETCDF4,NF90_NOCLOBBER), fid)
          msg = "ncdf_open: nf90_create("//trim(fname)//",..."
          call nc_error_handler(fid, ncstat, 2, msg=msg)

        case ("read", "r")
          !--- Open the netcdf file as read only
          imode = 3
          ncstat = nf90_open(fname, NF90_NOWRITE, fid)
          msg = "ncdf_open: nf90_open("//trim(fname)//",..."
          call nc_error_handler(fid, ncstat, 3, msg=msg)

        case ("write", "w")
          !--- Open an existing netcdf file in data mode and ready to write
          imode = 4
          ncstat = nf90_open(fname, NF90_WRITE, fid)
          msg = "ncdf_open: nf90_open("//trim(fname)//",..."
          call nc_error_handler(fid, ncstat, 4, msg=msg)

        case default
          write(6,*)"ncdf_open: Invalid mode = ",trim(mode)
          call xit("ncdf_open",-2)

      end select

      !--- Find an available file index for fid in ncdf_file_name
      idx = ncdf_file_index(fid, .false.)

      !--- Add this file name to ncdf_file_name
      ncdf_file_name(idx)%fid  = fid
      ncdf_file_name(idx)%name = trim(fname)
      ncdf_file_name(idx)%mode = imode

    end subroutine ncdf_open

    !***************************************************************************
    !--- Close an open netcdf file
    !***************************************************************************
    subroutine ncdf_close(fid)
      !--- The netcdf file ID
      integer(kind=4) :: fid

      !--- Local
      integer(kind=4) :: ncstat, idx

      !--- Find an existing file index for fid in ncdf_file_name
      idx = ncdf_file_index(fid, .true.)

      ncstat = nf90_close(fid)
      if (ncstat /= NF90_NOERR) then
        !--- Issue an error message and abort
        write(6,'(a)') trim( nf90_strerror(ncstat) )
        write(6,*)"ncdf_close: Error on close. file = ",trim(ncdf_file_name(idx)%name)
        call xit("ncdf_close",-1)
      else
        if ( ncdf_file_name(idx)%mode == 1 .or. ncdf_file_name(idx)%mode == 2 ) then
          write(6,*)"Created ",trim(ncdf_file_name(idx)%name)
        endif
        !--- Recycle this element of ncdf_file_name
        ncdf_file_name(idx)%fid  = -1
        ncdf_file_name(idx)%name = " "
        ncdf_file_name(idx)%mode = -1
      endif

    end subroutine ncdf_close

    !***************************************************************************
    !--- Determine the index in ncdf_file_name corresponding to a netcdf file ID
    !***************************************************************************
    function ncdf_file_index(fid, existing) result(idx_file)
      !--- The netcdf file ID
      integer(kind=4) :: fid

      logical :: existing

      !--- The index in ncdf_file_name corresponding fid is returned
      integer :: idx_file

      !--- Local
      integer :: idx1

      idx_file = 0
      if ( existing ) then
        !--- Find the index of an existing file ID
        do idx1=1,size(ncdf_file_name)
          if ( fid == ncdf_file_name(idx1)%fid ) then
            idx_file = idx1
            exit
          endif
        enddo
      else
        !--- Find the first undefined index in ncdf_file_name
        do idx1=1,size(ncdf_file_name)
          if ( ncdf_file_name(idx1)%fid < 0 ) then
            idx_file = idx1
            exit
          endif
        enddo
      endif

      if ( idx_file == 0 ) then
        write(6,*)"ncdf_file_index: Unable to find file index."
        write(6,*)"ncdf_file_index: fid=",fid,"  existing=",existing
        call xit("ncdf_file_index",-1)
      endif

    end function ncdf_file_index

    !***************************************************************************
    !--- Put a netcdf file into define mode
    !--- In define mode dimensions, variables, and attributes can be added or
    !--- renamed and attributes can be deleted
    !***************************************************************************
    subroutine ncdf_define_mode(fid)
      !--- The netcdf file ID
      integer(kind=4) :: fid

      !--- Local
      integer(kind=4) :: ncstat
      integer :: verbose = 0
      character(64) :: msg

      ncstat = nf90_redef(fid)
      if (ncstat /= NF90_NOERR) then
        if (ncstat == NF90_EINDEFINE) then
          !--- Ignore the "Operation not allowed in define mode" error (NF90_EINDEFINE)
          !--- because this means that we were already in define mode before the call
          if ( verbose > 0 ) then
            !--- Optionally issue a warning
            write(6,'(a)') trim( nf90_strerror(ncstat) )
            write(6,*)"ncdf_define_mode: Ignored ncstat=", ncstat
          endif
        else
          !--- Abort on any other error from nf90_redef
          msg = " "
          write(msg,'(a,i6)')"ncdf_define_mode: ncstat = ",ncstat
          call nc_error_handler(fid, ncstat, 1, msg=msg)
        endif
      endif

    end subroutine ncdf_define_mode

    !***************************************************************************
    !--- Put a netcdf file into data mode
    !--- In data mode existing data values can be accessed and changed,
    !--- existing attributes can be changed as long as they do not grow,
    !--- but nothing can be added
    !***************************************************************************
    subroutine ncdf_data_mode(fid)
      !--- The netcdf file ID
      integer(kind=4) :: fid

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      ncstat = nf90_enddef(fid)
      msg = " "
      write(msg,'(a,i6)')"ncdf_data_mode: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

    end subroutine ncdf_data_mode

    !***************************************************************************
    !--- Return info about dimensions found in an open file
    !***************************************************************************
    subroutine ncdf_inquire_dims(fid, dim_names, dim_ids, dim_sizes, rec_dimid)
      !--- The netcdf file ID
      integer(kind=4) :: fid

      character(64),   pointer, intent(out), optional :: dim_names(:)
      integer(kind=4), pointer, intent(out), optional :: dim_ids(:)
      integer,         pointer, intent(out), optional :: dim_sizes(:)
      integer(kind=4),          intent(out), optional :: rec_dimid

      !--- Local
      integer(kind=4) :: ncstat
      character(64)   :: msg, dimname
      integer(kind=4) :: ndims, nvars, natts, rec_dim_id
      integer(kind=4) :: idx, dimlen, dimids(20)

      !--- Get the number of dims and the record dimension ID
      ncstat = nf90_Inquire(fid, nDimensions=ndims, unlimitedDimId=rec_dim_id)
      msg = "ncdf_inquire_dims: nf90_Inquire(fid, "
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(dim_names) ) then
        if ( associated(dim_names) ) deallocate(dim_names)
        allocate( dim_names(ndims) )
        dim_names(:) = " "
      endif
      if ( present(dim_ids) ) then
        if ( associated(dim_ids) ) deallocate(dim_ids)
        allocate( dim_ids(ndims) )
        dim_ids(:) = -1
      endif
      if ( present(dim_sizes) ) then
        if ( associated(dim_sizes) ) deallocate(dim_sizes)
        allocate( dim_sizes(ndims) )
        dim_sizes(:) = 0
      endif
      if ( present(rec_dimid) ) rec_dimid = rec_dim_id

      do idx=1,ndims
        ncstat = nf90_Inquire_Dimension(fid, idx, name=dimname, len=dimlen)
        msg = "ncdf_inquire_dims: nf90_Inquire_Dimension(fid, "
        call nc_error_handler(fid, ncstat, 1, msg=msg)
        if ( present(dim_names) ) dim_names(idx) = trim(dimname)
        if ( present(dim_ids) )   dim_ids(idx)   = idx
        if ( present(dim_sizes) ) dim_sizes(idx) = dimlen
      enddo

    end subroutine ncdf_inquire_dims

    !***************************************************************************
    !--- Return info about variables found in a file. The file is opened here.
    !***************************************************************************
    subroutine ncdf_inquire_vars_by_file_name(fname, var_names, var_ids, &
                                          rec_dimid, nvars)
      !--- The netcdf file name
      character(*), intent(in) :: fname

      character(64),   pointer, intent(inout), optional :: var_names(:)
      integer(kind=4), pointer, intent(inout), optional :: var_ids(:)
      integer(kind=4),          intent(inout), optional :: rec_dimid
      integer(kind=4),          intent(inout), optional :: nvars

      !--- Local
      integer(kind=4) :: fid
      character(64),   pointer :: vnames(:) => null()
      integer(kind=4), pointer :: vids(:) => null()
      integer(kind=4)          :: recdimid, n_vars

      !--- Open the netcdf file for reading
      call ncdf_open(fid, trim(fname), mode="read")

      !--- Obtain local copies of all output variables
      call ncdf_inquire_vars_by_ID(fid, &
           var_names=vnames, var_ids=vids, rec_dimid=recdimid, nvars=n_vars)

      if ( present(var_names) ) then
        if ( associated(var_names) ) deallocate(var_names)
        allocate( var_names(size(vnames)) )
        var_names = vnames
      endif

      if ( present(var_ids) ) then
        if ( associated(var_ids) ) deallocate(var_ids)
        allocate( var_ids(size(vids)) )
        var_ids = vids
      endif

      if ( present(rec_dimid) ) rec_dimid = recdimid
      if ( present(nvars)     ) nvars = n_vars

      !--- Close the file
      call ncdf_close(fid)

      !--- Clean up
      deallocate( vnames, vids )

    end subroutine ncdf_inquire_vars_by_file_name

    !***************************************************************************
    !--- Return info about variables found in an open file
    !***************************************************************************
    subroutine ncdf_inquire_vars_by_ID(fid, var_names, var_ids, &
                                       rec_dimid, nvars)
      !--- The netcdf file ID
      integer(kind=4) :: fid

      character(64),   pointer, intent(out), optional :: var_names(:)
      integer(kind=4), pointer, intent(out), optional :: var_ids(:)
      integer(kind=4),          intent(out), optional :: rec_dimid
      integer(kind=4),          intent(out), optional :: nvars

      !--- Local
      integer(kind=4) :: ncstat
      character(64)   :: msg, dimname, vname
      integer(kind=4) :: ndims, n_vars, rec_dim_id
      integer(kind=4) :: idx, dimlen, dimids(20)

      !--- Get the number of vars and the record dimension ID
      ncstat = nf90_Inquire(fid, nVariables=n_vars, unlimitedDimId=rec_dim_id)
      msg = "ncdf_inquire_vars_by_ID: nf90_Inquire(fid, "
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nvars) ) nvars = n_vars

      if ( present(var_names) ) then
        if ( associated(var_names) ) deallocate(var_names)
        allocate( var_names(n_vars) )
        var_names(:) = " "
      endif

      if ( present(var_ids) ) then
        if ( associated(var_ids) ) deallocate(var_ids)
        allocate( var_ids(n_vars) )
        var_ids(:) = -1
      endif

      if ( present(rec_dimid) ) rec_dimid = rec_dim_id

      do idx=1,n_vars
        ncstat = nf90_Inquire_Variable(fid, idx, name=vname, ndims=ndims, dimids=dimids)
        msg = "ncdf_inquire_vars_by_ID: nf90_Inquire_Variable(fid, "
        call nc_error_handler(fid, ncstat, 2, msg=msg)
        if ( present(var_names) ) var_names(idx) = trim(vname)
        if ( present(var_ids) )   var_ids(idx)   = idx
      enddo

    end subroutine ncdf_inquire_vars_by_ID

    !***************************************************************************
    !--- Given the name of a netcdf file and the name of a variable in that
    !--- file, return the variable ID associated with that variable
    !***************************************************************************
    function ncdf_get_varid(fname, vname) result(vid)
      integer(kind=4) :: vid

      !--- The name of a netcdf file
      character(*), intent(in) :: fname

      !--- The name of a variable in this netcdf file
      character(*), intent(in) :: vname

      !--- Local
      integer :: idx
      character(64),   pointer :: var_names(:) => null()
      integer(kind=4), pointer :: var_ids(:) => null()

      vid = -1
      call ncdf_inquire_vars_by_file_name(fname, var_names=var_names, var_ids=var_ids)
      do idx=1,size(var_names)
        if ( trim(adjustl(var_names(idx))) .eq. trim(adjustl(vname)) ) then
          vid = var_ids(idx)
          exit
        endif
      enddo

      if ( vid < 0 ) then
        write(6,*)"ncdf_get_varid: Unable to determine variable ID for ", &
                  trim(vname)," in file ",trim(fname)
        call xit("ncdf_get_varid",-1)
      endif

    end function ncdf_get_varid

    !***************************************************************************
    !--- Get a global or variable character attribute
    !***************************************************************************
    subroutine ncdf_get_att_char(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      character(*), intent(out):: val

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      ncstat = nf90_get_att( fid, varid, name, val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_get_att_char: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

    end subroutine ncdf_get_att_char

    !***************************************************************************
    !--- Get a global or variable scalar integer*4 attribute
    !***************************************************************************
    subroutine ncdf_get_att_si4(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      integer(kind=4), intent(out):: val

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      ncstat = nf90_get_att( fid, varid, name, val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_get_att_si4: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

    end subroutine ncdf_get_att_si4

    !***************************************************************************
    !--- Get a global or variable integer*4 array attribute
    !***************************************************************************
    subroutine ncdf_get_att_ai4(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      integer(kind=4), intent(out):: val(:)

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      ncstat = nf90_get_att( fid, varid, name, val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_get_att_ai4: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

    end subroutine ncdf_get_att_ai4

    !***************************************************************************
    !--- Get a global or variable scalar integer*8 attribute
    !***************************************************************************
    subroutine ncdf_get_att_si8(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      integer(kind=8), intent(out):: val

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      ncstat = nf90_get_att( fid, varid, name, val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_get_att_si8: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

    end subroutine ncdf_get_att_si8

    !***************************************************************************
    !--- Get a global or variable integer*8 array attribute
    !***************************************************************************
    subroutine ncdf_get_att_ai8(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      integer(kind=8), intent(out):: val(:)

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      ncstat = nf90_get_att( fid, varid, name, val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_get_att_ai8: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

    end subroutine ncdf_get_att_ai8

    !***************************************************************************
    !--- Get a global or variable scalar real*4 attribute
    !***************************************************************************
    subroutine ncdf_get_att_sr4(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      real(kind=4), intent(out):: val

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      ncstat = nf90_get_att( fid, varid, name, val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_get_att_sr4: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

    end subroutine ncdf_get_att_sr4

    !***************************************************************************
    !--- Get a global or variable real*4 array attribute
    !***************************************************************************
    subroutine ncdf_get_att_ar4(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      real(kind=4), intent(out):: val(:)

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      ncstat = nf90_get_att( fid, varid, name, val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_get_att_ar4: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

    end subroutine ncdf_get_att_ar4

    !***************************************************************************
    !--- Get a global or variable scalar real*8 attribute
    !***************************************************************************
    subroutine ncdf_get_att_sr8(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      real(kind=8), intent(out):: val

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      ncstat = nf90_get_att( fid, varid, name, val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_get_att_sr8: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

    end subroutine ncdf_get_att_sr8

    !***************************************************************************
    !--- Get a global or variable real*8 array attribute
    !***************************************************************************
    subroutine ncdf_get_att_ar8(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      real(kind=8), intent(out):: val(:)

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      ncstat = nf90_get_att( fid, varid, name, val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_get_att_ar8: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

    end subroutine ncdf_get_att_ar8

    !***************************************************************************
    !--- Add a global or variable character attribute
    !***************************************************************************
    subroutine ncdf_add_att_char(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      character(*), intent(in) :: val

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      ncstat = nf90_put_att( fid, varid, trim(name), trim(val) )
      msg = " "
      write(msg,'(a,i6)')"ncdf_add_att_char: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

    end subroutine ncdf_add_att_char

    !***************************************************************************
    !--- Add a global or variable scalar integer*4 attribute
    !***************************************************************************
    subroutine ncdf_add_att_si4(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      integer(kind=4), intent(in) :: val

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      ncstat = nf90_put_att( fid, varid, trim(name), val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_add_att_si4: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

    end subroutine ncdf_add_att_si4

    !***************************************************************************
    !--- Add a global or variable integer*4 array attribute
    !***************************************************************************
    subroutine ncdf_add_att_ai4(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      integer(kind=4), intent(in) :: val(:)

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      ncstat = nf90_put_att( fid, varid, trim(name), val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_add_att_ai4: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

    end subroutine ncdf_add_att_ai4

    !***************************************************************************
    !--- Add a global or variable scalar integer*8 attribute
    !***************************************************************************
    subroutine ncdf_add_att_si8(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      integer(kind=8), intent(in) :: val

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      ncstat = nf90_put_att( fid, varid, trim(name), val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_add_att_si8: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

    end subroutine ncdf_add_att_si8

    !***************************************************************************
    !--- Add a global or variable integer*8 array attribute
    !***************************************************************************
    subroutine ncdf_add_att_ai8(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      integer(kind=8), intent(in) :: val(:)

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      ncstat = nf90_put_att( fid, varid, trim(name), val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_add_att_ai8: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

    end subroutine ncdf_add_att_ai8

    !***************************************************************************
    !--- Add a global or variable scalar real*4 attribute
    !***************************************************************************
    subroutine ncdf_add_att_sr4(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      real(kind=4), intent(in) :: val

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      ncstat = nf90_put_att( fid, varid, trim(name), val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_add_att_sr4: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

    end subroutine ncdf_add_att_sr4

    !***************************************************************************
    !--- Add a global or variable real*4 array attribute
    !***************************************************************************
    subroutine ncdf_add_att_ar4(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      real(kind=4), intent(in) :: val(:)

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      ncstat = nf90_put_att( fid, varid, trim(name), val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_add_att_ar4: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

    end subroutine ncdf_add_att_ar4

    !***************************************************************************
    !--- Add a global or variable scalar real*8 attribute
    !***************************************************************************
    subroutine ncdf_add_att_sr8(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      real(kind=8), intent(in) :: val

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      ncstat = nf90_put_att( fid, varid, trim(name), val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_add_att_sr8: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

    end subroutine ncdf_add_att_sr8

    !***************************************************************************
    !--- Add a global or variable real*8 array attribute
    !***************************************************************************
    subroutine ncdf_add_att_ar8(fid, varid, name, val)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(in) :: varid

      !--- The name of the attribute
      character(*), intent(in) :: name

      !--- The value of the attribute
      real(kind=8), intent(in) :: val(:)

      !--- Local
      integer(kind=4) :: ncstat
      character(64) :: msg

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      ncstat = nf90_put_att( fid, varid, trim(name), val )
      msg = " "
      write(msg,'(a,i6)')"ncdf_add_att_ar8: ncstat = ",ncstat
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

    end subroutine ncdf_add_att_ar8

    !***************************************************************************
    !--- Copy all the attibute from one file to the others
    !***************************************************************************
    subroutine ncdf_copy_all_att(fid_in, varname_in, fid_out, varname_out)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid_in, fid_out
      !--- The name of the variable
      character(*), intent(in) :: varname_in, varname_out

      !--- Local
      integer(kind=4) :: ncstat, natt, attid,iatt, varid_in, varid_out
      character(64) :: msg, attname

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid_out)

      ! get the natt of the input file
      ncstat = nf90_inq_varid(fid_in, varname_in, varid_in)
      msg = " "
      write(msg,'(a,i6)')"ncdf_copy_all_att 1: ncstat = ",ncstat
      call nc_error_handler(fid_in, ncstat, 1, msg=msg)
      ncstat = nf90_inquire_variable(fid_in, varid_in, nAtts= natt)
      msg = " "
      write(msg,'(a,i6)')"ncdf_copy_all_att 2: ncstat = ",ncstat
      call nc_error_handler(fid_in, ncstat, 1, msg=msg)
      ncstat = nf90_inq_varid(fid_out, varname_out, varid_out)
      msg = " "
      write(msg,'(a,i6)')"ncdf_copy_all_att 3: ncstat = ",ncstat
      call nc_error_handler(fid_out, ncstat, 1, msg=msg)

      do iatt=1,natt
        ncstat = nf90_inq_attname(fid_in, varid_in, iatt, attname)
        msg = " "
        write(msg,'(a,a,a,i6)')"ncdf_copy_all_att ",trim(attname),": ncstat = ",ncstat
        if (trim(attname).eq."_FillValue") cycle
        call nc_error_handler(fid_in, ncstat, 1, msg=msg)
        ncstat = nf90_copy_att(fid_in, varid_in, attname, fid_out, varid_out)
        msg = " "
        write(msg,'(a,a,a,i6)')"ncdf_copy_all_att ",trim(attname),": ncstat = ",ncstat
        call nc_error_handler(fid_in, ncstat, 1, msg=msg)
      enddo

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid_out)

    end subroutine ncdf_copy_all_att


    !***************************************************************************
    !--- Define a dimension in an open netcdf file
    !***************************************************************************
    subroutine ncdf_add_dim(fid, dimid, name, dimlen, record_dim)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf dimension ID defined here
      integer(kind=4), intent(out) :: dimid

      !--- The name of the dimension
      character(*), intent(in) :: name

      !--- The length of the dimension
      !--- This is only optional when defining the record dimension
      integer, intent(in), optional :: dimlen

      !--- A logical flag used to determine if this coordinate
      !--- is associated with the record dimension
      logical, intent(in), optional :: record_dim

      !--- Local
      integer(kind=4) :: ncstat, dimlen4
      logical :: rec_dim
      character(64) :: msg

      if ( present(record_dim) ) then
        rec_dim = record_dim
      else
        rec_dim = .false.
      endif

      !--- local 4 byte version of dimlen
      if ( present(dimlen) ) then
        dimlen4 = dimlen
      else
        if ( .not. rec_dim ) then
          write(6,*)"ncdf_add_dim: dimlen must be present unless this is the record dimension"
          call xit("ncdf_add_dim",-1)
        endif
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      !--- Define the dimension ID
      if ( rec_dim ) then
        ncstat = nf90_def_dim(fid, trim(name), NF90_UNLIMITED, dimid)
        msg = "ncdf_add_dim: nf90_def_dim(fid, "//trim(name)//",NF90_UNLIMITED"
        call nc_error_handler(fid, ncstat, 1, msg=msg)
      else
        ncstat = nf90_def_dim(fid, trim(name), dimlen4, dimid)
        msg = "ncdf_add_dim: nf90_def_dim(fid, "//trim(name)//",..."
        call nc_error_handler(fid, ncstat, 2, msg=msg)
      endif

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

    end subroutine ncdf_add_dim

    !***************************************************************************
    !--- Define a coordinate variable in an open netcdf file
    !***************************************************************************
    subroutine ncdf_add_coord_any(fid, varid, dimid, cname, record_dim, xtype, &
                                  val_a1i4, val_a1i8, val_a1r4, val_a1r8, val_a1char)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID defined here
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension ID defined here
      integer(kind=4), intent(out) :: dimid

      !--- The name of the coordinate variable to be defined
      character(*), intent(in) :: cname

      !--- A logical flag used to determine if this coordinate
      !--- is associated with the record dimension
      logical, intent(in), optional :: record_dim

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- Coordinate values to be written to the file
      integer(4),   intent(in), optional :: val_a1i4(:)
      integer(8),   intent(in), optional :: val_a1i8(:)
      real(4),      intent(in), optional :: val_a1r4(:)
      real(8),      intent(in), optional :: val_a1r8(:)
      character(*), intent(in), optional :: val_a1char(:)

      !--- Local
      integer :: idx, vtype
      integer(kind=4) :: ncstat, n_dim, local_xtype
      logical :: rec_dim
      character(64) :: msg
      logical :: use_var(5)

      if ( present(record_dim) ) then
        rec_dim = record_dim
      else
        rec_dim = .false.
      endif

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_DOUBLE
      endif

      if ( present(val_a1char) ) then
        write(6,*)"ncdf_add_coord_any: Character type coordinates are not supported."
        call xit("ncdf_add_coord_any",-1)
        !--- The external type must always be NF90_CHAR in this case
        ! local_xtype = NF90_CHAR
      endif

      !--- Identify the incomming variable type
      use_var(:) = .false.
      if ( present(val_a1i4)   ) use_var( 1) = .true.
      if ( present(val_a1i8)   ) use_var( 2) = .true.
      if ( present(val_a1r4)   ) use_var( 3) = .true.
      if ( present(val_a1r8)   ) use_var( 4) = .true.
      if ( present(val_a1char) ) use_var( 5) = .true.

      !--- Ensure that there is exactly one incomming variable
      idx = count(use_var)
      if ( .not. rec_dim .and. idx == 0 ) then
        write(6,*)"ncdf_add_coord_any: User supplied coordinate data is missing."
        write(6,*)"ncdf_add_coord_any: Data is required unless this is the record dimension."
        call xit("ncdf_add_coord_any",-1)
      endif
      if ( idx > 1 ) then
        write(6,*)"ncdf_add_coord_any: Only one input variable array is allowed."
        call xit("ncdf_add_coord_any",-2)
      endif

      !--- Determine the index in use_var of the input data type
      vtype = 0
      do idx=1,size(use_var)
        if ( use_var(idx) ) then
          vtype = idx
          !--- There is only 1 variable allowed so we exit the loop
          exit
        endif
      enddo
      if ( .not. rec_dim .and. vtype == 0) then
        write(6,*)"ncdf_add_coord_any: Unable to determine input variable type and kind."
        call xit("ncdf_add_coord_any",-3)
      endif

      if (rec_dim) then
        n_dim = NF90_UNLIMITED
      else
        n_dim = 0
        select case(vtype)
          case (1) !--- 1D integer*4 input array
            n_dim = size(val_a1i4)
          case (2) !--- 1D integer*8 input array
            n_dim = size(val_a1i8)
          case (3) !--- 1D real*4 input array
            n_dim = size(val_a1r4)
          case (4) !--- 1D real*4 input array
            n_dim = size(val_a1r8)
        end select
        if ( n_dim == 0 ) then
          write(6,*)"ncdf_add_coord_any: No coordinate data is present."
          call xit("ncdf_add_coord_any",-4)
        endif
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      !--- Define the dimension ID
      ncstat = nf90_def_dim(fid, trim(cname), n_dim, dimid)
      msg = "ncdf_add_coord_any: nf90_def_dim(fid, "//trim(cname)//",..."
      call nc_error_handler(fid, ncstat, 2, msg=msg)

      !--- Define coordinate variable
      ncstat = nf90_def_var(fid, trim(cname), local_xtype, dimid, varid)
      msg = "ncdf_add_coord_any: nf90_def_var(fid, "//trim(cname)//",..."
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      if ( vtype > 0 ) then
        !--- Write coordinate values to the file
        select case(vtype)
          case (1) !--- 1D integer*4 input array
            ncstat = nf90_put_var(fid, varid, val_a1i4)
          case (2) !--- 1D integer*8 input array
            ncstat = nf90_put_var(fid, varid, val_a1i8)
          case (3) !--- 1D real*4 input array
            ncstat = nf90_put_var(fid, varid, val_a1r4)
          case (4) !--- 1D real*4 input array
            ncstat = nf90_put_var(fid, varid, val_a1r8)
        end select
        msg = "ncdf_add_coord_any: nf90_put_var(...)  coord="//trim(cname)
        call nc_error_handler(fid, ncstat, 4, msg=msg)
      endif

    end subroutine ncdf_add_coord_any

    !-----------------------------------------------------------------------------
    !--- Define a record dimension coordinate variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_add_coord_none(fid, varid, dimid, cname, record_dim, xtype)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension ID
      integer(kind=4), intent(out) :: dimid

      !--- The name of the coordinate variable to be defined
      character(*), intent(in) :: cname

      !--- A logical flag used to determine if this coordinate
      !--- is associated with the record dimension
      logical, intent(in), optional :: record_dim

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      if ( present(record_dim) .and. present(xtype) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, &
                                record_dim=record_dim, xtype=xtype)
      else if ( present(record_dim) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, record_dim=record_dim)
      else if ( present(xtype) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, xtype=xtype)
      else
        call ncdf_add_coord_any(fid, varid, dimid, cname)
      endif

    end subroutine ncdf_add_coord_none

    !-----------------------------------------------------------------------------
    !--- Add a 1D integer*4 coordinate variable to an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_add_coord_i4(fid, varid, dimid, cname, vals, record_dim, xtype)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension ID
      integer(kind=4), intent(out) :: dimid

      !--- The name of the coordinate variable to be defined
      character(*), intent(in) :: cname

      !--- The variable values to be written to the file
      integer(kind=4), intent(in) :: vals(:)

      !--- A logical flag used to determine if this coordinate
      !--- is associated with the record dimension
      logical, intent(in), optional :: record_dim

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      if ( present(record_dim) .and. present(xtype) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1i4=vals, &
                                record_dim=record_dim, xtype=xtype)
      else if ( present(record_dim) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1i4=vals, &
                                record_dim=record_dim)
      else if ( present(xtype) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1i4=vals, &
                                xtype=xtype)
      else
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1i4=vals)
      endif

    end subroutine ncdf_add_coord_i4

    !-----------------------------------------------------------------------------
    !--- Add a 1D integer*8 coordinate variable to an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_add_coord_i8(fid, varid, dimid, cname, vals, record_dim, xtype)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension ID
      integer(kind=4), intent(out) :: dimid

      !--- The name of the coordinate variable to be defined
      character(*), intent(in) :: cname

      !--- The variable values to be written to the file
      integer(kind=8), intent(in) :: vals(:)

      !--- A logical flag used to determine if this coordinate
      !--- is associated with the record dimension
      logical, intent(in), optional :: record_dim

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      if ( present(record_dim) .and. present(xtype) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1i8=vals, &
                                record_dim=record_dim, xtype=xtype)
      else if ( present(record_dim) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1i8=vals, &
                                record_dim=record_dim)
      else if ( present(xtype) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1i8=vals, &
                                xtype=xtype)
      else
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1i8=vals)
      endif

    end subroutine ncdf_add_coord_i8

    !-----------------------------------------------------------------------------
    !--- Add a 1D real*4 coordinate variable to an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_add_coord_r4(fid, varid, dimid, cname, vals, record_dim, xtype)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension ID
      integer(kind=4), intent(out) :: dimid

      !--- The name of the coordinate variable to be defined
      character(*), intent(in) :: cname

      !--- The variable values to be written to the file
      real(kind=4), intent(in) :: vals(:)

      !--- A logical flag used to determine if this coordinate
      !--- is associated with the record dimension
      logical, intent(in), optional :: record_dim

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      if ( present(record_dim) .and. present(xtype) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1r4=vals, &
                                record_dim=record_dim, xtype=xtype)
      else if ( present(record_dim) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1r4=vals, &
                                record_dim=record_dim)
      else if ( present(xtype) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1r4=vals, &
                                xtype=xtype)
      else
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1r4=vals)
      endif

    end subroutine ncdf_add_coord_r4

    !-----------------------------------------------------------------------------
    !--- Add a 1D real*8 coordinate variable to an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_add_coord_r8(fid, varid, dimid, cname, vals, record_dim, xtype)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension ID
      integer(kind=4), intent(out) :: dimid

      !--- The name of the coordinate variable to be defined
      character(*), intent(in) :: cname

      !--- The variable values to be written to the file
      real(kind=8), intent(in) :: vals(:)

      !--- A logical flag used to determine if this coordinate
      !--- is associated with the record dimension
      logical, intent(in), optional :: record_dim

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      if ( present(record_dim) .and. present(xtype) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1r8=vals, &
                                record_dim=record_dim, xtype=xtype)
      else if ( present(record_dim) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1r8=vals, &
                                record_dim=record_dim)
      else if ( present(xtype) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1r8=vals, &
                                xtype=xtype)
      else
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1r8=vals)
      endif

    end subroutine ncdf_add_coord_r8

    !-----------------------------------------------------------------------------
    !--- Add a 1D character coordinate variable to an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_add_coord_char(fid, varid, dimid, cname, vals, record_dim, xtype)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension ID
      integer(kind=4), intent(out) :: dimid

      !--- The name of the coordinate variable to be defined
      character(*), intent(in) :: cname

      !--- The variable values to be written to the file
      character(*), intent(in) :: vals(:)

      !--- A logical flag used to determine if this coordinate
      !--- is associated with the record dimension
      logical, intent(in), optional :: record_dim

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      if ( present(record_dim) .and. present(xtype) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1char=vals, &
                                record_dim=record_dim, xtype=xtype)
      else if ( present(record_dim) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1char=vals, &
                                record_dim=record_dim)
      else if ( present(xtype) ) then
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1char=vals, &
                                xtype=xtype)
      else
        call ncdf_add_coord_any(fid, varid, dimid, cname, val_a1char=vals)
      endif

    end subroutine ncdf_add_coord_char

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file
    !***************************************************************************
    subroutine ncdf_add_var_no_vals(fid, varid, vname, dimids, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_DOUBLE
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_no_vals: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_no_vals: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

    end subroutine ncdf_add_var_no_vals

    !***************************************************************************
    !--- Set the value of a single time coordinate
    !***************************************************************************
    subroutine ncdf_add_var_1d_i4(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- The variable values to be written to the file
      integer(kind=4), intent(in) :: vals(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_INT
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_1d_i4: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_1d_i4: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_1d_i4

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (1D integer*8)
    !***************************************************************************
    subroutine ncdf_add_var_1d_i8(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      integer(kind=8), intent(in) :: vals(:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_INT64
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_1d_i8: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_1d_i8: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_1d_i8

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (1D real*4)
    !***************************************************************************
    subroutine ncdf_add_var_1d_r4(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      real(kind=4), intent(in) :: vals(:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_FLOAT
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_1d_r4: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_1d_r4: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_1d_r4

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (1D real*8)
    !***************************************************************************
    subroutine ncdf_add_var_1d_r8(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      real(kind=8), intent(in) :: vals(:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_DOUBLE
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_1d_r8: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_1d_r8: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_1d_r8

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (2D integer*4)
    !***************************************************************************
    subroutine ncdf_add_var_2d_i4(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      integer(kind=4), intent(in) :: vals(:,:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_INT
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_2d_i4: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_2d_i4: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_2d_i4

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (2D integer*8)
    !***************************************************************************
    subroutine ncdf_add_var_2d_i8(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      integer(kind=8), intent(in) :: vals(:,:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_INT64
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_2d_i8: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_2d_i8: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_2d_i8

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (2D real*4)
    !***************************************************************************
    subroutine ncdf_add_var_2d_r4(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      real(kind=4), intent(in) :: vals(:,:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_FLOAT
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_2d_r4: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_2d_r4: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_2d_r4

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (2D real*8)
    !***************************************************************************
    subroutine ncdf_add_var_2d_r8(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      real(kind=8), intent(in) :: vals(:,:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_DOUBLE
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_2d_r8: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_2d_r8: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_2d_r8

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (3D integer*4)
    !***************************************************************************
    subroutine ncdf_add_var_3d_i4(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      integer(kind=4), intent(in) :: vals(:,:,:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_INT
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_3d_i4: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_3d_i4: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_3d_i4

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (3D integer*8)
    !***************************************************************************
    subroutine ncdf_add_var_3d_i8(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      integer(kind=8), intent(in) :: vals(:,:,:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_INT64
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_3d_i8: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_3d_i8: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_3d_i8

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (3D real*4)
    !***************************************************************************
    subroutine ncdf_add_var_3d_r4(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      real(kind=4), intent(in) :: vals(:,:,:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_FLOAT
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_3d_r4: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_3d_r4: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_3d_r4

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (3D real*8)
    !***************************************************************************
    subroutine ncdf_add_var_3d_r8(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      real(kind=8), intent(in) :: vals(:,:,:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_DOUBLE
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_3d_r8: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_3d_r8: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_3d_r8

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (4D integer*4)
    !***************************************************************************
    subroutine ncdf_add_var_4d_i4(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      integer(kind=4), intent(in) :: vals(:,:,:,:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_INT
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_4d_i4: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_4d_i4: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_4d_i4

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (4D integer*8)
    !***************************************************************************
    subroutine ncdf_add_var_4d_i8(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      integer(kind=8), intent(in) :: vals(:,:,:,:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_INT64
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_4d_i8: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_4d_i8: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_4d_i8

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (4D real*4)
    !***************************************************************************
    subroutine ncdf_add_var_4d_r4(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      real(kind=4), intent(in) :: vals(:,:,:,:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_FLOAT
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_4d_r4: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_4d_r4: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_4d_r4

    !***************************************************************************
    !--- Define a non-coordinate variable in an open netcdf file (4D real*8)
    !***************************************************************************
    subroutine ncdf_add_var_4d_r8(fid, varid, vname, dimids, vals, xtype, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The name of the variable to be defined
      character(*), intent(in) :: vname

      !--- Return the netcdf variable ID
      integer(kind=4), intent(out) :: varid

      !--- The netcdf dimension IDs for the current variable
      !--- The length of dimids must be the rank of vals
      integer(kind=4), intent(in) :: dimids(:)

      !--- External type for this variable
      integer(kind=4), intent(in), optional :: xtype

      !--- The variable values to be written to the file
      real(kind=8), intent(in) :: vals(:,:,:,:)

      !--- The index in the record dimension at which vals is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, local_xtype, idx, dimlen
      character(64) :: msg
      integer :: verbose=0

      if ( present(xtype) ) then
        local_xtype = xtype
      else
        local_xtype = NF90_DOUBLE
      endif

      !--- Put the netcdf file into define mode
      call ncdf_define_mode(fid)

      if ( verbose > 0 ) then
        write(6,*)"ncdf_add_var_4d_r8: dim name=",trim(vname),"  dimids=",dimids
      endif

      !--- Define this variable in the netcdf file, returning varid
      ncstat = nf90_def_var(fid, trim(vname), local_xtype, dimids, varid)
      msg = "ncdf_add_var_4d_r8: nf90_def_var(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      !--- Put the netcdf file back into data mode
      call ncdf_data_mode(fid)

      !--- Write coordinate values to the file
      if ( present(nt) ) then
        call ncdf_write_var(fid, varid, vals, nt=nt)
      else
        !--- Write to the first record index (if applicable) by default
        !--- since this variable has just been defined
        call ncdf_write_var(fid, varid, vals, nt=1)
      endif

    end subroutine ncdf_add_var_4d_r8

    !***************************************************************************
    !--- Read the entire record dimension coordinate variable (normally time)
    !***************************************************************************
    subroutine ncdf_get_time_by_ID(fid, time_vals, ntime, units)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The time coordinate values
      real(kind=8), pointer, intent(out) :: time_vals(:)

      !--- The number of time values found in the file
      integer, intent(out), optional :: ntime

      !--- The units associated with the time coordinate
      character(*), intent(out), optional :: units

      !--- Local
      integer(kind=4) :: ncstat, rec_dimid, time_varid, dimlen
      character(64)   :: msg, rec_name

      !--- Determine the record dimension ID
      ncstat = nf90_Inquire(fid, unlimitedDimId=rec_dimid)
      msg = "ncdf_get_time_by_ID: nf90_Inquire(fid, "
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Get the name and length of the record dimension
      ncstat = nf90_Inquire_Dimension(fid, rec_dimid, name=rec_name, len=dimlen)
      msg = "ncdf_get_time_by_ID: nf90_Inquire_Dimension(fid,"
      call nc_error_handler(fid, ncstat, 2, msg=msg)

      !--- Get the variable ID of the time coordinate from its name
      !--- This is the same name as the record (aka unlimited) dimension
      ncstat = nf90_inq_varid(fid, rec_name, time_varid)
      msg = "ncdf_get_time_by_ID: nf90_inq_varid(fid, "
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      if ( present(units) ) then
        !--- Read the units attribute from the time coordinate
        ncstat = nf90_get_att(fid, time_varid, "units", units)
        if ( ncstat /= NF90_NOERR ) then
          write(6,*)"ncdf_get_time_by_ID: The record dimension coordinate ",trim(rec_name), &
                    " does not contain a units attribute."
          units = " "
        endif
      endif

      !--- If there are no time slices written to the file then simply
      !--- return without comment
      !--- The caller is responsible for checking the value of ntime
      ntime = dimlen
      if ( ntime < 1 ) return

      if ( associated(time_vals) ) deallocate(time_vals)
      allocate( time_vals(dimlen) )

      write(6,*)"  ncdf_get_time_by_ID: Found ",ntime," time records"

      !--- Read this value to the time_vals array
      ncstat = nf90_get_var(fid, time_varid, time_vals)
      msg = "ncdf_get_time_by_ID: nf90_get_var(fid, "
      call nc_error_handler(fid, ncstat, 4, msg=msg)

    end subroutine ncdf_get_time_by_ID

    !***************************************************************************
    !--- Return info about variables found in a file. The file is opened here.
    !***************************************************************************
    subroutine ncdf_get_time_by_file_name(fname, time_vals, ntime, units)
      !--- The netcdf file name
      character(*), intent(in) :: fname

      !--- The time coordinate values
      real(kind=8), pointer, intent(inout) :: time_vals(:)

      !--- The number of time values found in the file
      integer, intent(inout), optional :: ntime

      !--- The units associated with the time coordinate
      character(*), intent(inout), optional :: units

      !--- Local
      integer(kind=4) :: fid

      !--- Open the netcdf file for reading
      call ncdf_open(fid, trim(fname), mode="read")

      if ( present(ntime) .and. present(units) ) then
        call ncdf_get_time_by_ID(fid, time_vals, ntime=ntime, units=units)
      else if ( present(ntime) ) then
        call ncdf_get_time_by_ID(fid, time_vals, ntime=ntime)
      else if ( present(units) ) then
        call ncdf_get_time_by_ID(fid, time_vals, units=units)
      else
        call ncdf_get_time_by_ID(fid, time_vals)
      endif

      !--- Close the file
      call ncdf_close(fid)

    end subroutine ncdf_get_time_by_file_name

    !***************************************************************************
    !--- Set a time variable in an open netcdf file
    !***************************************************************************
    subroutine ncdf_set_time(fid, time_val, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The time coordinate value
      real(kind=8), intent(in) :: time_val

      !--- The index in the record dimension at which time_val is written
      integer, intent(in), optional :: nt

      !--- Local
      integer(kind=4) :: ncstat, rec_dimid, time_varid, dimlen
      integer(kind=4) :: vstart(1), vcount(1)
      character(64)   :: msg, rec_name

      !--- Determine the record dimension ID
      ncstat = nf90_Inquire(fid, unlimitedDimId=rec_dimid)
      msg = "ncdf_set_time: nf90_Inquire(fid, "
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Get the name and length of the record dimension
      ncstat = nf90_Inquire_Dimension(fid, rec_dimid, name=rec_name, len=dimlen)
      msg = "ncdf_set_time: nf90_Inquire_Dimension(fid,"
      call nc_error_handler(fid, ncstat, 2, msg=msg)

      !--- Get the variable ID of the time coordinate from its name
      !--- This is the same name as the record (aka unlimited) dimension
      ncstat = nf90_inq_varid(fid, rec_name, time_varid)
      msg = "ncdf_set_time: nf90_inq_varid(fid, "
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      if ( present(nt) ) then
        !--- Write to a user supplied record index
        vstart(1) = nt
      else
        !--- Write to the largest record index
        if ( dimlen == 0 ) then
          vstart(1) = 1
        else
          vstart(1) = dimlen
        endif
      endif
      vcount(1) = 1

write(6,*)"ncdf_set_time: time_val=",time_val,"  vstart=",vstart,"  vcount=",vcount

      !--- Write this value to the time coordinate
      ncstat = nf90_put_var(fid, time_varid, (/time_val/), start=vstart, count=vcount)
      msg = "ncdf_set_time: nf90_put_var(fid, "
      call nc_error_handler(fid, ncstat, 4, msg=msg)

    end subroutine ncdf_set_time

    !****************************************************************************
    !--- Write data values to a known variable in an open netcdf file.
    !--- The input array containing the data to be written can be 1, 2, 3 or 4
    !--- dimensions. It will be written at a particular record dimesion index if
    !--- the target netcdf variable contains the record dimension.
    !--- If nt is present then nt is the record dimension index at which
    !--- the data is written. Otherwise the largest record index is used.
    !--- The input array must be of the same size and shape as the target netcdf
    !--- variable if the netcdf variable does not contain the record dimension.
    !--- If the netcdf variable does contain the record dimension then the input
    !--- array must be of the same size and shape of the netcdf variable
    !--- excluding the record dimension.
    !****************************************************************************
    subroutine ncdf_write_var_any(fid, varid, nt, &
                                  val_a1i4, val_a1i8, val_a1r4, val_a1r8, &
                                  val_a2i4, val_a2i8, val_a2r4, val_a2r8, &
                                  val_a3i4, val_a3i8, val_a3r4, val_a3r8, &
                                  val_a4i4, val_a4i8, val_a4r4, val_a4r8, &
                                  val_a1char)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- The variable values to be written to the file
      integer(4),   intent(in), optional :: val_a1i4(:)
      integer(8),   intent(in), optional :: val_a1i8(:)
      real(4),      intent(in), optional :: val_a1r4(:)
      real(8),      intent(in), optional :: val_a1r8(:)
      integer(4),   intent(in), optional :: val_a2i4(:,:)
      integer(8),   intent(in), optional :: val_a2i8(:,:)
      real(4),      intent(in), optional :: val_a2r4(:,:)
      real(8),      intent(in), optional :: val_a2r8(:,:)
      integer(4),   intent(in), optional :: val_a3i4(:,:,:)
      integer(8),   intent(in), optional :: val_a3i8(:,:,:)
      real(4),      intent(in), optional :: val_a3r4(:,:,:)
      real(8),      intent(in), optional :: val_a3r8(:,:,:)
      integer(4),   intent(in), optional :: val_a4i4(:,:,:,:)
      integer(8),   intent(in), optional :: val_a4i8(:,:,:,:)
      real(4),      intent(in), optional :: val_a4r4(:,:,:,:)
      real(8),      intent(in), optional :: val_a4r8(:,:,:,:)
      character(*), intent(in), optional :: val_a1char(:)

      !--- Local
      integer :: vtype, verbose=0
      integer(kind=4) :: ncstat, idx, dimlen, dimlen_in, rec_dimid
      integer(kind=4) :: ndims, dimids(7), vstart(7), vcount(7), rec_count
      real(8) :: time_val
      character(64) :: msg, vname, dim_name
      logical :: use_var(17), has_rec_dim, problem

      !--- Identify the incomming variable type
      use_var(:) = .false.
      if ( present(val_a1i4)   ) use_var( 1) = .true.
      if ( present(val_a1i8)   ) use_var( 2) = .true.
      if ( present(val_a1r4)   ) use_var( 3) = .true.
      if ( present(val_a1r8)   ) use_var( 4) = .true.
      if ( present(val_a2i4)   ) use_var( 5) = .true.
      if ( present(val_a2i8)   ) use_var( 6) = .true.
      if ( present(val_a2r4)   ) use_var( 7) = .true.
      if ( present(val_a2r8)   ) use_var( 8) = .true.
      if ( present(val_a3i4)   ) use_var( 9) = .true.
      if ( present(val_a3i8)   ) use_var(10) = .true.
      if ( present(val_a3r4)   ) use_var(11) = .true.
      if ( present(val_a3r8)   ) use_var(12) = .true.
      if ( present(val_a4i4)   ) use_var(13) = .true.
      if ( present(val_a4i8)   ) use_var(14) = .true.
      if ( present(val_a4r4)   ) use_var(15) = .true.
      if ( present(val_a4r8)   ) use_var(16) = .true.
      if ( present(val_a1char) ) use_var(17) = .true.

      !--- Ensure that there is exactly one incomming variable
      idx = count(use_var)
      if ( idx == 0 ) then
        write(6,*)"ncdf_write_var_any: User supplied data is missing."
        call xit("ncdf_write_var_any",-1)
      endif
      if ( idx > 1 ) then
        write(6,*)"ncdf_write_var_any: Only one input variable array is allowed."
        call xit("ncdf_write_var_any",-2)
      endif

      !--- Determine the index in use_var of the input data type
      vtype = 0
      do idx=1,size(use_var)
        if ( use_var(idx) ) then
          vtype = idx
          !--- There is only 1 variable allowed so we exit the loop
          exit
        endif
      enddo
      if (vtype == 0) then
        write(6,*)"ncdf_write_var_any: Unable to determine input variable type and kind."
        call xit("ncdf_write_var_any",-3)
      endif

      !--- Determine the record dimension ID
      ncstat = nf90_Inquire(fid, unlimitedDimId=rec_dimid)
      msg = "ncdf_write_var_any: nf90_Inquire(fid, "
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Get name and dimension IDs for the current variable
      ncstat = nf90_Inquire_Variable(fid, varid, name=vname, ndims=ndims, dimids=dimids)
      msg = "ncdf_write_var_any: nf90_Inquire_Variable(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 2, msg=msg)

      !--- Does the netcdf variable contain the record dimension
      if ( dimids(ndims) == rec_dimid ) then
        !--- The record dimension must always be the slowest varying (ie dimids(ndims))
        has_rec_dim = .true.
      else
        has_rec_dim = .false.
      endif

      !--- rec_count is the number of record dimension elements that are to be written
      !--- This is normally 1 but may be larger in the special case of a 1D variable
      !--- associated with the record dimension (typically time) and it is reset below
      !--- in this special case
      rec_count = 1

      !--- Verify that the number of dimension in the input array is consistent with
      !--- the number of dimensions in the netcdf variable
      problem = .false.
      select case(vtype)
        case (1,2,3,4)     !--- 1D input array
          if ( (      has_rec_dim .and. ndims > 2) .or. &
               (.not. has_rec_dim .and. ndims /= 1) ) then
            problem = .true.
          endif
          if ( has_rec_dim .and. ndims == 1 ) then
            !--- This is a special case
            !--- The input array must contain values for the coordinate variable
            !--- associated with the record dimension (normaly time)
            !--- Set a logical flag that will be used below
            select case(vtype)
              case (1) !--- 1D integer*4 input array
                rec_count = size(val_a1i4)
              case (2) !--- 1D integer*8 input array
                rec_count = size(val_a1i8)
              case (3) !--- 1D real*4 input array
                rec_count = size(val_a1r4)
              case (4) !--- 1D real*8 input array
                rec_count = size(val_a1r8)
            end select
          endif
        case (5,6,7,8)     !--- 2D input array
          if ( (      has_rec_dim .and. ndims /= 3) .or. &
               (.not. has_rec_dim .and. ndims /= 2) ) then
            problem = .true.
          endif
        case (9,10,11,12)  !--- 3D input array
          if ( (      has_rec_dim .and. ndims /= 4) .or. &
               (.not. has_rec_dim .and. ndims /= 3) ) then
            problem = .true.
          endif
        case (13,14,15,16) !--- 4D input array
          if ( (      has_rec_dim .and. ndims /= 5) .or. &
               (.not. has_rec_dim .and. ndims /= 4) ) then
            problem = .true.
          endif
        case (17)          !--- 1D input character array
          if ( (      has_rec_dim .and. ndims > 3) .or. &
               (.not. has_rec_dim .and. ndims /= 2) ) then
            problem = .true.
          endif
      end select
      if ( problem ) then
        write(6,*)"ncdf_write_var_any: Number of input dimensions does not match variable"
        call xit("ncdf_write_var_any",-4)
      endif

      !--- Populate start and count arrays and verify dimension sizes
      do idx=1,ndims
        ncstat = nf90_inquire_dimension(fid, dimids(idx), name=dim_name, len=dimlen)
        msg = "ncdf_write_var_any: nf90_Inquire_dimension(fid, "//trim(vname)//" dim="//trim(dim_name)
        call nc_error_handler(fid, ncstat, 3, msg=msg)

        if ( verbose > 0 ) then
          write(6,*)"ncdf_write_var_any: idx, dimids(idx): ",idx,dimids(idx),"  dimlen=",dimlen, &
                    "  vname=",trim(vname),"  dim=",trim(dim_name)
        endif

        !--- Populate start and count arrays
        vstart(idx) = 1
        vcount(idx) = dimlen
        if ( dimids(idx) == rec_dimid ) then
          !--- This is the record dimension
          if ( present(nt) ) then
            !--- The user has supplied the record index
            vstart(idx) = nt
          else
            !--- Use a default record index
            if ( dimlen == 0 ) then
              !--- No time dependent data has be written
              !--- Write the first data for this variable
              vstart(idx) = 1
            else
              !--- Append data along the record dimension
              !--- The next unassigned record index is dimlen (0 based indexing)
              vstart(idx) = dimlen
            endif
          endif
          vcount(idx) = rec_count
          exit
        endif

        !--- Verify dimension sizes
        select case(vtype)
          case ( 1)  !--- 1D array integer*4
            dimlen_in = size(val_a1i4,dim=idx)
          case ( 2)  !--- 1D array integer*8
            dimlen_in = size(val_a1i8,dim=idx)
          case ( 3)  !--- 1D array real*4
            dimlen_in = size(val_a1r4,dim=idx)
          case ( 4)  !--- 1D array real*8
            dimlen_in = size(val_a1r8,dim=idx)
          case ( 5)  !--- 2D array integer*4
            dimlen_in = size(val_a2i4,dim=idx)
          case ( 6)  !--- 2D array integer*8
            dimlen_in = size(val_a2i8,dim=idx)
          case ( 7)  !--- 2D array real*4
            dimlen_in = size(val_a2r4,dim=idx)
          case ( 8)  !--- 2D array real*8
            dimlen_in = size(val_a2r8,dim=idx)
          case ( 9)  !--- 3D array integer*4
            dimlen_in = size(val_a3i4,dim=idx)
          case (10)  !--- 3D array integer*8
            dimlen_in = size(val_a3i8,dim=idx)
          case (11)  !--- 3D array real*4
            dimlen_in = size(val_a3r4,dim=idx)
          case (12)  !--- 3D array real*8
            dimlen_in = size(val_a3r8,dim=idx)
          case (13)  !--- 4D array integer*4
            dimlen_in = size(val_a4i4,dim=idx)
          case (14)  !--- 4D array integer*8
            dimlen_in = size(val_a4i8,dim=idx)
          case (15)  !--- 4D array real*4
            dimlen_in = size(val_a4r4,dim=idx)
          case (16)  !--- 4D array real*8
            dimlen_in = size(val_a4r8,dim=idx)
        end select
        if ( vtype == 17 ) then
          !--- 1D array character
          !--- Character variables will contain an extra dimension
          !--- Each character string is written as an 1D array of single characters
          !--- and the length of each string must be the fastest varying dimension
          if ( idx > 1 ) then
            dimlen_in = size(val_a1char,dim=idx-1)
            if ( dimlen /= dimlen_in ) then
              write(6,*)"ncdf_write_var_any: Dimension ",idx, &
                        " length is inconsistent with the input data."
              write(6,*)"        variable name = ",trim(vname)
              write(6,*)"       dimension name = ",trim(dim_name)
              write(6,*)"    netcdf dim length = ",dimlen
              write(6,*)"input data dim length = ",dimlen_in
              call xit("ncdf_write_var_any",-1)
            endif
          endif
        else
          if ( dimlen /= dimlen_in ) then
            write(6,*)"ncdf_write_var_any: Dimension ",idx, &
                      " length is inconsistent with the input data."
            write(6,*)"        variable name = ",trim(vname)
            write(6,*)"       dimension name = ",trim(dim_name)
            write(6,*)"    netcdf dim length = ",dimlen
            write(6,*)"input data dim length = ",dimlen_in
            call xit("ncdf_write_var_any",-1)
          endif
        endif
      enddo

      if ( verbose > 0 ) then
        write(6,*)"ncdf_write_var_any: start = ",vstart(1:ndims)
        write(6,*)"ncdf_write_var_any: count = ",vcount(1:ndims)
      endif

      !--- Write variable values to the file
      select case(vtype)
        case ( 1)  !--- 1D array integer*4
          ncstat = nf90_put_var(fid, varid, val_a1i4, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case ( 2)  !--- 1D array integer*8
          ncstat = nf90_put_var(fid, varid, val_a1i8, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case ( 3)  !--- 1D array real*4
          ncstat = nf90_put_var(fid, varid, val_a1r4, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case ( 4)  !--- 1D array real*8
          ncstat = nf90_put_var(fid, varid, val_a1r8, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case ( 5)  !--- 2D array integer*4
          ncstat = nf90_put_var(fid, varid, val_a2i4, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case ( 6)  !--- 2D array integer*8
          ncstat = nf90_put_var(fid, varid, val_a2i8, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case ( 7)  !--- 2D array real*4
          ncstat = nf90_put_var(fid, varid, val_a2r4, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case ( 8)  !--- 2D array real*8
          ncstat = nf90_put_var(fid, varid, val_a2r8, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case ( 9)  !--- 3D array integer*4
          ncstat = nf90_put_var(fid, varid, val_a3i4, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case (10)  !--- 3D array integer*8
          ncstat = nf90_put_var(fid, varid, val_a3i8, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case (11)  !--- 3D array real*4
          ncstat = nf90_put_var(fid, varid, val_a3r4, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case (12)  !--- 3D array real*8
          ncstat = nf90_put_var(fid, varid, val_a3r8, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case (13)  !--- 4D array integer*4
          ncstat = nf90_put_var(fid, varid, val_a4i4, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case (14)  !--- 4D array integer*8
          ncstat = nf90_put_var(fid, varid, val_a4i8, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case (15)  !--- 4D array real*4
          ncstat = nf90_put_var(fid, varid, val_a4r4, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case (16)  !--- 4D array real*8
          ncstat = nf90_put_var(fid, varid, val_a4r8, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
        case (17)  !--- 1D array character
          ncstat = nf90_put_var(fid, varid, val_a1char, &
                                start=vstart(1:ndims), count=vcount(1:ndims))
      end select
      msg = "ncdf_write_var_any: nf90_put_var(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 4, msg=msg)

    end subroutine ncdf_write_var_any

    !-----------------------------------------------------------------------------
    !--- Write 1D character values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_1d_char(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      character(*), intent(in)    :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a1char=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a1char=vals)
      endif

    end subroutine ncdf_write_var_1d_char

    !-----------------------------------------------------------------------------
    !--- Write 1D integer*4 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_1d_i4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      integer(kind=4), intent(in) :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a1i4=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a1i4=vals)
      endif

    end subroutine ncdf_write_var_1d_i4

    !-----------------------------------------------------------------------------
    !--- Write 1D integer*8 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_1d_i8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      integer(kind=8), intent(in) :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a1i8=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a1i8=vals)
      endif

    end subroutine ncdf_write_var_1d_i8

    !-----------------------------------------------------------------------------
    !--- Write 1D real*4 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_1d_r4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      real(kind=4),    intent(in) :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a1r4=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a1r4=vals)
      endif

    end subroutine ncdf_write_var_1d_r4

    !-----------------------------------------------------------------------------
    !--- Write 1D real*8 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_1d_r8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      real(kind=8),    intent(in) :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a1r8=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a1r8=vals)
      endif

    end subroutine ncdf_write_var_1d_r8

    !-----------------------------------------------------------------------------
    !--- Write 2D integer*4 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_2d_i4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      integer(kind=4), intent(in) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a2i4=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a2i4=vals)
      endif

    end subroutine ncdf_write_var_2d_i4

    !-----------------------------------------------------------------------------
    !--- Write 2D integer*8 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_2d_i8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      integer(kind=8), intent(in) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a2i8=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a2i8=vals)
      endif

    end subroutine ncdf_write_var_2d_i8

    !-----------------------------------------------------------------------------
    !--- Write 2D real*4 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_2d_r4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      real(kind=4),    intent(in) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a2r4=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a2r4=vals)
      endif

    end subroutine ncdf_write_var_2d_r4

    !-----------------------------------------------------------------------------
    !--- Write 2D real*8 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_2d_r8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      real(kind=8),    intent(in) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a2r8=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a2r8=vals)
      endif

    end subroutine ncdf_write_var_2d_r8

    !-----------------------------------------------------------------------------
    !--- Write 3D integer*4 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_3d_i4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      integer(kind=4), intent(in) :: vals(:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a3i4=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a3i4=vals)
      endif

    end subroutine ncdf_write_var_3d_i4

    !-----------------------------------------------------------------------------
    !--- Write 3D integer*8 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_3d_i8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      integer(kind=8), intent(in) :: vals(:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a3i8=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a3i8=vals)
      endif

    end subroutine ncdf_write_var_3d_i8

    !-----------------------------------------------------------------------------
    !--- Write 3D real*4 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_3d_r4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      real(kind=4),    intent(in) :: vals(:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a3r4=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a3r4=vals)
      endif

    end subroutine ncdf_write_var_3d_r4

    !-----------------------------------------------------------------------------
    !--- Write 3D real*8 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_3d_r8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      real(kind=8),    intent(in) :: vals(:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a3r8=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a3r8=vals)
      endif

    end subroutine ncdf_write_var_3d_r8

    !-----------------------------------------------------------------------------
    !--- Write 4D integer*4 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_4d_i4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      integer(kind=4), intent(in) :: vals(:,:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a4i4=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a4i4=vals)
      endif

    end subroutine ncdf_write_var_4d_i4

    !-----------------------------------------------------------------------------
    !--- Write 4D integer*8 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_4d_i8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      integer(kind=8), intent(in) :: vals(:,:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a4i8=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a4i8=vals)
      endif

    end subroutine ncdf_write_var_4d_i8

    !-----------------------------------------------------------------------------
    !--- Write 4D real*4 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_4d_r4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      real(kind=4),    intent(in) :: vals(:,:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a4r4=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a4r4=vals)
      endif

    end subroutine ncdf_write_var_4d_r4

    !-----------------------------------------------------------------------------
    !--- Write 4D real*8 values to a given variable in an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_write_var_4d_r8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be written to the file
      real(kind=8),    intent(in) :: vals(:,:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_write_var_any(fid, varid, nt=nt, val_a4r8=vals)
      else
        call ncdf_write_var_any(fid, varid, val_a4r8=vals)
      endif

    end subroutine ncdf_write_var_4d_r8

    !=============================================================================
    ! Read a variable from an open netcdf file
    !=============================================================================
    subroutine ncdf_read_var_any(fid, varid, nt, &
                                  val_a1i4, val_a1i8, val_a1r4, val_a1r8, &
                                  val_a2i4, val_a2i8, val_a2r4, val_a2r8, &
                                  val_a3i4, val_a3i8, val_a3r4, val_a3r8, &
                                  val_a4i4, val_a4i8, val_a4r4, val_a4r8, &
                                  val_a1char)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- The variable values to be read from the file
      integer(4),   pointer, intent(inout), optional :: val_a1i4(:)
      integer(8),   pointer, intent(inout), optional :: val_a1i8(:)
      real(4),      pointer, intent(inout), optional :: val_a1r4(:)
      real(8),      pointer, intent(inout), optional :: val_a1r8(:)
      integer(4),   pointer, intent(inout), optional :: val_a2i4(:,:)
      integer(8),   pointer, intent(inout), optional :: val_a2i8(:,:)
      real(4),      pointer, intent(inout), optional :: val_a2r4(:,:)
      real(8),      pointer, intent(inout), optional :: val_a2r8(:,:)
      integer(4),   pointer, intent(inout), optional :: val_a3i4(:,:,:)
      integer(8),   pointer, intent(inout), optional :: val_a3i8(:,:,:)
      real(4),      pointer, intent(inout), optional :: val_a3r4(:,:,:)
      real(8),      pointer, intent(inout), optional :: val_a3r8(:,:,:)
      integer(4),   pointer, intent(inout), optional :: val_a4i4(:,:,:,:)
      integer(8),   pointer, intent(inout), optional :: val_a4i8(:,:,:,:)
      real(4),      pointer, intent(inout), optional :: val_a4r4(:,:,:,:)
      real(8),      pointer, intent(inout), optional :: val_a4r8(:,:,:,:)
      character(*), pointer, intent(inout), optional :: val_a1char(:)

      !--- Local
      integer :: vtype, verbose=0
      integer(kind=4) :: ncstat, idx, dimlen, dimlen_in, rec_dimid
      integer(kind=4) :: ndims, dimids(7), read_ndims
      integer(kind=4), pointer :: vshape(:) => null(), vstart(:) => null(), vcount(:) => null()
      real(8) :: time_val
      character(256) :: msg = " "
      character(64)  :: vname, dim_name, dname_strng
      character(64)  :: vshape_strng, vstart_strng, vcount_strng, vtype_strng
      logical :: use_var(17), has_rec_dim, problem

      !--- Identify the incomming variable type
      use_var(:) = .false.
      if ( present(val_a1i4)   ) use_var( 1) = .true.
      if ( present(val_a1i8)   ) use_var( 2) = .true.
      if ( present(val_a1r4)   ) use_var( 3) = .true.
      if ( present(val_a1r8)   ) use_var( 4) = .true.
      if ( present(val_a2i4)   ) use_var( 5) = .true.
      if ( present(val_a2i8)   ) use_var( 6) = .true.
      if ( present(val_a2r4)   ) use_var( 7) = .true.
      if ( present(val_a2r8)   ) use_var( 8) = .true.
      if ( present(val_a3i4)   ) use_var( 9) = .true.
      if ( present(val_a3i8)   ) use_var(10) = .true.
      if ( present(val_a3r4)   ) use_var(11) = .true.
      if ( present(val_a3r8)   ) use_var(12) = .true.
      if ( present(val_a4i4)   ) use_var(13) = .true.
      if ( present(val_a4i8)   ) use_var(14) = .true.
      if ( present(val_a4r4)   ) use_var(15) = .true.
      if ( present(val_a4r8)   ) use_var(16) = .true.
      if ( present(val_a1char) ) use_var(17) = .true.

      !--- Ensure that there is exactly one incomming variable
      idx = count(use_var)
      if ( idx == 0 ) then
        write(6,*)"ncdf_read_var_any: User supplied data is missing."
        call xit("ncdf_read_var_any",-1)
      endif
      if ( idx > 1 ) then
        write(6,*)"ncdf_read_var_any: Only one input variable array is allowed."
        call xit("ncdf_read_var_any",-2)
      endif

      !--- Determine the index in use_var of the input data type
      vtype = 0
      do idx=1,size(use_var)
        if ( use_var(idx) ) then
          vtype = idx
          !--- There is only 1 variable allowed so we exit the loop
          exit
        endif
      enddo
      if (vtype == 0) then
        write(6,*)"ncdf_read_var_any: Unable to determine input variable type and kind."
        call xit("ncdf_read_var_any",-3)
      endif

      !--- Determine the record dimension ID
      ncstat = nf90_Inquire(fid, unlimitedDimId=rec_dimid)
      msg = "ncdf_read_var_any: nf90_Inquire(fid, "
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      !--- Get name and dimension IDs for the current variable
      ncstat = nf90_Inquire_Variable(fid, varid, name=vname, ndims=ndims, dimids=dimids)
      msg = "ncdf_read_var_any: nf90_Inquire_Variable(fid, "//trim(vname)
      call nc_error_handler(fid, ncstat, 2, msg=msg)

      !--- Does the netcdf variable contain the record dimension
      if ( dimids(ndims) == rec_dimid ) then
        !--- The record dimension must always be the slowest varying (ie dimids(ndims))
        has_rec_dim = .true.
      else
        has_rec_dim = .false.
      endif

      !--- Allocate and assign the shape and slice arrays
      dname_strng = "dims="
      allocate( vshape(ndims), vstart(ndims), vcount(ndims) )
      do idx=1,ndims
        ncstat = nf90_Inquire_Dimension(fid, dimids(idx), name=dim_name, len=vshape(idx))
        msg = "ncdf_read_var_any: nf90_Inquire_Dimension(fid, "//trim(vname)
        call nc_error_handler(fid, ncstat, 1, msg=msg)
        if ( len_trim(dname_strng) < 1 ) then
          dname_strng = trim(dim_name)
        else
          dname_strng = trim(dname_strng)//" "//trim(dim_name)
        endif
      enddo
      vstart = 1
      vcount = vshape

      read_ndims = ndims
      if ( has_rec_dim ) then
        if ( present(nt) ) then
          !--- Read a single time slice at the time index supplied by the user
          if ( nt < 1 .or. nt > vshape(ndims) ) then
            write(6,*)"ncdf_read_var_any: nt=",nt," is out of range."
            write(6,*)"    valid range is 1 to ",vshape(ndims)
            call xit("ncdf_read_var_any",-4)
          endif
          vstart(ndims) = nt
          vcount(ndims) = 1
          if ( ndims > 1 ) then
            !--- The variable returned to the caller will not contain the record dimension
            read_ndims = ndims - 1
          endif
        else
          if ( ndims == 1 ) then
            !--- Time is the only dimension, in this case read all values
            vstart(ndims) = 1
            vcount(ndims) = vshape(ndims)
          else
            !--- Read the last time step
            vstart(ndims) = vshape(ndims)
            vcount(ndims) = 1
            !--- The variable returned to the caller will not contain the record dimension
            read_ndims = ndims - 1
          endif
        endif
      endif

      !--- Assign char variables containing info to be used for output below
      vtype_strng = " "
      vshape_strng = " "
      vstart_strng = " "
      vcount_strng = " "
      write(vshape_strng,'("shape=",7i5)')vshape
      write(vstart_strng,'("start=",7i5)')vstart
      write(vcount_strng,'("count=",7i5)')vcount

      !--- Read the variable
      select case(vtype)
        case ( 1)  !--- 1D array integer*4
          vtype_strng = "type=1D_i4"
          if ( read_ndims /= 1 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    1D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-70)
          endif
          if ( associated(val_a1i4) ) deallocate(val_a1i4)
          allocate( val_a1i4( vshape(1) ) )
          ncstat = nf90_get_var(fid, varid, val_a1i4, start=vstart, count=vcount)
        case ( 2)  !--- 1D array integer*8
          vtype_strng = "type=1D_i8"
          if ( read_ndims /= 1 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    1D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-71)
          endif
          if ( associated(val_a1i8) ) deallocate(val_a1i8)
          allocate( val_a1i8( vshape(1) ) )
          ncstat = nf90_get_var(fid, varid, val_a1i8, start=vstart, count=vcount)
        case ( 3)  !--- 1D array real*4
          vtype_strng = "type=1D_r4"
          if ( read_ndims /= 1 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    1D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-72)
          endif
          if ( associated(val_a1r4) ) deallocate(val_a1r4)
          allocate( val_a1r4( vshape(1) ) )
          ncstat = nf90_get_var(fid, varid, val_a1r4, start=vstart, count=vcount)
        case ( 4)  !--- 1D array real*8
          vtype_strng = "type=1D_r8"
          if ( read_ndims /= 1 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    1D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-73)
          endif
          if ( associated(val_a1r8) ) deallocate(val_a1r8)
          allocate( val_a1r8( vshape(1) ) )
          ncstat = nf90_get_var(fid, varid, val_a1r8, start=vstart, count=vcount)
        case ( 5)  !--- 2D array integer*4
          vtype_strng = "type=2D_i4"
          if ( read_ndims /= 2 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    2D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-74)
          endif
          if ( associated(val_a2i4) ) deallocate(val_a2i4)
          allocate( val_a2i4( vshape(1), vshape(2) ) )
          ncstat = nf90_get_var(fid, varid, val_a2i4, start=vstart, count=vcount)
        case ( 6)  !--- 2D array integer*8
          vtype_strng = "type=2D_i8"
          if ( read_ndims /= 2 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    2D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-75)
          endif
          if ( associated(val_a2i8) ) deallocate(val_a2i8)
          allocate( val_a2i8( vshape(1), vshape(2) ) )
          ncstat = nf90_get_var(fid, varid, val_a2i8, start=vstart, count=vcount)
        case ( 7)  !--- 2D array real*4
          vtype_strng = "type=2D_r4"
          if ( read_ndims /= 2 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    2D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-76)
          endif
          if ( associated(val_a2r4) ) deallocate(val_a2r4)
          allocate( val_a2r4( vshape(1), vshape(2) ) )
          ncstat = nf90_get_var(fid, varid, val_a2r4, start=vstart, count=vcount)
        case ( 8)  !--- 2D array real*8
          vtype_strng = "type=2D_r8"
          if ( read_ndims /= 2 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    2D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-77)
          endif
          if ( associated(val_a2r8) ) deallocate(val_a2r8)
          allocate( val_a2r8( vshape(1), vshape(2) ) )
          ncstat = nf90_get_var(fid, varid, val_a2r8, start=vstart, count=vcount)
        case ( 9)  !--- 3D array integer*4
          vtype_strng = "type=3D_i4"
          if ( read_ndims /= 3 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    3D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-78)
          endif
          if ( associated(val_a3i4) ) deallocate(val_a3i4)
          allocate( val_a3i4( vshape(1), vshape(2), vshape(3) ) )
          ncstat = nf90_get_var(fid, varid, val_a3i4, start=vstart, count=vcount)
        case (10)  !--- 3D array integer*8
          vtype_strng = "type=3D_i8"
          if ( read_ndims /= 3 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    3D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-79)
          endif
          if ( associated(val_a3i8) ) deallocate(val_a3i8)
          allocate( val_a3i8( vshape(1), vshape(2), vshape(3) ) )
          ncstat = nf90_get_var(fid, varid, val_a3i8, start=vstart, count=vcount)
        case (11)  !--- 3D array real*4
          vtype_strng = "type=3D_r4"
          if ( read_ndims /= 3 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    3D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-80)
          endif
          if ( associated(val_a3r4) ) deallocate(val_a3r4)
          allocate( val_a3r4( vshape(1), vshape(2), vshape(3) ) )
          ncstat = nf90_get_var(fid, varid, val_a3r4, start=vstart, count=vcount)
        case (12)  !--- 3D array real*8
          vtype_strng = "type=3D_r8"
          if ( read_ndims /= 3 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    3D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-81)
          endif
          if ( associated(val_a3r8) ) deallocate(val_a3r8)
          allocate( val_a3r8( vshape(1), vshape(2), vshape(3) ) )
          ncstat = nf90_get_var(fid, varid, val_a3r8, start=vstart, count=vcount)
        case (13)  !--- 4D array integer*4
          vtype_strng = "type=4D_i4"
          if ( read_ndims /= 4 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    4D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-82)
          endif
          if ( associated(val_a4i4) ) deallocate(val_a4i4)
          allocate( val_a4i4( vshape(1), vshape(2), vshape(3), vshape(4) ) )
          ncstat = nf90_get_var(fid, varid, val_a4i4, start=vstart, count=vcount)
        case (14)  !--- 4D array integer*8
          vtype_strng = "type=4D_i8"
          if ( read_ndims /= 4 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    4D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-83)
          endif
          if ( associated(val_a4i8) ) deallocate(val_a4i8)
          allocate( val_a4i8( vshape(1), vshape(2), vshape(3), vshape(4) ) )
          ncstat = nf90_get_var(fid, varid, val_a4i8, start=vstart, count=vcount)
        case (15)  !--- 4D array real*4
          vtype_strng = "type=4D_r4"
          if ( read_ndims /= 4 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    4D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-84)
          endif
          if ( associated(val_a4r4) ) deallocate(val_a4r4)
          allocate( val_a4r4( vshape(1), vshape(2), vshape(3), vshape(4) ) )
          ncstat = nf90_get_var(fid, varid, val_a4r4, start=vstart, count=vcount)
        case (16)  !--- 4D array real*8
          vtype_strng = "type=4D_r8"
          if ( read_ndims /= 4 ) then
            write(6,*)"ncdf_read_var_any: Dimension size mismatch."
            write(6,*)"    4D input variable but found ",read_ndims," dims in file"
            call xit("ncdf_read_var_any",-85)
          endif
          if ( associated(val_a4r8) ) deallocate(val_a4r8)
          allocate( val_a4r8( vshape(1), vshape(2), vshape(3), vshape(4) ) )
          ncstat = nf90_get_var(fid, varid, val_a4r8, start=vstart, count=vcount)
        case (17)  !--- 1D array character
          vtype_strng = "type=1D_ch"
          write(6,*)"ncdf_read_var_any: Reading a 1D character array is not yet supported."
          call xit("ncdf_read_var_any",-90)
          ncstat = nf90_get_var(fid, varid, val_a1char, start=vstart, count=vcount)
      end select
      msg = "ncdf_read_var_any: nf90_get_var(...)  var="//trim(vname)//" "//trim(vtype_strng) &
            //" "//trim(vshape_strng)//" "//trim(vstart_strng)//" "//trim(vcount_strng) &
            //" "//trim(dname_strng)
      call nc_error_handler(fid, ncstat, 3, msg=msg)

      if ( verbose > 1 ) write(6,'(a)')trim(msg)

      !--- Clean up
      deallocate( vshape, vstart, vcount )

    end subroutine ncdf_read_var_any

    !-----------------------------------------------------------------------------
    !--- Read 1D integer*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_1d_i4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      integer(kind=4), pointer, intent(inout) :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a1i4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a1i4=vals)
      endif

    end subroutine ncdf_read_var_1d_i4

    !-----------------------------------------------------------------------------
    !--- Read 1D integer*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_1d_i4(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      integer(kind=4), pointer, intent(inout) :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_1d_i4: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a1i4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a1i4=vals)
      endif

    end subroutine ncdf_read_var_byname_1d_i4

    !-----------------------------------------------------------------------------
    !--- Read 1D real*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_1d_r4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      real(kind=4), pointer, intent(inout) :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a1r4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a1r4=vals)
      endif

    end subroutine ncdf_read_var_1d_r4

    !-----------------------------------------------------------------------------
    !--- Read 1D real*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_1d_r4(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      real(kind=4), pointer, intent(inout) :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_1d_r4: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a1r4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a1r4=vals)
      endif

    end subroutine ncdf_read_var_byname_1d_r4

    !-----------------------------------------------------------------------------
    !--- Read 1D real*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_1d_r8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      real(kind=8), pointer, intent(inout) :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a1r8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a1r8=vals)
      endif

    end subroutine ncdf_read_var_1d_r8

    !-----------------------------------------------------------------------------
    !--- Read 1D real*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_1d_r8(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      real(kind=8), pointer, intent(inout) :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_1d_r8: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a1r8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a1r8=vals)
      endif

    end subroutine ncdf_read_var_byname_1d_r8

    !-----------------------------------------------------------------------------
    !--- Read 1D integer*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_1d_i8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      integer(kind=8), pointer, intent(inout) :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a1i8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a1i8=vals)
      endif

    end subroutine ncdf_read_var_1d_i8

    !-----------------------------------------------------------------------------
    !--- Read 1D integer*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_1d_i8(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      integer(kind=8), pointer, intent(inout) :: vals(:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_1d_i8: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a1i8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a1i8=vals)
      endif

    end subroutine ncdf_read_var_byname_1d_i8

    !-----------------------------------------------------------------------------
    !--- Read 2D integer*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_2d_i4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      integer(kind=4), pointer, intent(inout) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a2i4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a2i4=vals)
      endif

    end subroutine ncdf_read_var_2d_i4

    !-----------------------------------------------------------------------------
    !--- Read 2D integer*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_2d_i4(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      integer(kind=4), pointer, intent(inout) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_2d_i4: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a2i4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a2i4=vals)
      endif

    end subroutine ncdf_read_var_byname_2d_i4

    !-----------------------------------------------------------------------------
    !--- Read 2D real*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_2d_r4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      real(kind=4), pointer, intent(inout) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a2r4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a2r4=vals)
      endif

    end subroutine ncdf_read_var_2d_r4

    !-----------------------------------------------------------------------------
    !--- Read 2D real*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_2d_r4(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      real(kind=4), pointer, intent(inout) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_2d_r4: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a2r4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a2r4=vals)
      endif

    end subroutine ncdf_read_var_byname_2d_r4

    !-----------------------------------------------------------------------------
    !--- Read 2D integer*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_2d_i8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      integer(kind=8), pointer, intent(inout) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a2i8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a2i8=vals)
      endif

    end subroutine ncdf_read_var_2d_i8

    !-----------------------------------------------------------------------------
    !--- Read 2D integer*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_2d_i8(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      integer(kind=8), pointer, intent(inout) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_2d_i8: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a2i8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a2i8=vals)
      endif

    end subroutine ncdf_read_var_byname_2d_i8

    !-----------------------------------------------------------------------------
    !--- Read 2D real*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_2d_r8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      real(kind=8), pointer, intent(inout) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a2r8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a2r8=vals)
      endif

    end subroutine ncdf_read_var_2d_r8

    !-----------------------------------------------------------------------------
    !--- Read 2D real*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_2d_r8(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      real(kind=8), pointer, intent(inout) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_2d_r8: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a2r8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a2r8=vals)
      endif

    end subroutine ncdf_read_var_byname_2d_r8

    !-----------------------------------------------------------------------------
    !--- Read 2D real*8 values from a netcdf file
    !--- This routine will open the file, read the data, close the file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_by_file_and_name_2d_r8(fname, vname, vals, nt)
      !--- The netcdf file name
      character(*), intent(in) :: fname

      !--- The netcdf variable name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      real(kind=8), pointer, intent(inout) :: vals(:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: fid
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Open the file in read only mode
      call ncdf_open(fid, trim(fname), mode="read")

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_by_file_and_name_2d_r8: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a2r8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a2r8=vals)
      endif

      !--- Close the file before return
      call ncdf_close(fid)

    end subroutine ncdf_read_var_by_file_and_name_2d_r8

    !-----------------------------------------------------------------------------
    !--- Read 3D integer*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_3d_i4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      integer(kind=4), pointer, intent(inout) :: vals(:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a3i4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a3i4=vals)
      endif

    end subroutine ncdf_read_var_3d_i4

    !-----------------------------------------------------------------------------
    !--- Read 3D integer*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_3d_i4(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      integer(kind=4), pointer, intent(inout) :: vals(:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_3d_i4: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a3i4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a3i4=vals)
      endif

    end subroutine ncdf_read_var_byname_3d_i4

    !-----------------------------------------------------------------------------
    !--- Read 3D real*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_3d_r4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      real(kind=4), pointer, intent(inout) :: vals(:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a3r4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a3r4=vals)
      endif

    end subroutine ncdf_read_var_3d_r4

    !-----------------------------------------------------------------------------
    !--- Read 3D real*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_3d_r4(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      real(kind=4), pointer, intent(inout) :: vals(:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_3d_r4: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a3r4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a3r4=vals)
      endif

    end subroutine ncdf_read_var_byname_3d_r4

    !-----------------------------------------------------------------------------
    !--- Read 3D integer*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_3d_i8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      integer(kind=8), pointer, intent(inout) :: vals(:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a3i8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a3i8=vals)
      endif

    end subroutine ncdf_read_var_3d_i8

    !-----------------------------------------------------------------------------
    !--- Read 3D integer*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_3d_i8(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      integer(kind=8), pointer, intent(inout) :: vals(:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_3d_i8: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a3i8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a3i8=vals)
      endif

    end subroutine ncdf_read_var_byname_3d_i8

    !-----------------------------------------------------------------------------
    !--- Read 3D real*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_3d_r8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      real(kind=8), pointer, intent(inout) :: vals(:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a3r8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a3r8=vals)
      endif

    end subroutine ncdf_read_var_3d_r8

    !-----------------------------------------------------------------------------
    !--- Read 3D real*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_3d_r8(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      real(kind=8), pointer, intent(inout) :: vals(:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_3d_r8: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a3r8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a3r8=vals)
      endif

    end subroutine ncdf_read_var_byname_3d_r8

    !-----------------------------------------------------------------------------
    !--- Read 4D integer*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_4d_i4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      integer(kind=4), pointer, intent(inout) :: vals(:,:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a4i4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a4i4=vals)
      endif

    end subroutine ncdf_read_var_4d_i4

    !-----------------------------------------------------------------------------
    !--- Read 4D integer*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_4d_i4(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      integer(kind=4), pointer, intent(inout) :: vals(:,:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_4d_i4: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a4i4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a4i4=vals)
      endif

    end subroutine ncdf_read_var_byname_4d_i4

    !-----------------------------------------------------------------------------
    !--- Read 4D real*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_4d_r4(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      real(kind=4), pointer, intent(inout) :: vals(:,:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a4r4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a4r4=vals)
      endif

    end subroutine ncdf_read_var_4d_r4

    !-----------------------------------------------------------------------------
    !--- Read 4D real*4 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_4d_r4(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      real(kind=4), pointer, intent(inout) :: vals(:,:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_4d_r4: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a4r4=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a4r4=vals)
      endif

    end subroutine ncdf_read_var_byname_4d_r4

    !-----------------------------------------------------------------------------
    !--- Read 4D integer*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_4d_i8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      integer(kind=8), pointer, intent(inout) :: vals(:,:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a4i8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a4i8=vals)
      endif

    end subroutine ncdf_read_var_4d_i8

    !-----------------------------------------------------------------------------
    !--- Read 4D integer*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_4d_i8(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      integer(kind=8), pointer, intent(inout) :: vals(:,:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_4d_i8: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a4i8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a4i8=vals)
      endif

    end subroutine ncdf_read_var_byname_4d_i8

    !-----------------------------------------------------------------------------
    !--- Read 4D real*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_4d_r8(fid, varid, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf variable ID
      integer(kind=4), intent(in) :: varid

      !--- The variable values to be read
      real(kind=8), pointer, intent(inout) :: vals(:,:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a4r8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a4r8=vals)
      endif

    end subroutine ncdf_read_var_4d_r8

    !-----------------------------------------------------------------------------
    !--- Read 4D real*8 values from an open netcdf file
    !-----------------------------------------------------------------------------
    subroutine ncdf_read_var_byname_4d_r8(fid, vname, vals, nt)
      !--- The netcdf file ID
      integer(kind=4), intent(in) :: fid

      !--- The netcdf file name
      character(*), intent(in) :: vname

      !--- The variable values to be read
      real(kind=8), pointer, intent(inout) :: vals(:,:,:,:)

      !--- The index in the record dimension at which the data is written
      integer, intent(in), optional :: nt

      !--- Local
      !--- The netcdf variable ID
      integer(kind=4) :: varid, ncstat
      character(64) :: msg

      !--- Get the variable ID for the current name
      ncstat = nf90_inq_varid(fid, vname, varid)
      msg = "ncdf_read_var_byname_4d_r8: nf90_inq_varid(...)  var="//trim(vname)
      call nc_error_handler(fid, ncstat, 1, msg=msg)

      if ( present(nt) ) then
        call ncdf_read_var_any(fid, varid, nt=nt, val_a4r8=vals)
      else
        call ncdf_read_var_any(fid, varid, val_a4r8=vals)
      endif

    end subroutine ncdf_read_var_byname_4d_r8

end module ncdf
