! This program is for testing ESMF based remapping of NEMO4-type fields
! 
! Based on the remaping exercices made by NCS 2022-12-02

program remap
    use ESMF
    use ncdf
    use iso_fortran_env, only:  int8, int16, int32, int64, real32, real64
    use esmf_utils, only : define_esmf_grid, define_esmf_field, define_yg_grid_from_file, define_orca_grid_from_file, define_core_grid_from_file
    use file_readers, only: read_yg, read_NEMO_mesh_mask_file
    use spreading, only: runoff_push
    implicit none
    character(64) :: esmf_error_txt
    integer, parameter :: sp        = real32
    integer, parameter :: dp        = real64
    integer :: rc, localrc
    type(ESMF_grid)  :: src_grid, dst_grid
    type(ESMF_field) :: srcfield_data, fld_lfno_dst
    type(ESMF_field) :: area_field, dstfield_data
    type(ESMF_routeHandle) erh_dst_to_src, erh_src_to_dst
    type(ESMF_RegridMethod_Flag) regrid_methods_opt
    type(ESMF_NormType_Flag) normType_opt

    ! variables to be read from file
    real(ESMF_KIND_R8), dimension(:, :), pointer   :: srcdata_array=> null()
    real(ESMF_KIND_R8), dimension(:, :), pointer   :: srcdata_area=> null()
    integer(ESMF_KIND_I4), dimension(:, :), pointer   :: srcdata_mask=> null()
    real(kind=8), dimension(: ), pointer   :: srctime_vals=> null()
    integer srctime_nt
    character(200) tmp_char
    ! remaped grid variables
    real(ESMF_KIND_R8), pointer   ::  dstdata_array(:,:) => null()
    real(ESMF_KIND_R8), pointer   ::  dstdata_lat(:,:) => null()
    real(ESMF_KIND_R8), pointer   ::  dstdata_lon(:,:) => null()
    real(ESMF_KIND_R8), allocatable, dimension(:,:) ::  dstdata_array_out
    real(ESMF_KIND_R8), pointer   ::  dstdata_area(:,:)=> null()
    integer(ESMF_KIND_I4), pointer   ::  dstdata_mask(:,:)=> null()
    integer(ESMF_KIND_I4), allocatable, dimension(:,:) ::  mask_valid(:,:)
    integer(ESMF_KIND_I4), dimension(:,:), pointer :: rivers_mask => null()

    ! For NC output
    integer(kind=4) :: i, dst_fid, xid, yid, xdim, ydim, tid, tdim, &
    varid_dst, &
    varid_dst_lon,  & 
    varid_dst_lat, &
    varid_dstfrc_orca_dst, &
    ncstat,varid_src,  & 
    srcgrid_dims(2),dstgrid_dims(2)
    integer(kind=4), ALLOCATABLE, DIMENSION(:) :: xi, yi, xo,yo

    ! Namelist inputs
    character(300) :: srcgrd_file,dstgrd_file,srcdata_file,dstdata_file
    character(30)  :: methods,norm_type,srcgrd_name,srcgrd_type,dstgrd_name,dstgrd_type,srcdata_field,dstdata_field
    logical push_runoff_to_sea, use_landmask,srcgrd_mask,dstgrd_mask
    namelist /nl_remap/ methods,norm_type,push_runoff_to_sea
    namelist /nl_src_grid/ srcgrd_file,srcgrd_name,srcgrd_type,srcgrd_mask
    namelist /nl_dst_grid/ dstgrd_file,dstgrd_name,dstgrd_type,dstgrd_mask
    namelist /nl_src_fields/ srcdata_file,srcdata_field
    namelist /nl_dst_fields/ dstdata_file,dstdata_field
    integer, parameter :: iunit=87
    integer nt,src_fid
    logical :: first_step = .true.

    call ESMF_Initialize(rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! get the remaping options from the namelist 
    write(6,*)"Remmaping methods:"
    open(iunit, file="namelist_remap", status='old', form='formatted')
    read(iunit, nml=nl_remap) 
    close(iunit)
    select case (methods)
    case ("conservative");      regrid_methods_opt=ESMF_REGRIDMETHOD_CONSERVE
    case ("conservative_2nd");  regrid_methods_opt=ESMF_REGRIDMETHOD_CONSERVE_2ND
    case ("bilinear");          regrid_methods_opt=ESMF_REGRIDMETHOD_BILINEAR
    case ("patch");             regrid_methods_opt=ESMF_REGRIDMETHOD_PATCH
    case ("nearest_stop");      regrid_methods_opt=ESMF_REGRIDMETHOD_NEAREST_STOD
    case ("nearest_dtos");      regrid_methods_opt=ESMF_REGRIDMETHOD_NEAREST_DTOS
    case default; 
        write(*,*) methods,' is not an opion (conservative,conservative_2nd,bilinear,patch,nearest_stop,nearest_dtos)'
        stop
    end select
    select case (norm_type)
    case ("fracarea");  normType_opt=ESMF_NORMTYPE_FRACAREA
    case ("dstarea");   normType_opt=ESMF_NORMTYPE_DSTAREA
    case default; 
        write(*,*) norm_type,' is not an opion (fracarea,dstarea)'
        stop
    end select
    write(6,*)"  ",methods,',',norm_type
    

    ! get src grid info 
    write(6,*)"Source grid info:"
    open(iunit, file="namelist_remap", status='old', form='formatted')
    read(iunit, nml=nl_src_grid) 
    close(iunit)
    use_landmask=srcgrd_mask

    select case (srcgrd_type)
    case ("core")
        ! define CanNEMO like grid
        call define_core_grid_from_file(srcgrd_file, srcgrd_name, src_grid, grid_dims=srcgrid_dims,use_mask=use_landmask)
    case ("nemo")
        ! define CanNEMO like grid
        call define_orca_grid_from_file(srcgrd_file, srcgrd_name, src_grid, grid_dims=srcgrid_dims,use_mask=use_landmask)
    case ("yg")
        ! Define the yg grid
        call define_yg_grid_from_file(srcgrd_file, srcgrd_name, esmf_yggrid=dst_grid)
    case default
        write(*,*)  srcgrd_type," is not a valid option (nemo or yg)"
        stop
    end select

    ! get dst grid info 
    write(6,*)"Destination grid info:"
    open(iunit, file="namelist_remap", status='old', form='formatted')
    read(iunit, nml=nl_dst_grid) 
    close(iunit)
    use_landmask=dstgrd_mask

    select case (dstgrd_type)
    case ("core")
        ! define CanNEMO like grid
        call define_core_grid_from_file(srcgrd_file, srcgrd_name, src_grid, grid_dims=srcgrid_dims,use_mask=use_landmask)
    case ("nemo")
        ! define CanNEMO like grid
        call define_orca_grid_from_file(dstgrd_file, dstgrd_name, dst_grid, grid_dims=dstgrid_dims,use_mask=use_landmask)
    case ("yg")
        ! Define the yg grid
        call define_yg_grid_from_file(dstgrd_file, dstgrd_name, esmf_yggrid=dst_grid)
    case default
        write(*,*)  dstgrd_type," is not a valid option (nemo or yg)"
        stop
    end select
   
!!==========================================================================
    ! Setup fields and remap   
    ! get dst grid info 
    write(6,*)"Inputs fields info:"
    open(iunit, file="namelist_remap", status='old', form='formatted')
    read(iunit, nml=nl_src_fields) 
    close(iunit)
    write(6,*)"  ",trim(srcdata_field),',',trim(srcdata_file)

    ! open, get the time axis (remaping done on every time step)
    call ncdf_open(src_fid, srcdata_file, "read")
    call ncdf_get_time(srcdata_file,time_vals=srctime_vals,ntime=srctime_nt)

    ! create the outputs file
    write(6,*)"Outputs fields info:"
    open(iunit, file="namelist_remap", status='old', form='formatted')
    read(iunit, nml=nl_dst_fields) 
    close(iunit)
    write(6,*)"  ",trim(dstdata_field),',',trim(dstdata_file)
    call ncdf_open(dst_fid, dstdata_file, 'c')
    allocate(yo(dstgrid_dims(2)),xo(dstgrid_dims(1)))
    do i=1,dstgrid_dims(2)
       yo(i) = i
    end do
    do i=1,dstgrid_dims(1)
        xo(i) = i
    end do
    call ncdf_add_coord(dst_fid, xid, xdim, "x", vals=xo)
    call ncdf_add_coord(dst_fid, yid, ydim, "y", vals=yo)
    call ncdf_add_coord(dst_fid, tid, tdim, "time",record_dim=.true. )
    
    ! add some attribute to the time variables
    call ncdf_copy_all_att(src_fid, "time", dst_fid, "time" )
    call ncdf_write_var(dst_fid, tid, srctime_vals)

    ! put lat/lon in file
    call read_NEMO_mesh_mask_file(dstgrd_file, mid_lon=dstdata_lon, mid_lat=dstdata_lat )
    call ncdf_add_var(dst_fid, varid_dst, 'lon', (/xdim,ydim/))
    call ncdf_write_var(dst_fid, varid_dst, dstdata_lon)

    call ncdf_add_var(dst_fid, varid_dst, 'lat',(/xdim,ydim/))
    call ncdf_write_var(dst_fid, varid_dst, dstdata_lat)

    call ncdf_add_var(dst_fid, varid_dst, dstdata_field, &
                      (/xdim,ydim,tdim/), xtype=NF90_FLOAT)
    call ncdf_copy_all_att(src_fid,  srcdata_field, dst_fid, dstdata_field )

    print*,'Data processing...'
    do nt=1,srctime_nt 

      ! read the iteration of the scr_field
      call ncdf_read_var(src_fid,srcdata_field,srcdata_array,nt)

      ! define fiels for input and ourtputs 
      call define_esmf_field(src_grid, srcdata_array, srcfield_data)
      call define_esmf_field(dst_grid, field=dstfield_data)

      ! Precompute a Field regridding operation and return a RouteHandle (first time only)
      if (first_step) then
        call ESMF_FieldRegridStore(srcField=srcfield_data, &
                            dstField=dstfield_data, &
                            dstMaskValues=(/0/), &
                            regridmethod=regrid_methods_opt, &
                            normType=normType_opt, &
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            ignoreDegenerate=.true., &
                            routeHandle=erh_dst_to_src,  &
                            rc=rc)
        
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
        first_step=.false.
        
        ! get areas for the sum calculation ( only the first iteration)
        call ESMF_GridGetItem(src_grid,ESMF_GRIDITEM_AREA,farrayPtr=srcdata_area,rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
        call ESMF_GridGetItem(src_grid,ESMF_GRIDITEM_MASK,farrayPtr=srcdata_mask,rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
        call ESMF_GridGetItem(dst_grid,ESMF_GRIDITEM_AREA,farrayPtr=dstdata_area,rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
        call ESMF_GridGetItem(dst_grid,ESMF_GRIDITEM_MASK,farrayPtr=dstdata_mask,rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
        allocate(mask_valid(dstgrid_dims(1),dstgrid_dims(2)))
        if (push_runoff_to_sea) then
          ! only get the landmask  
          call read_NEMO_mesh_mask_file(dstgrd_file,  imask=rivers_mask )
          ! replace the dstdata_mask (that are all ones) by the rivers_mask
          mask_valid=rivers_mask
          allocate(dstdata_array_out(dstgrid_dims(1),dstgrid_dims(2)))
        else
          mask_valid=dstdata_mask
        endif

        write(*,*) 'Global sum calculation (*1e6)'
        write(*,'(A5,A15,A15,A15)'),'nt','source','destination','Difference(%)'

      endif

      ! Do the regrid 
      call ESMF_FieldRegrid(srcfield_data, dstfield_data, &
          routeHandle=erh_dst_to_src, zeroregion=ESMF_REGION_TOTAL,rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      call ESMF_FieldGet(field=dstfield_data, farrayPtr=dstdata_array, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      if (push_runoff_to_sea) then
    
        ! push runoff to the mask 
        dstdata_array_out=0.
        !call read_NEMO_mesh_mask_file(dstgrd_file,  imask=rivers_mask )
        call runoff_push(dstdata_array, dstdata_lon, dstdata_lat, rivers_mask, dstdata_area, dstdata_array_out)
        dstdata_array=dstdata_array_out
        mask_valid=rivers_mask
      else
        mask_valid=dstdata_mask
      endif

      ! write the iteration 
      call ncdf_write_var(dst_fid, varid_dst, dstdata_array, nt=nt)

      write(*,'(I5,F15.3,F15.3,F15.3)'),nt,sum(srcdata_array*srcdata_area,srcdata_mask==1)/1e6,sum(dstdata_array*dstdata_area,mask_valid==1)/1e6, &
                                        (sum(dstdata_array*dstdata_area,mask_valid==1)-sum(srcdata_array*srcdata_area,srcdata_mask==1))/ & 
                                         sum(srcdata_array*srcdata_area,srcdata_mask==1)*100
    enddo

    call ncdf_close(dst_fid)  
    call ncdf_close(src_fid)  


         !==========================================================================================================
       
    !call ESMF_Finalize()
    print*,'DONE'
    
end program remap

