! Some utility subroutines for easily definining ESMF trypes, includes grids and fields
!
! NCS 2022-12-02
module esmf_utils
    use ESMF
    use ncdf
    use iso_fortran_env, only:  int32, int64, real32, real64
    use file_readers, only: read_yg, read_NEMO_mesh_mask_file, read_CORE_mesh_mask_file
    private
    public ::  define_esmf_grid, define_esmf_field, define_yg_grid_from_file, define_orca_grid_from_file, define_core_grid_from_file
    
    !implicit none
    integer, parameter :: dp = real64 
    integer, parameter :: ip = int32 
    contains

    subroutine define_esmf_grid(gridtype, gridname, lons, lats, lon_bnds, lat_bnds, mask, area, my_esmf_grid)
        implicit none
        integer, parameter :: sp = real32
        integer, parameter :: dp = real64  
        integer, parameter :: ip = int32 
        character(*), intent(in) :: gridtype         ! Type ESM grid: either '1peri' or 'noperi'  
        character(*), intent(in) :: gridname          ! Name of the grid in the ESMF object
        type(ESMF_grid), intent(out) :: my_esmf_grid   ! ESMF grid object
        real(dp), intent(in), dimension(:,:) :: lons  ! incoming lon at grid center
        real(dp), intent(in), dimension(:,:) :: lats  ! incoming lats at grid center
        real(dp), intent(in), dimension(:,:) :: lon_bnds  ! incoming lon at grid corners
        real(dp), intent(in), dimension(:,:) :: lat_bnds  ! incoming lats at grid corners
        real(dp), intent(in), dimension(:,:), optional :: area  ! optional incoming area 
        integer(ESMF_KIND_I4), intent(in), dimension(:,:), optional :: mask  ! optional incoming (integer!) mask 

        real(ESMF_KIND_R8), pointer    :: ptr_grid_mid_lon(:,:) => null()
        real(ESMF_KIND_R8), pointer    :: ptr_grid_mid_lat(:,:) => null()
        real(ESMF_KIND_R8), pointer    :: ptr_grid_corn_lon(:,:) => null()
        real(ESMF_KIND_R8), pointer    :: ptr_grid_corn_lat(:,:) => null()
        real(ESMF_KIND_R8), pointer    :: ptr_grid_area(:,:) => null()
        integer(ESMF_KIND_I4), pointer    :: ptr_grid_mask(:,:) => null()

        integer :: rc
        integer :: nx, ny, gdims(2)

        gdims = shape(lons) ! lons should be a 2D array, lon x lat
        nx = gdims(1)
        ny = gdims(2)
        write(*,*) "  ",gridtype, nx, ny
        select case(gridtype)
            case('noperi')
                ! Create a non perdioc grid object
                my_esmf_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/nx, ny/), &
                        name=gridname, coordSys=ESMF_COORDSYS_SPH_DEG, &
                        rc=rc)  
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
            case('1peri')
                my_esmf_grid = ESMF_GridCreate1PeriDim(maxIndex=(/nx, ny/), &
                    coordSys=ESMF_COORDSYS_SPH_DEG, &
                    indexflag=ESMF_INDEX_GLOBAL,    &
                    name=gridname, rc=rc)
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
            case Default
                call xit("define_esmf_grid: gridtype must be one of '1peri' or 'noperi'",-1)
        end select
        ! Add center and corner coordinate references
        ! Add grid centers
        
        ! Add center storage
        call ESMF_GridAddCoord(my_esmf_grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
        ! Get a pointer array for cell center lons
        call ESMF_GridGetCoord(my_esmf_grid, coordDim=1, localDE=0, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                farrayPtr=ptr_grid_mid_lon, rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
        ptr_grid_mid_lon(:,:) = lons(:,:)
        
        !Get a pointer array for cell center lat 
        call ESMF_GridGetCoord(my_esmf_grid, coordDim=2, localDE=0, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, &
                                farrayPtr=ptr_grid_mid_lat, rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
        ptr_grid_mid_lon = lats(:,:)
        
        ! Add grid corners
        ! Add corner storage
        call ESMF_GridAddCoord(my_esmf_grid, staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
        !--- Retrieve a pointer to the cell corner lon (dim 1) coordinate values
        !--- add assign these values using input data
        call ESMF_GridGetCoord(my_esmf_grid, coordDim=1, localDE=0, &
                               staggerloc=ESMF_STAGGERLOC_CORNER, &
                               farrayPtr=ptr_grid_corn_lon, rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
        ptr_grid_corn_lon(:,:) = lon_bnds(:,:)

        !--- Retrieve a pointer to the cell corner lat (dim 2) coordinate values
        !--- add assign these values using data read from the file
        call ESMF_GridGetCoord(my_esmf_grid, coordDim=2, localDE=0, &
                                staggerloc=ESMF_STAGGERLOC_CORNER, &
                                farrayPtr=ptr_grid_corn_lat, rc=rc)
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
        ptr_grid_corn_lat = lat_bnds(:,:)   
    
        ! Add a mask reference, if it was passed
        if ( present(mask) ) then
            write(*,*) '  Adding mask to ', gridname, minval(mask), maxval(mask)
            call ESMF_GridAddItem(my_esmf_grid, itemflag=ESMF_GRIDITEM_MASK, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
            call ESMF_GridGetItem(my_esmf_grid, &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, &
                                    itemflag=ESMF_GRIDITEM_MASK, &
                                    farrayPtr=ptr_grid_mask, rc=rc)
            ptr_grid_mask(:,:) = mask(:,:)                           
        endif
        ! Add an area reference, if it was passed
        if ( present(area) ) then
            write(*,*) '  Adding area to ', gridname
            call ESMF_GridAddItem(my_esmf_grid, itemflag=ESMF_GRIDITEM_AREA, &
                                staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
            call ESMF_GridGetItem(my_esmf_grid, &
                                    staggerloc=ESMF_STAGGERLOC_CENTER, &
                                    itemflag=ESMF_GRIDITEM_AREA, &
                                    farrayPtr=ptr_grid_area, rc=rc)
            ptr_grid_area(:,:) = area(:,:)                           
        endif
    end subroutine define_esmf_grid

    subroutine define_esmf_field(grid, inarray, field)
        implicit none
        type(ESMF_grid), intent(in) :: grid
        type(ESMF_field), intent(out) :: field
        real(ESMF_KIND_R8), dimension(:,:), intent(in), optional :: inarray 
        real(ESMF_KIND_R8), pointer :: ptr_esmf_var(:,:) => null()
        integer rc

        field = ESMF_FieldCreate(grid,              &             
                staggerloc=ESMF_STAGGERLOC_CENTER,  &
                typekind=ESMF_TYPEKIND_R8,          & 
                rc=rc) 
        if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)      
        if (present(inarray)) then
            call ESMF_FieldGet(field=field, farrayPtr=ptr_esmf_var, & 
                              rc=rc)
            if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
            ptr_esmf_var(:,:) = inarray(:,:)
        endif
    end subroutine define_esmf_field

    subroutine define_yg_grid_from_file(filename, gridname, with_flnd, esmf_yggrid)
        implicit none
        character(*), intent(in) :: filename  
        character(*), intent(in) :: gridname          ! Name of the grid in the ESMF object
        type(ESMF_grid), intent(out) :: esmf_yggrid   ! ESMF grid object
        real(dp), dimension(:,:), pointer  :: lons, lats
        real(dp), dimension(:,:,:), pointer  :: lon_vertices, lat_vertices
        real(dp), dimension(:,:), pointer  :: flnd
        integer(ip), dimension(:,:), allocatable :: imask
        real(dp), dimension(:,:), allocatable :: lon_bnds, lat_bnds
        logical, optional :: with_flnd
        logical :: local_with_flnd
        integer :: nx, ny

        if (present(with_flnd)) then
            local_with_flnd = with_flnd
        else
            local_with_flnd = .false.
        endif

        if(local_with_flnd) then
            ! Read the elements from an input nc file. Normally this file
            ! has been created by fst2nc, including the bounds
             call read_yg(filename, lons, lats, lon_vertices, lat_vertices, flnd)
        else
            call read_yg(filename, lons, lats, lon_vertices, lat_vertices)
        endif
        

        nx = size(lons, dim=1)
        ny = size(lons, dim=2)
        allocate(imask(nx,ny), lon_bnds(nx+1,ny+1), lat_bnds(nx+1,ny+1) )

        ! Build up n+1 corner arrays, from the input vertices.
        lon_bnds(1:nx,1:ny) = lon_vertices(1,:,:) ! inside
        lon_bnds(nx+1,1:ny) = lon_vertices(2,nx,:) ! right column
        lon_bnds(1:nx,ny+1) = lon_vertices(4,:,ny)  ! top row, save right
        lon_bnds(nx+1,ny+1) = lon_vertices(3,nx,ny) ! top right point
        lat_bnds(1:nx,1:ny) = lat_vertices(1,:,:) ! inside
        lat_bnds(nx+1,1:ny) = lat_vertices(2,nx,:) ! right column
        lat_bnds(1:nx,ny+1) = lat_vertices(4,:,ny)  ! top row, save right
        lat_bnds(nx+1,ny+1) = lat_vertices(3,nx,ny) ! top right point

        ! Define a binary/integer mask wherever there is any fraction of land
        if(local_with_flnd) then
            where (flnd>0)
                imask = 1
            elsewhere
                imask = 0
            endwhere
        else
            imask = 1
        endif
        ! Use these inputs to get an ESMF grid object and return it
        call define_esmf_grid('noperi', gridname, lons, lats, lon_bnds, lat_bnds, mask=imask, my_esmf_grid=esmf_yggrid)
    end subroutine define_yg_grid_from_file

    subroutine define_orca_grid_from_file(filename, gridname, esmf_orcagrid,grid_dims,use_mask)
        character(*), intent(in) :: filename  
        character(*), intent(in) :: gridname          ! Name of the grid in the ESMF object
        type(ESMF_grid), intent(out) :: esmf_orcagrid   ! ESMF grid object
        !--- The shape of the grid ie (nx,ny)
        integer, intent(out), optional :: grid_dims(2)
        logical, intent(in), optional :: use_mask
        real(dp), dimension(:,:), pointer :: mid_lon => null()
        real(dp), dimension(:,:), pointer :: mid_lat => null()
        real(dp), dimension(:,:), pointer :: corn_lon => null()
        real(dp), dimension(:,:), pointer :: corn_lat => null()
        real(dp), dimension(:,:), pointer :: tmask => null()
        real(dp), dimension(:,:), pointer :: orca_area => null()
        integer :: orca_grid_dims(2), onx, ony
        integer(ESMF_KIND_I4), dimension(:,:), allocatable :: tmaski
        logical luse_mask
        
        call read_NEMO_mesh_mask_file(filename, grid_dims=orca_grid_dims)
        onx=orca_grid_dims(1)
        ony=orca_grid_dims(2)
        allocate(tmaski(onx, ony))!mid_lon(onx, ony), mid_lat(onx, ony), corn_lon(onx+1, ony+1), corn_lat(onx+1, ony+1), tmask(onx, ony), orca_area(onx, ony), &
                 

        call read_NEMO_mesh_mask_file(filename,  mask=tmask, &
            mid_lon=mid_lon, mid_lat=mid_lat,corn_lon=corn_lon, corn_lat=corn_lat, &
            area=orca_area)
        where (mid_lon<0) mid_lon=mid_lon+360
        where (corn_lon<0) corn_lon=corn_lon+360

        tmaski = int(tmask)
        if ( present(use_mask).and..not.use_mask ) then
          ! set all the mask valule to one if user decide to not use mask. 
          tmaski=1
        endif
        call define_esmf_grid('1peri', gridname, mid_lon, mid_lat, corn_lon, corn_lat, mask=tmaski, area=orca_area, my_esmf_grid=esmf_orcagrid)

        if ( present(grid_dims) ) then
        grid_dims = orca_grid_dims
        endif

    end subroutine define_orca_grid_from_file

    subroutine define_core_grid_from_file(filename, gridname, esmf_coregrid,grid_dims,use_mask)
        character(*), intent(in) :: filename  
        character(*), intent(in) :: gridname          ! Name of the grid in the ESMF object
        type(ESMF_grid), intent(out) :: esmf_coregrid   ! ESMF grid object
        !--- The shape of the grid ie (nx,ny)
        integer, intent(out), optional :: grid_dims(2)
        logical, intent(in), optional :: use_mask
        real(dp), dimension(:,:), pointer :: mid_lon => null()
        real(dp), dimension(:,:), pointer :: mid_lat => null()
        real(dp), dimension(:,:), pointer :: corn_lon => null()
        real(dp), dimension(:,:), pointer :: corn_lat => null()
        real(dp), dimension(:,:), pointer :: tmask => null()
        real(dp), dimension(:,:), pointer :: core_area => null()
        integer :: core_grid_dims(2), onx, ony
        integer(ESMF_KIND_I4), dimension(:,:), allocatable :: tmaski
        logical luse_mask
        
        call read_CORE_mesh_mask_file(filename, grid_dims=core_grid_dims)
        onx=core_grid_dims(1)
        ony=core_grid_dims(2)
        allocate(tmaski(onx, ony))!mid_lon(onx, ony), mid_lat(onx, ony), corn_lon(onx+1, ony+1), corn_lat(onx+1, ony+1), tmask(onx, ony), core_area(onx, ony), &
                 

        call read_CORE_mesh_mask_file(filename,  mask=tmask, &
            mid_lon=mid_lon, mid_lat=mid_lat,corn_lon=corn_lon, corn_lat=corn_lat, &
            area=core_area)
        where (mid_lon<0) mid_lon=mid_lon+360
        where (corn_lon<0) corn_lon=corn_lon+360

        tmaski = int(tmask)
        if ( present(use_mask).and..not.use_mask ) then
          ! set all the mask valule to one if user decide to not use mask. 
          tmaski=1
        endif
        call define_esmf_grid('1peri', gridname, mid_lon, mid_lat, corn_lon, corn_lat, mask=tmaski, area=core_area, my_esmf_grid=esmf_coregrid)

        if ( present(grid_dims) ) then
        grid_dims = core_grid_dims
        endif

    end subroutine define_core_grid_from_file
end module esmf_utils
