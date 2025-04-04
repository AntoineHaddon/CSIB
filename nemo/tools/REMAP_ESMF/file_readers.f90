module file_readers
    use ESMF
    use ncdf
    use iso_fortran_env, only:  int32, int64, real32, real64
    !use esmf_utils, only : define_esmf_grid, define_esmf_field 
    private
    public :: read_yg, read_NEMO_mesh_mask_file, read_CORE_mesh_mask_file

    contains

    subroutine read_yg(file_name, lon, lat, lon_bnds, lat_bnds, flnd, bego)
        implicit none

        !--- The name of the NEMO mesh_mask file
        logical :: exists
        character(*), intent(in)  :: file_name
        integer(kind=4) :: fid
        real(kind=8), pointer, dimension(:, :), intent(out) :: lon
        real(kind=8), pointer, dimension(:, :), intent(out) :: lat
        real(kind=8), pointer, dimension(:, :)  :: rlon=>null() , rlat=>null()
        real(kind=8), pointer, dimension(:, :), intent(out), optional :: bego, flnd
        real(kind=8), pointer, dimension(:, :, :), intent(out) :: lon_bnds, lat_bnds
        
        inquire(file=trim(file_name), exist=exists)
        if ( .not. exists ) then
            write(6,*)"read_yin: File does not exist ",trim(file_name)
            call flush(6)
            call xit("read yin",-2)
        endif

        !--- Open the yin file for reading
        call ncdf_open(fid, trim(file_name), mode="read")
        call ncdf_read_var(fid, "lon", lon, nt=1)
        call ncdf_read_var(fid, "lat", lat, nt=1)

        if (present(flnd) ) then
            call ncdf_read_var(fid, "FLND", flnd, nt=1)
        endif
        if( present(bego)) then
            call ncdf_read_var(fid, "BEGO", bego, nt=1)
        endif
        call ncdf_read_var(fid, "lon_bnds", lon_bnds, nt=1)
        call ncdf_read_var(fid, "lat_bnds", lat_bnds, nt=1)
        call ncdf_close(fid)
    end subroutine read_yg

    ! Several routines in this module adapted are from CanCPL existing routines
    subroutine mm_corner_lon(clon, offset, corn_lon)
        !--- clon contains corner longitudes as read from the mesh_mask file
        real(kind=8), intent(in),  pointer :: clon(:,:)
  
        !--- offset is what needs to be added to (i,j) to get to the correct elements
        !--- in the clon arrays (read from the mesh_mask file) that correspond to the
        !--- northeast corner of each cell
        integer, intent(in) :: offset(2)
  
        !--- corn_lon is output and will contain the corner longitudes
        !--- suitable for use in ESMF grid definitions
        real(kind=8), intent(out), pointer :: corn_lon(:,:)
  
        !--- Local
        real(kind=8), parameter :: circle = 360.0_8 ! assumes corn_lon is in degrees
        integer :: nx, ny, i, j, ic, jc, imid
  
        if ( .not. associated(clon) ) then
          write(6,*)"mm_corner_lon: Input clon array is not allocated"
          call flush(6)
          call xit("mm_corner_lon",-1)
        endif
  
        !--- Get the dimension size from the input clon array
        nx = size(clon,dim=1)
        ny = size(clon,dim=2)
  
        !--- Always destroy then recreate the corn_lon output array
        if ( associated(corn_lon) ) deallocate(corn_lon)
        allocate ( corn_lon(nx+1,ny+1) )
        corn_lon = 0.0_8
  
        !-----------------------------------------------------------------------
        ! Corners must be ordered as follows:
        !     4 3
        !     1 2
        !
        ! Create a pair of arrays for the corners where corn_lon(1,1) corresponds
        ! to the south west corner of a box centered at (glam(1,1),gphi(1,1))
  
        !--- Assign the internal points
        do j = 2,ny
          do i = 2,nx
            ic = i + offset(1) - 1
            jc = j + offset(2) - 1
            corn_lon(i,j) = clon(ic,jc)
          enddo
        enddo
  
        !--- Then the boundary points
        !--- Assume that in the data read from the coordinates file there are 2 overlap
        !--- longitudes so that the following conditions hold (glam == lon, gphi == lat)
        !---     glam(1,:) == glam(nx-1,:)     gphi(1,:) == gphi(nx-1,:)
        !---     glam(2,:) == glam(nx,:)       gphi(2,:) == gphi(nx,:)
        !--- and
        !---     clon(1,:) == clon(nx-1,:)   clat(1,:) == clat(nx-1,:)
        !---     clon(2,:) == clon(nx,:)     clat(2,:) == clat(nx,:)
        do i=2,nx
          ic = i + offset(1) - 1
          !--- First (south/bottom) row
          corn_lon(i,1) = clon(ic,1)
  
          !--- Last (north/top) row
          corn_lon(i,ny+1) = clon(ic,ny)
        enddo
  
        do j=2,ny
          jc = j + offset(2) - 1
          !--- First (west/left) column
          corn_lon(1,j) = clon(nx-2,jc)
  
          !--- Last (east/right) column
          corn_lon(nx+1,j) = clon(2,jc)
        enddo
  
        !--- Then the corner points
        ic = nx-2
        corn_lon(1,1) = clon(ic,1)
  
        ic = 2
        corn_lon(nx+1,1) = clon(ic,1)
  
        imid = (nx-1)/2 + 1
        ic = 2*imid - nx + 1
        corn_lon(1,ny+1) = clon(ic,ny)
  
        ic = 2*imid
        corn_lon(nx+1,ny+1) = clon(ic,ny)
  
        !--- Zero out small values
        where (abs(corn_lon)<1e-20) corn_lon=0.0_8
  
      end subroutine mm_corner_lon
  
      subroutine mm_corner_lat(clat, offset, corn_lat)
        !--- clat contains corner latitudes as read from the mesh_mask file
        real(kind=8), intent(in),  pointer :: clat(:,:)
  
        !--- offset is what needs to be added to (i,j) to get to the correct elements
        !--- in the clat arrays (read from the mesh_mask file) that correspond to the
        !--- northeast corner of each cell
        integer, intent(in) :: offset(2)
  
        !--- corn_lat is output and will contain the corner latitudes
        !--- suitable for use in ESMF grid definitions
        real(kind=8), intent(out), pointer :: corn_lat(:,:)
  
        !--- Local
        integer :: nx, ny, i, j, ic, jc, imid
  
        if ( .not. associated(clat) ) then
          write(6,*)"mm_corner_lat: Input clat array is not allocated"
          call flush(6)
          call xit("mm_corner_lat",-1)
        endif
  
        !--- Get the dimension size from the input clat array
        nx = size(clat,dim=1)
        ny = size(clat,dim=2)
  
        !--- Always destroy then recreate the corn_lat output array
        if ( associated(corn_lat) ) deallocate(corn_lat)
        allocate ( corn_lat(nx+1,ny+1) )
        corn_lat = 0.0_8
  
        !-----------------------------------------------------------------------
        ! Corners must be ordered as follows:
        !     4 3
        !     1 2
        !
        ! Create a pair of arrays for the corners where corn_lat(1,1) corresponds
        ! to the south west corner of a box centered at (glam(1,1),gphi(1,1))
  
        !--- Assign the internal points
        do j = 2,ny
          do i = 2,nx
            ic = i + offset(1) - 1
            jc = j + offset(2) - 1
            corn_lat(i,j) = clat(ic,jc)
          enddo
        enddo
  
        !--- Then the boundary points
        !--- Assume that in the data read from the coordinates file there are 2 overlap
        !--- longitudes so that the following conditions hold (glam == lon, gphi == lat)
        !---     glam(1,:) == glam(nx-1,:)     gphi(1,:) == gphi(nx-1,:)
        !---     glam(2,:) == glam(nx,:)       gphi(2,:) == gphi(nx,:)
        !--- and
        !---     clon(1,:) == clon(nx-1,:)   clat(1,:) == clat(nx-1,:)
        !---     clon(2,:) == clon(nx,:)     clat(2,:) == clat(nx,:)
        do i=2,nx
          ic = i + offset(1) - 1
          !--- First (south/bottom) row
          corn_lat(i,1) = clat(ic,1) - (clat(ic,2)-clat(ic,1))
  
          !--- Last (north/top) row
          corn_lat(i,ny+1) = clat(ic,ny)
        enddo
  
        do j=2,ny
          jc = j + offset(2) - 1
          !--- First (west/left) column
          corn_lat(1,j) = clat(nx-2,jc)
  
          !--- Last (east/right) column
          corn_lat(nx+1,j) = clat(2,jc)
        enddo
  
        !--- Then the corner points
        ic = nx-2
        corn_lat(1,1) = clat(ic,1) - (clat(ic,2)-clat(ic,1))
  
        ic = 2
        corn_lat(nx+1,1) = clat(ic,1) - (clat(ic,2)-clat(ic,1))
  
        imid = (nx-1)/2 + 1
        ic = 2*imid - nx + 1
        corn_lat(1,ny+1) = clat(ic,ny)
  
        ic = 2*imid
        corn_lat(nx+1,ny+1) = clat(ic,ny)
  
        !--- Zero out small values
        where (abs(corn_lat)<1e-20) corn_lat=0.0_8
  
      end subroutine mm_corner_lat

    subroutine read_CORE_mesh_mask_file(file_name, grid, overlap, &
        mid_lon, mid_lat, corn_lon, corn_lat, &
        mask, imask, area, &
        title, grid_size, grid_dims, grid_edges)

        implicit none

        !--- The name of the NEMO mesh_mask file
        character(*), intent(in)  :: file_name

        !--- The NEMO grid type to extract (T,U,V,F) default is T grid
        character(*), intent(in), optional :: grid

        !--- If overlap = F then do not include overlaping grid cells in any
        !--- dimension in the output grid description
        logical, intent(in), optional :: overlap

        !--- lon/lat of grid cell centers
        real(kind=8), intent(out), pointer, optional ::  mid_lon(:,:)
        real(kind=8), intent(out), pointer, optional ::  mid_lat(:,:)

        !--- lon/lat of grid cell corners
        !--- It is assumed that there are 4 corners per cell
        real(kind=8), intent(out), pointer, optional :: corn_lon(:,:)
        real(kind=8), intent(out), pointer, optional :: corn_lat(:,:)

        !--- Land mask at level 1 of the T grid (real or integer)
        real(kind=8), intent(out), pointer, optional ::  mask(:,:)
        integer,      intent(out), pointer, optional :: imask(:,:)

        !--- The grid cell area (m2)
        real(kind=8), intent(out), pointer, optional :: area(:,:)

        !--- The grid title
        character(*), intent(out), optional :: title

        !--- The total number of points in the grid
        integer, intent(out), optional :: grid_size

        !--- The shape of the grid ie (nx,ny)
        integer, intent(out), optional :: grid_dims(2)

        !--- The number of overlap grids (lower and upper) in each dimension
        integer, intent(out), optional :: grid_edges(2,2)

        !--- Local
        logical :: exists, resize
        integer(kind=4) :: fid
        character(32)  :: lon_name="xc", lat_name="yc"
        character(32)  :: lon_corn_name="xv", lat_corn_name="yv"
        character(32)  :: area_name="area", mask_name="mask"
        real(kind=8), pointer, dimension(:,:,:)   :: clon=>null(), clat=>null()
        real(kind=8), pointer, dimension(:,:) :: xmask=>null()
        integer     , pointer, dimension(:,:) :: ximask=>null()
        integer :: nx, ny, idx
        integer :: offset(2)
        character(64),   pointer :: dim_names(:)
        integer,         pointer :: dim_sizes(:)

        if ( len_trim(file_name) < 1 ) then
        write(6,*)"read_CORE_mesh_mask_file: file_name is blank"
        call flush(6)
        call xit("read_CORE_mesh_mask_file",-1)
        endif

        inquire(file=trim(file_name), exist=exists)
        if ( .not. exists ) then
        write(6,*)"read_CORE_mesh_mask_file: File does not exist ",trim(file_name)
        call flush(6)
        call xit("read_CORE_mesh_mask_file",-2)
        endif

        if ( present(overlap) ) then
        !--- If overlap = T then keep the overlap longitude in output arrays
        !--- If overlap = F then remove the overlap longitude before return
        resize = .not. overlap
        else
        !--- Remove the overlap longitude by default
        resize = .true.
        endif

        !--- Open the mesh_mask file for reading
        call ncdf_open(fid, trim(file_name), mode="read")

        !--- Read dimension sizes from the file
        call ncdf_inquire_dims(fid, dim_names=dim_names, dim_sizes=dim_sizes)
        nx = -1
        ny = -1
        do idx=1,size(dim_names)
        if ( trim(adjustl(dim_names(idx))) .eq. "ni") nx = dim_sizes(idx)
        if ( trim(adjustl(dim_names(idx))) .eq. "nj") ny = dim_sizes(idx)
        enddo
        if ( nx < 0 .or. ny < 0 ) then
        write(6,*)"read_CORE_mesh_mask_file: Unable to determine nx/ny from ",trim(file_name)
        call flush(6)
        call xit("read_CORE_mesh_mask_file",-4)
        endif

        write(6,*)"  read_CORE_mesh_mask_file: ",trim(file_name),"   nx=",nx,"   ny=",ny
        call flush(6)

        if ( present(title) ) then
        title = " "
        if ( nx > 999 .or. ny > 999 ) then
        write(title,'(3a,i4,a,i4)')"CORE "," -- ",nx,"x",ny
        else
        write(title,'(3a,i3.3,a,i3.3)')"CORE "," -- ",nx,"x",ny
        endif
        write(6,*)"read_CORE_mesh_mask_file: title = ",trim(title)
        call flush(6)
        endif

        if ( present(grid_edges) ) then
        grid_edges = 0
        !--- Assume 2 overlap logitudes
        grid_edges(2,1) = 2
        !write(6,*)"read_CORE_mesh_mask_file: grid_edges = ",grid_edges
        call flush(6)
        endif

        if ( present(mid_lon) ) then
        !--- Read grid cell center longitudes from the file
        !--- This call will destroy then reallocate mid_lon
        call ncdf_read_var(fid, trim(lon_name), mid_lon)
        
        endif

        if ( present(mid_lat) ) then
        !--- Read grid cell center latitudes from the file
        !--- This call will destroy then reallocate mid_lat
        call ncdf_read_var(fid, trim(lat_name), mid_lat)
        
        
        endif

        if ( present(corn_lon) ) then
        !--- Read grid cell corner longitudes from the file
        !--- This call will destroy then reallocate clon
        call ncdf_read_var(fid, trim(lon_corn_name), clon )
        !--- Reformat these corner lon values to add missing boundary points
        !--- corn_lon is returned with 1 more row and 1 more column than clon
        call mm_corner_lon(clon(:,:,3), offset, corn_lon)
        
        endif

        if ( present(corn_lat) ) then
        !--- Read grid cell corner latitudes from the file
        !--- This call will destroy then reallocate clat
        call ncdf_read_var(fid, trim(lat_corn_name), clat)
        
        !--- Reformat these corner lat values to add missing boundary points
        !--- corn_lat is returned with 1 more row and 1 more column than clat
        call mm_corner_lat(clat(:,:,3), offset, corn_lat)
        
        endif

        if ( present(mask) ) then
        !--- Read T grid land mask from the file
        call ncdf_read_var(fid, trim(mask_name), xmask)
        
        !--- Assign mask with surface (level 1) mask values from
        !--- the 3D mask that was read from the file
        if ( associated(mask) ) deallocate(mask)
        allocate( mask(size(xmask,dim=1),size(xmask,dim=2)) )
        mask = xmask
        deallocate(xmask)
        endif

        if ( present(imask) ) then
        !--- Read T grid land mask from the file
        call ncdf_read_var(fid, trim(mask_name), ximask)
        
        !--- Assign mask with surface (level 1) mask values from
        !--- the 3D mask that was read from the file
        if ( associated(imask) ) deallocate(imask)
        allocate( imask(size(ximask,dim=1),size(ximask,dim=2)) )
        imask = ximask
        deallocate(ximask)
        
        endif

        if ( present(area) ) then
        !--- Read scale lengths from the file
        
        call ncdf_read_var(fid, trim(area_name), area)
        
        endif

        !--- Define grid_dims and grid_size here so that any changes to nx,ny due to
        !--- removing overlap longitudes will be reflected in these output values
        if ( present(grid_dims) ) then
        grid_dims(1) = nx
        grid_dims(2) = ny
        endif

        if ( present(grid_size) ) then
        grid_size = nx*ny
        endif

        !--- Close the file
        call ncdf_close(fid)

    end subroutine read_CORE_mesh_mask_file

    subroutine read_NEMO_mesh_mask_file(file_name, grid, overlap, &
        mid_lon, mid_lat, corn_lon, corn_lat, &
        mask, imask, area, &
        title, grid_size, grid_dims, grid_edges)

        implicit none

        !--- The name of the NEMO mesh_mask file
        character(*), intent(in)  :: file_name

        !--- The NEMO grid type to extract (T,U,V,F) default is T grid
        character(*), intent(in), optional :: grid

        !--- If overlap = F then do not include overlaping grid cells in any
        !--- dimension in the output grid description
        logical, intent(in), optional :: overlap

        !--- lon/lat of grid cell centers
        real(kind=8), intent(out), pointer, optional ::  mid_lon(:,:)
        real(kind=8), intent(out), pointer, optional ::  mid_lat(:,:)

        !--- lon/lat of grid cell corners
        !--- It is assumed that there are 4 corners per cell
        real(kind=8), intent(out), pointer, optional :: corn_lon(:,:)
        real(kind=8), intent(out), pointer, optional :: corn_lat(:,:)

        !--- Land mask at level 1 of the T grid (real or integer)
        real(kind=8), intent(out), pointer, optional ::  mask(:,:)
        integer,      intent(out), pointer, optional :: imask(:,:)

        !--- The grid cell area (m2)
        real(kind=8), intent(out), pointer, optional :: area(:,:)

        !--- The grid title
        character(*), intent(out), optional :: title

        !--- The total number of points in the grid
        integer, intent(out), optional :: grid_size

        !--- The shape of the grid ie (nx,ny)
        integer, intent(out), optional :: grid_dims(2)

        !--- The number of overlap grids (lower and upper) in each dimension
        integer, intent(out), optional :: grid_edges(2,2)

        !--- Local
        logical :: exists, resize
        integer(kind=4) :: fid
        character(32)  :: lon_name="glamt", lat_name="gphit"
        character(32)  :: lon_corn_name="glamf", lat_corn_name="gphif"
        character(32)  :: e1_name="e1t", e2_name="e2t", mask_name="tmask"
        character(1)   :: gtype="T"
        real(kind=8), pointer, dimension(:,:)   :: clon=>null(), clat=>null()
        real(kind=8), pointer, dimension(:,:)   :: e1=>null(), e2=>null()
        real(kind=8), pointer, dimension(:,:,:) :: xmask=>null()
        integer     , pointer, dimension(:,:,:) :: ximask=>null()
        integer :: nx, ny, idx
        integer :: offset(2)
        character(64),   pointer :: dim_names(:)
        integer,         pointer :: dim_sizes(:)

        if ( len_trim(file_name) < 1 ) then
        write(6,*)"read_NEMO_mesh_mask_file: file_name is blank"
        call flush(6)
        call xit("read_NEMO_mesh_mask_file",-1)
        endif

        inquire(file=trim(file_name), exist=exists)
        if ( .not. exists ) then
        write(6,*)"read_NEMO_mesh_mask_file: File does not exist ",trim(file_name)
        call flush(6)
        call xit("read_NEMO_mesh_mask_file",-2)
        endif

        if ( present(overlap) ) then
        !--- If overlap = T then keep the overlap longitude in output arrays
        !--- If overlap = F then remove the overlap longitude before return
        resize = .not. overlap
        else
        !--- Remove the overlap longitude by default
        resize = .true.
        endif

        if ( present(grid) ) then
        gtype = grid(1:1)
        !--- Define names of variables found in the mesh mask file that correspond to
        !--- grid cell centers, corners and scale factors
        !--- offset is what needs to be added to (i,j) to get to the correct elements
        !--- in the corner arrays (as read from the file) that correspond to the
        !--- northeast corner of each cell
        select case(gtype)
        case ("t", "T")
        gtype = 'T'
        lon_name = "glamt"
        lat_name = "gphit"
        lon_corn_name = "glamf"
        lat_corn_name = "gphif"
        mask_name = "tmask"
        e1_name = "e1t"
        e2_name = "e2t"
        offset = (/0,0/)
        case ("u", "U")
        gtype = 'U'
        lon_name = "glamu"
        lat_name = "gphiu"
        lon_corn_name = "glamv"
        lat_corn_name = "gphiv"
        mask_name = "umask"
        e1_name = "e1u"
        e2_name = "e2u"
        offset = (/1,0/)
        case ("v", "V")
        gtype = 'V'
        lon_name = "glamv"
        lat_name = "gphiv"
        lon_corn_name = "glamu"
        lat_corn_name = "gphiu"
        mask_name = "vmask"
        e1_name = "e1v"
        e2_name = "e2v"
        offset = (/0,1/)
        case ("f", "F")
        !--- The values for lon_corn_name, lat_corn_name and the values for the offset
        !--- need to be verified for the F grid
        write(6,*)"read_NEMO_mesh_mask_file: Grid type F not currently supported"
        call flush(6)
        call xit("read_NEMO_mesh_mask_file",-9)
        gtype = 'F'
        lon_name = "glamf"
        lat_name = "gphif"
        lon_corn_name = "glamt"
        lat_corn_name = "gphit"
        mask_name = "fmask"
        e1_name = "e1f"
        e2_name = "e2f"
        offset = (/-1,0/)
        case default
        write(6,*)"read_NEMO_mesh_mask_file: Invalid grid type ",gtype
        call flush(6)
        call xit("read_NEMO_mesh_mask_file",-3)
        end select
        else
        !--- Default is the T grid
        gtype = 'T'
        lon_name = "glamt"
        lat_name = "gphit"
        lon_corn_name = "glamf"
        lat_corn_name = "gphif"
        mask_name = "tmask"
        e1_name = "e1t"
        e2_name = "e2t"
        offset = (/0,0/)
        endif

        !--- Open the mesh_mask file for reading
        call ncdf_open(fid, trim(file_name), mode="read")

        !--- Read dimension sizes from the file
        call ncdf_inquire_dims(fid, dim_names=dim_names, dim_sizes=dim_sizes)
        nx = -1
        ny = -1
        do idx=1,size(dim_names)
        if ( trim(adjustl(dim_names(idx))) .eq. "x") nx = dim_sizes(idx)
        if ( trim(adjustl(dim_names(idx))) .eq. "y") ny = dim_sizes(idx)
        enddo
        if ( nx < 0 .or. ny < 0 ) then
        write(6,*)"read_NEMO_mesh_mask_file: Unable to determine nx/ny from ",trim(file_name)
        call flush(6)
        call xit("read_NEMO_mesh_mask_file",-4)
        endif

        write(6,*)"  read_NEMO_mesh_mask_file: ",trim(file_name),"   nx=",nx,"   ny=",ny
        call flush(6)

        if ( present(title) ) then
        title = " "
        if ( nx > 999 .or. ny > 999 ) then
        write(title,'(3a,i4,a,i4)')"ORCA ",gtype," grid -- ",nx,"x",ny
        else
        write(title,'(3a,i3.3,a,i3.3)')"ORCA ",gtype," grid -- ",nx,"x",ny
        endif
        write(6,*)"read_NEMO_mesh_mask_file: title = ",trim(title)
        call flush(6)
        endif

        if ( present(grid_edges) ) then
        grid_edges = 0
        !--- Assume 2 overlap logitudes
        grid_edges(2,1) = 2
        !write(6,*)"read_NEMO_mesh_mask_file: grid_edges = ",grid_edges
        call flush(6)
        endif

        if ( present(mid_lon) ) then
        !--- The file should contain
        !---     float glamt(t, y, x)
        !---     float glamu(t, y, x)
        !---     float glamv(t, y, x)
        !---     float glamf(t, y, x)
        !--- Read grid cell center longitudes from the file
        !--- This call will destroy then reallocate mid_lon
        call ncdf_read_var(fid, trim(lon_name), mid_lon, nt=1)
        
        endif

        if ( present(mid_lat) ) then
        !--- The file should contain
        !---     float gphit(t, y, x)
        !---     float gphiu(t, y, x)
        !---     float gphiv(t, y, x)
        !---     float gphif(t, y, x)
        !--- Read grid cell center latitudes from the file
        !--- This call will destroy then reallocate mid_lat
        call ncdf_read_var(fid, trim(lat_name), mid_lat, nt=1)
        
        
        endif

        if ( present(corn_lon) ) then
        !--- The file should contain
        !---     float glamt(t, y, x)
        !---     float glamu(t, y, x)
        !---     float glamv(t, y, x)
        !---     float glamf(t, y, x)
        !--- Read grid cell corner longitudes from the file
        !--- This call will destroy then reallocate clon
        call ncdf_read_var(fid, trim(lon_corn_name), clon, nt=1)
        
        !--- Reformat these corner lon values to add missing boundary points
        !--- corn_lon is returned with 1 more row and 1 more column than clon
        call mm_corner_lon(clon, offset, corn_lon)
        
        endif

        if ( present(corn_lat) ) then
        !--- The file should contain
        !---     float gphit(t, y, x)
        !---     float gphiu(t, y, x)
        !---     float gphiv(t, y, x)
        !---     float gphif(t, y, x)
        !--- Read grid cell corner latitudes from the file
        !--- This call will destroy then reallocate clat
        call ncdf_read_var(fid, trim(lat_corn_name), clat, nt=1)
        
        !--- Reformat these corner lat values to add missing boundary points
        !--- corn_lat is returned with 1 more row and 1 more column than clat
        call mm_corner_lat(clat, offset, corn_lat)
        
        endif

        if ( present(mask) ) then
        !--- TODO --- This should return tmask_i
        !write(6,*)"read_NEMO_mesh_mask_file: Optional mask argument not yet supported"
        !call xit("read_NEMO_mesh_mask_file",-98)
        !--- The file should contain
        !---     byte tmask(t, z, y, x)
        !---     byte umask(t, z, y, x)
        !---     byte vmask(t, z, y, x)
        !---     byte fmask(t, z, y, x)
        !--- Read T grid land mask from the file
        call ncdf_read_var(fid, trim(mask_name), xmask, nt=1)
        
        !--- Assign mask with surface (level 1) mask values from
        !--- the 3D mask that was read from the file
        if ( associated(mask) ) deallocate(mask)
        allocate( mask(size(xmask,dim=1),size(xmask,dim=2)) )
        mask = xmask(:,:,1)
        deallocate(xmask)
        endif

        if ( present(imask) ) then
        !--- TODO --- This should return tmask_i
        !write(6,*)"read_NEMO_mesh_mask_file: Optional imask argument not yet supported"
        !call xit("read_NEMO_mesh_mask_file",-99)
        !--- The file should contain
        !---     byte tmask(t, z, y, x)
        !---     byte umask(t, z, y, x)
        !---     byte vmask(t, z, y, x)
        !---     byte fmask(t, z, y, x)
        !--- Read T grid land mask from the file
        call ncdf_read_var(fid, trim(mask_name), ximask, nt=1)
        
        !--- Assign mask with surface (level 1) mask values from
        !--- the 3D mask that was read from the file
        if ( associated(imask) ) deallocate(imask)
        allocate( imask(size(ximask,dim=1),size(ximask,dim=2)) )
        imask = ximask(:,:,1)
        deallocate(ximask)
        
        endif

        if ( present(area) ) then
        !--- The file should contain
        !---     double e1t(t, y, x)
        !---     double e1u(t, y, x)
        !---     double e1v(t, y, x)
        !---     double e1f(t, y, x)
        !---
        !---     double e2t(t, y, x)
        !---     double e2u(t, y, x)
        !---     double e2v(t, y, x)
        !---     double e2f(t, y, x)
        !--- Read scale lengths from the file
        
        call ncdf_read_var(fid, trim(e1_name), e1, nt=1)
        call ncdf_read_var(fid, trim(e2_name), e2, nt=1)
        
        !--- Grid cell area is the product of these two scale lengths
        if ( associated(area) ) deallocate(area)
        allocate( area(size(e1,dim=1),size(e1,dim=2)) )
        area = e1 * e2
        deallocate( e1, e2 )
        endif

        !--- Define grid_dims and grid_size here so that any changes to nx,ny due to
        !--- removing overlap longitudes will be reflected in these output values
        if ( present(grid_dims) ) then
        grid_dims(1) = nx
        grid_dims(2) = ny
        endif

        if ( present(grid_size) ) then
        grid_size = nx*ny
        endif

        !--- Close the file
        call ncdf_close(fid)

    end subroutine read_NEMO_mesh_mask_file
end module file_readers
