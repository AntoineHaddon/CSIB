MODULE spreading

    private
    public :: runoff_push

    integer, parameter :: dp=8
    Integer, allocatable :: ii_cloast_routing(:,:),jj_cloast_routing(:,:)
    LOGICAL, save :: first = .true.

CONTAINS

  subroutine runoff_push(runoff_in, runoff_lon, runoff_lat, riv_mask, ocean_areas, runoff_modified)
      !======================================================================
      ! Redirect river runoff to riv_mask (closest point)
      !
      ! Nicolas Lambert, 11/2022
      !
      !=======================================================================
      IMPLICIT NONE
      REAL(kind=8), intent (in), DIMENSION(:,:)  :: runoff_in, ocean_areas
      INTEGER(kind=4), intent (in), DIMENSION(:,:)  ::  riv_mask
      REAL(kind=8), intent (in), DIMENSION(:,:)  :: runoff_lon,runoff_lat
      REAL(kind=8), intent (out), DIMENSION(:,:) :: runoff_modified
      REAL dist, dist_min
      INTEGER idx1, idy1, idx2, idy2, nx, ny
      REAL(dp) glob_runoff_in, glob_runoff_out                                ! Globally integrated input and output runoff

      nx = size(runoff_in,1)
      ny = size(runoff_in,2)
      if (first) then
        ! Do the initalization of the routing arrays at the first iteration
        allocate(ii_cloast_routing(nx,ny),jj_cloast_routing(nx,ny))
        ii_cloast_routing(:,:)=0
        jj_cloast_routing(:,:)=0
        first =.false.

      endif

      ! runoff_in in the valid grid cell of runoff_modified
      runoff_modified(:,:)=0.
      where (riv_mask.eq.1) runoff_modified = runoff_in
      ! re-route all the runoff values to the closest cell on riv_mask
        do idx1=1,nx
        do idy1=1,ny
            if (runoff_in(idx1,idy1).eq.0.) cycle ! skip re-routing if no runoff
            if (riv_mask(idx1,idy1).eq.1)   cycle  ! skip re-routing if already a oean point

            if (ii_cloast_routing(idx1,idy1)==0.or.jj_cloast_routing(idx1,idy1)==0) then
                ! only look for re-routing indices if not done yet (should all be done at he first iteration)
                dist_min = 1e20
                ii_cloast_routing(idx1,idy1)=idx1 ! default, will be replaced
                jj_cloast_routing(idx1,idy1)=idy1
                do idx2=1,nx
                do idy2=1,ny
                    if (riv_mask(idx2,idy2).eq.0) cycle
                    dist = distance(runoff_lon(idx1,idy1),runoff_lat(idx1,idy1),runoff_lon(idx2,idy2),runoff_lat(idx2,idy2))
                    if (dist.lt.dist_min) then
                        dist_min = dist
                        ii_cloast_routing(idx1,idy1)=idx2
                        jj_cloast_routing(idx1,idy1)=idy2
                    endif
                enddo
                enddo
            endif

            ! re-use the indices variable now that the loop is over
            idx2=ii_cloast_routing(idx1,idy1)
            idy2=jj_cloast_routing(idx1,idy1)

            runoff_modified(idx2,idy2)=runoff_modified(idx2,idy2)+( runoff_in(idx1,idy1)*ocean_areas(idx1,idy1) ) / ocean_areas(idx2,idy2)

        enddo
        enddo
        ! Compute and print the global runoff integrals before and after spreading.
        glob_runoff_in = SUM(runoff_in * ocean_areas )
        glob_runoff_out = SUM(runoff_modified * ocean_areas )
        !glob_runoff_out = SUM(runoff_modified * ocean_areas * riv_mask)
        !write(6,*)'runoff_cloast: Input runoff (1e9 kg/s):', glob_runoff_in/1e9_dp
        !write(6,*)'runoff_cloast: Out runoff (1e9 kg/s)  :', glob_runoff_out/1e9_dp
        !deallocate(ii_cloast_routing,jj_cloast_routing)

  end subroutine


      real function distance(lon1,lat1,lon2,lat2)
!     giv the distance in kilometers between two coordinates
      real(kind=8) lon1,lat1,lon2,lat2
      real(kind=8) radius, pi, rad
      parameter (pi = 3.14159)
      parameter (radius = 6371.009)

      !haversine formula formula for
      rad=2.*asind(sqrt(sind(abs(lat2-lat1)/2.)**2.+cosd(lat1)*cosd(lat2)*sind(abs(lon2-lon1)/2.)**2.))

      distance=rad*pi/180.*radius


      end


     real(kind=8) function asind(x)
     real(kind=8) x
     asind=asin(x/180.0*3.141592)
     return
     end

     real(kind=8) function cosd(x)
     real(kind=8) x
     cosd=cos(x/180.0*3.141592)
     return
     end

     real(kind=8) function sind(x)
     real(kind=8) x
     sind=sin(x/180.0*3.141592)
     return
     end



END MODULE spreading
