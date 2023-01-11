!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains necessary routines for performing an 
!     nearest interpolation. For now, it will produce 4 weight points
!     like "bilinear" but only the first wil have a weight of 1 and 
!     the 3 others will have zero. This make it easier in NEMO...
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: remap_nearest.f,v 1.6 2001/08/22 18:20:40 pwjones Exp $
!
!     Copyright (c) 1997, 1998 the Regents of the University of 
!       California.
!
!     This software and ancillary information (herein called software) 
!     called SCRIP is made available under the terms described here.  
!     The software has been approved for release with associated 
!     LA-CC Number 98-45.
!
!     Unless otherwise indicated, this software has been authored
!     by an employee or employees of the University of California,
!     operator of the Los Alamos National Laboratory under Contract
!     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
!     Government has rights to use, reproduce, and distribute this
!     software.  The public may copy and use this software without
!     charge, provided that this Notice and any statement of authorship
!     are reproduced on all copies.  Neither the Government nor the
!     University makes any warranty, express or implied, or assumes
!     any liability or responsibility for the use of this software.
!
!     If software is modified to produce derivative works, such modified
!     software should be clearly marked, so as not to confuse it with 
!     the version available from Los Alamos National Laboratory.
!
!***********************************************************************

      module remap_nearest

!-----------------------------------------------------------------------

      use kinds_mod     ! defines common data types
      use constants     ! defines common constants
      use grids         ! module containing grid info
      use remap_vars    ! module containing remap info

      implicit none

!-----------------------------------------------------------------------

      integer (kind=int_kind), parameter :: &
          max_iter = 100   ! max iteration count for i,j iteration

      real (kind=dbl_kind), parameter :: &
           converge = 1.e-10_dbl_kind  ! convergence criterion

!***********************************************************************

      contains

!***********************************************************************

      subroutine remap_near

!-----------------------------------------------------------------------
!
!     this routine computes the weights for a nearest interpolation.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n,icount, &
           dst_add,         & ! destination address
           iter,            & ! iteration counter
           nmap            ! index of current map being computed

      integer (kind=int_kind), dimension(4) ::  &
           src_add         ! address for the four source points

      real (kind=dbl_kind), dimension(4)  :: &
           src_lats,        & ! latitudes  of four nearest corners
           src_lons,        & ! longitudes of four nearest corners
           wgts            ! nearest weights for four corners

      real (kind=dbl_kind) :: &
           plat, plon,        & ! lat/lon coords of destination point
           iguess, jguess,    & ! current guess for nearest coordinate
           thguess, phguess,  & ! current guess for lat/lon coordinate
           deli, delj,        & ! corrections to i,j
           dth1, dth2, dth3,  & ! some latitude  differences
           dph1, dph2, dph3,  & ! some longitude differences
           dthp, dphp,        & ! difference between point and sw corner
           mat1, mat2, mat3, mat4,  & ! matrix elements
           determinant,       & ! matrix determinant
           sum_wgts          ! sum of weights for normalization

!-----------------------------------------------------------------------
!
!     compute mappings from grid1 to grid2
!
!-----------------------------------------------------------------------

      nmap = 1
      if (grid1_rank /= 2) then
        stop 'Can not do nearest interpolation when grid_rank /= 2'
      endif

      !***
      !*** loop over destination grid 
      !***

      grid_loop1: do dst_add = 1, grid2_size

        src_add(:)=0
        wgts(:)=0

        if ( grid2_mask(dst_add)==1 ) then

          plat = grid2_center_lat(dst_add)
          plon = grid2_center_lon(dst_add)

          !***
          !*** find nearest square of grid points on source grid
          !***

          call grid_search_nearest(src_add, src_lats, src_lons,  &
                               plat, plon,  &
                               grid1_center_lat, grid1_center_lon, &
                               grid1_mask, bin_addr1)

          ! for the nearest, wgts(1) = 1.
          wgts(1) = 1.
          wgts(2) = 0.
          wgts(3) = 0.
          wgts(4) = 0.
         
        endif
        call store_link_nearest(dst_add, src_add, wgts, nmap)


      end do grid_loop1

!-----------------------------------------------------------------------
!
!     compute mappings from grid2 to grid1 if necessary
!
!-----------------------------------------------------------------------

      if (num_maps > 1) then

      nmap = 2
      if (grid2_rank /= 2) then
        stop 'Can not do nearest interpolation when grid_rank /= 2'
      endif

      !***
      !*** loop over destination grid 
      !***

      grid_loop2: do dst_add = 1, grid1_size

        src_add(:)=0
        wgts(:)=0

        if ( grid1_mask(dst_add)==1 ) then

          plat = grid1_center_lat(dst_add)
          plon = grid1_center_lon(dst_add)

          !***
          !*** find nearest square of grid points on source grid
          !***

          call grid_search_nearest(src_add, src_lats, src_lons,  &
                               plat, plon,  &
                               grid2_center_lat, grid2_center_lon, &
                               grid2_mask, bin_addr2)

          ! for the nearest, wgts(1) = 1.
          wgts(1) = 1.
          wgts(2) = 0.
          wgts(3) = 0.
          wgts(4) = 0.
         
        endif
        call store_link_nearest(dst_add, src_add, wgts, nmap)

      end do grid_loop2

      endif ! nmap=2

!-----------------------------------------------------------------------

      end subroutine remap_near

!***********************************************************************

      subroutine grid_search_nearest(src_add, src_lats, src_lons,  &
                                   plat, plon,  &
                                   src_center_lat, src_center_lon, &
                                   src_mask, &
                                   src_bin_add)

!-----------------------------------------------------------------------
!
!     this routine finds the location of the search point plat, plon
!     in the source grid and returns the corners needed for a bilinear
!     interpolation.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), dimension(4), intent(out) :: &
              src_add  ! address of each corner point enclosing P

      real (kind=dbl_kind), dimension(4), intent(out) :: &
              src_lats,  & ! latitudes  of the four corner points
              src_lons  ! longitudes of the four corner points

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), intent(in) :: &
              plat,    & ! latitude  of the search point
              plon    ! longitude of the search point

      real (kind=dbl_kind), dimension(:), intent(in) :: &
              src_center_lat,  & ! latitude  of each src grid center 
              src_center_lon  ! longitude of each src grid center

      integer (1), dimension(:),  intent(in) :: &
              src_mask ! mask of src 

      integer (kind=int_kind), dimension(:,:), intent(in) :: &
              src_bin_add ! latitude bins for restricting

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n, next_n, srch_add,    & ! dummy indices
          min_add, max_add ! addresses for restricting search

      real (kind=dbl_kind) ::   & 
            dist_min, distance ! for computing dist

!-----------------------------------------------------------------------
!
!     restrict search first using bins
!
!-----------------------------------------------------------------------

      src_add = 0

      min_add = size(src_center_lat)
      max_add = 1
      do n=1,num_srch_bins
        if (plat >= bin_lats(1,n) .and. plat <= bin_lats(2,n) .and. &
            plon >= bin_lons(1,n) .and. plon <= bin_lons(2,n)) then
          min_add = min(min_add, src_bin_add(1,n))
          max_add = max(max_add, src_bin_add(2,n))
        endif
      end do
      if (max_add.eq.1) then 
         min_add = 1
         max_add = size(src_center_lat)
      endif 
 
!-----------------------------------------------------------------------
!
!     now perform a more detailed search 
!
!-----------------------------------------------------------------------

      dist_min = 1e20 ! bug number to be replaced at first loop 

      srch_loop: do srch_add = min_add,max_add

        !*** first check bounding box
        if (src_mask(srch_add)==0) cycle srch_loop

        distance = acos(cos(plat)*cos(src_center_lat(srch_add))* &
                       (cos(plon)*cos(src_center_lon(srch_add)) + &
                        sin(plon)*sin(src_center_lon(srch_add)))+ &
                        sin(plat)*sin(src_center_lat(srch_add)))
        if (distance <= dist_min) then
           dist_min = distance 

          !***
          !*** we are at a closer point, replace previous values
          !***

          src_add(1) = srch_add
          src_add(2) = srch_add
          src_add(3) = srch_add
          src_add(4) = srch_add

          src_lats(1) = src_center_lat(srch_add)
          src_lats(2) = src_center_lat(srch_add)
          src_lats(3) = src_center_lat(srch_add)
          src_lats(4) = src_center_lat(srch_add)

          src_lons(1) = src_center_lon(srch_add)
          src_lons(2) = src_center_lon(srch_add)
          src_lons(3) = src_center_lon(srch_add)
          src_lons(4) = src_center_lon(srch_add)
!
          !***
          !*** otherwise move on to next cell
          !***

        endif !bounding box check
      end do srch_loop

!-----------------------------------------------------------------------

      end subroutine grid_search_nearest 

!***********************************************************************

      subroutine store_link_nearest(dst_add, src_add, weights, nmap)

!-----------------------------------------------------------------------
!
!     this routine stores the address and weight for four links 
!     associated with one destination point in the appropriate address 
!     and weight arrays and resizes those arrays if necessary.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: &
              dst_add,   & ! address on destination grid
              nmap      ! identifies which direction for mapping

      integer (kind=int_kind), dimension(4), intent(in) :: &
              src_add   ! addresses on source grid

      real (kind=dbl_kind), dimension(4), intent(in) :: &
              weights ! array of remapping weights for these links

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n,  & ! dummy index
             num_links_old          ! placeholder for old link number

!-----------------------------------------------------------------------
!
!     increment number of links and check to see if remap arrays need
!     to be increased to accomodate the new link.  then store the
!     link.
!
!-----------------------------------------------------------------------

      select case (nmap)
      case(1)

        num_links_old  = num_links_map1
        num_links_map1 = num_links_old + 4

        if (num_links_map1 > max_links_map1)  &
           call resize_remap_vars(1,resize_increment)

        do n=1,4
          grid1_add_map1(num_links_old+n) = src_add(n)
          grid2_add_map1(num_links_old+n) = dst_add
          wts_map1    (1,num_links_old+n) = weights(n)
        end do

      case(2)

        num_links_old  = num_links_map2
        num_links_map2 = num_links_old + 4

        if (num_links_map2 > max_links_map2)  &
           call resize_remap_vars(2,resize_increment)

        do n=1,4
          grid1_add_map2(num_links_old+n) = dst_add
          grid2_add_map2(num_links_old+n) = src_add(n)
          wts_map2    (1,num_links_old+n) = weights(n)
        end do

      end select

!-----------------------------------------------------------------------

      end subroutine store_link_nearest

!***********************************************************************

      end module remap_nearest

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
