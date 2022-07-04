!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     this module contains necessary routines for performing an 
!     fountain interpolation. This look for non-zero value in the
!     data to re-direct them at the closest point in NEMO. The
!     weight can be one-to-one or according the the erea (conservative)
!
!-----------------------------------------------------------------------
!
!     CVS:$Id: remap_fountain.f,v 1.6 2001/08/22 18:20:40 pwjones Exp $
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

      module remap_fountainhead

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

      subroutine remap_fountain

!-----------------------------------------------------------------------
!
!     this routine computes the weights for a fountain interpolation.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n,icount, &
           src_add,         & ! destination address
           iter,            & ! iteration counter
           nmap            ! index of current map being computed

      integer (kind=int_kind)  ::  &
           dst_add         ! address for the four source points

      real (kind=dbl_kind)   :: &
           wgts            ! closest weights for four corners

      real (kind=dbl_kind) :: &
           plat, plon,        & ! lat/lon coords of destination point
           iguess, jguess,    & ! current guess for closest coordinate
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
        stop 'Can not do fountain interpolation when grid_rank /= 2'
      endif
      print*,'Fountain remap'
      print*,'============================'

      !***
      !*** loop over destination grid 
      !***

      grid_loop1: do src_add = 1, grid1_size

        dst_add=0
        wgts=0

        if ( grid1_mask(src_add)==1 ) then

          plat = grid1_center_lat(src_add)
          plon = grid1_center_lon(src_add)

          !***
          !*** find closest square of grid points on source grid
          !***

          call grid_search_closest(dst_add, &
                               plat, plon,  &
                               grid2_center_lat, grid2_center_lon, &
                               grid2_mask, bin_addr2)
          if (dst_add.eq.0) cycle

          if (norm_opt .eq. norm_opt_frcarea) then
            if (.not.luse_grid1_area.or..not.luse_grid2_area) & 
               stop 'luse_grid1_area and luse_grid2_area need to be true for  norm_opt = frcarea '
                    !'
            wgts = grid1_area_in(src_add) / grid2_area_in(dst_add)
            ! to conserve if runoff is in kg/m2/s  (area)

          elseif (norm_opt .eq. norm_opt_none) then
            ! simply put one in tha case
            wgts = 1.
          else ; stop 'unknown norm_opt for fouintain interpolation'
          endif
         
          call store_link_fountain(dst_add, src_add, wgts, nmap)
        endif

      end do grid_loop1
      where (grid1_add_map1.eq.0) grid1_add_map1 = grid1_size ! NEMO need src .ne. to zeros. Put it at grid1_size put no worries, wgt()=0.

!-----------------------------------------------------------------------
!
!     compute mappings from grid2 to grid1 if necessary
!
!-----------------------------------------------------------------------

      if (num_maps > 1) then

      nmap = 2
      if (grid2_rank /= 2) then
        stop 'Can not do fountain interpolation when grid_rank /= 2'
      endif

      !***
      !*** loop over destination grid 
      !***

      grid_loop2: do src_add = 1, grid2_size

        dst_add=0
        wgts=0

        if ( grid2_mask(src_add) ==1) then

          plat = grid2_center_lat(src_add)
          plon = grid2_center_lon(src_add)

          !***
          !*** find closest square of grid points on source grid
          !***

          call grid_search_closest(src_add, &
                               plat, plon,  &
                               grid1_center_lat, grid1_center_lon, &
                               grid1_mask, bin_addr1)
          if (dst_add.eq.0) cycle

          if (norm_opt .eq. norm_opt_frcarea) then
            if (.not.luse_grid1_area.or..not.luse_grid2_area) & 
               stop 'luse_grid1_area and luse_grid2_area need to be true for  norm_opt = frcarea '
                    !'
            wgts = grid1_area_in(src_add) / grid2_area_in(dst_add)

          elseif (norm_opt .eq. norm_opt_none) then
            ! simply put one in tha case
            wgts = 1.
          else ; stop 'unknown norm_opt for fouintain interpolation'
          endif
         
          call store_link_fountain(dst_add, src_add, wgts, nmap)
         
        endif

      end do grid_loop2
      where (grid1_add_map2.eq.0) grid1_add_map2 = grid2_size ! NEMO need src .ne. to zeros. Put it at grid2_size put no worries, wgt()=0.

      endif ! nmap=2

!-----------------------------------------------------------------------

      end subroutine remap_fountain

!***********************************************************************

      subroutine grid_search_closest(src_add, &
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

      integer (kind=int_kind), intent(out) :: &
              src_add  ! address of each corner point enclosing P


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
            dist_min, dist_max,distance ! for computing dist

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

          src_add = srch_add

          !***
          !*** otherwise move on to next cell
          !***

        endif !bounding box check
      end do srch_loop
      if (src_add.ne.0) then
        if(src_mask(src_add)==-1) src_add=0 !src_mask==-1   means on the out-of-bound catcher
      endif

!-----------------------------------------------------------------------

      end subroutine grid_search_closest 

!***********************************************************************

      subroutine store_link_fountain(dst_add, src_add, weights, nmap)

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

      integer (kind=int_kind),  intent(in) :: &
              src_add   ! addresses on source grid

      real (kind=dbl_kind),  intent(in) :: &
              weights ! array of remapping weights for these links

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n,  & ! dummy index
             num_links_add          ! placeholder for old link number

!-----------------------------------------------------------------------
!
!     increment number of links and check to see if remap arrays need
!     to be increased to accomodate the new link.  then store the
!     link.
!
!-----------------------------------------------------------------------

      select case (nmap)
      case(1)

          num_links_add = 3
          if (grid1_add_map1(dst_add*4-num_links_add).ne.0) num_links_add = 2
          if (grid1_add_map1(dst_add*4-num_links_add).ne.0) num_links_add = 1
          if (grid1_add_map1(dst_add*4-num_links_add).ne.0) num_links_add = 0
          if (grid1_add_map1(dst_add*4-num_links_add).ne.0) &
                stop 'More than 4 fountain on one cell' 
          num_links_map1 = max_links_map1 
          grid1_add_map1(dst_add*4-num_links_add)= src_add
          grid2_add_map1(dst_add*4-num_links_add)= dst_add
          wts_map1    (1,dst_add*4-num_links_add)=  weights

      case(2)

          num_links_add = 3
          if (grid1_add_map2(dst_add*4-num_links_add).ne.0) num_links_add = 2
          if (grid1_add_map2(dst_add*4-num_links_add).ne.0) num_links_add = 1
          if (grid1_add_map2(dst_add*4-num_links_add).ne.0) num_links_add = 0
          if (grid1_add_map2(dst_add*4-num_links_add).ne.0) &
                stop 'More than 4 fountain on one cell' 
          num_links_map2 = max_links_map2 
          grid1_add_map2(dst_add*4-num_links_add)= src_add
          grid2_add_map2(dst_add*4-num_links_add)= dst_add
          wts_map2    (1,dst_add*4-num_links_add)=  weights


      end select

!-----------------------------------------------------------------------

      end subroutine store_link_fountain

!***********************************************************************

      end module remap_fountainhead

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
