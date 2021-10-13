
MODULE traspp
   !!==============================================================================
   !!                       ***  MODULE  trasbc  ***
   !! Ocean active tracers: Parameterization of salt plumes due to brine rejection
   !!==============================================================================
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_sbc       : update the salinity trend in the ocean surface boundary layer
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE sbc_oce        ! surface boundary condition: ocean
   USE dom_oce        ! ocean space domain variables
   USE phycst         ! physical constant
   USE sbcmod         ! ln_rnf
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE iom            ! xIOS server
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE timing         ! Timing
   USE zdfmxl, only : nmln, hmlp, zdf_mxl

   REAL(wp) :: rn_spp_rho_c = 0.02_wp ! Density criterion to determine mixed layer depth
   INTEGER  :: nn_power     = 5       ! Affects the shape of the power law. 0: uniform distribution

CONTAINS

   SUBROUTINE tra_spp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdfmxl  ***
      !!
      !! ** Purpose : Distribute the salt flux due to sea-ice formation
      !!              within the surface boundary layer
      !!
      !!
      !! ** Method  : As a parameterization of the salt plumes that form
      !!              due to brine rejection during sea-ice formation
      !!              spread the salt flux vertically following a power
      !!              power law distribution after
      !!              Nguyen, A. T., D. Menemenlis, and R. Kwok (2009),
      !!              Improved modeling of the Arctic halocline with a subgrid-scale brine
      !!                 rejection parameterization, JGR
      !!
      !! ** Action  : tsa(:,:,:,jp_sal)

      real(wp) :: power_z
      real(wp) :: dimension(jpk)

      call zdf_mxl( kt, rn_spp_rho_c )



      call zdf_mxl( kt )


   END SUBROUTINE tra_spp

END MODULE traspp