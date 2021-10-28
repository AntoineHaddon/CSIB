
MODULE sbcspp
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
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE iom            ! xIOS server
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC tra_spp
   PUBLIC sbc_spp_div

   LOGICAL,  PUBLIC :: ln_vertspp = .false.   ! If true, use the salt plume parameterization
   LOGICAL,  PUBLIC :: ln_spp_c_grad = .true. ! If true, the density criterion is a local gradient as opposed
                                              ! to a density change from the surface
   REAL(wp), PUBLIC :: rn_spp_rho_c = 0.02_wp ! Density criterion to determine mixed layer depth
   INTEGER , PUBLIC :: nn_power     = 5       ! Affects the shape of the power law. 0: uniform distribution
   REAL(wp), PUBLIC :: rn_spp_z_max = 100.   ! Maximum depth of the salt plume

#  include "vectopt_loop_substitute.h90"

CONTAINS

   SUBROUTINE sbc_spp_div( phdivn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_spp_div***
      !!
      !! ** Purpose : Apply the salt flux associated with sea-ice formation as
      !!              to the salinity trend
      !!
      !! ** Method  : As a parameterization of the salt plumes that form
      !!              due to brine rejection during sea-ice formation
      !!              spread the mass flux vertically following a power
      !!              power law distribution after
      !!              Nguyen, A. T., D. Menemenlis, and R. Kwok (2009),
      !!              Improved modeling of the Arctic halocline with a subgrid-scale brine
      !!                 rejection parameterization, JGR
      !!
      !! ** Action  : phdivn
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT( INOUT ) ::   phdivn   ! horizontal divergence
      REAL(wp)                                    ::   zfact    ! Timestepping factor
      INTEGER :: kl_salt_plume
      REAL :: wt

      REAL(wp), DIMENSION(jpk) :: z_power
      INTEGER :: ji, jj, jk

      zfact = 0.5
      DO jj = 2, jpj
         DO ji = fs_2, fs_jpim1
            ! Only distribute salt flux if flux is positive in this or the previous time step. Note that
            ! this could lead to a freshening at depth if sfx + sfx_b < 0., but is necessary to ensure
            ! symmetry in the leap frog timestepping
            IF (fmmflx_b(ji,jj) > 0. .or. fmmflx(ji,jj) > 0.) THEN
               CALL spp_coeffs_col( ji, jj, kl_salt_plume, z_power )
               wt = (zfact*r1_rau0)*(fmmflx_b(ji,jj)+fmmflx(ji,jj))

               ! Distribute divergence in the vertical
               DO jk = 1,kl_salt_plume
                  phdivn(ji,jj,jk) = phdivn(ji,jj,jk) + (wt*z_power(jk))/e3t_n(ji,jj,jk)
               END DO
            ELSE
               phdivn(ji,jj,1) = phdivn(ji,jj,1) + &
                               & ( fmmflx(ji,jj) + fmmflx_b(ji,jj) ) * zfact * r1_rau0 / e3t_n(ji,jj,1)
            ENDIF
         END DO
      END DO

   END SUBROUTINE sbc_spp_div

   SUBROUTINE spp_coeffs_col( ji, jj, kl_salt_plume, z_power )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_coeffs_col***
      !!
      !! ** Purpose : Calculate the weighting factors consistent with the
      !!              power law distribution proposed in Nguyen et al.
      !!
      !! ** Method  : As a parameterization of the salt plumes that form
      !!              due to brine rejection during sea-ice formation
      !!              spread the mass flux vertically following a power
      !!              power law distribution after
      !!              Nguyen, A. T., D. Menemenlis, and R. Kwok (2009),
      !!              Improved modeling of the Arctic halocline with a subgrid-scale brine
      !!                 rejection parameterization, JGR
      !!
      !! ** Action  : kl_salt_plume, z_power
      !!----------------------------------------------------------------------

      INTEGER, INTENT(IN   ) :: ji,jj         ! Indices of the column
      INTEGER, INTENT(  OUT) :: kl_salt_plume ! Index of last layer within the salt plume
      REAL   , DIMENSION(jpk), INTENT(  OUT) :: z_power       ! Weighting factor used to distribute flux

      REAL(wp) :: density_criterion, h_salt_plume, n2_crit, z_crit
      REAL(wp) :: z_power_sum

      INTEGER :: jk, ki_salt_plume

      z_crit = MIN(rn_spp_z_max,gdepw_n(ji,jj,mbkt(ji,jj)))
      ! Determine the depth of the salt plume based on either a local gradient density criterion
      ! or density difference from the surface

      n2_crit = grav*rn_spp_rho_c*r1_rau0
      IF (ln_spp_c_grad) THEN
         DO jk = 2,jpk-1
            IF ( rn2b(ji,jj,jk) >= n2_crit .or. gdepw_n(ji,jj,jk) >= z_crit) THEN
               ki_salt_plume = jk
               exit
            ENDIF
         ENDDO
      ELSE
         density_criterion = 0.
         DO jk=2,jpk-1
            density_criterion = density_criterion + MAX(rn2b(ji,jj,jk), 0.)*e3w_n(ji,jj,jk)
            IF ( density_criterion >= n2_crit .or. gdepw_n(ji,jj,jk) >= z_crit ) THEN
               ki_salt_plume = jk
               exit
            ENDIF
         ENDDO
      ENDIF
      kl_salt_plume = ki_salt_plume - 1

      ! Calculate coefficient used in the (Eq. 9) by discretizing the constraint in Eq. 10
      z_power_sum = 0.
      DO jk=1,kl_salt_plume
        z_power(jk) = gdept_n(ji,jj,jk)**nn_power
        z_power_sum = z_power_sum + z_power(jk)
      ENDDO

      ! Normalize the weighting factor
      DO jk=1,kl_salt_plume
         z_power(jk) = z_power(jk)/z_power_sum
      ENDDO

   END SUBROUTINE spp_coeffs_col

END MODULE sbcspp