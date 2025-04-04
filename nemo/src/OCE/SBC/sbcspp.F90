
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

   PUBLIC sbc_spp_div

   LOGICAL,  PUBLIC :: ln_vertspp = .false.   ! If true, use the salt plume parameterization
   LOGICAL,  PUBLIC :: ln_spp_c_grad = .true. ! If true, the density criterion is a local gradient as opposed
                                              ! to a density change from the surface
   REAL(wp), PUBLIC :: rn_spp_rho_c = 0.02_wp ! Density criterion to determine mixed layer depth
   INTEGER , PUBLIC :: nn_power     = 5       ! Affects the shape of the power law. 0: uniform distribution
   REAL(wp), PUBLIC :: rn_spp_z_max = 100.   ! Maximum depth of the salt plume
   REAL(wp), PUBLIC :: rn_spp_z_min = 10.   ! Maximum depth of the salt plume

!#  include "vectopt_loop_substitute.h90"
#  include "domzgr_substitute.h90"

CONTAINS

   SUBROUTINE sbc_spp_div( Kmm, phdivn )
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
      INTEGER, INTENT(IN   ) :: Kmm           ! Time
      REAL(wp), DIMENSION(:,:,:), INTENT( INOUT ) ::   phdivn   ! horizontal divergence
      REAL(wp)                                    ::   zfact    ! Timestepping factor
      INTEGER :: kl_salt_plume
      REAL :: wt

      REAL(wp), DIMENSION(jpk) :: z_power
      INTEGER :: ji, jj, jk

      zfact = 0.5
      DO jj = 2, jpj
         DO ji = 2, jpi
            ! Only distribute salt flux if flux is positive in this or the previous time step. Note that
            ! this could lead to a freshening at depth if sfx + sfx_b < 0., but is necessary to ensure
            ! symmetry in the leap frog timestepping
            IF (fmmflx(ji,jj) > 0.) THEN
               CALL spp_coeffs_col( ji, jj,Kmm, kl_salt_plume, z_power )
               wt = r1_rho0*fmmflx(ji,jj)

               ! Distribute divergence in the vertical
               DO jk = 1,kl_salt_plume
                  phdivn(ji,jj,jk) = phdivn(ji,jj,jk) + (wt*z_power(jk))/e3t(ji,jj,jk,Kmm)
               END DO
            ELSE
               phdivn(ji,jj,1) = phdivn(ji,jj,1) + &
                               & fmmflx(ji,jj) * r1_rho0 / e3t(ji,jj,1,Kmm)
            ENDIF
         END DO
      END DO

   END SUBROUTINE sbc_spp_div

   SUBROUTINE spp_coeffs_col( ji, jj,Kmm, kl_salt_plume, z_power )
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
      INTEGER, INTENT(IN   ) :: Kmm           ! Time
      INTEGER, INTENT(  OUT) :: kl_salt_plume ! Index of last layer within the salt plume
      REAL   , DIMENSION(jpk), INTENT(  OUT) :: z_power       ! Weighting factor used to distribute flux

      REAL(wp) :: rhoc, h_salt_plume, n2_crit, z_crit_min, z_crit_max
      REAL(wp) :: z_power_sum

      INTEGER :: jk, ki_salt_plume, ki_z_min, ki_z_max

      ! Set a maximum bound to the depth of the salt plume
      z_crit_max = MIN(rn_spp_z_max, gdepw(ji,jj,mbkt(ji,jj),Kmm))
      ! Set a minimum bound on the depth of the salt plume
      z_crit_min = MIN(rn_spp_z_min, gdepw(ji,jj,mbkt(ji,jj),Kmm))

      n2_crit = grav*rn_spp_rho_c*r1_rho0
      ki_salt_plume = 2
      IF (ln_spp_c_grad) THEN
         DO jk=2,jpk
            IF ( gdepw(ji,jj,jk,Kmm) < z_crit_min ) THEN
               CYCLE
            ELSEIF ( rn2b(ji,jj,jk) >= n2_crit .or. gdepw(ji,jj,jk+1,Kmm) >= z_crit_max ) THEN
               ki_salt_plume = jk
               EXIT
            ENDIF
         ENDDO
      ELSE
         rhoc = 0.
         DO jk=2,jpk
            rhoc = rhoc + MAX(rn2b(ji,jj,jk), 0.)*e3w(ji,jj,jk,Kmm)
            IF ( gdepw(ji,jj,jk,Kmm) < z_crit_min ) THEN
               CYCLE
            ELSEIF ( rhoc >= n2_crit           .or. gdepw(ji,jj,jk+1,Kmm) >= z_crit_max ) THEN
               ki_salt_plume = jk
               EXIT
            ENDIF
         ENDDO
      ENDIF
      kl_salt_plume = ki_salt_plume - 1

      ! Calculate coefficient used in the (Eq. 9) by discretizing the constraint in Eq. 10
      z_power_sum = 0.
      DO jk=1,kl_salt_plume
        z_power(jk) = gdept(ji,jj,jk,Kmm)**nn_power
        z_power_sum = z_power_sum + z_power(jk)
      ENDDO

      ! Normalize the weighting factor
      DO jk=1,kl_salt_plume
         z_power(jk) = z_power(jk)/z_power_sum
      ENDDO

   END SUBROUTINE spp_coeffs_col

END MODULE sbcspp
