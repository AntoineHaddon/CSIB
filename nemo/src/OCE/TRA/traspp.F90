
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

   LOGICAL,  PUBLIC :: ln_vertspp = .false.   ! If true, use the salt plume parameterization
   LOGICAL,  PUBLIC :: ln_spp_c_grad = .true. ! If true, the density criterion is a local gradient as opposed
                                              ! to a density change from the surface
   REAL(wp), PUBLIC :: rn_spp_rho_c = 0.02_wp ! Density criterion to determine mixed layer depth
   INTEGER , PUBLIC :: nn_power     = 5       ! Affects the shape of the power law. 0: uniform distribution
   REAL(wp), PUBLIC :: rn_spp_z_max = 1000.   ! Maximum depth of the salt plume

#  include "vectopt_loop_substitute.h90"

CONTAINS

   SUBROUTINE tra_spp( kt, zfact, nn_sfx_sign )
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
      INTEGER,  INTENT(IN) :: kt
      REAL(wp), INTENT(IN) :: zfact
      INTEGER,  INTENT(IN) :: nn_sfx_sign

      REAL(wp) :: wt, density_criterion, h_salt_plume, n2_crit, z_crit
      REAL(wp), DIMENSION(jpk) :: z_power, tend_col
      REAL(wp) :: z_power_sum
      real(wp), DIMENSION(jpi,jpj,jpk) :: spp_tend_3d
      real(wp), DIMENSION(jpi,jpj)     :: spp_thick

      INTEGER :: ji, jj, jk
      INTEGER :: ki_salt_plume ! Interface index of the salt plume
      INTEGER :: kl_salt_plume ! Index of last layer within the salt plume

      ! Convert density criterion to an equivalent N2 criterion
      n2_crit = grav*rn_spp_rho_c*r1_rau0
      IF (iom_use("spp_tend")) THEN
         spp_tend_3d(:,:,:) = 0.
      ENDIF
      IF (iom_use("spp_thick")) THEN
         spp_thick(:,:) = 0.
      ENDIF
      DO jj = 2, jpj
         DO ji = fs_2, fs_jpim1
            ! Only distribute salt flux if ice is being formed in this or the previous time step. Note that
            ! this could lead to a freshening at depth if sfx + sfx_b < 0., but is necessary to ensure
            ! symmetry in the leap frog timestepping

            IF (nn_sfx_sign*sfx(ji,jj) > 0. .or. nn_sfx_sign*sfx_b(ji,jj) > 0.) THEN
               z_crit = MIN(rn_spp_z_max,gdepw_n(ji,jj,mbkt(ji,jj)))
               ! Determine the depth of the salt plume based on either a local gradient density criterion
               ! or density difference from the surface
               tend_col(:) = 0.

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
               IF (iom_use("spp_thick")) THEN
                  spp_thick(ji,jj) = gdepw_n(ji,jj,ki_salt_plume)
               ENDIF
               kl_salt_plume = ki_salt_plume - 1

               ! Calculate coefficient used in the (Eq. 9) by discretizing the constraint in Eq. 10
               z_power_sum = 0.
               do jk=1,kl_salt_plume
                 z_power(jk) = gdept_n(ji,jj,jk)**nn_power
                 z_power_sum = z_power_sum + z_power(jk)
               enddo
               wt = nn_sfx_sign*zfact*(sfx_b(ji,jj)+sfx(ji,jj))
               wt = wt/z_power_sum
               wt = wt*r1_rau0

               ! Distribute tendencies in the vertical
               DO jk = 1,kl_salt_plume
                  tend_col(jk) = wt*z_power(jk)
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) + tend_col(jk)/e3t_n(ji,jj,jk)
               END DO

               IF (iom_use("spp_tend")) THEN
                  spp_tend_3d(ji,jj,:) = tend_col(:)
               ENDIF
            ENDIF
         END DO
      END DO

      CALL iom_put("spp_tend" , spp_tend_3d)
      CALL iom_put("spp_thick", spp_thick)

   END SUBROUTINE tra_spp

END MODULE traspp