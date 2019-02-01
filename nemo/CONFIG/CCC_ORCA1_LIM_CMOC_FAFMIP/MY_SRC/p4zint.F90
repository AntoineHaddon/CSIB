MODULE p4zint
   !!======================================================================
   !!                         ***  MODULE p4zint  ***
   !! TOP :   PISCES interpolation and computation of various accessory fields
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
#if defined key_pisces
   USE oce_trc       !  shared variables between ocean and passive tracers
   USE trc           !  passive tracers common variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   update_salt_avg_2d
   PUBLIC   glob_avg_area_wt 

   CONTAINS

   !> Abiotic tracers rely on the estimate of globally and annually averaged surface salinity from the previous year (or
   !! average of the initial salinity for the first year). To reduce the number of global sums, a 2D array of the
   !! thickness-weighted column averaged salinity is stored. The time-averaging is done using the online calculation
   !! described by West (1979) so that salt_annual_2d always is a reasonable value. This algorithm should also be
   !! robust to non-uniform timestep lengths
   SUBROUTINE update_salt_avg_2d( salt_surf, dt, salt_avg, dt_sum )
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(IN   ) :: salt_surf !< Current value of surface salinity
      REAL(wp),                         INTENT(IN   ) :: dt        !< Timestep length used to time weight 'salt'
      REAL(wp), DIMENSION(jpi,jpj),     INTENT(INOUT) :: salt_avg  !< Array containing running value of salinity
      REAL(wp),                         INTENT(INOUT) :: dt_sum    !< Time elapsed so far

      INTEGER  :: ji, jj, jk
      REAL(wp) :: wt
      
      ! Update the stream mean
      dt_sum = dt_sum + dt
      wt = dt / dt_sum
      DO jj = 1,jpj
         DO ji = 1,jpi
            salt_avg(ji,jj) = salt_avg(ji,jj) + wt*(salt_surf(ji,jj) - salt_avg(ji,jj))
         ENDDO
      ENDDO
      
   END SUBROUTINE update_salt_avg_2d

   !> Calculate the global average from the 2D column average fields
   REAL(wp) FUNCTION glob_avg_area_wt( arr2d )
      REAL(wp), DIMENSION(jpi,jpj) :: arr2d 

      glob_avg_area_wt = glob_sum( arr2d(:,:)*e1e2t(:,:) ) / glob_sum(e1e2t(:,:)) 

   END FUNCTION glob_avg_area_wt 
#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_int                   ! Empty routine
      WRITE(*,*) 'p4z_int: You should not have seen this print! error?'
   END SUBROUTINE p4z_int
#endif 

   !!======================================================================
END MODULE  p4zint
