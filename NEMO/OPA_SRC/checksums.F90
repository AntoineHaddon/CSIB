MODULE checksum
   !!======================================================================
   !!                       ***  MODULE checksums  ***
   !! Utility functions: Calculate checksums of arrays. Robust across PE
   !!                    count and domain decomposition
   !!
   !!======================================================================
   !! History :  3.4  ! Initial implementation  2018-04  (A.Shao)
   !!----------------------------------------------------------------------
   !!   chksum   :  Calculates a checksum of a given 2d/3d array
   !!----------------------------------------------------------------------
   USE dom_oce       : tmask, umask, vmask
   USE libmpp,  only : mpp_min, mpp_sum
   USE oce,     only : un, vn, wn, tsn, ua, va, tsa, ub, vb, tsb
   USE par_oce, only : jpi, jpj, jpim1, jpjm1, jpk, jp_tem
   IMPLICIT NONE
   PRIVATE

   PUBLIC    chksum
   PUBLIC    before_state_chksum, now_state_chksum, after_state_chksum
   PUBLIC    now_ts_chksum, after_ts_chksum
   INTERFACE chksum
      MODULE PROCEDURE chksum_2d, chksum_3d
   END INTERFACE

   LOGICAL, PUBLIC    :: nn_chksum = .FALSE. !< If true, calculate and write checksums

   ! Module variables
   INTEGER, PARAMETER :: bitlen = 1000000000 !< Length of the checksum
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: module_example 2737 2011-04-11 10:30:51Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   !> Does a bit count and a reproducing sum on a 2D array to keep track of the global array
   !! as NEMO integrates. Primarily intended for use as a debugging tool.
   SUBROUTINE chksum_2d( msg, array, mask, istart, iend, jstart, jend )
      CHARACTER(LEN=*)                       :: msg      !< Prefix phrase for print message, usually
                                                         !! 'varname before/after process'
      REAL(wp), DIMENSION(jpi,jpj)           :: array    !< Array to be checksummed
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL :: mask     !< Array to be checksummed
      INTEGER, OPTIONAL                      :: istart   !< First index to use along the i-axis
      INTEGER, OPTIONAL                      :: iend     !< Last index to use along the i-axis
      INTEGER, OPTIONAL                      :: jstart   !< First index to use along the i-axis
      INTEGER, OPTIONAL                      :: jend     !< Last index to use along the i-axis
      ! Local variables
      INTEGER  :: ji, jj         ! Loop variables
      INTEGER  :: is, ie, js, je ! Beginning and end start indices
      INTEGER  :: bc             ! bitcount of array
      REAL(wp) :: minarray, maxarray

      ! By default set indices to only the "inner" part of the array, not the halo
      is = 2 ; ie = jpim1 ; js = 2  ; je = jpjm1
      IF (PRESENT(istart)) is = istart
      IF (PRESENT(jstart)) js = jstart
      IF (PRESENT(iend))   ie = iend
      IF (PRESENT(jend))   je = jend

      ! Set the initial min/max values to the first value to be checked in the array
      minarray = array(is,js) ; maxarray = array(is,js)

      bc = 0
      IF (PRESENT(mask)) THEN
         DO ji = is, ie ; DO jj = js, je
            IF ( mask(ji,jj) > 0. ) THEN
               minarray = MIN( minarray, array(ji,jj) )
               maxarray = MAX( maxarray, array(ji,jj) )
               bc = bc + bitcount( array(ji,jj) )
            ENDIF
         ENDDO ; ENDDO
      ELSE
         DO ji = is, ie ; DO jj = js, je
            minarray = MIN) minarray, array(ji,jj) )
            maxarray = MAX( maxarray, array(ji,jj) )
            bc = bc + bitcount( array(ji,jj) )
         ENDDO ; ENDDO
      ENDIF

      CALL mpp_sum( bc )
      CALL mpp_min( minarray )
      CALL mpp_max( maxarray )

      bc = mod(bc, bitlen)

      IF (lwp) THEN
        WRITE(*,'(A,X,I10,X,A,ES25.16,X,A,ES25.16') &
              TRIM(msg), "chksum=", bc, "Global minimum=", minarray, "Global maximum=", maxarray
      ENDIF

   END SUBROUTINE chksum_2d

   !> Does a bit count and a reproducing sum on a 3D array to keep track of the global array
   !! as NEMO integrates. Primarily intended for use as a debugging tool.
   SUBROUTINE chksum_3d( msg, array, mask, istart, iend, jstart, jend, kstart, kend )
      CHARACTER(LEN=*)                           :: msg      !< Name of the array to be checksummed
      REAL(wp), DIMENSION(jpi,jpj,jpk)           :: array    !< Array to be checksummed
      REAL(wp), DIMENSION(jpi,jpj,jpk), OPTIONAL :: mask     !< Array to be checksummed
      INTEGER, OPTIONAL                          :: istart   !< First index to use along the i-axis
      INTEGER, OPTIONAL                          :: iend     !< Last index to use along the i-axis
      INTEGER, OPTIONAL                          :: jstart   !< First index to use along the i-axis
      INTEGER, OPTIONAL                          :: jend     !< Last index to use along the i-axis
      INTEGER, OPTIONAL                          :: kstart   !< First index to use along the i-axis
      INTEGER, OPTIONAL                          :: kend     !< Last index to use along the i-axis
      ! Local variables
      INTEGER  :: ji, jj, jk             ! Loop variables
      INTEGER  :: is, ie, js, je, ks, ke ! Beginning and end start indices
      INTEGER  :: bc                     ! bitcount of array
      REAL(wp) :: minarray, maxarray

      ! By default set indices to only the "inner" part of the array, not the halo
      is = 2 ; ie = jpim1 ; js = 2  ; je = jpjm1 ; ks = 1 ; ke = jpk
      IF (PRESENT(istart)) is = istart
      IF (PRESENT(jstart)) js = jstart
      IF (PRESENT(kstart)) js = kstart
      IF (PRESENT(iend))   ie = iend
      IF (PRESENT(jend))   je = jend
      IF (PRESENT(kend))   ke = kend

      ! Set the initial min/max values to the first value to be checked in the array
      minarray = array(is,js,ks) ; maxarray = array(is,js,ks)
      bc = 0
      IF (PRESENT(mask)) THEN
        DO jk = ks, ke ; DO ji = is, ie ; DO jj = js, je
           IF (mask(ji,jj,k)>0.) THEN
              minarray = MIN(minarray, array(ji,jj,jk))
              maxarray = MAX(maxarray, array(ji,jj,jk))
              bc = bc + bitcount( array(ji,jj,jk) )
           ENDIF
        ENDDO ; ENDDO ; ENDDO
      ELSE
        DO jk = ks, ke ; DO ji = is, ie ; DO jj = js, je
           minarray = MIN(minarray, array(ji,jj,jk))
           maxarray = MAX(maxarray, array(ji,jj,jk))
           bc = bc + bitcount( array(ji,jj,jk) )
        ENDDO ; ENDDO ; ENDDO
      ENDIF

      CALL mpp_sum( bc )
      CALL mpp_min( minarray )
      CALL mpp_max( maxarray )

      bc = mod(bc, bitlen)

      IF (lwp) THEN
        WRITE(*,'(A,X,I10,X,A,E25.16,X,A,E25.16') &
              TRIM(msg), "chksum=", bc, "Global minimum=", minarray, "Global maximum=", maxarray
      ENDIF

   END SUBROUTINE chksum_3d
   !!!! Full state checksums
   !> Convenience routine to do a chksum of current state of the model 'after' arrays of u,v,w,T,S
   SUBROUTINE before_state_chksum(msg)
      CHARACTER(LEN=*) :: msg !< The point of the algorithm that the checksum is being done

      CALL chksum( "u before array "//TRIM(msg), ub, umask)
      CALL chksum( "v before array "//TRIM(msg), vb, vmask)
      CALL chksum( "T before array "//TRIM(msg), tsb(:,:,:,jp_tem), tmask)
      CALL chksum( "S before array "//TRIM(msg), tsb(:,:,:,jp_sal), tmask)

   END SUBROUTINE before_state_chksum

   !> Convenience routine to do a chksum of current state of the model 'now' arays of u,v,w,T,S
   SUBROUTINE now_state_chksum(msg)
      CHARACTER(LEN=*) :: msg !< The point of the algorithm that the checksum is being done

      CALL chksum( "u now array "//TRIM(msg), un, umask)
      CALL chksum( "v now array "//TRIM(msg), vn, vmask)
      CALL chksum( "w now array "//TRIM(msg), wn, tmask)
      CALL chksum( "T now array "//TRIM(msg), tsn(:,:,:,jp_tem), tmask)
      CALL chksum( "S now array "//TRIM(msg), tsn(:,:,:,jp_sal), tmask)

   END SUBROUTINE now_state_chksum

   !> Convenience routine to do a chksum of current state of the model 'after' arrays of u,v,w,T,S
   SUBROUTINE after_state_chksum(msg)
      CHARACTER(LEN=*) :: msg !< The point of the algorithm that the checksum is being done

      CALL chksum( "u after array "//TRIM(msg), ua, umask)
      CALL chksum( "v after array "//TRIM(msg), va, vmask)
      CALL chksum( "T after array "//TRIM(msg), tsa(:,:,:,jp_tem), tmask)
      CALL chksum( "S after array "//TRIM(msg), tsa(:,:,:,jp_sal), tmask)

   END SUBROUTINE after_state_chksum
   !!!! T/S cchecksums
   !> Convenience routine to do a chksum of current state of the model 'now' arays of T,S
   SUBROUTINE now_ts_chksum(msg)
      CHARACTER(LEN=*) :: msg !< The point of the algorithm that the checksum is being done

      CALL chksum( "T now array "//TRIM(msg), tsn(:,:,:,jp_tem), tmask)
      CALL chksum( "S now array "//TRIM(msg), tsn(:,:,:,jp_sal), tmask)

   END SUBROUTINE now_ts_chksum

   !> Convenience routine to do a chksum of current state of the model 'after' arrays of T,S
   SUBROUTINE after_ts_chksum(msg)
      CHARACTER(LEN=*) :: msg !< The point of the algorithm that the checksum is being done

      CALL chksum( "T after array "//TRIM(msg), tsa(:,:,:,jp_tem), tmask)
      CALL chksum( "S after array "//TRIM(msg), tsa(:,:,:,jp_sal), tmask)

   END SUBROUTINE after_ts_chksum
   !!======================================================================
END MODULE checksums
