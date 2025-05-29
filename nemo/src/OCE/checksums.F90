MODULE checksums
   !!======================================================================
   !!                       ***  MODULE checksums  ***
   !! Utility functions: Calculate checksums of arrays. Robust across PE
   !!                    count and domain decomposition
   !!
   !!======================================================================
   !! History :  3.4  ! Initial implementation  2018-04  (A.Shao)
   !!            4.2  ! Copied (and shortened) from 3.4 implementation 
   !!                                               2025-05 (J. Izett)
   !!----------------------------------------------------------------------
   !!   chksum   :  Calculates a checksum of a given 2d/3d array
   !!----------------------------------------------------------------------
   USE par_kind,       only : wp
   USE dom_oce,        only : tmask, umask, vmask, narea
   USE lib_mpp,        only : mpp_min, mpp_max, mpp_sum
   USE oce,            only : uu, vv, ww, ts
   USE par_oce,        only : jpi, jpj, jpk, jp_tem, jp_sal
   USE in_out_manager, only : lwp
   IMPLICIT NONE
   PRIVATE

   PUBLIC    chksum
   PUBLIC    now_state_chksum
   INTERFACE chksum
      MODULE PROCEDURE chksum_2d, chksum_3d
   END INTERFACE

   LOGICAL, PUBLIC    :: ln_chksum = .FALSE. !< If true, calculate and write checksums

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
   SUBROUTINE chksum_2d( array, mask, istart, iend, jstart, jend, alt_unit, msg, bc_out )
      REAL(wp), DIMENSION(jpi,jpj)          , INTENT(IN   ) :: array    !< Array to be checksummed
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(IN   ) :: mask     !< Array to be checksummed
      INTEGER, OPTIONAL                     , INTENT(IN   ) :: istart   !< First index to use along the i-axis
      INTEGER, OPTIONAL                     , INTENT(IN   ) :: iend     !< Last index to use along the i-axis
      INTEGER, OPTIONAL                     , INTENT(IN   ) :: jstart   !< First index to use along the i-axis
      INTEGER, OPTIONAL                     , INTENT(IN   ) :: jend     !< Last index to use along the i-axis
      INTEGER, OPTIONAL                     , INTENT(IN   ) :: alt_unit !< Where to write the message 
      CHARACTER(LEN=*), OPTIONAL            , INTENT(IN   ) :: msg      !< Prefix phrase for print message
      INTEGER, OPTIONAL                     , INTENT(  OUT) :: bc_out   !< Return the bitcount if requested 
      ! Local variables
      INTEGER  :: ji, jj         ! Loop variables
      INTEGER  :: is, ie, js, je ! Beginning and end start indices
      INTEGER  :: bc             ! bitcount of array
      INTEGER  :: write_unit
      REAL(wp) :: minarray, maxarray

      ! By default set indices to only the "inner" part of the array, not the halo
      is = jpi+1 ; ie = jpi-1
      js = jpj+1 ; je = jpj-1
      IF (PRESENT(istart)) is = istart
      IF (PRESENT(jstart)) js = jstart
      IF (PRESENT(iend))   ie = iend
      IF (PRESENT(jend))   je = jend
      ! By default write ot the screen
      write_unit = 6 ; IF( PRESENT(alt_unit) ) write_unit = alt_unit
      
      ! Set the initial min/max values to be ridiculous values
      minarray = HUGE(minarray) ; maxarray = -HUGE(maxarray)

      bc = 0
      IF (PRESENT(mask)) THEN
         DO jj = js, je ; DO ji = is, ie
            IF ( mask(ji,jj) > 0. ) THEN
               minarray = MIN( minarray, array(ji,jj) )
               maxarray = MAX( maxarray, array(ji,jj) )
               bc = bc + bitcount( array(ji,jj) )
            ENDIF
         ENDDO ; ENDDO
      ELSE
         DO jj = js, je ; DO ji = is, ie
            minarray = MIN( minarray, array(ji,jj) )
            maxarray = MAX( maxarray, array(ji,jj) )
            bc = bc + bitcount( array(ji,jj) )
         ENDDO ; ENDDO
      ENDIF

      CALL mpp_sum( 'checksums', bc )
      CALL mpp_min( 'checksums', minarray )
      CALL mpp_max( 'checksums', maxarray )

      bc = mod(bc, bitlen)

      IF (lwp .AND. PRESENT(msg)) THEN
        WRITE(write_unit,'(A,X,A,I10.10,X,A,ES25.16,X,A,ES25.16)') &
              TRIM(msg), "chksum=", bc, "Global min=", minarray, "Global max=", maxarray
      ENDIF
      IF( PRESENT(bc_out) ) bc_out = bc
   END SUBROUTINE chksum_2d

   !> Does a bit count and a reproducing sum on a 3D array to keep track of the global array
   !! as NEMO integrates. Primarily intended for use as a debugging tool.
   SUBROUTINE chksum_3d( array, mask, istart, iend, jstart, jend, kstart, kend, alt_unit, msg, bc_out )
      REAL(wp), DIMENSION(jpi,jpj,jpk)          , INTENT(IN   ) :: array    !< Array to be checksummed
      REAL(wp), DIMENSION(jpi,jpj,jpk), OPTIONAL, INTENT(IN   ) :: mask     !< Array to be checksummed
      INTEGER, OPTIONAL                         , INTENT(IN   ) :: istart   !< First index to use along the i-axis
      INTEGER, OPTIONAL                         , INTENT(IN   ) :: iend     !< Last index to use along the i-axis
      INTEGER, OPTIONAL                         , INTENT(IN   ) :: jstart   !< First index to use along the i-axis
      INTEGER, OPTIONAL                         , INTENT(IN   ) :: jend     !< Last index to use along the i-axis
      INTEGER, OPTIONAL                         , INTENT(IN   ) :: kstart   !< First index to use along the i-axis
      INTEGER, OPTIONAL                         , INTENT(IN   ) :: kend     !< Last index to use along the i-axis
      INTEGER, OPTIONAL                         , INTENT(IN   ) :: alt_unit !< Where to write the message 
      CHARACTER(LEN=*), OPTIONAL                , INTENT(IN   ) :: msg      !< Name of the array to be checksummed
      INTEGER, OPTIONAL                         , INTENT(  OUT) :: bc_out   !< Return the bitcount if requested 
      ! Local variables
      INTEGER  :: ji, jj, jk             ! Loop variables
      INTEGER  :: is, ie, js, je, ks, ke ! Beginning and end start indices
      INTEGER  :: bc                     ! bitcount of array
      INTEGER  :: write_unit
      REAL(wp) :: minarray, maxarray

      ! By default set indices to only the "inner" part of the array, not the halo
      is = 2 ; ie = jpi-1
      js = 2 ; je = jpj-1
      ks = 1 ; ke = jpk
      IF (PRESENT(istart)) is = istart
      IF (PRESENT(jstart)) js = jstart
      IF (PRESENT(kstart)) ks = kstart
      IF (PRESENT(iend))   ie = iend
      IF (PRESENT(jend))   je = jend
      IF (PRESENT(kend))   ke = kend
      ! By default write to the screen
      write_unit = 6 ; IF( PRESENT(alt_unit) ) write_unit = alt_unit

      ! Set the initial min/max values to be ridiculous values
      minarray = HUGE(minarray) ; maxarray = -HUGE(maxarray)
      bc = 0
      IF (PRESENT(mask)) THEN
        DO jk = ks, ke ; DO jj = js, je ; DO ji = is, ie
           IF (mask(ji,jj,jk)>0.) THEN
              minarray = MIN(minarray, array(ji,jj,jk))
              maxarray = MAX(maxarray, array(ji,jj,jk))
              bc = bc + bitcount( array(ji,jj,jk) )
           ENDIF
        ENDDO ; ENDDO ; ENDDO
      ELSE
        DO jk = ks, ke ; DO jj = js, je ; DO ji = is, ie
           minarray = MIN(minarray, array(ji,jj,jk))
           maxarray = MAX(maxarray, array(ji,jj,jk))
           bc = bc + bitcount( array(ji,jj,jk) )
        ENDDO ; ENDDO ; ENDDO
      ENDIF

      CALL mpp_sum( 'checksums', bc )
      CALL mpp_min( 'checksums', minarray )
      CALL mpp_max( 'checksums', maxarray )

      bc = mod(bc, bitlen)

      IF (lwp .AND. PRESENT(msg)) THEN
        WRITE(write_unit,'(A,X,A,I10.10,X,A,E25.16,X,A,E25.16)') &
              TRIM(msg), "chksum=", bc, "Global min=", minarray, "Global max=", maxarray
      ENDIF
      IF( PRESENT(bc_out) ) bc_out = bc

   END SUBROUTINE chksum_3d
   !!!! Full state checksums
   !> Convenience routine to do a chksum of current state of the model 'now' arays of u,v,w,T,S
   SUBROUTINE now_state_chksum(msg, Nstep, alt_unit, state_bc)
      CHARACTER(LEN=*),  INTENT(IN   ) :: msg !< The point of the algorithm that the checksum is being done
      INTEGER,           INTENT(IN   ) :: Nstep ! The time level
      INTEGER, OPTIONAL, INTENT(IN   ) :: alt_unit
      INTEGER, OPTIONAL, INTENT(  OUT) :: state_bc
      INTEGER :: write_unit, bc, bc_sum

      write_unit = 6; IF( PRESENT(alt_unit) ) write_unit = alt_unit
      bc_sum = 0
      CALL chksum( uu(:,:,:,Nstep), mask = umask, alt_unit=write_unit, bc_out = bc                 , msg = "u now array "//TRIM(msg))
      bc_sum = bc_sum + bc
      CALL chksum( vv(:,:,:,Nstep), mask = vmask, alt_unit=write_unit, bc_out = bc                 , msg = "v now array "//TRIM(msg))
      bc_sum = bc_sum + bc
      CALL chksum( ww(:,:,:), mask = tmask, alt_unit=write_unit, bc_out = bc                 , msg = "w now array "//TRIM(msg))
      bc_sum = bc_sum + bc
      CALL chksum( ts(:,:,:,jp_tem,Nstep), mask = tmask, alt_unit = write_unit, bc_out = bc, msg = "T now array "//TRIM(msg))
      bc_sum = bc_sum + bc
      CALL chksum( ts(:,:,:,jp_sal,Nstep), mask = tmask, alt_unit = write_unit, bc_out = bc, msg = "S now array "//TRIM(msg))
      bc_sum = bc_sum + bc

      IF( PRESENT( state_bc )) state_bc = bc_sum
   END SUBROUTINE now_state_chksum

   !> Calculates the bitcount of a real number by transferring its memory representation to an
   !! integer of the same byte-size and then using BTEST to check it bit by bit
   INTEGER FUNCTION bitcount( scalar )
      REAL(wp)    :: scalar     !< Scalar to be bit-counted
      INTEGER(wp) :: scalar_int !< Integer with memory representation of 'scalar'

      INTEGER :: bit
      bitcount = 0
      scalar_int = TRANSFER(scalar, scalar_int)
      DO bit=0,BIT_SIZE(scalar_int)-1
        IF ( BTEST(scalar_int, bit) ) bitcount = bitcount + 1
      ENDDO

   END FUNCTION bitcount
   !!======================================================================
END MODULE checksums
