MODULE checksum
   !!======================================================================
   !!                       ***  MODULE checksums  ***
   !! Utility functions: Calculate checksums of arrays. Robust across PE
   !!                    count and domain decomposition
   !!
   !!======================================================================
   !! History :  3.4  ! Initial implementation  2018-04  (A.Shao)
   !!----------------------------------------------------------------------
#if defined key_example
   !!----------------------------------------------------------------------
   !!   'key_example'  :                brief description of the key option
   !!----------------------------------------------------------------------
   !!   exa_mpl       : liste of module subroutine (caution, never use the
   !!   exa_mpl_init  : name of the module for a routine)
   !!   exa_mpl_stp   : Please try to use 3 letter block for routine names
   !!----------------------------------------------------------------------

   USE libmpp,  only :: mpp_min, mpp_sum
   USE par_oce, only :: jpi, jpj, jpim1, jpjm1, jpk

   IMPLICIT NONE
   PRIVATE

   PUBLIC  :: chksum
   INTERFACE chksum
      MODULE PROCEDURE chksum_2d, chksum_3d
   END INTERFACE

   INTEGER, PARAMETER :: bitlen = 1000000000 !< Length of the checksum

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: module_example 2737 2011-04-11 10:30:51Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
  !> Does a bit count and a reproducing sum on a 2D array to keep track of the global array
  !! as NEMO integrates. Primarily intended for use as a debugging tool.
  SUBROUTINE chksum_2d( name, array, mask, istart, iend, jstart, jend )
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
       WRITE(*,'(A,X,I10,X,A,E25.16,X,A,E25.16') &
             name, "chksum=", bc, "Global minimum=", minarray, "Global maximum=", maxarray
     ENDIF

  END SUBROUTINE ckhsum_2d

  !> Does a bit count and a reproducing sum on a 3D array to keep track of the global array
  !! as NEMO integrates. Primarily intended for use as a debugging tool.
  SUBROUTINE chksum_3d( name, array, mask, istart, iend, jstart, jend, kstart, kend )
     CHARACTER(LEN=*)                           :: name     !< Name of the array to be checksummed
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
     minarray = array(is,js) ; maxarray = array(is,js)

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
             name, "chksum=", bc, "Global minimum=", minarray, "Global maximum=", maxarray
     ENDIF

  END SUBROUTINE ckhsum_3d

   !!======================================================================
END MODULE checksums
