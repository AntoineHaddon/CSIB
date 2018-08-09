MODULE trc_util
   !!======================================================================
   !!                      ***  MODULE  trc_util  ***
   !! Passive tracers   : Common utilities for passive tracers
   !!======================================================================
   !! History :
   !!   CanNEMO      3.4.1  ! Add netCDF routines for reading netcdf variables
   !!----------------------------------------------------------------------
   USE lib_mpp,   only : ctl_stop
   USE obs_utils, only : chkerr
   USE par_kind,  only : wp
   USE netcdf

   IMPLICIT NONE
   PRIVATE

   PUBLIC read_var1d
   PUBLIC read_var2d

   CONTAINS
   !> read_var1d: Read a 1d variable from a netcdf file. Allocate the output
   !! array, fill with values, and return
   SUBROUTINE read_var1d(ncid, varname, varout, dimlen)
      INTEGER                            , INTENT(IN   ) :: ncid    !< netcdf file id
      CHARACTER(LEN=*)                   , INTENT(IN   ) :: varname !< name of the variable to be read
      REAL(wp), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: varout  !< Output array
      INTEGER,                             INTENT(  OUT), OPTIONAL :: dimlen  !< Length of the dimension

      INTEGER :: varid, ndim
      INTEGER :: dimlen_loc
      INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: dimid_loc

      ! Get the variable id name and other information about the variable
      CALL chkerr( nf90_inq_varid(ncid, varname, varid), 'trc_util', 0 )
      CALL chkerr( nf90_inquire_variable(ncid, varid, ndims=ndim, dimids=dimid_loc), 'trc_util', 0)
      ! Consistency checks
      IF (ndim /= 1) &
          CALL ctl_stop( 'STOP', 'read_var1d: read_var1d was called for a variable with more than 1 dimension')

      ! Allocate if necessary
      IF (.NOT. ALLOCATED(varout)) THEN
         CALL chkerr(nf90_inquire_dimension(ncid, dimid_loc(1), len = dimlen_loc), 'trc_util', 0)
         ALLOCATE( varout(dimlen_loc) )
      ENDIF

      varout(:) = 0.
      CALL chkerr( nf90_get_var(ncid, varid, varout), 'trc_util', 0 )

      IF (PRESENT(dimlen)) dimlen = dimlen_loc

   END SUBROUTINE read_var1d
   !> read_var2d: Read a 2d variable from a netcdf file. Allocate the output
   !! array, fill with values, and return
   SUBROUTINE read_var2d(ncid, varname, varout, dimlen1, dimlen2)
      INTEGER                              , INTENT(IN   ) :: ncid    !< netcdf file id
      CHARACTER(LEN=*)                     , INTENT(IN   ) :: varname !< name of the variable to be read
      REAL(wp), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: varout  !< Output array
      INTEGER, DIMENSION(2), INTENT(  OUT), OPTIONAL :: dimlen1  !< Length of dimension 1
      INTEGER, DIMENSION(2), INTENT(  OUT), OPTIONAL :: dimlen2  !< Length of dimension 2

      INTEGER :: varid, ndim
      INTEGER :: dimlen_loc1, dimlen_loc2
      INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: dimid_loc

      ! Get the variable id name and other information about the variable
      CALL chkerr( nf90_inq_varid(ncid, varname, varid), 'trc_util', 0 )
      CALL chkerr( nf90_inquire_variable(ncid, varid, ndims=ndim, dimids=dimid_loc), 'trc_util', 0)
      ! Consistency checks
      IF (ndim /= 2) &
          CALL ctl_stop( 'STOP', 'read_var2d: read_var2d was called for a variable with more than 1 dimension')

      ! Allocate if necessary
      IF (.NOT. ALLOCATED(varout)) THEN
         CALL chkerr(nf90_inquire_dimension(ncid, dimid_loc(1), len = dimlen_loc1), 'trc_util', 0)
         CALL chkerr(nf90_inquire_dimension(ncid, dimid_loc(2), len = dimlen_loc2), 'trc_util', 0)
         ALLOCATE( varout(dimlen_loc1,dimlen_loc2) )
      ENDIF

      varout(:,:) = 0.
      CALL chkerr( nf90_get_var(ncid, varid, varout), 'trc_util', 0 )

      IF (PRESENT(dimlen1)) dimlen1 = dimlen_loc1
      IF (PRESENT(dimlen2)) dimlen2 = dimlen_loc2

   END SUBROUTINE read_var2d

   !!======================================================================
END MODULE trc_util
