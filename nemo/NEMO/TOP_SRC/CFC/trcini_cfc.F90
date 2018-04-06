

MODULE trcini_cfc
   !!======================================================================
   !! *** MODULE trcini_cfc ***
   !! TOP : initialisation of the CFC tracers
   !!======================================================================
   !! History : 2.0 ! 2007-12 (C. Ethe, G. Madec)
   !!----------------------------------------------------------------------
#if defined key_cfc
   !!----------------------------------------------------------------------
   !! 'key_cfc' CFC tracers
   !!----------------------------------------------------------------------
   !! trc_ini_cfc : CFC model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc ! Ocean variables
   USE par_trc ! TOP parameters
   USE trc ! TOP variables
   USE trcsms_cfc ! CFC sms trends
   USE trcnam_cfc, ONLY : cfc_nc_file
   USE obs_utils, ONLY : chkerr
   USE netcdf
   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_ini_cfc ! called by trcini.F90 module

   CHARACTER (len=34) :: clname = 'cfc1112.atm' ! ???

   INTEGER :: inum ! unit number
   REAL(wp) :: ylats = -10. ! 10 degrees south
   REAL(wp) :: ylatn = 10. ! 10 degrees north

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcini_cfc.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_cfc
      !!----------------------------------------------------------------------
      !! *** trc_ini_cfc ***
      !!
      !! ** Purpose : initialization for cfc model
      !!
      !! ** Method : - Initializes cfc values to 0 (cold start) or to zero
      !! -
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj, jn, jl, jm, js, io, ierr
      REAL :: zyd
      LOGICAL :: read_nc ! If true, then read atmospheric values from a netcdf file
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_cfc: initialisation of CFC chemical model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~'
      IF( trc_sms_cfc_alloc() /= 0 ) CALL ctl_stop( 'STOP', 'trc_ini_cfc: unable to allocate CFC arrays' )

      ! If cfc_nc_file is not blank then read from a netcdf file
      read_nc = ( TRIM(cfc_nc_file) /= '' )

      ! Initialization of boundaries conditions
      ! ---------------------------------------
      xphem (:,:) = 0._wp
      p_cfc(:,:,:) = 0._wp

      ! Initialization of qint in case of no restart
      !----------------------------------------------
      qtr_cfc(:,:,:) = 0._wp
      IF( .NOT. ln_rsttr ) THEN
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'Initialization de qint ; No restart : qint equal zero '
         ENDIF
         qint_cfc(:,:,:) = 0._wp
         DO jl = 1, jp_cfc
            jn = jp_cfc0 + jl - 1
            trn(:,:,:,jn) = 0._wp
         END DO
      ENDIF

      ! Read in atmospheric surface values either from 1) formatted text file 2) netcdf file
      IF ( READ_NC ) THEN
        CALL read_from_netcdf( )
      ELSE
        CALL read_from_formatted_input( )
      ENDIF

      ! Interpolation factor of atmospheric partial pressure
      ! Linear interpolation between 2 hemispheric function of latitude between ylats and ylatn
      !---------------------------------------------------------------------------------------
      zyd = ylatn - ylats
      DO jj = 1 , jpj
         DO ji = 1 , jpi
            IF( gphit(ji,jj) >= ylatn ) THEN ; xphem(ji,jj) = 1.e0
            ELSEIF( gphit(ji,jj) <= ylats ) THEN ; xphem(ji,jj) = 0.e0
            ELSE ; xphem(ji,jj) = ( gphit(ji,jj) - ylats) / zyd
            ENDIF
         END DO
      END DO
      !
      IF(lwp) WRITE(numout,*) 'Initialization of CFC tracers done'
      IF(lwp) WRITE(numout,*) ' '
      !
   END SUBROUTINE trc_ini_cfc

   !> This is the original NEMO 3.4 behavior which opens and reads in CFC-11 and CFC-12 values from
   !! a formatted text file whose provenance is currently unknown?
   !! This routine should be substituted for something more general so that we do not rely on explicit
   !! array indexing to keep track of the time. As of now, the atmospheric history of CFCs/SF6 has the following format
   !! 6 lines in the header
   !! Atmospheric values starting in 1931 formatted with the following columns:
   !! YEAR CFC-11[North] CFC-12[North] SF6[North] CFC-11[South] CFC-12[South] SF6[South]
   SUBROUTINE read_from_formatted_input( )
      INTEGER :: ji, jj, jn, jl, jm, js, io, ierr
      INTEGER :: iskip = 6 ! number of 1st descriptor lines
      REAL(wp) :: zyy, zyd

      IF(lwp) WRITE(numout,*) 'read of formatted file cfc1112atm'

      CALL ctl_opn( inum, clname, 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      REWIND(inum)

      ! compute the number of year in the file
      ! file starts in 1931 do jn represent the year in the century
      jn = 31
      DO
        READ(inum,'(1x)',END=100)
        jn = jn + 1
      END DO
 100 jpyear = jn - 1 - iskip
      IF ( lwp) WRITE(numout,*) '    ', jpyear ,' years read'
      ! ! Allocate CFC arrays

      ALLOCATE( p_cfc(jpyear,jphem,jp_cfc), STAT=ierr )
      IF( ierr > 0 ) THEN
         CALL ctl_stop( 'trc_ini_cfc: unable to allocate p_cfc array' ) ; RETURN
      ENDIF

      REWIND(inum)

      DO jm = 1, iskip ! Skip over 1st six descriptor lines
         READ(inum,'(1x)')
      END DO
      ! file starts in 1931 do jn represent the year in the century.jhh
      ! Read file till the end
      jn = 31
      DO
        ! File is assumed to have 7 columns: Year, CFC-11 North, CFC-12 North, SF6 North
        ! CFC-11 South, CFC-12 South, SF6 South
        READ(inum,*, IOSTAT=io) zyy, p_cfc(jn,1,1), p_cfc(jn,1,2), p_cfc(jn,1,3), &
                                     p_cfc(jn,2,1), p_cfc(jn,2,2), p_cfc(jn,2,3)
        IF( io < 0 ) exit
        jn = jn + 1
      END DO

      p_cfc(32,1:2,1) = 5.e-4 ! modify the values of the first years
      p_cfc(33,1:2,1) = 8.e-4
      p_cfc(34,1:2,1) = 1.e-6
      p_cfc(35,1:2,1) = 2.e-3
      p_cfc(36,1:2,1) = 4.e-3
      p_cfc(37,1:2,1) = 6.e-3
      p_cfc(38,1:2,1) = 8.e-3
      p_cfc(39,1:2,1) = 1.e-2

      IF(lwp) THEN ! Control print
         WRITE(numout,*)
         WRITE(numout,*) ' Year   p11HN    p11HS    p12HN    p12HS   sf6HN   sf6HS'
         DO jn = 30, jpyear
            WRITE(numout, '( 1I4, 4F9.2)') jn, p_cfc(jn,1,1), p_cfc(jn,2,1), &
                                               p_cfc(jn,1,2), p_cfc(jn,2,2), &
                                               p_cfc(jn,1,3), p_cfc(jn,2,3)
         END DO
      ENDIF

   END SUBROUTINE read_from_formatted_input

   !> Reads the atmospheric history from a netcdf file
   SUBROUTINE read_from_netcdf( )
      INTEGER :: ncid, varid, dimlen, dimid
      CHARACTER(LEN=255) :: dimname

      ! Open netcdf and get the dimension
      CALL chkerr(nf90_open( cfc_nc_file, NF90_NOWRITE, ncid ), 'trcini_cfc', 0)
      CALL chkerr(nf90_inq_dimid(ncid, "index", dimid), 'trcini_cfc', 0)
      CALL chkerr(nf90_inquire_dimension(ncid, dimid, dimname, len = dimlen), 'trcini_cfc', 0)

      ! Allocate arrays now that we know how many years are in the file
      ALLOCATE(p_cfc_year(dimlen)) ; p_cfc_year(:) = 0.
      ALLOCATE(p_cfc(dimlen, jphem, jp_cfc)) ; p_cfc(:,:,:) = 0.

      ! Read all the necessary fields
      CALL read_var1d( ncid, "Year", p_cfc_year(:) )
      CALL read_var1d( ncid, "CFC11NH", p_cfc(:,1,1) )
      CALL read_var1d( ncid, "CFC11SH", p_cfc(:,2,1) )
      CALL read_var1d( ncid, "CFC12NH", p_cfc(:,1,2) )
      CALL read_var1d( ncid, "CFC12SH", p_cfc(:,2,2) )
      CALL read_var1d( ncid, "SF6NH", p_cfc(:,1,3) )
      CALL read_var1d( ncid, "SF6SH", p_cfc(:,2,3) )

      CALL chkerr(nf90_close( ncid ), 'trcini_cfc',0)

   END SUBROUTINE read_from_netcdf

   !> Read a vector variable given its name
   SUBROUTINE read_var1d( ncid, varname, varout )
      INTEGER , INTENT(IN ) :: ncid !< File ID for an already opened netcdf file
      CHARACTER(LEN=*) , INTENT(IN ) :: varname !< Name of variable to be read
      REAL(wp), DIMENSION(:), INTENT(INOUT) :: varout !< Variable data

      INTEGER :: varid
      CALL chkerr( nf90_inq_varid(ncid, varname, varid), 'trcini_cfc', 0 )
      CALL chkerr( nf90_get_var(ncid, varid, varout), 'trcini_cfc', 0 )

   END SUBROUTINE read_var1d

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                         No CFC tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ini_cfc             ! Empty routine
   END SUBROUTINE trc_ini_cfc
#endif

END MODULE trcini_cfc
