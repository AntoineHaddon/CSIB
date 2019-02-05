

MODULE trcnam_cfc
   !!======================================================================
   !! *** MODULE trcnam_cfc ***
   !! TOP : initialisation of some run parameters for CFC chemical model
   !!======================================================================
   !! History : 2.0 ! 2007-12 (C. Ethe, G. Madec) from trcnam.cfc.h90
   !!----------------------------------------------------------------------
#if defined key_cfc
   !!----------------------------------------------------------------------
   !! 'key_cfc' CFC tracers
   !!----------------------------------------------------------------------
   !! trc_nam_cfc : CFC model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc ! Ocean variables
   USE par_trc ! TOP parameters
   USE trc ! TOP variables
   USE trcsms_cfc ! CFC specific variable
   USE iom ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_nam_cfc ! called by trcnam.F90 module
   CHARACTER(LEN=255), PUBLIC :: cfc_nc_file = '' ! Name of the netcdf file containing atmospheric history of CFCs

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_cfc.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_cfc
      !!-------------------------------------------------------------------
      !! *** ROUTINE trc_nam_cfc ***
      !!
      !! ** Purpose : Definition some run parameter for CFC model
      !!
      !! ** Method : Read the namcfc namelist and check the parameter
      !! values called at the first timestep (nittrc000)
      !!
      !! ** input : Namelist namcfc
      !!----------------------------------------------------------------------
      INTEGER :: numnatc
      INTEGER :: jl, jn
      TYPE(DIAG), DIMENSION(jp_cfc_2d) :: cfcdia2d
      TYPE(DIAG), DIMENSION(jp_cfc_3d) :: cfcdia3d
      !!
      NAMELIST/namcfcparam/ offset_cfc_year, cfc_nc_file, ln_reset_cfc, nn_reset_cfc
      NAMELIST/namcfcdia/ cfcdia2d ! additional diagnostics
      !!-------------------------------------------------------------------

      ! If this is blank, then the atmospheric values are expected to be
     ! to be read in from a formatted text file
      cfc_nc_file = ''

      ! ! Open namelist file
      CALL ctl_opn( numnatc, 'namelist_cfc', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )

      READ( numnatc , namcfcparam ) ! read namelist

      IF(lwp) THEN ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' trc_nam: Read namdates, namelist for CFC chemical model'
         WRITE(numout,*) ' ~~~~~~~'
         WRITE(numout,*) '    offset from the model year                offset_cfc_year = ' , offset_cfc_year
         WRITE(numout,*) '    netcdf file with surface values           cfc_nc_file     = ' , cfc_nc_file
         WRITE(numout,*) '    Set CFCs to 0 until a specified year      ln_reset_cfc    = ' , ln_reset_cfc 
         WRITE(numout,*) '    Year until which CFCs will be set to 0    nn_cfc_reset    = ' , nn_reset_cfc
      ENDIF

      IF( .NOT.lk_iomput .AND. ln_diatrc ) THEN
         !
         ! Namelist namcfcdia
         ! -------------------
         ! Surface fluxes
         DO jl = 1, jp_cfc_2d
            WRITE(cfcdia2d(jl)%sname,'("2D_",I1)') jl ! short name
            WRITE(cfcdia2d(jl)%lname,'("2D DIAGNOSTIC NUMBER ",I2)') jl ! long name
            cfcdia2d(jl)%units = ' ' ! units
         END DO

         REWIND( numnatc ) ! read natrtd
         READ ( numnatc, namcfcdia )

         DO jl = 1, jp_cfc_2d
            jn = jp_cfc0_2d + jl - 1
            ctrc2d(jn) = TRIM( cfcdia2d(jl)%sname )
            ctrc2l(jn) = TRIM( cfcdia2d(jl)%lname )
            ctrc2u(jn) = TRIM( cfcdia2d(jl)%units )
         END DO

         IF(lwp) THEN ! control print
            WRITE(numout,*)
            WRITE(numout,*) ' Namelist : namcfcdia'
            DO jl = 1, jp_cfc_2d
               jn = jp_cfc0_2d + jl - 1
               WRITE(numout,*) '  2d diag nb : ', jn, '    short name : ', ctrc2d(jn), &
                 & '  long name  : ', ctrc2l(jn), '   unit : ', ctrc2u(jn)
            END DO
            WRITE(numout,*) ' '
         ENDIF
         !
      ENDIF

   END SUBROUTINE trc_nam_cfc
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                                No CFC
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_cfc                      ! Empty routine
   END  SUBROUTINE  trc_nam_cfc
#endif

   !!======================================================================
END MODULE trcnam_cfc
