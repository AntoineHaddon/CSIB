MODULE sbcfaf
   !!======================================================================
   !!                       ***  MODULE  sbcanom  ***
   !! Ocean surface flux anomalies
   !!=====================================================================
   !! History :  3.4  !  2019-01 Andrew shao
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   namanom   : flux formulation namlist
   !!   sbc_anom  : flux formulation as ocean surface boundary condition (forced mode, fluxes read in NetCDF files)
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean fields
   USE phycst          ! physical constants
   USE fldread         ! read input fields
   USE iom             ! IOM library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_fafmip            ! routine called by step.F90
   PUBLIC sbc_fafmip_init       ! routine called by step.
   LOGICAL, PUBLIC :: ln_faftau  = .false.  ! If true, read in and apply windstress anomalies
   LOGICAL, PUBLIC :: ln_fafemp  = .false.  ! If true, read in and apply freshwater anomalies
   LOGICAL, PUBLIC :: ln_fafhflx = .false.  ! If true, read in heat flux anomalies which MAY be used later
   LOGICAL, PUBLIC :: ln_fafheat = .false.  ! If true, use the FAFMIP temperature variable for use in flux calculations

   INTEGER , PARAMETER ::   jpfld   = 4   ! maximum number of files to read 
   INTEGER , PUBLIC,  PARAMETER ::   jp_utau = 1   ! index of wind stress (i-component) file
   INTEGER , PUBLIC,  PARAMETER ::   jp_vtau = 2   ! index of wind stress (j-component) file
   INTEGER , PUBLIC,  PARAMETER ::   jp_emp  = 3   ! index of evaporation-precipation file
   INTEGER , PUBLIC,  PARAMETER ::   jp_hflx = 4   ! index of heat flux

   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:)       :: sf_fafmip    ! structure of input fields (file informations, fields read)
   REAL,      PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)   :: Tr_sbc       ! Surface flux for redistributed heat tracer


   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO-consortium (2010) 
   !! $Id: sbcflx.F90 2715 2011-03-30 15:58:35Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_fafmip( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_flx  ***
      !!                   
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!                (momentum, heat, freshwater and runoff) 
      !!
      !! ** Method  : - READ each fluxes in NetCDF files:
      !!                   i-component of the stress              utau  (N/m2)
      !!                   j-component of the stress              vtau  (N/m2)
      !!                   net downward heat flux                 qtot  (watt/m2)
      !!                   net downward radiative flux            qsr   (watt/m2)
      !!                   net upward freshwater (evapo - precip) emp   (kg/m2/s)
      !!
      !!      CAUTION :  - never mask the surface stress fields
      !!                 - the stress is assumed to be in the mesh referential
      !!                   i.e. the (i,j) referential
      !!
      !! ** Action  :   update at each time-step
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        10m wind module at T-point
      !!              - qns, qsr    non-slor and solar heat flux
      !!              - emp, emps   evaporation minus precipitation
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj, jf            ! dummy indices
      INTEGER  ::   ierror                ! return error code
      REAL(wp) ::   zfact                 ! temporary scalar
      REAL(wp) ::   zrhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3       ! drag coefficient
      REAL(wp) ::   ztx, zty, zmod, zcoef ! temporary variables
      !!
      !!---------------------------------------------------------------------
      !
      ! Check to see if anomalies should be added by each flux type
      IF ( ln_faftau ) THEN
         CALL fld_read( kt, nn_fsbc, sf_fafmip(jp_utau:jp_utau) )
         CALL fld_read( kt, nn_fsbc, sf_fafmip(jp_vtau:jp_vtau) )
      ENDIF
      IF ( ln_fafemp ) THEN
         CALL fld_read( kt, nn_fsbc, sf_fafmip(jp_emp:jp_emp) )
      ENDIF
      IF ( ln_fafhflx ) THEN
         ! Even though this is read in here, this is not applied until traqsr.
         CALL fld_read( kt, nn_fsbc, sf_fafmip(jp_hflx:jp_hflx) )
      ENDIF
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN                        ! update ocean fluxes at each SBC frequency
         ! Treat the perturbation to wind stress first
         IF ( ln_faftau ) THEN
!CDIR COLLAPSE
            DO jj = 1, jpj                                           ! set the ocean fluxes from read fields
               DO ji = 1, jpi
                  utau(ji,jj) = utau(ji,jj) + sf_fafmip(jp_utau)%fnow(ji,jj,1)
                  vtau(ji,jj) = vtau(ji,jj) + sf_fafmip(jp_vtau)%fnow(ji,jj,1)
               END DO
            END DO
            !                                                        ! module of wind stress and wind speed at T-point
            zcoef = 1. / ( zrhoa * zcdrag )
!CDIR NOVERRCHK
            DO jj = 2, jpjm1
!CDIR NOVERRCHK
               DO ji = fs_2, fs_jpim1   ! vect. opt.
                  ztx = utau(ji-1,jj  ) + utau(ji,jj) 
                  zty = vtau(ji  ,jj-1) + vtau(ji,jj) 
                  zmod = 0.5 * SQRT( ztx * ztx + zty * zty )
                  taum(ji,jj) = zmod
                  wndm(ji,jj) = SQRT( zmod * zcoef )
               END DO
            END DO
            CALL lbc_lnk( taum(:,:), 'T', 1. )   ;   CALL lbc_lnk( wndm(:,:), 'T', 1. )
         ENDIF
        ! Now apply freshwater perturbations 
         IF ( ln_fafemp ) THEN
            DO jj = 1, jpj                                           ! set the ocean fluxes from read fields
               DO ji = 1, jpi
                  ! FAFMIP provides Water flux into the ocean
                  emp (ji,jj) = emp(ji,jj) + -1.*sf_fafmip(jp_emp)%fnow(ji,jj,1)
               END DO
            END DO
            emps(:,:) = emp (:,:)                                    ! Initialization of emps (needed when no ice model)

            CALL iom_put("fafemp_x_sss",sf_fafmip(jp_emp)%fnow(:,:,1)*tsn(:,:,1,jp_sal))
         ENDIF

      ENDIF
   END SUBROUTINE sbc_fafmip

   SUBROUTINE sbc_fafmip_init( )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_init ***
      !!
      !! ** Purpose :   Initialisation of the FAFMIP anomaly configurations
      !!
      !! ** Method  :   Read the namsbc_fafmip namelist and set derived parameters
      !!
      !! ** Action  : - read namsbc parameters
      !!----------------------------------------------------------------------
      CHARACTER(len=100) ::  cn_dir                               ! Root directory for location of flx files
      TYPE(FLD_N) ::   sn_utau, sn_vtau, sn_qtot, sn_emp          ! information about the fields to be read
      NAMELIST/namsbc_fafmip/ ln_faftau, ln_fafhflx, ln_fafemp, ln_fafheat, cn_dir, sn_utau, sn_vtau, sn_qtot, sn_emp 
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i                    ! array of namelist information structures
      INTEGER :: ji, ierror

      ! set file information
      cn_dir = './'        ! directory in which the model is executed
      ! ... default values (NB: frequency positive => hours, negative => months)
      !                   !  file   ! frequency !  variable  ! time intep !  clim   ! 'yearly' or ! weights  ! rotation  !
      !                   !  name   !  (hours)  !   name     !   (T/F)    !  (T/F)  !  'monthly'  ! filename ! pairs     !
      sn_utau = FLD_N(  'utau' ,    -1     ,  'utau'    ,  .true.    , .true.  ,   'yearly'  , ''       , ''        )
      sn_vtau = FLD_N(  'vtau' ,    -1     ,  'vtau'    ,  .true.    , .true.  ,   'yearly'  , ''       , ''        )
      sn_qtot = FLD_N(  'qtot' ,    -1     ,  'qtot'    ,  .true.    , .true.  ,   'yearly'  , ''       , ''        )
      sn_emp  = FLD_N(  'emp'  ,    -1     ,  'emp'     ,  .true.    , .true.  ,   'yearly'  , ''       , ''        )
      !
      REWIND ( numnam )                         ! read in namlist namflx
      READ   ( numnam, namsbc_fafmip ) 
      
      IF(lwp) THEN               ! Control print
         WRITE(numout,*) '        Namelist namsbc_fafmip                                    '
         WRITE(numout,*) '          Read in and apply wind stress anomalies                 ln_faftau  = ', ln_faftau
         WRITE(numout,*) '          Read in and apply water flux anomalies                  ln_fafemp  = ', ln_fafemp
         WRITE(numout,*) '          Read in heat flux anomalies which MAY be used later     ln_fafhflx = ', ln_fafhflx
         WRITE(numout,*) '          Read in and apply heat flux anomalies                   ln_fafheat = ', ln_fafheat
      ENDIF
      !
      !                                         ! store namelist information in an array
      slf_i(jp_utau) = sn_utau   ;   slf_i(jp_vtau) = sn_vtau
      slf_i(jp_hflx) = sn_qtot
      slf_i(jp_emp ) = sn_emp
      !
      ALLOCATE( sf_fafmip(jpfld), STAT=ierror )        ! set sf structure
      IF( ierror > 0 ) THEN   
         CALL ctl_stop( 'sbc_flx: unable to allocate sf structure' )   ;   RETURN  
      ENDIF
      DO ji= 1, jpfld
         ALLOCATE( sf_fafmip(ji)%fnow(jpi,jpj,1) )
         IF( slf_i(ji)%ln_tint ) ALLOCATE( sf_fafmip(ji)%fdta(jpi,jpj,1,2) )
      END DO
      !                                         ! fill sf with slf_i and control print
      CALL fld_fill( sf_fafmip, slf_i, cn_dir, 'sbc_anom', 'flux anomalies for ocean surface boundary condition', 'namsbc_fafmip' )
      
      ! Allocate the "redistributed heat" tracer flux array
      ALLOCATE( Tr_sbc(jpi,jpj,jpk ) )
      Tr_sbc(:,:,:) = 0.

   END SUBROUTINE sbc_fafmip_init
   !!======================================================================
END MODULE sbcfaf
