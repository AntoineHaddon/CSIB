MODULE trcsrc
   !!======================================================================
   !!                         ***  MODULE trcsrc  ***
   !! TOP :   TOP external inputs of nutrients
   !!======================================================================
   !! History :   0.0  !  2022-07 (O. Riche) test
   !!----------------------------------------------------------------------
   !!   trc_src        :  Read and interpolate time-varying nutrients fields
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 

   USE iom             !  I/O manager
   USE fldread         !  time interpolation

   USE prtctl_trc      !  print control for debugging
   USE in_out_manager  ! I/O manager
   USE dom_oce         ! ocean space and time domain 
   USE timing          ! Timing
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   USE sms_top         ! access index/array definitions for ext. sources
   USE trc_closeabgc   ! tmask_bgc_closea
      
   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_src_init
   PUBLIC trc_src3d
   PUBLIC trc_src2d
   PUBLIC trc_src_fedep
   PUBLIC trc_src_fesed

   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)    ::  sf_src3d   ! structure of input 3D fields (file informations, fields read)
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)    ::  sf_src2d   ! structure of input 2D fields (file informations, fields read)

   INTEGER, SAVE, PUBLIC :: nb_src3d
   INTEGER, SAVE, PUBLIC :: nb_src2d

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)  ::   src3d_dta       !: 3d source arrays
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,  :)  ::   src2d_dta       !: 2d source arrays
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)   ::  irondep_src
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)   ::  ironsed_src
   
   REAL(wp), SAVE, PUBLIC :: dustsolub0   = 0.014_wp      !: dust0 solubility      (fraction?)
   REAL(wp), SAVE, PUBLIC :: wdust0       = 2.0_wp        !: dust0 sinking speed   (m s^-1)
   REAL(wp), SAVE, PUBLIC :: sedfeinput0  = 1000._wp      !: coastal iron release (?)

   ! O Riche Aug 25th 2022
   ! Should we have this?
   ! code source of trc_srcfe
   ! contains a refernce to substituted indices
   ! !!* Substitution
!#  include "top_substitute.h90"
! O Riche Aug 25th 2022
! Only needed are fs_2/fs_jpim1
#   include "vectopt_loop_substitute.h90"

CONTAINS

   SUBROUTINE trc_src_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_src_init  ***
      !!
      !! ** Purpose :   Open files that will be queried in trc_src
	    !! 
      !! ** Method  :   fld_fill
      !!
      !!----------------------------------------------------------------------
      ! 
      INTEGER :: ios, numnml            ! error handling and namelist unit number  
      INTEGER :: ierr0, ierr1, ierr2    ! more error handling
      !
      ! Assign values to ext. source indices
      CHARACTER(len=34)  :: cltra
      INTEGER            :: jn
      !
      CHARACTER(len=100) :: cn_dir   ! Root directory for location of external source files
      ! set max number of files to a set value (<=100) but will very unlikely reach that high of a number
      TYPE(FLD_N), DIMENSION(25)        ::   sn_src3d    ! informations about the 3D external sources
      REAL(wp)   , DIMENSION(25)        ::   rn_src3d    ! scaling factors
      TYPE(FLD_N), DIMENSION(25)        ::   sn_src2d    ! informations about the 2D external sources
      REAL(wp)   , DIMENSION(25)        ::   rn_src2d    ! scaling factors
      !!
      !!----------------------------------------------------------------------
      !
      IF (lwp) THEN
        WRITE(numout,*)
        WRITE(numout,*) '   trc_src_init: initialization subroutine'
        WRITE(numout,*) '   trc_src_init: opening files'
        WRITE(numout,*) 
        CALL FLUSH(numout)
      ENDIF
      !
      NAMELIST/namtrcsrcfe/ sedfeinput0, dustsolub0, wdust0
      NAMELIST/namtrc_src3d/ cn_dir, nb_src3d, sn_src3d, rn_src3d
      NAMELIST/namtrc_src2d/ cn_dir, nb_src2d, sn_src2d, rn_src2d
      !
      ios = 0 ; ierr0 = 0  ;  ierr1 = 0  ;  ierr2 = 0 
      !
      IF( ln_timing )   CALL timing_start('trc_src_init') 
      !
      ! All the info are stored in this section of namelist_top_cfg
      CALL ctl_opn( numnml, 'namelist_top_cfg'   , 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      !
      !!!!!!!!!! Read iron dust and iron sedimentation parameters
      ! Check if the namelist has updated values of the parameters
      !
      REWIND( numnml )
      READ  ( numnml, namtrcsrcfe, IOSTAT = ios, ERR = 801)
801   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtrcsrcfe in reference namelist' )
      !       
      !!!!!!!!!! Read namelist info about external sources
      !     
      REWIND( numnml )              ! Namelist namtrc_dta in configuration namelist
      READ  ( numnml, namtrc_src3d , IOSTAT = ios, ERR = 901 )
901   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtrc_src in configuration namelist_top_cfg' )
      !
      REWIND( numnml )              ! Namelist namtrc_dta in configuration namelist
      READ  ( numnml, namtrc_src2d , IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtrc_src in configuration namelist_top_cfg' )
      !
      !
      ALLOCATE( sf_src3d(nb_src3d), sf_src2d(nb_src2d), STAT=ierr0 )
      IF( ierr0 > 0 ) THEN
        CALL ctl_stop( 'trc_src_init: unable to allocate sf_src3d and sf_src2d structures' )   ;   RETURN
      ENDIF
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_src_init: reading information on external source files:'
         WRITE(numout,*) ' number of external sources to be read:', nb_src2d+nb_src3d
         WRITE(numout,*) ' 3D fields: ', nb_src3d
         WRITE(numout,*) ' 2D fields: ', nb_src2d
      ENDIF
       !
      CALL fld_fill( sf_src3d, (/ sn_src3d /), cn_dir, 'trc_src_init', 'external 3D fields', 'namtrc_src' )
      CALL fld_fill( sf_src2d, (/ sn_src2d /), cn_dir, 'trc_src_init', 'external 2D fields', 'namtrc_src' )
      !
      ! Allocate space for each 3D-field structure associated with a file
      DO jn = 1, nb_src3d
        ALLOCATE( sf_src3d(jn)%fnow(jpi,jpj,jpk) , STAT=ierr1 )
        IF( sn_src3d(jn)%ln_tint )  ALLOCATE( sf_src3d(jn)%fdta(jpi,jpj,jpk,2) , STAT=ierr2)
        IF( ierr1 + ierr2 > 0 ) THEN
         CALL ctl_stop( 'trc_src_init: unable to allocate sf_src3d arrays' )   ;   RETURN
        ENDIF
      END DO
      
      ! Allocate space for each 2D-field structure associated with a file
      DO jn = 1, nb_src2d
        ALLOCATE( sf_src2d(jn)%fnow(jpi,jpj,1) , STAT=ierr1 )
        IF( sn_src2d(jn)%ln_tint )  ALLOCATE( sf_src2d(jn)%fdta(jpi,jpj,1,2) , STAT=ierr2)
        IF( ierr1 + ierr2 > 0 ) THEN
         CALL ctl_stop( 'trc_src_init: unable to allocate sf_src2d arrays' )   ;   RETURN
        ENDIF
      END DO
      
      ! Assign values to external source indices
      DO jn = 1, nb_src3d 
       if (lwp)  write(numout,*) sn_src3d(jn)%clvar
       cltra = TRIM( sn_src3d(jn)%clvar )
       IF( cltra == 'NO3'      )   js3d_no3 = jn      !: NO3
       IF( cltra == 'Si'       )   js3d_si  = jn      !: silicic acid
       IF( cltra == 'PO4'      )   js3d_po4 = jn      !: phosphate   
       IF( cltra == 'DOC'      )   js3d_doc = jn      !: diss. org. C
       IF( cltra == 'Fer'      )   js3d_fe  = jn      !: diss. iron  
       IF( cltra == 'epsdb'    )   js3d_hyfe= jn      !: diss. iron  
      END DO
      DO jn = 1, nb_src2d 
       if (lwp)  write(numout,*) sn_src2d(jn)%clvar
       cltra = TRIM( sn_src2d(jn)%clvar )
       IF( cltra == 'CHLA'        )   js2d_chla     = jn      !: chlorophyll-a
       IF( cltra == 'dust'        )   js2d_dust     = jn      !: dust
       IF( cltra == 'fr_par'      )   js2d_par      = jn      !: PAR
       IF( cltra == 'femask'      )   js2d_femask   = jn      !: Iron mask
       IF( cltra == 'ndep'        )   js2d_ndep     = jn      !: N deposition
       IF( cltra == 'riverdic'    )   js2d_rdic     = jn     
       IF( cltra == 'riverdoc'    )   js2d_rdoc     = jn     
       IF( cltra == 'riverpoc'    )   js2d_rpoc     = jn     
       IF( cltra == 'solubility1' )   js2d_fsol1    = jn     
       IF( cltra == 'solubility2' )   js2d_fsol2    = jn     
      END DO
      !
      ! Write info about the namelist section in output.namelist.top
      IF(lwm) WRITE ( numont, namtrcsrcfe )
      IF(lwm) WRITE ( numont, namtrc_src3d )
      IF(lwm) WRITE ( numont, namtrc_src2d )
      !
      ! These are used only to store values of the iron sources. See trc_src_fe
      ALLOCATE( irondep_src(jpi,jpj,jpk),ironsed_src(jpi,jpj,jpk), STAT=ierr0 )
      IF( ierr0 /= 0 )   CALL ctl_stop( 'STOP', 'trc_src_init: failed to allocate trc_src_fe arrays for trc_src' ) 
      !
      ! Now allocate space for the 3D and 2D ext. source arrays
      ALLOCATE( src3d_dta(jpi,jpj,jpk,nb_src3d), STAT=ierr0 )
      IF( ierr0 /= 0 )   CALL ctl_stop( 'STOP', 'trc_src_init: failed to allocate 3d arrays for trc_src' ) 
      !
      ALLOCATE( src2d_dta(jpi,jpj,nb_src2d),   STAT=ierr0 )
      IF( ierr0 /= 0 )   CALL ctl_stop( 'STOP', 'trc_src_init: failed to allocate 2d arrays for trc_src' ) 
      !
      IF( ln_timing )   CALL timing_stop('trc_src_init')
      !
     
   END SUBROUTINE trc_src_init

   !!======================================================================

   SUBROUTINE trc_src3d( kt , jn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_src  ***
      !!
      !! ** Purpose :   Read and interpolate the other nutrient IC fields
	    !!                that are not state variables but still useful for param.
      !! ** Method  :   Read the files and interpolate the appropriate variables
      !!
      !!----------------------------------------------------------------------
      ! 
      INTEGER, INTENT(in)   ::   kt          ! ocean time-step
      INTEGER, INTENT(in)   ::   jn          ! array index
      TYPE(FLD)             ::  sf_src       ! dummy variable to avoid exception to the shape matching run of passing inputs

      IF( ln_timing )   CALL timing_start('trc_src3d')
      !
      IF (lwp) THEN
        WRITE(numout,*)
        WRITE(numout,*) 'trc_src3d, reading data: ', sf_src3d(jn)%clvar
        WRITE(numout,*)
      ENDIF
      CALL fld_read( kt, 1, sf_src3d )
      !
      src3d_dta(:,:,:,jn) = sf_src3d(jn)%fnow(:,:,:)*tmask(:,:,:)
      !
      !
      IF( ln_timing )   CALL timing_stop('trc_src3d')
      !
     
   END SUBROUTINE trc_src3d
   !!======================================================================

   SUBROUTINE trc_src2d( kt, jn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_src  ***
      !!
      !! ** Purpose :   Read and interpolate the other nutrient IC fields
	    !!                that are not state variables but still useful for param.
      !! ** Method  :   Read the files and interpolate the appropriate variables
      !!
      !!----------------------------------------------------------------------
      ! 
      INTEGER, INTENT(in)   ::   kt          ! ocean time-step
      INTEGER, INTENT(in)   ::   jn          ! array index
      !
      IF( ln_timing )   CALL timing_start('trc_src2d')
      !
      IF (lwp) THEN
        WRITE(numout,*)
        WRITE(numout,*) 'trc_src2d, reading data: ', sf_src2d(jn)%clvar
        WRITE(numout,*)
      ENDIF
      CALL fld_read( kt, 1, sf_src2d )
      !
      src2d_dta(:,:,jn) = sf_src2d(jn)%fnow(:,:,1)*tmask(:,:,1)
      !
      !
      IF( ln_timing )   CALL timing_stop('trc_src2d')
      !
     
  END SUBROUTINE trc_src2d
   !!======================================================================

  SUBROUTINE trc_src_fedep( kt )
      ! compute iron sources: surface deposition from the atm. 
      !                       based on CanESM5/CanOE code.
      !
      INTEGER  :: jk                          !: loop variables
      INTEGER  :: ierr, ios, kt               !: working variables
      REAL(wp) :: ryyss                       !: number of seconds per year
      REAL(wp) :: rmtss                       !: number of seconds per month
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zirondep, zafe
      !
      IF( ln_timing )   CALL timing_start('trc_src_fedep')
      !
      IF (lwp) THEN
        WRITE(numout,*)
        WRITE(numout,*) 'trc_src: calling trc_src_fedep'
        WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*)
      ENDIF
      !      
      ALLOCATE( zirondep(jpi,jpj,jpk), zafe(jpi,jpj,jpk), STAT=ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'trc_src_fedep: failed to allocate 3d arrays for trc_src_fedep' )
      !
      ryyss    = REAL(nyear_len(1), wp ) * rday ! nyear_len is integer.
      rmtss    = ryyss / raamo
      ! 
      IF(lwp) WRITE(numout,*) 'computation of iron aeolian depostion'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)
      IF(lwp) CALL FLUSH(numout) 
      !
      !
      IF(lwp) THEN
         WRITE(numout,*) 'Namelist : namtrcsrcfe '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   coastal iron release sedfeinput0 = ', sedfeinput0
         WRITE(numout,*) '   dust solubility      dustsolub0  = ', dustsolub0
         WRITE(numout,*) '   dust sinking speed   wdust0      = ', wdust0
         WRITE(numout,*)
         CALL FLUSH(numout)
      ENDIF
      
      ! iron aeolian depostion
      CALL trc_src2d( kt, js2d_dust )
      !
      zirondep(:,:,:) = 0.e0          ! Initialisation of variables USEd to compute deposition

      ! Iron deposition at the surface
      ! -------------------------------------
      ! dust0 is in kgFe m^-2 month^-1; zirondep is in nmolFe m^-3 s^-1
      zirondep(:,:,1) = dustsolub0 * src2d_dta(:,:,js2d_dust) / ( 55.85 * rmtss * e3t_n(:,:,1) ) * 1.E+12

      ! Iron solubilization of particles in the water column
      ! ----------------------------------------------------
      DO jk = 2, jpkm1
         zirondep(:,:,jk) = src2d_dta(:,:,js2d_dust) / ( wdust0 * 55.85 * rmtss ) * 1.e-4 * EXP( -gdept_n(:,:,jk) / 1000. ) * 1.E+12
      END DO

      ! Diagnostics
      IF( lk_iomput ) THEN
        zafe(:,:,:) = zirondep(:,:,:) * 1.E-9 * tmask_bgc_closea(:,:,:)      ! zirondep and ironsed are in nmol m^-3 s^-1
        CALL iom_put( "Irondep", zafe  )  ! surface downward net flux of iron
      ENDIF  
      irondep_src(:,:,:) = zirondep(:,:,:)  

      IF( ln_timing )   CALL timing_stop('trc_src_fedep')
  
  END SUBROUTINE trc_src_fedep

  SUBROUTINE trc_src_fesed
      ! compute iron sources: bottom flux from sediments
      !                       based on CanESM5/CanOE code.
      INTEGER  :: ji, jj, jk                  !: loop variables
      INTEGER  :: ierr, inum, ios             !: working variables
      !
      REAL(wp) :: zexpide, zdenitide, zmaskt  ! iron sed working variables
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zironsed, zbfe, zcmask
      ! 
      IF( ln_timing )   CALL timing_start('trc_src_fesed')
      !
      IF (lwp) THEN
        WRITE(numout,*)
        WRITE(numout,*) 'trc_src: calling trc_src_fesed'
        WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*)
      ENDIF
      !      
      ALLOCATE( zironsed(jpi,jpj,jpk), zbfe(jpi,jpj,jpk), zcmask(jpi,jpj,jpk), STAT=ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'trc_src_fesed: failed to allocate 3d arrays for trc_src_fesed' )
      !

      ! Iron coastal flux
      ! -------------------------
      ! coastal and island masks
      ! ------------------------
      IF(lwp) WRITE(numout,*) 'computation of an island mask to enhance coastal supply of iron'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)
      IF(lwp) CALL FLUSH(numout)
      !
      IF(lwp) THEN
         WRITE(numout,*) 'Namelist : namtrcsrcfe '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   coastal iron release sedfeinput0 = ', sedfeinput0
         WRITE(numout,*) '   dust solubility      dustsolub0  = ', dustsolub0
         WRITE(numout,*) '   dust sinking speed   wdust0      = ', wdust0
         WRITE(numout,*)
         CALL FLUSH(numout)
      ENDIF
      !     
      CALL iom_open('bathy.orca.nc', inum)
      CALL iom_get(inum, jpdom_data,'bathy',zcmask(:,:,:), lrowattr=ln_use_jattr)
      CALL iom_close(inum)
      !
      DO jk = 1, 5
        DO jj = 2, jpjm1
           DO ji = fs_2, fs_jpim1     ! These if required are added with the include statement just above the CONTAINS statement
              IF( tmask_bgc_closea(ji,jj,jk) /= 0. ) THEN
                 zmaskt = tmask_bgc_closea(ji+1,jj,jk) * tmask_bgc_closea(ji-1,jj,jk) & 
                    &   * tmask_bgc_closea(ji,jj+1,jk) * tmask_bgc_closea(ji,jj-1,jk) &
                    &   * tmask_bgc_closea(ji,jj,jk+1)
                 IF( zmaskt == 0. )  zcmask(ji,jj,jk ) = MAX( 0.1, zcmask(ji,jj,jk) ) 
              END IF
           END DO
        END DO
      END DO
      CALL lbc_lnk('trc_src_fesed', zcmask(:,:,:) , 'T', 1. )      ! lateral boundary conditions on cmask   (sign unchanged)
      DO jk = 1, jpk
        DO jj = 1, jpj
           DO ji = 1, jpi
              zexpide   = MIN( 8.,( gdept_n(ji,jj,jk) / 500. )**(-1.5) )
              zdenitide = -0.9543 + 0.7662 * LOG( zexpide ) - 0.235 * LOG( zexpide )**2
              zcmask(ji,jj,jk) = zcmask(ji,jj,jk) * MIN( 1., EXP( zdenitide ) / 0.5 )
           END DO
        END DO
      END DO
      
      ! Coastal supply of iron
      ! -------------------------
      zironsed(:,:,jpk) = 0._wp
      DO jk = 1, jpkm1
        zironsed(:,:,jk) = sedfeinput0 * zcmask(:,:,jk) / ( e3t_n(:,:,jk) * rday )
      END DO

      ! Diagnostics
      IF( lk_iomput ) THEN
        zbfe(:,:,:) = zironsed(:,:,:) * 1.E-9 * tmask_bgc_closea(:,:,:)
        CALL iom_put( "Ironsed", zbfe  )  ! iron from sediments        
      ENDIF  
      ironsed_src(:,:,:) = zironsed(:,:,:)
      !
      DEALLOCATE( zcmask )

      IF( ln_timing )   CALL timing_stop('trc_src_fesed')
      
  END SUBROUTINE trc_src_fesed
 

END MODULE trcsrc        
