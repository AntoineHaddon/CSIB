MODULE trcsrc_canbgc
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

   USE prtctl      !  print control for debugging
   USE in_out_manager  ! I/O manager
   USE dom_oce         ! ocean space and time domain 
   USE timing          ! Timing
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   USE sms_cmoc          ! shared variables
   USE sms_top_canbgc    ! access index/array definitions for ext. sources
   USE trc_closea_canbgc ! tmask_bgc_closea
   
   IMPLICIT NONE
   PRIVATE

   ! General external source subroutines
   PUBLIC trc_src_init
   PUBLIC trc_src3d
   PUBLIC trc_src2d
   ! CMOC specific source subroutines
   PUBLIC trc_src_fedep
   PUBLIC trc_src_fesed
   PUBLIC trc_src_criver
   PUBLIC trc_bott_cmoc
   PUBLIC trc_n2fx_denit_cmoc
   PUBLIC trc_n2fx_init_cmoc

   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)    ::  sf_src3d   ! structure of input 3D fields (file informations, fields read)
   TYPE(FLD), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:)    ::  sf_src2d   ! structure of input 2D fields (file informations, fields read)

   INTEGER, SAVE, PUBLIC :: nb_src3d
   INTEGER, SAVE, PUBLIC :: nb_src2d

   ! External CMOC iron sources
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)   ::  irondep_cmoc
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)   ::  ironsed_cmoc

   ! External CMOC river DIC/DOC sources
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:  )   ::  cotdep_cmoc
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:  )   ::  rivinp_cmoc   

   ! specific CMOC RHS terms derived from cotdep_cmoc and rivinp_cmoc
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:  )   ::  no3river_cmoc
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:  )   ::  dicriver_cmoc   
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:  )   ::  talriver_cmoc
   ! specific CMOC RHS terms needed for bottom instantaneous remineralization
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:  )   ::  dicbott_cmoc
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:  )   ::  talbott_cmoc
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:  )   ::  no3bott_cmoc
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:  )   ::  oxybott_cmoc
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:  )   ::  pocbott_cmoc

   ! Iron dust parameters
   REAL(wp), SAVE, PUBLIC :: dustsolub0   = 0.014_wp      !: dust0 solubility      (fraction?)
   REAL(wp), SAVE, PUBLIC :: wdust0       = 2.0_wp        !: dust0 sinking speed   (m s^-1)
   REAL(wp), SAVE, PUBLIC :: sedfeinput0  = 1000._wp      !: coastal iron release (?)

   ! External source switches
   LOGICAL, SAVE, PUBLIC  :: ln_dust0  = .false. 
   LOGICAL, SAVE, PUBLIC  :: ln_river0 = .false. 
   LOGICAL, SAVE, PUBLIC  :: ln_ndepo0 = .false.  
   !
   ! Conversion coefficients
   REAL(wp), SAVE, PUBLIC :: ryyssb   !: number of seconds per year
   REAL(wp), SAVE, PUBLIC :: rmtssb  !: number of seconds per month
   !
   ! specific CMOC RHS terms needed for N2 fixation and
   ! denitrification
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: n2fix_cmoc
   REAL(wp), SAVE, PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: denit_cmoc
   !
   ! O Riche Aug 25th 2022
   ! Should we have this?
   ! code source of trc_srcfe
   ! contains a refernce to substituted indices
   ! !!* Substitution
!#  include "top_substitute.h90"
! O Riche Aug 25th 2022
! Only needed are fs_2/fs_jpim1
!#   include "vectopt_loop_substitute.h90"
#  include "domzgr_substitute.h90"

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
      TYPE(FLD_N), DIMENSION(25)  ::   sn_src3d    ! informations about the 3D external sources
      REAL(wp)   , DIMENSION(25)  ::   rn_src3d    ! scaling factors
      TYPE(FLD_N), DIMENSION(25)  ::   sn_src2d    ! informations about the 2D external sources
      REAL(wp)   , DIMENSION(25)  ::   rn_src2d    ! scaling factors
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
      ! Conversion factors used later in trc_src_fedep and trc_src_criver
      ryyssb = REAL(nyear_len(1), wp ) * rday ! nyear_len is integer.
      rmtssb = ryyssb / raamo                 ! raamo, number of months in a year, from USE phycst (called in oce_trc) 
      !
      ! Read namelists
      NAMELIST/namtrcsrcfe/ sedfeinput0, dustsolub0, wdust0
      NAMELIST/namtrcsrclog/ ln_dust0, ln_river0, ln_ndepo0
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
      !       !
      REWIND( numnml )
      READ  ( numnml, namtrcsrclog, IOSTAT = ios, ERR = 802)
802   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtrcsrclog in reference namelist' )
      !       !
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
      ALLOCATE( irondep_cmoc(jpi,jpj,jpk),ironsed_cmoc(jpi,jpj,jpk), STAT=ierr0 )
      IF( ierr0 /= 0 )   CALL ctl_stop( 'STOP', 'trc_src_init: failed to allocate trc_src_fe arrays for trc_src' ) 
      !
      ! These are used only to store values of the rivers sources. See trc_src_criver
      ALLOCATE( cotdep_cmoc(jpi,jpj),rivinp_cmoc(jpi,jpj), STAT=ierr0 )
      IF( ierr0 /= 0 )   CALL ctl_stop( 'STOP', 'trc_src_init: failed to allocate trc_src_criver arrays for trc_src' ) 
      ! These are used only to store values of N2 fixation and denitrification. See trc_n2fx_denit_cmoc
      ALLOCATE( n2fix_cmoc(jpi,jpj,jpk),denit_cmoc(jpi,jpj,jpk), STAT=ierr0 )
      IF( ierr0 /= 0 )   CALL ctl_stop( 'STOP', 'trc_src_init: failed to allocate trc_n2fx_denit_cmoc arrays for trc_src' ) 
      !
      ALLOCATE( no3river_cmoc(jpi,jpj), dicriver_cmoc(jpi,jpj), &
              & talriver_cmoc(jpi,jpj),          STAT=ierr0 )
      IF( ierr0 /= 0 )   CALL ctl_stop( 'STOP', 'trc_src_init: failed to allocate CMOC RHS trc_src_criver arrays for trc_src' ) 
      ALLOCATE( no3bott_cmoc(jpi,jpj), dicbott_cmoc(jpi,jpj), &
              & talbott_cmoc(jpi,jpj), oxybott_cmoc(jpi,jpj), &
              & pocbott_cmoc(jpi,jpj),                        STAT=ierr0 )
      IF( ierr0 /= 0 )   CALL ctl_stop( 'STOP', 'trc_src_init: failed to allocate CMOC RHS trc_bott_cmoc arrays for trc_src' ) 
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

  SUBROUTINE trc_src_fedep( kt, Kbb, Kmm, Krhs )
      ! compute iron sources: surface deposition from the atm. 
      !                       based on CanESM5/CanOE code.
      !
      INTEGER, INTENT(in) :: kt, Kbb, Kmm, Krhs 
      !
      INTEGER  :: jk                          !: loop variables
      INTEGER  :: ierr, ios                   !: working variables
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
      
      ! iron aeolian deposition
      CALL trc_src2d( kt, js2d_dust )
      !
      zirondep(:,:,:) = 0.e0          ! Initialisation of variables USEd to compute deposition

      ! Iron deposition at the surface
      ! -------------------------------------
      ! dust0 is in kgFe m^-2 month^-1; zirondep is in nmolFe m^-3 s^-1
      zirondep(:,:,1) = dustsolub0 * src2d_dta(:,:,js2d_dust) / ( 55.85 * rmtssb * e3t(:,:,1,Kmm) ) * 1.E+12

      ! Iron solubilization of particles in the water column
      ! ----------------------------------------------------
      DO jk = 2, jpkm1
         zirondep(:,:,jk) = src2d_dta(:,:,js2d_dust) / ( wdust0 * 55.85 * rmtssb ) * 1.e-4 * EXP( -gdept(:,:,jk,Kmm) / 1000. ) * 1.E+12
      END DO

      DO jk = 1, jpkm1
         tr(:,:,jk,jrfer, Krhs) = tr(:,:,jk,jrfer, Krhs) + zirondep(:,:,jk) * qfact2
      END DO

      ! Diagnostics
      IF( lk_iomput ) THEN
        zafe(:,:,:) = zirondep(:,:,:) * 1.E-9 * tmask_bgc_closea(:,:,:)      ! zirondep and ironsed are in nmol m^-3 s^-1
        CALL iom_put( "Irondep", zafe  )  ! surface downward net flux of iron
      ENDIF  
      irondep_cmoc(:,:,:) = zirondep(:,:,:)  

      IF( ln_timing )   CALL timing_stop('trc_src_fedep')
  
  END SUBROUTINE trc_src_fedep

  SUBROUTINE trc_src_fesed ( Kbb, Kmm, Krhs ) 
      ! compute iron sources: bottom flux from sediments
      !                       based on CanESM5/CanOE code.
      INTEGER, INTENT(in) :: Kbb, Kmm, Krhs
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
!      CALL iom_open('bathy.orca.nc', inum)
!      CALL iom_get(inum, jpdom_data,'bathy',zcmask(:,:,:), lrowattr=ln_use_jattr)
!      CALL iom_close(inum)
      !
!      DO jk = 1, 5
!        DO jj = 2, jpjm1
!           DO ji = fs_2, fs_jpim1     ! These if required are added with the include statement just above the CONTAINS statement
!              IF( tmask_bgc_closea(ji,jj,jk) /= 0. ) THEN
!                 zmaskt = tmask_bgc_closea(ji+1,jj,jk) * tmask_bgc_closea(ji-1,jj,jk) & 
!                    &   * tmask_bgc_closea(ji,jj+1,jk) * tmask_bgc_closea(ji,jj-1,jk) &
!                    &   * tmask_bgc_closea(ji,jj,jk+1)
!                 IF( zmaskt == 0. )  zcmask(ji,jj,jk ) = MAX( 0.1, zcmask(ji,jj,jk) ) 
!              END IF
!           END DO
!        END DO
!      END DO
!      CALL lbc_lnk('trc_src_fesed', zcmask(:,:,:) , 'T', 1. )      ! lateral boundary conditions on cmask   (sign unchanged)
      zcmask(:,:,:) = 1.
      DO jk = 1, jpk
        DO jj = 1, jpj
           DO ji = 1, jpi
              zexpide   = MIN( 8.,( gdept(ji,jj,jk,Kmm) / 500. )**(-1.5) )
              zdenitide = -0.9543 + 0.7662 * LOG( zexpide ) - 0.235 * LOG( zexpide )**2
              zcmask(ji,jj,jk) = MIN( 1., EXP( zdenitide ) / 0.5 )
              !zcmask(ji,jj,jk) = zcmask(ji,jj,jk) * MIN( 1., EXP( zdenitide ) / 0.5 )
           END DO
        END DO
      END DO
      
      ! Coastal supply of iron
      ! -------------------------
      zironsed(:,:,jpk) = 0._wp
      DO jk = 1, jpkm1
        zironsed(:,:,jk) = sedfeinput0 * zcmask(:,:,jk) / ( e3t(:,:,jk,Kmm) * rday )
      END DO

      DO jj = 1, jpj
       DO ji = 1, jpi
        jk  = mbkt(ji,jj)
        tr(ji,jj,jk,jrfer, Krhs) = tr(ji,jj,jk,jrfer, Krhs) + zironsed(ji,jj,jk) * qfact2
       END DO
      END DO
      
      ! Diagnostics
      IF( lk_iomput ) THEN
        zbfe(:,:,:) = zironsed(:,:,:) * 1.E-9 * tmask_bgc_closea(:,:,:)
        CALL iom_put( "Ironsed", zbfe  )  ! iron from sediments        
      ENDIF  
      ironsed_cmoc(:,:,:) = zironsed(:,:,:)
      !
      DEALLOCATE( zironsed, zbfe, zcmask )

      IF( ln_timing )   CALL timing_stop('trc_src_fesed')
      
  END SUBROUTINE trc_src_fesed
 
  SUBROUTINE trc_src_criver( kt, Krhs, write_rhs_flag )
      ! compute dic and doc sources from rivers
      !                       based on CanESM5/CMOC code.
      INTEGER, INTENT(in)  :: kt, Krhs
      INTEGER              :: ji, jj
      !
      LOGICAL, OPTIONAL, INTENT(in) :: write_rhs_flag   ! 
      LOGICAL                       :: write_rhs_flag0  ! 
      !
      REAL(wp), DIMENSION(jpi,jpj) :: zcoef
      !
      IF( ln_timing )   CALL timing_start('trc_src_criver')
      !
      IF (lwp) THEN
        WRITE(numout,*)
        WRITE(numout,*) 'trc_src: calling trc_src_criver'
        WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*)
      ENDIF
      !
      CALL trc_src2d( kt , js2d_rdic )
      CALL trc_src2d( kt , js2d_rdoc )
      !
      ! O Riche Nov 24th 2022
      ! Test if any external source data point is too large
      WRITE(numout,*)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF( lwp .AND. ABS(src2d_dta(ji,jj,js2d_rdic)) > HUGE(1._wp) ) WRITE(numout,*) 'trc_src_criver: js2d_rdic has reached a huge value at ji = ', ji, ' jj =', jj
          IF( lwp .AND. ABS(src2d_dta(ji,jj,js2d_rdoc)) > HUGE(1._wp) ) WRITE(numout,*) 'trc_src_criver: js2d_rdoc has reached a huge value at ji = ', ji, ' jj =', jj
          CALL FLUSH(numout)
        END DO
      END DO
      WRITE(numout,*)
      !
      zcoef(:,:)       =   ryyssb * cvol(:,:,1)
      cotdep_cmoc(:,:) =   src2d_dta(:,:,js2d_rdic)                              * 1.e9 / (12.  * zcoef(:,:) + rtrn )
      rivinp_cmoc(:,:) = ( src2d_dta(:,:,js2d_rdic) + src2d_dta(:,:,js2d_rdoc) ) * 1.e9 / (31.6 * zcoef(:,:) + rtrn )
      ! RHS terms derived from above external sources
      no3river_cmoc(:,:) = qfact2 *   rivinp_cmoc(:,:)
      dicriver_cmoc(:,:) = qfact2 *   rivinp_cmoc(:,:) * 2.631
      talriver_cmoc(:,:) = qfact2 * ( cotdep_cmoc(:,:) - rivinp_cmoc(:,:) * ncrr_cmoc)
      !
      IF ( .NOT. PRESENT(write_rhs_flag) ) THEN
        write_rhs_flag0 = .true.
      ELSE
        write_rhs_flag0 = write_rhs_flag
      ENDIF
      !
      IF( write_rhs_flag0 ) THEN
        tr(:,:,1,jqno3, Krhs) = tr(:,:,1,jqno3, Krhs) + no3river_cmoc(:,:)
        tr(:,:,1,jqdic, Krhs) = tr(:,:,1,jqdic, Krhs) + dicriver_cmoc(:,:)
        tr(:,:,1,jqtal, Krhs) = tr(:,:,1,jqtal, Krhs) + talriver_cmoc(:,:)
      END IF
      IF( ln_timing )   CALL timing_stop('trc_src_criver')
      !    
  END SUBROUTINE trc_src_criver

  SUBROUTINE trc_bott_cmoc( Kmm, Krhs, write_rhs_flag )
      ! Fate of POC reaching the ocean floor: complete remineralization
      ! into DIC, DIN and sink of O2 and TALK 
      !
      INTEGER, INTENT(in) ::    Kmm, Krhs  ! time level indices
      LOGICAL, OPTIONAL, INTENT(in) :: write_rhs_flag   ! 
      LOGICAL                       :: write_rhs_flag0  ! 
      !      
      INTEGER  :: ji, jj, jk, ikt             !: loop variables
      INTEGER  :: ierr                        !: working variables
      !
      REAL(wp) :: zwsbio32, zwsmax, zdep                   !: working variables
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwsbio3   !: working variables
      !
      IF( ln_timing )   CALL timing_start('trc_bott_cmoc')
      !
      IF (lwp) THEN
        WRITE(numout,*)
        WRITE(numout,*) 'trc_src: calling trc_bott_cmoc'
        WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*)
      ENDIF
      !      
      ALLOCATE( zwsbio3(jpi,jpj,jpk), STAT=ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'trc_bott_cmoc: failed to allocate 3d array' )
      !      
      ! limit the values of the sinking speeds to avoid numerical instabilities
      zwsbio3(:,:,:) = ws_cmoc
      !
      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zwsmax = 0.8 * e3t(ji,jj,jk,Kmm) / xstepb
               zwsbio3(ji,jj,jk) = MIN( zwsbio3(ji,jj,jk), zwsmax )
            END DO
         END DO
      END DO
      !
      IF ( .NOT. PRESENT(write_rhs_flag) ) THEN
        write_rhs_flag0 = .true.
      ELSE
        write_rhs_flag0 = write_rhs_flag
      ENDIF
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstepb / e3t(ji,jj,ikt,Kmm)
            zwsbio32 = zwsbio3(ji,jj,ikt) * zdep
            dicbott_cmoc(ji,jj) =  tr(ji,jj,ikt,jqpoc, Kmm) * zwsbio32 
            talbott_cmoc(ji,jj) = -tr(ji,jj,ikt,jqpoc, Kmm) * zwsbio32 * ncrr_cmoc
            no3bott_cmoc(ji,jj) =  tr(ji,jj,ikt,jqpoc, Kmm) * zwsbio32 
            oxybott_cmoc(ji,jj) = -tr(ji,jj,ikt,jqpoc, Kmm) * zwsbio32 
            pocbott_cmoc(ji,jj) = -tr(ji,jj,ikt,jqpoc, Kmm) * zwsbio32 
            !IF( write_rhs_flag0 ) THEN
              tr(ji,jj,ikt,jqdic, Krhs) = tr(ji,jj,ikt,jqdic, Krhs) + dicbott_cmoc(ji,jj)
              tr(ji,jj,ikt,jqtal, Krhs) = tr(ji,jj,ikt,jqtal, Krhs) + talbott_cmoc(ji,jj)
              tr(ji,jj,ikt,jqno3, Krhs) = tr(ji,jj,ikt,jqno3, Krhs) + no3bott_cmoc(ji,jj)
              tr(ji,jj,ikt,jqoxy, Krhs) = tr(ji,jj,ikt,jqoxy, Krhs) + oxybott_cmoc(ji,jj)
              tr(ji,jj,ikt,jqpoc, Krhs) = tr(ji,jj,ikt,jqpoc, Krhs) + pocbott_cmoc(ji,jj)
            !END IF      
            !
         END DO
      END DO
      !
      DEALLOCATE( zwsbio3 )
      !
      IF( ln_timing )   CALL timing_stop('trc_cmoc_bott')
      !
  END SUBROUTINE trc_bott_cmoc


  SUBROUTINE trc_n2fx_denit_cmoc( zpar, kt, jnt , Kbb, Kmm, Krhs,write_rhs_flag )
      ! compute N2 fixation and denitrification
      ! as prescribed in CanESM5/CMOC
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) :: zpar  ! any PAR array
      INTEGER, INTENT(in) ::    Kbb, Kmm, Krhs  ! time level indices
      !
      LOGICAL, OPTIONAL, INTENT(in) :: write_rhs_flag   ! 
      LOGICAL                       :: write_rhs_flag0  ! 
      !
      INTEGER                       :: ji, jj, jk      ! nested loop indices
      INTEGER, INTENT(in)           :: kt, jnt ! ocean time step
      ! <CMOC code OR 10/15/2015> arrays for total water column remineralisation, 
      ! total euphotic zone nitrogen fixation, temporary array for DNF diagnostics, 
      ! pon flux (euphotic zone bottom) for PIC burial diagnostics, PIC flux at the 
      ! bottom, bottom POC
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zn2fixtot, zwork
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zdenittot
      ! <CMOC code OR 10/15/2015> arrays for depth-dependent rates, zJNd is used to 
      !compute the balance between denitrification and nitrogen fixation
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zn2fix,   zJNd
      REAL(wp)   :: zrtn
      CHARACTER (len=25) :: charout
      !
      ALLOCATE( zn2fix(jpi, jpj, jpk), zJNd (jpi, jpj, jpk) )
      ALLOCATE( zn2fixtot(jpi, jpj  ), zdenittot(jpi, jpj ), zwork(jpi, jpj ) )
      ! Nitrogen fixation and denitrification
      ! ----------------------------------------------------------

      ! <CMOC code OR 10/15/2015> Initialization of CMOC arrays
      zn2fix   (:,:,:) = 0._wp
      zn2fixtot(:,:)   = 0._wp
      ! <CMOC code OR 12/11/2015> Total denitrification diagnostics
      zdenittot(:,:)   = 0._wp
      zJNd     (:,:,:) = 0._wp
      zwork    (:,:)   = 0._wp
      !
      DO jk = 1, jk_eud_cmoc
         DO jj = 1, jpj
            DO ji = 1, jpi
                   zn2fix(ji,jj,jk) = pnf_cmoc * cnrr_cmoc * 1e-12_wp / 3600._wp * qfact2         & ! reference rate
                   !
                   &                 * kn_cmoc * 1e-6_wp / ( kn_cmoc * 1e-6_wp                    &
                   &                                         + tr(ji,jj,jk,jqno3, Kmm) + rtrn)        & ! N inhibition
                   !
                   &                 * zpar(ji,jj,jk) / inf_cmoc                                  & ! ligh sensitivity
                   !
                   &                 * ( max(ts(ji,jj,jk,jp_tem,Kmm), tnfmi_cmoc ) - tnfmi_cmoc )    &
                   &                 / ( tnfMa_cmoc - tnfmi_cmoc ) &                               ! temperature dependence
                   !
                   &                 * ( phinf_cmoc * exp( 1._wp ) * anf_cmoc * gdept(ji,jj,jk,Kmm) &
                   &                 * exp ( -anf_cmoc * gdept(ji,jj,jk,Kmm) ) + phi0_cmoc )        & ! diazotroph abundance dependence
                   &                 * oomask(ji,jj) * tmask_bgc_closea(ji,jj,jk)                             ! open ocean / land mask
                   !
                   ! total nitrogen fixation on the current 1/4 time step, is this still true, depends on qnrdttrc
                   zn2fixtot(ji,jj) = zn2fixtot(ji,jj) + zn2fix(ji,jj,jk) * e3t(ji,jj,jk,Kmm)
                   zJNd(ji,jj,jk)   = zn2fix(ji,jj,jk)
               END DO
          END DO
      END DO
      !
      ! Store 3D N2 fixation rate
      n2fix_cmoc(:,:,:) = zJNd(:,:,:)
      !
      ! Denitrication prescribed by N2 fixation and remineralization rates
      DO jk = jk_eud_cmoc+1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                  zJNd(ji,jj,jk)  =  -zn2fixtot(ji,jj) *                                &
                   &                 ( redet(ji,jj,jk) / (redettot(ji,jj) + rtrn) )     & 
                   &                                   * tmask_bgc_closea(ji,jj,jk) * oomask(ji,jj)

                  zdenittot(ji,jj) = zdenittot(ji,jj) + zJNd(ji,jj,jk) * e3t(ji,jj,jk,Kmm)
               END DO
          END DO
      END DO
      !
      ! Store 3D denitrication rate
      denit_cmoc(:,:,:) = zJNd(:,:,:) - n2fix_cmoc(:,:,:)
      !
       WRITE(numout,*) 'DNF sum:', SUM(zn2fixtot(:,:)) , SUM(zdenittot(:,:))

      !     --------------------------------------------------------------------
      !     Update the arrays TRA which contain the biological sources and sinks
      !     --------------------------------------------------------------------
      !
      IF ( .NOT. PRESENT(write_rhs_flag) ) THEN
        write_rhs_flag0 = .true.
      ELSE
        write_rhs_flag0 = write_rhs_flag
      ENDIF
      !
      DO jk = 1, jpkm1
        IF( write_rhs_flag0 ) THEN  
          tr(:,:,jk,jqno3, Kbb) = tr(:,:,jk,jqno3, Kbb) +  zJNd(:,:,jk)
        END IF
      END DO
      !
      ! print mean trends (used for debugging)
      IF( sn_cfctl%l_prttrc )   THEN
         WRITE(charout, FMT="('rem6')")
         CALL prt_ctl_info(charout, cdcomp = 'top')
         CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask_bgc_closea, clinfo=ctrcnm)
      ENDIF
      !
      IF( lk_iomput ) THEN
         IF( jnt == qnrdttrc ) THEN
            ! <CMOC code OR 10/15/2015> 1.e+3_wp is to convert from L^-1 to m^-3
            !  (left in the sum line #119); the diagnostics has to be rescaled 
            ! to per second by dividing by rfact2.
            zwork(:,:)  =  zn2fixtot(:,:) * ncrr_cmoc * 1.e+3_wp * qfact2r * tmask_bgc_closea(:,:,1)
            ! nitrogen fixation in molN m^-2 s^-1 
            CALL iom_put( "Nfix"   , zwork )
            ! <CMOC code OR 12/11/2015> 1.e+3_wp is to convert from L^-1 to 
            ! m^-3 (left in the sum line #119); the diagnostics has to be 
            ! rescaled to per second by dividing by rfact2; NOTE: land mask 
            ! already taken into account
            zwork(:,:)  = -zdenittot(:,:) * ncrr_cmoc * 1.e+3_wp * qfact2r
            CALL iom_put( "Denit"  , zwork ) ! denitrification in molN m^-2 s^-1 
       ENDIF
      ENDIF
      !
      DEALLOCATE(zJNd, zdenittot, zn2fix, zn2fixtot, zwork)
      !
  END SUBROUTINE trc_n2fx_denit_cmoc


  SUBROUTINE trc_n2fx_init_cmoc
      !
      INTEGER ::   ios  
      ! 
      NAMELIST/namcmocnfx/ phinf_cmoc, phi0_cmoc, anf_cmoc, pnf_cmoc, inf_cmoc, tnfMa_cmoc, tnfmi_cmoc
      REWIND( numnatp_refb )              ! Namelist namcmocnfx in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcmocnfx, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcmocnfx in reference namelist_cmoc' )
      REWIND( numnatp_cfgb )              ! Namelist namcmocnfx in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcmocnfx, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcmocnfx in configuration namelist_cmoc' )
      !
      IF(lwm) WRITE( numonpb, namcmocnfx )      
      !
      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for dinitrogen fix. , namcmocnfx'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Maximum ref. diazotroph concentration    phinf_cmoc =',  phinf_cmoc
         WRITE(numout,*) '    Surface ref. diazotroph concentration     phi0_cmoc =',   phi0_cmoc
         WRITE(numout,*) '    Inverse depth of diazotroph conc.max.      anf_cmoc =',    anf_cmoc
         WRITE(numout,*) '    Maximum ref. rate of dinitrogen fix.       pnf_cmoc =',    pnf_cmoc
         WRITE(numout,*) '    Maximum ref. dinitrogen fix surf. irr.     inf_cmoc =',    inf_cmoc
         WRITE(numout,*) '    Maximum ref. dinitrogen fix SST          tnfMa_cmoc =',  tnfMa_cmoc
         WRITE(numout,*) '    Minimum ref. dinitrogen fix SST          tnfmi_cmoc =',  tnfmi_cmoc
         WRITE(numout,*) ' '
      END IF   
      !
  END SUBROUTINE trc_n2fx_init_cmoc

END MODULE trcsrc_canbgc

