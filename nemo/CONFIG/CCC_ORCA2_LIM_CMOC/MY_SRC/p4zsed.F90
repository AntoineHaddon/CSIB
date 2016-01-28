MODULE p4zsed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :   PISCES Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) USE of fldread
   !!           CMOC1  !  203-15   (O. Riche) river forcing, bottom POC remineralization, and dissolved iron mask
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_sed        :  Compute loss of organic matter in the sediments
   !!   p4z_sbc        :  Read and interpolate time-varying nutrients fluxes
   !!   p4z_sed_init   :  Initialization of p4z_sed
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE iom             !  I/O manager
   USE fldread         !  time interpolation
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sed   
   PUBLIC   p4z_sed_init   
   PUBLIC   p4z_sed_alloc

   !! * Shared module variables
   LOGICAL  :: ln_river    = .FALSE.    !: boolean for river input of nutrients

   !! * Module variables
   REAL(wp) :: ryyss                  !: number of seconds per year 
   LOGICAL  :: ll_sbc

   !! * File structures
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_riverdic  ! structure of input riverdic
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_riverdoc  ! structure of input riverdoc
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_fmsk      ! structure of input iron limitation mask ! <CMOC code OR 10/22/2015>

   INTEGER , PARAMETER :: nbtimes = 365  !: maximum number of times record in a file
   INTEGER  :: ntimes_riv

   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:) :: rivinp, cotdep    !: river input fields

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Header:$ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_sed( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed  ***
      !!
      !! ** Purpose :   Compute loss of organic matter in the sediments. This
      !!              is by no way a sediment model. The loss is simply 
      !!              computed to balance the inout from rivers and dust
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER  ::   ji, jj, jk, ikt
      REAL(wp) ::   zwsbio3, zdep
      CHARACTER (len=25) :: charout
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sed')
      !

      IF( jnt == 1 .AND. ll_sbc ) CALL p4z_sbc( kt )

      ! Add the external input of nutrients, carbon and alkalinity
      ! ----------------------------------------------------------
      trn(:,:,1,jpno3) = trn(:,:,1,jpno3) + rivinp(:,:) * rfact2
      trn(:,:,1,jpdic) = trn(:,:,1,jpdic) + rivinp(:,:) * 2.631 * rfact2
      trn(:,:,1,jptal) = trn(:,:,1,jptal) + (cotdep(:,:) - rno3*rivinp(:,:) ) * rfact2

      ! Fate of POC reaching the ocean floor: complete remineralization into DIC, DIN
      ! and sink of O2 and TALK
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstep / fse3t(ji,jj,ikt)
            zwsbio3 = wsbio3(ji,jj,ikt) * zdep

            trn(ji,jj,ikt,jpdic) = trn(ji,jj,ikt,jpdic)                       &
               &                             + trn(ji,jj,ikt,jppoc) * zwsbio3 
            trn(ji,jj,ikt,jptal) = trn(ji,jj,ikt,jptal)                       &
               &                             - trn(ji,jj,ikt,jppoc) * zwsbio3 
            trn(ji,jj,ikt,jpno3) = trn(ji,jj,ikt,jpno3)                       &
               &                             + trn(ji,jj,ikt,jppoc) * zwsbio3 
            trn(ji,jj,ikt,jpoxy) = trn(ji,jj,ikt,jpoxy)                       &
               &                             - trn(ji,jj,ikt,jppoc) * zwsbio3 
            trn(ji,jj,ikt,jppoc) = trn(ji,jj,ikt,jppoc)                       &
                                             - trn(ji,jj,ikt,jppoc) * zwsbio3 
         END DO
      END DO

      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sed')
      !
   END SUBROUTINE p4z_sed

   SUBROUTINE p4z_sbc( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_sbc  ***
      !!
      !! ** purpose :   read and interpolate the external sources of 
      !!                nutrients
      !!
      !! ** method  :   read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * local declarations
      INTEGER  :: ji,jj 
      REAL(wp) :: zcoef
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sbc')
      !
      ! N release due to coastal rivers
      ! Compute river at nit000 or only if there is more than 1 time record in river file
      ! -----------------------------------------
      IF( ln_river ) THEN
         IF( kt == nit000 .OR. ( kt /= nit000 .AND. ntimes_riv > 1 ) ) THEN
            CALL fld_read( kt, 1, sf_riverdic )
            CALL fld_read( kt, 1, sf_riverdoc )
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zcoef = ryyss * cvol(ji,jj,1) 
                  cotdep(ji,jj) =   sf_riverdic(1)%fnow(ji,jj,1)                                  * 1E9 / ( 12. * zcoef + rtrn )
                  rivinp(ji,jj) = ( sf_riverdic(1)%fnow(ji,jj,1) + sf_riverdoc(1)%fnow(ji,jj,1) ) * 1E9 / ( 31.6* zcoef + rtrn )
               END DO
            END DO
         ENDIF
      ENDIF

      IF( nn_timing == 1 )  CALL timing_stop('p4z_sbc')
!      !
   END SUBROUTINE p4z_sbc

   SUBROUTINE p4z_sed_init

      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_sed_init  ***
      !!
      !! ** purpose :   initialization of the external sources of nutrients
      !!
      !! ** method  :   read the files and compute the budget
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER  :: ji, jj, jk, jm
      INTEGER  :: numfmsk, numriv                            ! <CMOC code OR 10/22/2015> add numfmsk (iron limitation mask)
      INTEGER  :: ierr, ierr1, ierr2
      REAL(wp), DIMENSION(nbtimes) :: zsteps                 ! times records
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zriverdic, zriverdoc !
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE :: zfmask      ! <CMOC code OR 10/22/2015> add zfmask for iron limitation mask
      !
      CHARACTER(len=100) ::  cn_dir                          ! Root directory for location of ssr files
      TYPE(FLD_N) ::   sn_fmsk, sn_riverdoc, sn_riverdic     ! <CMOC code OR 10/22/2015> add sn_fmsk (iron limitation mask)        ! informations about the fields to be read

      NAMELIST/nampissed/cn_dir, sn_riverdic, sn_riverdoc,                         &
        &                ln_river,                                                 &
        &                sn_fmsk
        
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sed_init')
      !
      !                                    ! number of seconds per year and per month
        ryyss    = nyear_len(1) * rday
      !                            !* set file information
      cn_dir  = './'            ! directory in which the model is executed
      ! ... default values (NB: frequency positive => hours, negative => months)
      !                  !   file       ! frequency !  variable   ! time intep !  clim   ! 'yearly' or ! weights  ! rotation !
      !                  !   name       !  (hours)  !   name      !   (T/F)    !  (T/F)  !  'monthly'  ! filename ! pairs    !
      sn_riverdic = FLD_N( 'river'      ,   -12     ,  'riverdic' ,  .false.   , .true.  ,   'yearly'  , ''       , '' )
      sn_riverdoc = FLD_N( 'river'      ,   -12     ,  'riverdoc' ,  .false.   , .true.  ,   'yearly'  , ''       , '' )
      sn_fmsk     = FLD_N( 'fermask'    ,   -12     ,  'femask'   ,  .false.   , .true.  ,   'yearly'  , ''       , '' )   ! <CMOC code OR 10/22/2015> add the iron limitation mask

      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampissed )

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' namelist : nampissed '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~ '
         WRITE(numout,*) '    river input of nutrients                 ln_river    = ', ln_river
       END IF

      IF( ln_river ) THEN
           ll_sbc = .TRUE.
      ELSE
           ll_sbc = .FALSE.
      ENDIF

      ! iron mask/iron limitation for the ocean ! <CMOC code OR 10/22/2015>
      ! ---------------------------------------
         IF(lwp) WRITE(numout,*) '    initialize CMOC iron limitation'
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         !
         ierr = 0 
         ALLOCATE( sf_fmsk(1), STAT=ierr )           !* allocate and fill sf_fmsk (forcing structure) with sn_fmsk
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_init: unable to allocate sf_fmsk structure' )
         ierr = 0
         
         ! I think this section copy the file structure in memory
         CALL fld_fill( sf_fmsk, (/ sn_fmsk /), cn_dir, 'p4z_sed_init', 'Iron limitation mask', 'nampissed' )
                                   ALLOCATE( sf_fmsk(1)%fnow(jpi,jpj,1), STAT=ierr )  ! fnow current values based on interpolation (OR)?
                                   IF( ierr > 0 ) THEN
                                            CALL ctl_stop( 'p4zsed: iron limitation mask, unable to allocate iron limitation array' )     ;    RETURN 
                                   ENDIF 
         !
         ! Open the the channel numfmsk associated with  file 'sn_fmsk%clname'
         CALL iom_open (  TRIM( sn_fmsk%clname ) , numfmsk )
         IF(lwp) WRITE(numout,*) '    iron mask file opened' 
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         CALL flush(numout)
         ! allocate dynamically mem space for the mask and get the data with iom_get
         ! jpdom_data must hold some information on the grid tile
         ALLOCATE(zfmask(jpi,jpj))
         CALL iom_get  ( numfmsk, jpdom_data, TRIM( sn_fmsk%clvar ), zfmask(:,:), 1 ) ! <CMOC code OR 10/22/2015> get the data from the iron mask file without worrying about the time step
         ! because data time scale is yearly
         IF(lwp) WRITE(numout,*) '    iron mask data read' 
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~'
         CALL flush(numout) 
         CALL iom_close( numfmsk )
           DO ji = 1, jpi
              DO jj = 1, jpj
                    xlimnfecmoc(ji,jj) = zfmask(ji,jj)
              END DO
           END DO
         IF(lwp) WRITE(numout,*) '    iron mask data save in memory' 
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         CALL flush(numout)
         DEALLOCATE(zfmask) 
         IF(lwp) WRITE(numout,*) '    iron limitation initialization: done'
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         CALL flush(numout)  

      ! nutrient input from rivers
      ! --------------------------
      IF( ln_river ) THEN
         ALLOCATE( sf_riverdic(1), STAT=ierr1 )           !* allocate and fill sf_sst (forcing structure) with sn_sst
         ALLOCATE( sf_riverdoc(1), STAT=ierr2 )           !* allocate and fill sf_sst (forcing structure) with sn_sst
         IF( ierr1 + ierr2 > 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_init: unable to allocate sf_apr structure' )
         !
         CALL fld_fill( sf_riverdic, (/ sn_riverdic /), cn_dir, 'p4z_sed_init', 'Input DOC from river ', 'nampissed' )
         CALL fld_fill( sf_riverdoc, (/ sn_riverdoc /), cn_dir, 'p4z_sed_init', 'Input DOC from river ', 'nampissed' )
                                   ALLOCATE( sf_riverdic(1)%fnow(jpi,jpj,1)   )
                                   ALLOCATE( sf_riverdoc(1)%fnow(jpi,jpj,1)   )
         IF( sn_riverdic%ln_tint ) ALLOCATE( sf_riverdic(1)%fdta(jpi,jpj,1,2) )
         IF( sn_riverdoc%ln_tint ) ALLOCATE( sf_riverdoc(1)%fdta(jpi,jpj,1,2) )
 
	 ! Get total input rivers ; need to compute total river supply in a year
         CALL iom_open ( TRIM( sn_riverdic%clname ), numriv )
         CALL iom_gettime( numriv, zsteps, kntime=ntimes_riv)
         ALLOCATE( zriverdic(jpi,jpj,ntimes_riv) )   ;     ALLOCATE( zriverdoc(jpi,jpj,ntimes_riv) )
         DO jm = 1, ntimes_riv
            CALL iom_get( numriv, jpdom_data, TRIM( sn_riverdic%clvar ), zriverdic(:,:,jm), jm )
            CALL iom_get( numriv, jpdom_data, TRIM( sn_riverdoc%clvar ), zriverdoc(:,:,jm), jm )
         END DO
         CALL iom_close( numriv )

      ELSE
         rivinp(:,:) = 0._wp
         cotdep(:,:) = 0._wp

      END IF
!      !
      IF( ll_sbc ) CALL p4z_sbc( nit000 ) 
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sed_init')
      !
   END SUBROUTINE p4z_sed_init

     INTEGER FUNCTION p4z_sed_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed_alloc  ***
      !!----------------------------------------------------------------------

      ALLOCATE( rivinp(jpi,jpj)     , cotdep(jpi,jpj), STAT=p4z_sed_alloc )  

      IF( p4z_sed_alloc /= 0 ) CALL ctl_warn('p4z_sed_alloc : failed to allocate arrays.')

     END FUNCTION p4z_sed_alloc
#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sed                         ! Empty routine
   END SUBROUTINE p4z_sed
#endif 

   !!======================================================================
END MODULE  p4zsed
