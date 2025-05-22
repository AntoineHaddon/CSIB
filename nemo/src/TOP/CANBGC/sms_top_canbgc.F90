MODULE sms_top_canbgc
        !!----------------------------------------------------------------------
        !!                     ***  sms_bgcm.F90  ***  
        !! TOP :   All common variables to all BGCMs
        !!----------------------------------------------------------------------
        !! History :   1.0  !  2000-02 (O. Aumont) original code
        !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
        !!                  !  2022-06 (O. Riche) formerly sms_pisces.F90
        !!                  !                     necessary and can't be plugged
        !!                     !                  in trcsms.F90 without
        !!                     !                  circular dependency err. in make
        !!                  !  2022-11 (J. Christian) 2D and 3D carbon chem modes
        !!----------------------------------------------------------------------
        USE par_oce
        USE par_trc
        USE oce_trc                       ! O Riche June 29th 2022: to be able to use numout file ID

        IMPLICIT NONE
        PUBLIC

        !!*  Time variables
        INTEGER,  SAVE ::   qnrdttrc          !: BGCM time-step multiplier, i.e. ocean time step * qnrdttrc = BGCM time step
        INTEGER,  SAVE ::   qndayflxtr        !: use to check for new day when updating carbon chemistry state in CANOE    
        REAL(wp), SAVE ::   qfact , qfactr    !: BGCM time-step qfact, and its inverse, both used if Euler scheme is in use
        REAL(wp), SAVE ::   qfact2, qfact2r   !: BGCM time-step qfact2, and inverse,    both used if Leap-Frog scheme is in use   

  !!* Mass conservation
  ! LOGICAL  ::  ln_check_mass_canoe  !: Flag to check mass conservation

        !!* Variable for chemistry of the CO2 cycle
        !! 3D carbon chem
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qakb3       !: K_B
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qak13       !: K1
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qak23       !: K_2
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qaksp       !: K_sp
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qakw3       !: K_w
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qakp13      !: K_PO4_1
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qakp23      !: K_PO4_2
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qakp33      !: K_PO4_3
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qaksi3      !: K_si
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qasi3       !: [Si(OH)4-]
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qborat3     !: Borate conc.
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qhi         !: [H+] to compute pH and alkalinity
        !! 2D carbon chem
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qakb2  
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qak12 
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qak22
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qakw2
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qakp12
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qakp22
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qakp32
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qaksi2 
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qasi2
        REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qborat2 

         ! External sources, fluxes, and quantities for BGCMs
         ! 3d array indices
        INTEGER,  PUBLIC, SAVE :: js3d_no3
        INTEGER,  PUBLIC, SAVE :: js3d_si
        INTEGER,  PUBLIC, SAVE :: js3d_po4
        INTEGER,  PUBLIC, SAVE :: js3d_doc
        INTEGER,  PUBLIC, SAVE :: js3d_fe
        INTEGER,  PUBLIC, SAVE :: js3d_hyfe  ! hydrovent iron
        ! 2d array indices
        INTEGER,  PUBLIC, SAVE :: js2d_chla
        INTEGER,  PUBLIC, SAVE :: js2d_dust
        INTEGER,  PUBLIC, SAVE :: js2d_par
        INTEGER,  PUBLIC, SAVE :: js2d_femask ! CMOC iron mask
        INTEGER,  PUBLIC, SAVE :: js2d_ndep   ! atm. N deposition
        INTEGER,  PUBLIC, SAVE :: js2d_rdic
        INTEGER,  PUBLIC, SAVE :: js2d_rdoc 
        INTEGER,  PUBLIC, SAVE :: js2d_rpoc
        INTEGER,  PUBLIC, SAVE :: js2d_fsol1
        INTEGER,  PUBLIC, SAVE :: js2d_fsol2
        
        REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)  ::   src3d_dta       !: 3d source arrays
        REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,  :)  ::   src2d_dta       !: 2d source arrays

        ! CMOC-specific arrays
        REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)      ::   xlimnfecmoc     !: iron mask

!========================================================
! sms_bgcm section
!========================================================

        INTEGER ::   numnatp_refb = -1           !! Logical unit for the ref namelist for the parameters of one of the CANBGC models
        INTEGER ::   numnatp_cfgb = -1           !! Logical unit for the cfg namelist for the parameters of one of the CANBGC models
        INTEGER ::   numonpb      = -1           !! Logical unit for the above ref/cfg namelists output

        !!*  Time variables
        REAL(wp) ::   xstepb             !: Time step duration for biology
        REAL(wp) ::   xfactb             !: possibly useful as 1 thousand divided by biology/BGC time step
        REAL(wp) ::   ryyssb             !: number of seconds per year 
        REAL(wp) ::   r1_ryyssb          !: inverse number of seconds per year 

        PUBLIC sms_top_alloc
        PUBLIC trc_xnegtr

   CONTAINS

    SUBROUTINE trc_xnegtr( jptra0 , jptra1, Kbb, Kmm, Krhs , qnegtr )
      ! 
      ! Check the effect of the trend on the current array
      ! and if any tracer goes beyond zero reduce the time step
      ! inside the whole trn array
      ! Assume a Leapfrog scheme, but trb will be set to trn
      ! already if this is the 1st time step in the calling
      ! subroutine, i.e. trcsms_cmoc or trcsms_canoe.
      !
      USE trc, ONLY: tr
      !
      INTEGER,                          INTENT(in)    ::  jptra0, jptra1   !: tracer indices
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) :: qnegtr
      INTEGER             ::  jn, ji, jj, jk   !: dummy loop indices
      REAL(wp)            ::  ztra
      !
      IF( lwp ) THEN
        WRITE(numout,*)
        WRITE(numout,*), 'trc_xnegtr: correct offshooting tra trend array'
        WRITE(numout,*), 'by reducing the time step for all the tracers. '
        WRITE(numout,*), '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*)
        CALL FLUSH(numout)
      END IF
      !
      qnegtr(:,:,:) = 1.e0      
      !
      DO jn = jptra0, jptra1
        DO jk = 1, jpk
           DO jj = 1, jpj
              DO ji = 1, jpi
                 IF( ( tr(ji,jj,jk,jn, Kbb) + tr(ji,jj,jk,jn, Krhs) ) < 0.e0 ) THEN
                    ztra             = ABS( ( tr(ji,jj,jk,jn, Kbb) - rtrn ) & 
                    &                     / ( tr(ji,jj,jk,jn, Krhs) + rtrn ) )
                    qnegtr(ji,jj,jk) = MIN( qnegtr(ji,jj,jk),  ztra )
                 ENDIF
             END DO
           END DO
        END DO
      END DO    
      !                                ! where at least 1 tracer concentration becomes negative
      !                                ! and by tracer we mean only the CMOC or shared BGC tracer.
      !
    END SUBROUTINE trc_xnegtr

    INTEGER FUNCTION sms_top_alloc()
      !!----------------------------------------------------------------------
      !!        *** ROUTINE sms_pisces_alloc ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_stop
      USE iom
    !
      INTEGER ::   ierr(22)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !* Variables for CO2 chemistry 
      ALLOCATE(         qak13(jpi,jpj,jpk)      ,  qakb3(jpi,jpj,jpk)   ,     &
                &       qak23(jpi,jpj,jpk)      ,  qaksp(jpi,jpj,jpk)   ,     &
                &       qhi  (jpi,jpj,jpk)      ,  qakw3(jpi,jpj,jpk)   ,     &
                &      qakp13(jpi,jpj,jpk)      , qakp23(jpi,jpj,jpk)   ,     &
                &      qakp33(jpi,jpj,jpk)      , qaksi3(jpi,jpj,jpk)   ,     &
                &     qborat3(jpi,jpj,jpk)      ,  qasi3(jpi,jpj,jpk)   ,     &
                &       qak12(jpi,jpj)          ,  qakb2(jpi,jpj)       ,     &
                &       qak22(jpi,jpj)          ,  qakw2(jpi,jpj)       ,     &
                &      qakp12(jpi,jpj)          , qakp22(jpi,jpj)       ,     &
                &      qakp32(jpi,jpj)          , qaksi2(jpi,jpj)       ,     &
                &     qborat2(jpi,jpj)          ,  qasi2(jpi,jpj)       ,     & 
                &     STAT=ierr(1) )
    !* CMOC PP terms/factors
      ALLOCATE( xlimnfecmoc(jpi,jpj)     , STAT=ierr(2))

    !ALLOCATE( etotb(jpi,jpj,jpk), nelnb(jpi,jpj), heupb(jpi,jpj),    &
    ! &       heup_01b(jpi,jpj) , STAT=ierr(3) )
    ! &       heup_01b(jpi,jpj) , xksib(jpi,jpj)               ,  STAT=ierr(2) )
  
      sms_top_alloc = MAXVAL( ierr )
                !
      IF( sms_top_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sms_top_alloc: failed to allocate arrays' ) 
                !
    ! OR Jan 19th 2023
    ! default concentration (useful to flag the other layers when only trc_che_2D is used).
      qhi(:,:,:) = 1.e-9_wp
    !
    END FUNCTION sms_top_alloc

END MODULE sms_top_canbgc

