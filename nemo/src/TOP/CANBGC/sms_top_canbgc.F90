MODULE sms_top_canbgc
	!!----------------------------------------------------------------------
	!!                     ***  sms_bgcm.F90  ***  
	!! TOP :   All common variables to all BGCMs
	!!----------------------------------------------------------------------
	!! History :   1.0  !  2000-02 (O. Aumont) original code
	!!             3.2  !  2009-04 (C. Ethe & NEMO team) style
	!!                  !  2022-06 (O. Riche) formerly sms_pisces.F90
	!!                  !                     necessary and can't be plugged
	!! 				   !                     in trcsms.F90 without
	!! 				   !                     circurlar dependency err. in make
	!!----------------------------------------------------------------------
	USE par_oce
	USE par_trc
  USE oce_trc                       ! O Riche June 29th 2022: to be able to use numout file ID

	IMPLICIT NONE
	PUBLIC

	!!*  Time variables
	INTEGER  ::   nrdttrc           !: BGCM time-step multiplier, i.e. ocean time step * nrdttrc = BGCM time step
	INTEGER  ::   ndayflxtr         !: use to check for new day when updating carbon chemistry state in CANOE    
	REAL(wp) ::   rfact , rfactr    !: BGCM time-step rfact, and its inverse, both used if Euler scheme is in use
	REAL(wp) ::   rfact2, rfact2r   !: BGCM time-step rfact2, and inverse,    both used if Leap-Frog scheme is in use   

  !!* Mass conservation
  ! LOGICAL  ::  ln_check_mass_canoe  !: Flag to check mass conservation

	!!* Variable for chemistry of the CO2 cycle
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akb3       !: ???
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak13       !: ???
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak23       !: ???
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksp       !: ???
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akw3       !: ???
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akp13
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akp23
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akp33
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksi3      !: [Si(OH)4-]
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   asi3		    !: [Si(OH)4-]
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   borat      !: Borate conc.
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hi         !: [H+] to compute pH and alkalinity

	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   excess     !: needed in sms_pisces/CanOE, excess in phyto as zoop food wrt Redfield ratios
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aphscale   !: absolute pH scale   

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
  ! PP
  REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)      ::   xlimnfecmoc     !: iron mask

!========================================================
! sms_bgcm section
!========================================================

   INTEGER ::   numnatp_refb = -1           !! Logical units for namelist top
   INTEGER ::   numnatp_cfgb = -1           !! Logical units for namelist top
   INTEGER ::   numonpb      = -1           !! Logical unit for namelist top output

   !!* Model used


   !!*  Time variables
   ! INTEGER  ::   nrdttrc           !: ???
   ! REAL(wp) ::   rfact , rfactr    !: ???
   ! REAL(wp) ::   rfact2, rfact2r   !: ???
   REAL(wp) ::   xstepb             !: Time step duration for biology
   REAL(wp) ::   ryyssb             !: number of seconds per year 
   REAL(wp) ::   r1_ryyssb          !: inverse number of seconds per year 


   ! !!*  Biological parameters 


   ! !!*  diagnostic parameters 

   ! !!* restoring

   ! !!* Mass conservation

   ! !!*  PAR variables and diagnostic variables
   ! INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  nelnb  !: number of T-levels + 1 in the euphotic layer
   ! REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heupb  !: euphotic layer depth
   ! REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etotb  !: par (photosynthetic available radiation)
   ! REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etot_ndcyb      !: PAR over 24h in case of diurnal cycle
   ! REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  emoyb           !: averaged PAR in the mixed layer
   ! REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heup_01b !: Absolute euphotic layer depth
   ! ! REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  xksib  !:  LOBSTER : zooplakton closure

   ! !!*  Biological fluxes for primary production

   ! !!*  Sinking speed

   ! !!*  SMS for the organic matter

   !!* Variable for chemistry of the CO2 cycle

   ! !!* Temperature dependancy of SMS terms


	CONTAINS

		INTEGER FUNCTION sms_top_alloc()
		!!----------------------------------------------------------------------
		!!        *** ROUTINE sms_pisces_alloc ***
		!!----------------------------------------------------------------------
		USE lib_mpp , ONLY: ctl_stop
    !
    USE iom
    !
		INTEGER ::   ierr(10)        ! Local variables
    !
		!!----------------------------------------------------------------------
		ierr(:) = 0
		!* Variable for chemistry of the CO2 cycle
		ALLOCATE( 	ak13(jpi,jpj,jpk)      ,  akb3(jpi,jpj,jpk)   ,     &
				&       ak23(jpi,jpj,jpk)      ,  aksp(jpi,jpj,jpk)   ,     &
				&       hi  (jpi,jpj,jpk)      ,  akw3(jpi,jpj,jpk)   ,     &
				&      akp13(jpi,jpj,jpk)      , akp23(jpi,jpj,jpk)   ,     &
				&      akp33(jpi,jpj,jpk)      , aksi3(jpi,jpj,jpk)   ,     &
				&      borat(jpi,jpj,jpk)      ,  asi3(jpi,jpj,jpk)   ,     &
				&   aphscale(jpi,jpj,jpk)      , excess(jpi,jpj,jpk)  ,     &  
        &  	    STAT=ierr(1) )

    !* CMOC PP terms/factors
    ALLOCATE( xlimnfecmoc(jpi,jpj)     , STAT=ierr(2))

      ! ALLOCATE( etotb(jpi,jpj,jpk), nelnb(jpi,jpj), heupb(jpi,jpj),    &
        ! &       heup_01b(jpi,jpj) , STAT=ierr(3) )
        ! ! &       heup_01b(jpi,jpj) , xksib(jpi,jpj)               ,  STAT=ierr(2) )
      ! !
  
      ! !*  Biological fluxes for light 
      ! ALLOCATE(  etot_ndcyb(jpi,jpj,jpk), emoyb(jpi,jpj,jpk)  ,  STAT=ierr(4) ) 

		!
		sms_top_alloc = MAXVAL( ierr )
		!
		IF( sms_top_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sms_top_alloc: failed to allocate arrays' ) 
		!

		
		END FUNCTION sms_top_alloc


END MODULE sms_top_canbgc
