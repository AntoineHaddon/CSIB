MODULE sms_top
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
        USE closea, ONLY : closea_mask,   & ! O Riche June 29th 2022: keep open the option to use mask
            &                       ln_closea,     & ! as it was necessary in NEMO3.4 in ORCA1 config.
            &                       l_sbc_clo
        USE oce_trc                       ! O Riche June 29th 2022: to be able to use numout file ID

	IMPLICIT NONE
	PUBLIC

	!!*  Time variables
	INTEGER  ::   nrdttrc           !: BGCM time-step multiplier, i.e. ocean time step * nrdttrc = BGCM time step
        INTEGER  ::   ndayflxtr         !: use to check for new day when updating carbon chemistry state in CANOE
	REAL(wp) ::   qfact , qfactr    !: BGCM time-step qfact, and its inverse, both used if Euler scheme is in use
	REAL(wp) ::   qfact2, qfact2r   !: BGCM time-step qfact2, and inverse,    both used if Leap-Frog scheme is in use   

	!!* Variable for chemistry of the CO2 cycle
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

	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qakb2  
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qak12 
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qak22
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qakw2
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qakp12
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qakp22
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qakp32
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qaksi2 
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qasi2
	REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   qborat2 

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

        !REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)  ::   src3d_dta       !: 3d source arrays
        !REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,  :)  ::   src2d_dta       !: 2d source arrays

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


	CONTAINS

		INTEGER FUNCTION sms_top_alloc()
		!!----------------------------------------------------------------------
		!!        *** ROUTINE sms_top_alloc ***
		!!----------------------------------------------------------------------
		USE lib_mpp , ONLY: ctl_stop
		INTEGER ::   ierr(22)        ! Local variables
		!!----------------------------------------------------------------------
		ierr(:) = 0
		!* Variable for chemistry of the CO2 cycle
		ALLOCATE( 	qak13(jpi,jpj,jpk)      ,  qakb3(jpi,jpj,jpk)   ,     &
				&       qak23(jpi,jpj,jpk)      ,  qaksp(jpi,jpj,jpk)   ,     &
				&       qhi  (jpi,jpj,jpk)      ,  qakw3(jpi,jpj,jpk)   ,     &
				&      qakp13(jpi,jpj,jpk)      , qakp23(jpi,jpj,jpk)   ,     &
				&      qakp33(jpi,jpj,jpk)      , qaksi3(jpi,jpj,jpk)   ,     &
				&      qborat3(jpi,jpj,jpk)     ,  qasi3(jpi,jpj,jpk)   ,     &
		                &      qak12(jpi,jpj)       ,  qakb2(jpi,jpj)   ,     &
				&       qak22(jpi,jpj)      ,  qakw2(jpi,jpj)   ,     &
				&      qakp12(jpi,jpj)      , qakp22(jpi,jpj)   ,     &
				&      qakp32(jpi,jpj)      , qaksi2(jpi,jpj)   ,     &
				&      qborat2(jpi,jpj)      ,  qasi2(jpi,jpj)  ,     &
                &  	    STAT=ierr(1) )
		!
		sms_top_alloc = MAXVAL( ierr )
                write(numout,*) ierr
                write(numout,*) sms_top_alloc
		!
		IF( sms_top_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sms_top_alloc: failed to allocate arrays' ) 

		END FUNCTION sms_top_alloc

END MODULE sms_top
