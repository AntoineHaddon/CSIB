MODULE sms_pisces   
   !!----------------------------------------------------------------------
   !!                     ***  sms_pisces.F90  ***  
   !! TOP :   PISCES Source Minus Sink variables
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-02 (O. Aumont) original code
   !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                         PISCES model
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc

   IMPLICIT NONE
   PUBLIC

   INTEGER ::   numnatp
   INTEGER ::   numcmoc  ! <CMOC OR 03/08/2014> CMOC namelist unit number

   !!*  Time variables
   INTEGER  ::   nrdttrc           !: ???
   INTEGER  ::   ndayflxtr         !: ???
   REAL(wp) ::   rfact , rfactr    !: ???
   REAL(wp) ::   rfact2, rfact2r   !: ???
   REAL(wp) ::   xstep             !: Time step duration for biology

   !!*  Biological parameters 
   REAL(wp) ::   rno3              !: ???
   REAL(wp) ::   o2ut              !: ???
   REAL(wp) ::   po4r              !: ???
   REAL(wp) ::   rdenit            !: ???
   REAL(wp) ::   rdenita           !: ???
   REAL(wp) ::   o2nit             !: ???
   REAL(wp) ::   wsbio, wsbio2     !: ???
   REAL(wp) ::   xkmort            !: ???
   ! <CMOC OR 06/17/2014> Code trimming !  REAL(wp) ::   ferat3            !: ???

   !!* Damping 
   LOGICAL  ::   ln_pisdmp         !: relaxation or not of nutrients to a mean value
   INTEGER  ::   nn_pisdmp         !: frequency of relaxation or not of nutrients to a mean value
   LOGICAL  ::   ln_pisclo         !: Restoring or not of nutrients to initial value
                                   !: on close seas

   !!*  Biological fluxes for light
   INTEGER , ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::  neln       !: number of T-levels + 1 in the euphotic layer
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::  heup       !: euphotic layer depth

   !!*  Biological fluxes for primary production
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::   xksi       !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::   xksimax    !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanono3   !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatno3   !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanonh4   !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatnh4   !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimphy    !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdia    !: ???
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::   xlimnfecmoc!: <CMOC OR 01/22/2014> iron limitation mask
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::   xrcico     !: <CMOC OR 03/13/2014> reverse to 2D xrcico <CMOC OR 02/19/2014> rain ratio 
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::   xn2fixdia  !: <CMOC OR 04/28/2014> fix for the N2-fixation diagnostics
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concdfe    !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concnfe    !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimnfe    !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdfe    !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimsi     !: ???

   ! <CMOC OR 03/08/2014> CMOC block start
   !!*  CMOC model parameters
   !  Phytoplankton Growth
   !REAL(wp)    :: apar_cmoc  
   !REAL(wp)    :: kw_cmoc    
   !REAL(wp)    :: kchl_cmoc  
   REAL(wp)    :: achl_cmoc  
   REAL(wp)    :: thm_cmoc   
   REAL(wp)    :: tau_cmoc   
   REAL(wp)    :: itau_cmoc
   REAL(wp)    :: ep_cmoc    
   REAL(wp)    :: tvm_cmoc   
   REAL(wp)    :: vm_cmoc    
   REAL(wp)    :: kn_cmoc    
   !  Mortality Growth
   REAL(wp)    :: mpd_cmoc   
   REAL(wp)    :: mpd2_cmoc 
   !  Zooplankton
   REAL(wp)    :: rm_cmoc    
   REAL(wp)    :: kp_cmoc    
   REAL(wp)    :: ga_cmoc    
   REAL(wp)    :: mzn_cmoc   
   REAL(wp)    :: mzd_cmoc   
   REAL(wp)    :: mz2_cmoc   
   !   Detritus and remineralization
   REAL(wp)    :: ed_cmoc    
   REAL(wp)    :: ws_cmoc    
   REAL(wp)    :: reref_cmoc 
   !   Calcite parametrization
   REAL(wp)    :: rmcico_cmoc 
   REAL(wp)    :: trcico_cmoc
   REAL(wp)    :: aci_cmoc   
   REAL(wp)    :: dci_cmoc   
   !   Dinitrogen fixation
   REAL(wp)    :: phinf_cmoc 
   REAL(wp)    :: phi0_cmoc  
   REAL(wp)    :: anf_cmoc   
   REAL(wp)    :: pnf_cmoc   
   REAL(wp)    :: inf_cmoc   
   REAL(wp)    :: tnfMa_cmoc 
   REAL(wp)    :: tnfmi_cmoc 
   !   Redfield ratio and euphotic zone
   REAL(wp)    :: cnrr_cmoc  
   REAL(wp)    :: ncrr_cmoc  
   REAL(wp)    :: deup_cmoc  
   REAL(wp)    :: ideup_cmoc  
   ! <CMOC OR 03/08/2014> CMOC block end


   !!*  SMS for the organic matter
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xfracal    !: ??
   ! <CMOC OR 06/17/2014> Code trimming !  REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   nitrfac    !: ??
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xlimbac    !: ??
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xdiss      !: ??
! <CMOC OR 06/17/2014> Code trimming !      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodcal    !: Calcite production
! <CMOC OR 06/13/2014> Code trimming !      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   grazing    !: Total zooplankton grazing

   !!* Variable for chemistry of the CO2 cycle
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akb3       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak13       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak23       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksp       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akw3       !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   borat      !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hi         !: ???
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   excess     !: ???

   !!* Temperature dependancy of SMS terms
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc    !: Temp. dependancy of various biological rates
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc2   !: Temp. dependancy of mesozooplankton rates

   !!* Array used to indicate negative tracer values
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xnegtr     !: ???

! <CMOC OR 06/13/2014> Code trimming !  #if defined key_kriest
   !!*  Kriest parameter for aggregation
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp) ::   xkr_eta                            !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp) ::   xkr_zeta                           !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp) ::   xkr_massp                          !: ???
! <CMOC OR 06/13/2014> Code trimming !     REAL(wp) ::   xkr_mass_min, xkr_mass_max         !: ???
! <CMOC OR 06/13/2014> Code trimming !  #endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: sms_pisces.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sms_pisces_alloc()
      !!----------------------------------------------------------------------
      !!        *** ROUTINE sms_pisces_alloc ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_warn
      INTEGER ::   ierr(7)        ! <CMOC OR 05/21/2014> add another element to the vector because of xn2fixdia line #196, error handling ! ierr(6)            ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !*  Biological fluxes for light
      ALLOCATE( neln(jpi,jpj), heup(jpi,jpj),                   STAT=ierr(1) )
      !
      !*  Biological fluxes for primary production
! <CMOC OR 06/13/2014> Code trimming !        ALLOCATE( xksimax(jpi,jpj)     , xksi(jpi,jpj)        ,       &
! <CMOC OR 06/13/2014> Code trimming !           &      xnanono3(jpi,jpj,jpk), xdiatno3(jpi,jpj,jpk),       &
! <CMOC OR 06/13/2014> Code trimming !           &      xnanonh4(jpi,jpj,jpk), xdiatnh4(jpi,jpj,jpk),       &
! <CMOC OR 06/13/2014> Code trimming !           &      xlimphy (jpi,jpj,jpk), xlimdia (jpi,jpj,jpk),       &
! <CMOC OR 06/13/2014> Code trimming !           &      xlimnfe (jpi,jpj,jpk), xlimdfe (jpi,jpj,jpk),       &
! <CMOC OR 06/13/2014> Code trimming !           &      xlimsi  (jpi,jpj,jpk), concdfe (jpi,jpj,jpk),       &
! <CMOC OR 06/13/2014> Code trimming !           &      concnfe (jpi,jpj,jpk), xlimnfecmoc(jpi,jpj),    STAT=ierr(2) ) ! <CMOC OR 01/22/2014> iron limitation mask
      ALLOCATE( xlimnfecmoc(jpi,jpj),    STAT=ierr(2) ) ! <CMOC OR 01/22/2014> iron limitation mask
         !
      !*  SMS for the organic matter
! <CMOC OR 06/13/2014> Code trimming !        ALLOCATE( xfracal (jpi,jpj,jpk), nitrfac(jpi,jpj,jpk),       &
! <CMOC OR 06/17/2014> Code trimming !        ALLOCATE( nitrfac (jpi,jpj,jpk),        &
        ALLOCATE(                           &
! <CMOC OR 06/13/2014> Code trimming !           &      prodcal(jpi,jpj,jpk) , grazing(jpi,jpj,jpk),       &
! <CMOC OR 06/17/2014> Code trimming !           &      prodcal(jpi,jpj,jpk) ,        &
! <CMOC OR 06/13/2014> Code trimming !           &      xlimbac (jpi,jpj,jpk), xdiss  (jpi,jpj,jpk),       &
         &      xdiss  (jpi,jpj,jpk),       &
         &      xrcico  (jpi,jpj),                             STAT=ierr(3) ) ! <CMOC OR 03/13/2014> reverse to 2D xrcico ! <CMOC OR 02/19/2014> rain ratio 
         !
      !* Variable for chemistry of the CO2 cycle
      ALLOCATE( akb3(jpi,jpj,jpk)    , ak13  (jpi,jpj,jpk) ,       &
         &      ak23(jpi,jpj,jpk)    , aksp  (jpi,jpj,jpk) ,       &
         &      akw3(jpi,jpj,jpk)    , borat (jpi,jpj,jpk) ,       &
         &      hi  (jpi,jpj,jpk)    , excess(jpi,jpj,jpk) ,   STAT=ierr(4) )
         !
      !* Temperature dependancy of SMS terms
! <CMOC OR 06/13/2014> Code trimming !        ALLOCATE( tgfunc(jpi,jpj,jpk)  , tgfunc2(jpi,jpj,jpk) ,   STAT=ierr(5) )
         !
      !* Array used to indicate negative tracer values  
      ALLOCATE( xnegtr(jpi,jpj,jpk)  ,                          STAT=ierr(6) )
      !
      ! <CMOC OR 04/28/2014> fix N2-fixation diagnostics !* Array used to indicate negative tracer values  
      ALLOCATE( xn2fixdia(jpi,jpj),                             STAT=ierr(7) )
      !
      sms_pisces_alloc = MAXVAL( ierr )
      !
      IF( sms_pisces_alloc /= 0 )   CALL ctl_warn('sms_pisces_alloc: failed to allocate arrays') 
      !
   END FUNCTION sms_pisces_alloc

#else
   !!----------------------------------------------------------------------   
   !!  Empty module :                                     NO PISCES model
   !!----------------------------------------------------------------------
#endif
   
   !!======================================================================   
END MODULE sms_pisces    
