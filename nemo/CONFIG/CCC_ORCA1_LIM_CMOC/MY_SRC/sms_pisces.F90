MODULE sms_pisces   
   !!----------------------------------------------------------------------
   !!                     ***  sms_pisces.F90  ***  
   !! TOP :   PISCES Source Minus Sink variables
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-02 (O. Aumont) original code
   !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
   !!            CMOC1 !  2013-15 (O. Riche) add CMOC global variables and namelist parameters
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
  ! <CMOC code OR 10/21/2015> CMOC namelist unit number
   INTEGER ::   numcmoc

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

   !!* Damping 
   LOGICAL  ::   ln_pisdmp         !: relaxation or not of nutrients to a mean value
   INTEGER  ::   nn_pisdmp         !: frequency of relaxation or not of nutrients to a mean value
   LOGICAL  ::   ln_pisclo         !: Restoring or not of nutrients to initial value
                                   !: on close seas
   !!*  Remineralization
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:  ) ::   oomask         ! Open ocean mask 
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::  redet          !: detritus remineralization
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)    ::  redettot       !: detritus remineralization integrated below the euphotic zone

   !!*  Biological fluxes for light
   INTEGER , ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::  neln       !: number of T-levels + 1 in the euphotic layer
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::  heup       !: euphotic layer depth

   !!*  Biological fluxes for primary production
   ! <CMOC code OR 10/21/2015>  iron limitation mask
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::   xlimnfecmoc!:
   ! <CMOC code OR 10/21/2015>  rain ratio 
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:)  ::   xrcico     !:

   ! <CMOC code OR 10/21/2015> CMOC block start
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
   INTEGER     :: jk_eud_cmoc   ! <CMOC code OR 01/23/2016> scale of the euphotic zone
   INTEGER     :: nk_bal_cmoc   ! <CMOC code OR 01/23/2016> open ocean criterion
   ! <CMOC code OR 10/21/2015> CMOC block end


   !!* Variable for chemistry of the CO2 cycle
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akb3       !: pH constant
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak13       !: ...
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak23       !: ...
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksp       !: ...
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akw3       !: ...
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akp13
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akp23
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akp33
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksi3
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   asi3
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   borat      !: borate constant
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hi         !: hydronium concentration
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hj         !: abiotic hydronium concentration
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hk         !: natural hydronium concentration
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hk         !: radiocarbon hydronium concentration
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   excess     !: calcite saturation(>0) / undersaturation(<0) state (p4zlys)

   !!* Array used to indicate negative tracer values
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xnegtr     !: time step correction

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
      ! <CMOC code OR 11/13/2015> removing user-defined DNF diagnostics, revert to PISCES diagnostics
      INTEGER ::   ierr(7)                    ! error handling ! ierr(6)            ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !*  Biological fluxes for light
      ALLOCATE( neln(jpi,jpj), heup(jpi,jpj),    STAT=ierr(1) )
      !
      !*  Biological fluxes for primary production
      ALLOCATE( xlimnfecmoc(jpi,jpj),            STAT=ierr(2) ) !  iron limitation mask
      ALLOCATE( xrcico  (jpi,jpj),               STAT=ierr(3) ) !  rain ratio 
         !
      !* Variable for chemistry of the CO2 cycle
      ALLOCATE( akb3(jpi,jpj,jpk)    , ak13  (jpi,jpj,jpk) ,       &
         &      ak23(jpi,jpj,jpk)    , aksp  (jpi,jpj,jpk) ,       &
         &      akw3(jpi,jpj,jpk)    , borat (jpi,jpj,jpk) ,       &
         &      akp13(jpi,jpj,jpk)   , akp23 (jpi,jpj,jpk) ,       &
         &      akp33(jpi,jpj,jpk)   , aksi3 (jpi,jpj,jpk) ,       &
         &      asi3 (jpi,jpj,jpk)   , hi    (jpi,jpj,jpk) ,       &
         &      hj  (jpi,jpj,jpk)    , hk  (jpi,jpj,jpk)   ,       &
         &      excess(jpi,jpj,jpk)  ,   STAT=ierr(4) )
         !
      !* Array used to indicate negative tracer values  
      ALLOCATE( xnegtr(jpi,jpj,jpk)  ,            STAT=ierr(6) )


      !*  Remineralization
      ALLOCATE(redet(jpi,jpj,jpk), redettot(jpi,jpj),             &
        &     oomask(jpi,jpj), STAT=ierr(7) ) 

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
