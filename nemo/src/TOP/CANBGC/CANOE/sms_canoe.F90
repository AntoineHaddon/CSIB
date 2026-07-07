MODULE sms_canoe   
   !!----------------------------------------------------------------------
   !!                     ***  sms_pisces.F90  ***  
   !! TOP :   PISCES Source Minus Sink variables
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-02 (O. Aumont) original code
   !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
   !!            CANOE2!  2022-23 (O. Riche & J Christian ) add CANOE global variables and namelist parameters
   !!----------------------------------------------------------------------

   USE par_kind   !: access wp kind

   IMPLICIT NONE
   PUBLIC
   
   ! Particles sinking speed
   REAL(wp), SAVE :: ws_canoe
   REAL(wp), SAVE :: ws_canoe2
   REAL(wp), SAVE :: ws_canoec  

   ! qnegtr block skipping switch
   LOGICAL, SAVE :: ln_canoenegtr

   ! Elemental ratios
   REAL(wp) ::   rr_c2n, rr_fe2n
   REAL(wp) ::   rr_n2c, rr_fe2c, rr_c2fe, rr_n2fe
  
   ! Temperature dependencies
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfuncp0  ! Temp. dependence of phytoplankton growth  (0 to distinguish for PISCES array)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfuncz0  ! Temp. dependence of microzooplankton resp (0 to distinguish for PISCES array)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfuncz20 ! Temp. dependence of mesozooplankton resp  (0 to distinguish for PISCES array)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfuncr0  ! Temp. dependence of remineralization      (0 to distinguish for PISCES array)

   ! Iron and nitrogen limitation function of small and large size primary production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xlimnfe0  ! small phyto iron limitation (0 to distinguish for PISCES array)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xlimdfe0  ! large phyto iron limitation (0 to distinguish for PISCES array)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xlimnn    ! small phyto nitrogen limitation
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xlimdn    ! large phyto nitrogen limitation

   ! Chlorophyll-a production arrays
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zprochln  ! nanophytoplankton  chl production 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zprochld  ! microphytoplankton chl production    

   ! Nitrogen fixation and denitrification
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zn2fix
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr

   ! Variables needed for dmsoce
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   grazing1   !: microzooplankton grazing
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   grazing2   !: mesozooplankton grazing on phytoplankton
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   grazing3   !: mesozooplankton grazing on microzooplankton
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zmortpn    !: nanophytoplankton mortality
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zmortpd    !: diatoms mortality
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zprocn     !: primary production by nanophyto
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zprocd     !: primary production by diatom

   
END MODULE sms_canoe
