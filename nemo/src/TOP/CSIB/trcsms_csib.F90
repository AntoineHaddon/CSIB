MODULE trcsms_csib
   !!======================================================================
   !!                         ***  MODULE trcsms_csib  ***
   !! TOP :   Main module of the CSIB tracers
   !!======================================================================
   !! History :      !  2025 (A. Haddon) Original code
   !!                !  2026 (T. Sou, A. Haddon) Sea ice DMS
   !!----------------------------------------------------------------------
   !! trc_sms_csib       : CSIB model main routine
   !! trc_sms_csib_alloc : allocate arrays specific to CSIB sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables   
   USE dom_oce         ! ocean domain
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc

   USE ice              ! ice variables
   USE phycst           ! physical constants: rhoi, rhos
   USE sbc_oce , ONLY : ssu_m, ssv_m         ! sea surface velocity, for computation of friction velocity
   USE sbc_oce , ONLY : tprecip, sprecip     ! total and solid precipitation
   USE sbc_oce , ONLY : sst_m                ! sea surface temperature (Celsius)
   USE sbc_ice , ONLY : qsr_ice              ! solar heat flux over ice (W m-2)
   USE zdfmxl  , ONLY : nmln                 ! level of mixed layer depth for dic/tak fluxes

   USE par_csib
   USE par_canoe        ! indices of CanOE model variables
   USE sms_canoe , ONLY : rr_c2n ! redfield ratio

   USE par_dms

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_csib       ! called by trcsms.F90 module
   PUBLIC   trc_sms_csib_alloc ! called by trcini_csib.F90 module

!! * Substitutions
#  include "do_loop_substitute.h90"

   !! Sea ice tracers are like sea ice model variables and have 2 versions: 
   !!    - one intensive ("equivalent") for biogeochemistry: in situ concentration 
   !!    - one extenisve ("global") for dynamics and ice transport: grid cell average 
   !! Both are related by sea ice concentration:
   !! grid cell average = in situ concentration           * sea ice concentration
   !!                   = mass / (sea ice area*thickness) * sea ice area / grid cell area  
   !!                   = mass / (grid cell area*thickness) 
   !! (thickness is the thickness of the considered sea ice layer, here only bottom ice)
   !!
   !! ** variables 
   !! name        | indice      |                                 | Unit       |
   !!-------------|-------------|---------------------------------|------------|
   !! icediac     | jridiac     | Ice algae C concentration     | mg C m-3   |
   !! icedian     | jridian     | Ice algae N concentration     | mmol N m-3 |
   !! icediach    | jridiach    | Ice algae Chl concentration   | mg Chl m-3 |
   !! iceno3      | jrino3      | Ice NO3 concentration           | mmol N m-3 |
   !! icenh4      | jrinh4      | Ice NH4 concentration           | mmol N m-3 |
   !!-------------|-------------|---------------------------------|------------|
   !! IF ln_dmsice = .true.
   !! icedmspd    | jridmspd    | Ice DMSPd concentration         | umol S m-3 |
   !! icedms      | jridms      | Ice DMS concentration           | umol S m-3 |
   !!-------------|-------------|---------------------------------|------------|
   !!

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)     :: icetra            !  Ice tracer per ice area (4d: 2d horizontal * ice category * ice tracers)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)     :: icetra_gca        !  Ice tracer grid cell average
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)       :: icetragca_2d      !  Ice tracer grid cell average, 2d version for ice model
   
   !! Biomass ratios (used to convert C fluxes to N or Chl fluxes)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)       :: qnidia             !  Ice algae N/C (gN gC-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)       :: qchidia            !  Ice algae Chl/C (gChl gC-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)       :: qnidiamax          ! Maximum ice algae N/C (gN gC-1)
   

   !!
   !! Sources and sinks
   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flushrate        !  Flushrate per ice category (m s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup            !  Flowrate of water uptake from bottom ice growth per ice category (m s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup            !  Flowrate of water uptake per ice area from lateral ice growth (m s-1)
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flush_dia        !  Loss rate of ice algae from flushing per ice category (mg C m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: slough_dia       !  Loss rate of ice algae from sloughing per ice category (mg C m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lamloss_dia      !  Loss rate from lateral melt per ice category (mg C m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: heatexp_dia      !  Loss rate from heating export per ice category (mg C m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: t_i_b            !  Mean sea ice temperature at previous time step (deg K)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dt_i             !  Mean sea ice temperature change (C s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup_dia        !  Ice algae uptake rate from bottom ice growth per ice category (mg C m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_dia        !  Ice algae uptake rate from lateral ice growth  per ice category (mg C m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nxsicedia        !  N excess export per ice category (mg N m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: cxsicedia        !  C excess export per ice category (mg C m-3 s-1)
      
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flush_no3        !  Loss rate of ice no3 from flushing per ice category (mmol m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flush_nh4        !  Loss rate of ice nh4 from flushing per ice category (mmol m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: slough_no3       !  Loss rate of ice no3 from sloughing per ice category (mmol m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: slough_nh4       !  Loss rate of ice nh4 from sloughing per ice category (mmol m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lamloss_no3      !  Loss rate from lateral melt per ice category (mmol m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lamloss_nh4      !  Loss rate from lateral melt per ice category (mmol m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_no3        !  NO3 uptake rate from lateral ice growth  per ice category (mmol m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_nh4        !  NH4 uptake rate from lateral ice growth  per ice category (mmol m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup_no3        !  No3 uptake rate from bottom ice growth per ice category (mmol m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup_nh4        !  Nh4 uptake rate from bottom ice growth per ice category (mmol m-3 s-1)
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flush_dmspd      ! Loss rate of ice dmspd from flushing per ice category (umolS m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: flush_dms        ! Loss rate of ice dms from flushing per ice category   (umolS m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: slough_dms       ! Loss rate of ice dms from sloughing per ice category (umolS m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: slough_dmspd     ! Loss rate of ice dmspd from sloughing per ice category (umolS m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lamloss_dmspd    ! Loss rate from lateral melt per ice category (umolS m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lamloss_dms      ! Loss rate from lateral melt per ice category (umolS m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_dmspd      ! DMSPd uptake rate from lateral ice growth  per ice category (umolS m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lagup_dms        ! DMS uptake rate from lateral ice growth  per ice category (umolS m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup_dmspd      ! DMSPd uptake rate from bottom ice growth per ice category (umolS m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: bogup_dms        ! DMS uptake rate from bottom ice growth per ice category (umolS m-3 s-1)

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fric_vel         !  Ice friction velocity (m s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: moldif_no3       !  Molecular diffusion rate at ice ocean interface for NO3 per ice category (mmol m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: moldif_nh4       !  Molecular diffusion rate at ice ocean interface for NH4 per ice category (mmol m-3 s-1)
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: par_bi_cat       !  Bottom ice PAR per ice category - recomputed here (W m-2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: phot_dia         !  Ice algae growth rate per ice category (mg C m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lim_PAR          !  Ice algae light limitation factor per ice category (-)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lim_nut          !  Ice algae nutrients limitation factor per ice category (-)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lim_ice          !  Ice growth limitation factor per ice category (-)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: diaupn           !  N uptake by ice algae per ice category (mg N m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: chlsyn           !  Ice algae Chl synthesis rate per ice category (mg Chl m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: mortlin_dia      !  Linear mortality rate of ice algae per ice category (mg C m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: mortquad_dia     !  Quadratic mortality rate of ice algae per ice category (mg C m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: remin_dia        !  N remineralization rate in ice per ice category (mmol m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nitri            !  Nitrification rate in ice per ice category (mmol m-3 s-1)
   
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dmsp_exud        ! DMPSd production from exudation (umolS m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dmsp_lysis       ! DMSPd production from lysis (umolS m-3 s-1)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: dms_phot         ! DMS photolysis (umolS m-3 s-1)

 

 ! model parameters
   REAL(wp), PUBLIC, SAVE ::   z_ia                 ! height of skeletal layer (m)
   ! ice algae
   REAL(wp), PUBLIC, SAVE ::   qnidiamin            ! Mininum ice algae N/C (mmolN mmolC-1)
   REAL(wp), PUBLIC, SAVE ::   cn_fct               ! Factor for target ice algae N/C as function of C:Chl (-)
   REAL(wp), PUBLIC, SAVE ::   cn_pow               ! Power exponent for target ice algae N/C as function of C:Chl (-)
   REAL(wp), PUBLIC, SAVE ::   qchidiaref           ! Reference ice algae Chl:C for uptake from ice growth (gChl gC-1)
   REAL(wp), PUBLIC, SAVE ::   ear=4498._wp         ! activation energy / R gas constant (deg K)
   REAL(wp), PUBLIC, SAVE ::   tempref=298.15_wp    ! Reference temperature for photosynthesis (deg K)
   REAL(wp), PUBLIC, SAVE ::   alrefidia            ! Reference initial slope of photosynthesis-irradiance curve (gC gChl-1 (W m-2)-1 d-1)) 
   REAL(wp), PUBLIC, SAVE ::   pcrefidia            ! Reference photosynthesis rate (d-1) (converted to s-1) 
   REAL(wp), PUBLIC, SAVE ::   betaidia             ! photo-inhibition for ice algae (gC gChl-1 (W m-2)-1 d-1)) (converted to s-1) 
   REAL(wp), PUBLIC, SAVE ::   cigr                 ! Critical ice growth rate for ice algal limitation  (m d-1) (converted to s-1) 
   REAL(wp), PUBLIC, SAVE ::   vnref                ! Reference N uptake rate (gN gC d-1) (converted to s-1)
   REAL(wp), PUBLIC, SAVE ::   knh4                 ! NH4 limitation half saturation constant (mmol N m-3)
   REAL(wp), PUBLIC, SAVE ::   kno3                 ! NO3 limitation half saturation constant (mmol N m-3)
   REAL(wp), PUBLIC, SAVE ::   etares               ! Respiratory cost of biosynthesis (gC gN-1)
   REAL(wp), PUBLIC, SAVE ::   ch2nmax              ! Maximum Chl synthesis rate to N uptake rate (gChl gN-1)
   REAL(wp), PUBLIC, SAVE ::   min_icedia           ! Mortality threshold for ice algae (mmol C m-3)
   REAL(wp), PUBLIC, SAVE ::   t_ia                 ! Temperature sensitivity coefficient for ice algae (C)-1
   REAL(wp), PUBLIC, SAVE ::   r_m1                 ! Linear Mortality rate for ice algae (d-1)  (converted to s-1)
   REAL(wp), PUBLIC, SAVE ::   r_m2                 ! Quadratic Mortality rate for ice algae (mmol C m-3 d-1)  (converted to s-1)
   REAL(wp), PUBLIC, SAVE ::   f_p2                 ! Seeding fraction (-)
   REAL(wp), PUBLIC, SAVE ::   f_flsh               ! Flushing fraction (-)
   REAL(wp), PUBLIC, SAVE ::   f_slgh               ! Sloughing fraction (-)
   REAL(wp), PUBLIC, SAVE ::   dt_mo                ! Heating export sea ice warming threshold (deg C d-1)(converted to s-1)
   REAL(wp), PUBLIC, SAVE ::   t_mo                 ! Heating export sea ice temp trheshold (deg C) + 273.15 = (deg K)
   REAL(wp), PUBLIC, SAVE ::   d_mo                 ! Heating export coeffecient ((deg C mg m-3)-1)
   ! ice N
   REAL(wp), PUBLIC, SAVE ::   f_rm                 ! Remineralization fraction (-)
   REAL(wp), PUBLIC, SAVE ::   r_ni                 ! Nitrification rate (d-1 W m-2)  (converted to s-1)
   REAL(wp), PUBLIC, SAVE ::   c_di                 ! Molecular diffusion coefficient for dissolved nutrients at the ice-water interface (m s-12)
   REAL(wp), PUBLIC, SAVE ::   c_nu                 ! Kinematic viscosity of seawater (m2 s-1)
   ! sea ice C pump
   LOGICAL , PUBLIC, SAVE ::   sicpump              ! Flag for activation of sea ice C pump
   REAL(wp), PUBLIC, SAVE ::   icedicref            ! Sea ice reference DIC (mmol C L-1)
   REAL(wp), PUBLIC, SAVE ::   icetalref            ! Sea ice reference TA (mmol C L-1)
   REAL(wp), PUBLIC, SAVE ::   f_dicsw              ! fraction of DIC rejected into seawater during growth (-)
   REAL(wp), PUBLIC, SAVE ::   f_dicsw_melt         ! fraction of DIC rejected into seawater during melt (-)
   ! DMS
   REAL(wp), PUBLIC, SAVE :: q_pi                   ! intracellular dmsp-to-chlorophyll ratio (umol S:mmol C)
   REAL(wp), PUBLIC, SAVE :: f_zi                   ! sloppy feeding fraction (-)
   REAL(wp), PUBLIC, SAVE :: f_ei                   ! exudation fraction (-)
   REAL(wp), PUBLIC, SAVE :: f_yieldi               ! bacterial conversion fraction (-)
   REAL(wp), PUBLIC, SAVE :: k_dmspdi               ! bacterial dmspd consumption rate constant (d-1) (converted to s-1)
   REAL(wp), PUBLIC, SAVE :: k_dmsi                 ! bacterial dms consumption rate constat (d-1) (converted to s-1)
   REAL(wp), PUBLIC, SAVE :: k_freei                ! free lyase rate constant (d-1) (converted to s-1)
   REAL(wp), PUBLIC, SAVE :: k_photoi               ! photlysis rate constant (d-1) (converted to s-1)
   REAL(wp), PUBLIC, SAVE :: h_ni                   ! half-saturation constant for monod equation
   ! optics
   LOGICAL , PUBLIC, SAVE :: ln_bipar               ! Flag to compute bottom ice PAR: false = use PAR from SI3 ; true = recompute PAR in CSIB
   REAL(wp), PUBLIC, SAVE :: i0_sdry                ! i0 fraction of light transmitted through SSL for dry snow (-)
   REAL(wp), PUBLIC, SAVE :: i0_swet                ! i0 fraction of light transmitted through SSL for wet snow (-)
   REAL(wp), PUBLIC, SAVE :: i0_ice                 ! i0 fraction of light transmitted through SSL for ice (-)
   REAL(wp), PUBLIC, SAVE :: sslh_sdry              ! surface scattering layer height for dry snow (m)
   REAL(wp), PUBLIC, SAVE :: sslh_swet              ! surface scattering layer height for wet snow (m)
   REAL(wp), PUBLIC, SAVE :: sslh_ice               ! surface scattering layer height for ice (m)
   REAL(wp), PUBLIC, SAVE :: parext_swet            ! PAR attenuation coefficent in wet snow (m-1) 
   REAL(wp), PUBLIC, SAVE :: parext_sdry            ! PAR attenuation coefficent in wet snow (m-1)
   REAL(wp), PUBLIC, SAVE :: parext_ice             ! PAR attenuation coefficent in ice (m-1)
  
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_csib( kt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_csib  ***
      !!
      !! ** Purpose :   main routine of CSIB model, sources and sinks from biogeochemistry + ice ocean exchanges
      !!
      !! ** Method  : - called after ice model has done transport of ice BGC (advection+ redistribution)
      !!              - Conversion from extensive (used for ice transport) to intensive (used for bgc) variables
      !!              - compute biogeochemical sources and sinks + sea ice ocean exchange processes
      !!              - timestepping of ice BGC variables with simple explicit Euler scheme (except for NO3 and NH4 diffusion implicit Euler)
      !!              - add to CanOE dynamics sources and sink at ocean surface from sea ice ocean exchange processes
      !!              - Conversion back to extensive variables for ice transport
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt               ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs   ! time level indices
      INTEGER  :: ji,jj,jl,jn          ! dummy loop index
      REAL(wp) :: zsicmin = 1.e-3_wp   ! threshold of ice concentration for ice BGC
      REAL(wp) :: mmc=12._wp           ! Molar mass of Carbon (g mol-1)
      REAL(wp) :: mmn=14._wp           ! Molar mass of Nitrogen (g mol-1)
      REAL(wp) :: zscale               ! scale factor between sea ice skeletal layer and ocean surface layer
      REAL(wp) :: if_below_min=0._wp   ! mortality switch
      REAL(wp) :: zdtDif               ! temp variable for molecular diff
      REAL(wp) :: zmaxia               ! for diagnostics/debug
      REAL(wp) :: zsigup_tot           ! total water uptake from sea ice growth per ice area
      REAL(wp) :: ztemp                ! Temperature factor
      REAL(wp) :: zcnredfield=5.6      ! Redfield C:N ratio (g g-1)
      REAL(wp) :: zcn                  ! Target C:N  (g g-1)
      REAL(wp) :: zpcmaxidia           ! temp variable for photosynthesis
      REAL(wp) :: zalphaidia           ! temp variable for photosynthesis
      REAL(wp) :: znut                 ! N switch
      REAL(wp) :: zlim_nh4, zlim_no3   ! NH4 and NO3 limitation factors
      REAL(wp) :: zrhoch               ! for Chl synthesis
      REAL(wp) :: zsimt                ! mean sea ice temperature
      REAL(wp) :: zsti                 ! sea ice temperature coeffecient for heating export
      REAL(wp) :: if_dti_higher        ! sea ice temperature change switch for heating export
      REAL(wp) :: zno3Old, znh4Old     ! for ice-ocean diffusion 
      REAL(wp) :: zphyn2c              ! phytoplankton N:C
      REAL(wp) :: ztotexp_icediac      ! total expected C export
      REAL(wp) :: zidmspd,zidms        ! temp var DMSPd, DMS 
      REAL(wp) :: zinpp,zlim_nut_dmspd ! for DMSPd 
      REAL(wp) :: zbipar               ! bottom ice PAR
      
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_csib')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_csib:  CSIB model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)

      IF (kt == nittrc000) THEN ! first iteration of run
         ! Check that CanOE is active (done here because CSIB init is in sea ice model init which is before CanOE init)
         IF (.NOT.  ln_canoe) THEN 
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'CSIB on (ln_csib=true in sea ice namelist) but CanOE off (ln_canoe=false in TOP namelist) - deactivating CSIB'
            IF(lwp) WRITE(numout,*) 
            ln_csib=.false.
            RETURN
         ENDIF
         ! For sea ice DMS check that ocean DMS is active
         IF (ln_dmsice .AND. .NOT.ln_dmsoce) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'Sea ice DMS on (ln_dmsice=true in sea ice namelist) but ocean DMS off (ln_dmsoce=true in TOP namelist) - deactivating sea ice DMS'
            IF(lwp) WRITE(numout,*) 
            ln_dmsice=.false.
            jp_csib=5
         ENDIF
      ENDIF

      ! Initiation from ocean surface concentrations (done here because CSIB init is in sea ice model init which is before CanOE init)
      IF ( (kt == nittrc000) .AND. ((.NOT. ln_rsttr) .OR. ln_ibgcspinup) ) THEN
         IF(lwp) WRITE(numout,*) '     Init CSIB from ocean surface'
         DO jl = 1, jpl ! loop ice categories
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( a_i(ji,jj,jl) > zsicmin ) THEN ! if ice
                     icetra(ji,jj,jl,jridiac) = tr(ji,jj,1,jrdia,Kmm) *mmc ! convert to mass units
                     icetra(ji,jj,jl,jridian) = icetra(ji,jj,jl,jridiac) /8._wp ! init at C:N=8
                     icetra(ji,jj,jl,jridiach) = icetra(ji,jj,jl,jridiac) /20._wp ! init at C:Chl=20
                     icetra(ji,jj,jl,jrino3) = tr(ji,jj,1,jqno3,Kmm)
                     icetra(ji,jj,jl,jrinh4) = tr(ji,jj,1,jrnh4,Kmm)
                     IF (ln_dmsice) THEN
                        icetra(ji,jj,jl,jridmspd) = tr(ji,jj,1,jrdmspd,Kmm)
                        icetra(ji,jj,jl,jridms)  = tr(ji,jj,1,jrdms,Kmm)
                     ENDIF
                  END IF
               ENDDO
            ENDDO
         ENDDO
         DO jn = 1,jp_csib
            icetra_gca(:,:,:,jn) = icetra(:,:,:,jn) * a_i(:,:,:)
         ENDDO
      END IF

      ! Conversion from extensive (grid cell average, used by ice transport model) to intensive (in situ concentrations, used by bgc model) variables
      ! grid cell average = in situ concentration * sea ice concentration
      ! Similar to conversion from "global" to "equivalent" varaibles by sea ice model SI3
      ! and set to 0 if low ice concentration 
      DO jl = 1, jpl ! loop ice categories
         DO jj = 1, jpj
            DO ji = 1, jpi
               
               IF( a_i(ji,jj,jl) > zsicmin ) THEN ! if ice
                  DO jn = 1,jp_csib
                     icetra(ji,jj,jl,jn) = icetra_gca(ji,jj,jl,jn) / a_i(ji,jj,jl)
                  ENDDO
                  
               ELSE ! low ice
                  ! set tracers to 0 
                  icetra(ji,jj,jl,:)=0._wp
                  ! send what was in ice to the ocean
                  zscale = z_ia / e3t_0(ji,jj,1) /rn_Dt ! dilution ratio and conversion to a rate
                  tr(ji,jj,1,jrdia,Krhs) = tr(ji,jj,1,jrdia,Krhs) + icetra_gca(ji,jj,jl,jridiac) * zscale /mmc ! convert to mmol
                  tr(ji,jj,1,jrdn ,Krhs) = tr(ji,jj,1,jrdn ,Krhs) + icetra_gca(ji,jj,jl,jridian) * zscale /mmn ! convert to mmol
                  tr(ji,jj,1,jrdch,Krhs) = tr(ji,jj,1,jrdch,Krhs) + icetra_gca(ji,jj,jl,jridiach) * zscale
                  tr(ji,jj,1,jqno3,Krhs) = tr(ji,jj,1,jqno3,Krhs) + icetra_gca(ji,jj,jl,jrino3) * zscale
                  tr(ji,jj,1,jrnh4,Krhs) = tr(ji,jj,1,jrnh4,Krhs) + icetra_gca(ji,jj,jl,jrinh4) * zscale
                  IF (ln_dmsice) THEN 
                     tr(ji,jj,1,jrdmspd,Krhs) = tr(ji,jj,1,jrdmspd,Krhs) + icetra_gca(ji,jj,jl,jridmspd) * zscale
                     tr(ji,jj,1,jrdms ,Krhs) =  tr(ji,jj,1,jrdms  ,Krhs) + icetra_gca(ji,jj,jl,jridms) * zscale
                  ENDIF

               ENDIF ! if ice
               
            ENDDO ! loop jpi
         ENDDO ! loop jpj
      ENDDO ! loop jpl ice categories
      
      ! reset rates: so that they are 0 where there is no ice
      flushrate(:,:,:) = 0._wp
      bogup(:,:,:) = 0._wp
      lagup(:,:,:) = 0._wp
      
      flush_dia(:,:,:) = 0._wp
      slough_dia(:,:,:) = 0._wp
      lamloss_dia(:,:,:) = 0._wp
      dt_i(:,:,:) = 0._wp
      heatexp_dia(:,:,:) = 0._wp
      bogup_dia(:,:,:) = 0._wp
      lagup_dia(:,:,:) = 0._wp
      nxsicedia(:,:,:) = 0._wp
      cxsicedia(:,:,:) = 0._wp

      flush_no3(:,:,:) = 0._wp
      slough_no3(:,:,:) = 0._wp
      lamloss_no3(:,:,:) = 0._wp
      moldif_no3(:,:,:) = 0._wp
      lagup_no3(:,:,:) = 0._wp
      bogup_no3(:,:,:) = 0._wp
      
      flush_nh4(:,:,:) = 0._wp
      slough_nh4(:,:,:) = 0._wp
      lamloss_nh4(:,:,:) = 0._wp
      moldif_nh4(:,:,:) = 0._wp
      lagup_nh4(:,:,:) = 0._wp
      bogup_nh4(:,:,:) = 0._wp
      
      phot_dia(:,:,:) = 0._wp
      lim_PAR(:,:,:) = 0._wp
      lim_nut(:,:,:) = 0._wp
      lim_ice(:,:,:) = 0._wp
      diaupn(:,:,:) = 0._wp
      chlsyn(:,:,:) = 0._wp
      mortlin_dia(:,:,:) = 0._wp
      mortquad_dia(:,:,:) = 0._wp
      remin_dia(:,:,:) = 0._wp
      nitri(:,:,:) = 0._wp
      
      qnidia(:,:,:) = 0._wp
      qchidia(:,:,:) = 0._wp
      qnidiamax(:,:,:) = 0._wp

      IF (ln_dmsice) THEN
         flush_dmspd(:,:,:) = 0._wp
         slough_dmspd(:,:,:) = 0._wp
         lamloss_dmspd(:,:,:) = 0._wp
         lagup_dmspd(:,:,:) = 0._wp
         bogup_dmspd(:,:,:) = 0._wp
         
         flush_dms(:,:,:)   =0._wp
         slough_dms(:,:,:) = 0._wp
         lamloss_dms(:,:,:) = 0._wp
         lagup_dms(:,:,:) = 0._wp
         bogup_dms(:,:,:) = 0._wp

         dmsp_exud(:,:,:) = 0._wp
         dmsp_lysis(:,:,:) = 0._wp
         dms_phot(:,:,:) = 0._wp
      ENDIF

      ! Compute friction velocity, for molecular diffusion at sea ice ocean interface
      CALL ice_friction_velocity
      
      if (ln_bipar) then
         ! Compute bottom ice PAR for photosynthesis
         CALL bottom_ice_PAR
      endif
      
      DO jj = 1, jpj
         DO ji = 1, jpi
            
            ! Bottom ice variables
            DO jl = 1, jpl ! loop ice cat
               
               IF( a_i(ji,jj,jl) > zsicmin ) THEN ! precence of ice

                  ! Ice algal ratios
                  qnidia(ji,jj,jl) = icetra(ji,jj,jl,jridian) /  (icetra(ji,jj,jl,jridiac)+rtrn )
                  qchidia(ji,jj,jl) = icetra(ji,jj,jl,jridiach) /  (icetra(ji,jj,jl,jridiac)+rtrn )
                  
               ! Loss of ice tracers from ice-ocean exchanges

                  ! flushrate: water flowrate per ice area (m s-1) 
                  flushrate(ji,jj,jl) = &
                           ! change of ice thickness from surface melt (m s-1) * ice density (kg m-3) / freshwater density (kg m-3)
                     &     ( dh_sum_cat(ji,jj,jl) ) * rhoi / rhow    &
                           ! flowrate of melt pond drainage volume per ice area (m s-1)
                     &     + dh_mpdrn_cat(ji,jj,jl)      &
                           ! SIC * (total precipation (Kg/m2 s-1) - solid precipitation (Kg/m2 s-1) ) / freshwater density (kg m-3)
                     &     + a_i(ji,jj,jl) * MAX( 0._wp, tprecip(ji,jj) - sprecip(ji,jj) ) / rhow

                  ! flushing of ice tracers: flushrate / skeletal layer height * ice tracers concentration * flushing fraction (-)
                  flush_dia(ji,jj,jl) = flushrate(ji,jj,jl)/z_ia  * icetra(ji,jj,jl,jridiac) * f_flsh  ! ice algae   
                  flush_no3(ji,jj,jl) = flushrate(ji,jj,jl)/z_ia  * icetra(ji,jj,jl,jrino3)  ! ice no3   
                  flush_nh4(ji,jj,jl) = flushrate(ji,jj,jl)/z_ia  * icetra(ji,jj,jl,jrinh4)  ! ice nh4 

                  ! sloughing of ice tracers: bottom ice melt rate / skeletal layer height * ice tracers concentration * flushing fraction (-)
                  slough_dia(ji,jj,jl) = dh_bom_cat(ji,jj,jl) * rhoi / rhow /z_ia  * icetra(ji,jj,jl,jridiac) * f_slgh ! ice algae   
                  slough_no3(ji,jj,jl) = dh_bom_cat(ji,jj,jl) * rhoi / rhow /z_ia  * icetra(ji,jj,jl,jrino3) ! ice no3   
                  slough_nh4(ji,jj,jl) = dh_bom_cat(ji,jj,jl) * rhoi / rhow /z_ia  * icetra(ji,jj,jl,jrinh4) ! ice nh4   
 
                  ! loss of ice tracers from lateral melt : fraction of ice concentration lost (s-1) * ice tracers concentration (mg m-3)
                  lamloss_dia(ji,jj,jl) = da_lam_cat(ji,jj,jl) * icetra(ji,jj,jl,jridiac)   ! ice algae
                  lamloss_no3(ji,jj,jl) = da_lam_cat(ji,jj,jl) * icetra(ji,jj,jl,jrino3)   ! ice no3
                  lamloss_nh4(ji,jj,jl) = da_lam_cat(ji,jj,jl) * icetra(ji,jj,jl,jrinh4)   ! ice nh4

                  ! loss of ice tracers from heating export : heating export coeff ((C mg m-3)-1) * ice temp rate (C s-1) * ice temp change switch (-) * ice temp coeff (-) * biomass^2 ((mg m-3)2)
                  zsimt = SUM(t_i(ji,jj,:,jl)) / nlay_i                                                     ! mean sea ice temperature (deg K)
                  dt_i(ji,jj,jl) = (zsimt - t_i_b(ji,jj,jl) ) / rn_Dt                                       ! mean sea ice temperature change (deg C s-1)
                  if_dti_higher = MAX( 0._wp , SIGN(1._wp, dt_i(ji,jj,jl) - dt_mo) )                        ! ice temp change switch (-)
                  zsti = MIN(1._wp, MAX(0._wp, ( (zsimt - t_mo) / (-1.8_wp+273.15_wp - t_mo) ))**0.2_wp )   ! ice temp coeff (-)
                  heatexp_dia(ji,jj,jl) = d_mo * dt_i(ji,jj,jl) * if_dti_higher * zsti * icetra(ji,jj,jl,jridiac)**2._wp

                  IF (ln_dmsice) THEN
                    ! flushrate / skeletal layer height * ice tracers concentration ( m s-1 /m * umolS m-3 = umolS m-3 s-1)
                    flush_dmspd(ji,jj,jl) = flushrate(ji,jj,jl)/z_ia  * icetra(ji,jj,jl,jridmspd) 
                    flush_dms(ji,jj,jl) =   flushrate(ji,jj,jl)/z_ia  * icetra(ji,jj,jl,jridms) 

                    ! sloughing of ice tracers: bottom ice melt rate   / skeletal layer height * ice tracers concentration * flushing  fraction (-) :: (m s-1 * 1/m * umolS m-3 = umolS m-3 s-1)
                    slough_dmspd(ji,jj,jl) = dh_bom_cat(ji,jj,jl) * rhoi / rhow /z_ia  *  icetra(ji,jj,jl,jridmspd)
                    slough_dms(ji,jj,jl) =   dh_bom_cat(ji,jj,jl) * rhoi / rhow /z_ia  *  icetra(ji,jj,jl,jridms)

                    ! loss of ice tracers from lateral melt : fraction of ice concentration lost  * ice tracers concentration :: (s-1 * umolS m-3 = umolS m-3 s-1))
                    lamloss_dmspd(ji,jj,jl) = da_lam_cat(ji,jj,jl) * icetra(ji,jj,jl,jridmspd) 
                    lamloss_dms(ji,jj,jl) =  da_lam_cat(ji,jj,jl) * icetra(ji,jj,jl,jridms)
                  ENDIF

                  
               ! Uptake of ice tracers from ice growth

                  ! bogup: water uptake flowrate per ice area from bottom ice growth (m s-1) 
                  ! = rate of ice thickness change (m s-1) * ice density (kg m-3) / freshwater density (kg m-3)
                  bogup(ji,jj,jl) = dh_bog_cat(ji,jj,jl) * rhoi / rhow

                  ! Uptake of ice tracers from bottom ice growth 
                  ! = flowrate per ice area (m s-1) * ocean surface concentration (mmol m-3)  s-1keletal layer (m) (*molar mass)
                  bogup_dia(ji,jj,jl) = bogup(ji,jj,jl) * tr(ji,jj,1,jrdia,Kbb) /z_ia *mmc
                  bogup_no3(ji,jj,jl) = bogup(ji,jj,jl) * tr(ji,jj,1,jqno3,Kbb) /z_ia
                  bogup_nh4(ji,jj,jl) = bogup(ji,jj,jl) * tr(ji,jj,1,jrnh4,Kbb) /z_ia

                  ! lagup: water uptake flowrate per ice area from lateral ice growth (m s-1) 
                  ! = (sic increase rate / sic) (1 s-1) * (height new ice) (s) * ice density (kg m-3) / freshwater density (kg m-3)
                  lagup(ji,jj,jl) = da_lag_cat(ji,jj,jl) * ht_i_new(ji,jl) * rhoi / rhow

                  ! Uptake of ice algae from lateral ice growth  
                  lagup_dia(ji,jj,jl) = &
                        ! water uptake flowrate per ice area * skeletal layer * ocean surface algae concentration *molar mass    
                        & lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrdia,Kbb) *mmc     &
                        ! - ice tracer conc * (sic increase rate / sic)
                        & - icetra(ji,jj,jl,jridiac) * da_lag_cat(ji,jj,jl)
                  ! Uptake of ice N from lateral ice growth  
                  lagup_no3(ji,jj,jl) = lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jqno3,Kbb) - icetra(ji,jj,jl,jrino3) * da_lag_cat(ji,jj,jl)
                  lagup_nh4(ji,jj,jl) = lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrnh4,Kbb) - icetra(ji,jj,jl,jrinh4) * da_lag_cat(ji,jj,jl)

                  IF (ln_dmsice) THEN
                    ! Uptake of ice tracers from ice growth :: (m s-1 *  umolS m-3 /m  =  umolS m-3 s-1)
                    bogup_dmspd(ji,jj,jl) = bogup(ji,jj,jl) * tr(ji,jj,1,jrdmspd,Kbb) /z_ia
                    bogup_dms(ji,jj,jl) =   bogup(ji,jj,jl) * tr(ji,jj,1,jrdms,Kbb) /z_ia

                    ! lateral ice growth
                    lagup_dmspd(ji,jj,jl) =  lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrdmspd,Kbb) - icetra(ji,jj,jl,jridmspd)  * da_lag_cat(ji,jj,jl)
                    lagup_dms(ji,jj,jl) =    lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrdms,Kbb) - icetra(ji,jj,jl,jridms)  * da_lag_cat(ji,jj,jl)
                  ENDIF
                  


               ! Biogeochemical processes
                  
                  ! Target C:N dependent on C:Chl 
                  zcn = (cn_fct * (1._wp/(qchidia(ji,jj,jl)+rtrn))**cn_pow + zcnredfield)/2._wp
                  ! Convert to optimal N:C, bounded by Redfield and C:N=20
                  qnidiamax(ji,jj,jl) = MIN(1.0_wp/zcnredfield, MAX(0.05_wp, 1._wp/zcn ))
                  ! N limitation factor (for photosynthesis) (-)
                  lim_nut(ji,jj,jl) = MIN( MAX(0._wp, (qnidia(ji,jj,jl) - qnidiamin) / (qnidiamax(ji,jj,jl) - qnidiamin) ) , 1._wp)

                  ! Temperature factor (-)
                  ztemp = EXP(-ear*( 1._wp/(sst_m(ji,jj)+273.15_wp) - 1._wp/tempref) ) 

                  ! Light limitation factor (-), with inhibition 
                  if (ln_bipar) then 
                     zbipar = par_bi_cat(ji,jj,jl) ! bottom ice PAR (W m-2) computed here
                  else 
                     zbipar = qtr_ice_bot(ji,jj,jl) ! shortwave radiation transmitted through ice (W m-2) computed by sea ice model
                  endif
                  zpcmaxidia = pcrefidia * ztemp * lim_nut(ji,jj,jl) / (qchidia(ji,jj,jl) + rtrn)
                  zalphaidia = alrefidia * qchidia(ji,jj,jl)
                  lim_PAR(ji,jj,jl) = (1._wp - EXP(-zalphaidia/(zpcmaxidia+rtrn) * zbipar ) ) &
                                    & * EXP(-betaidia/(zpcmaxidia+rtrn) * zbipar )   ! inhibition at high light
                  
                  ! Ice growth limitation
                  lim_ice(ji,jj,jl) = MAX( 0._wp , SIGN(1._wp,  1._wp - dh_bog_cat(ji,jj,jl)/cigr ) )
                  
                  ! Photosynthesis rate (mg C m-3 s-1)
                  phot_dia(ji,jj,jl) = pcrefidia * ztemp * lim_nut(ji,jj,jl) * lim_PAR(ji,jj,jl) * lim_ice(ji,jj,jl) * icetra(ji,jj,jl,jridiac)
                  
                  ! N uptake switch (-)
                  znut = ( MIN( MAX(0._wp, (qnidiamax(ji,jj,jl) - qnidia(ji,jj,jl)) / (qnidiamax(ji,jj,jl) - qnidiamin) ), 1._wp) )**0.05_wp

                  ! nh4 and no3 limitation factors (-)
                  zlim_nh4 = icetra(ji,jj,jl,jrinh4) / (icetra(ji,jj,jl,jrinh4) + knh4)
                  zlim_no3 = icetra(ji,jj,jl,jrino3) / (icetra(ji,jj,jl,jrino3) + kno3)

                  ! N uptake by ice algae (mg N m-3 s-1)
                  diaupn(ji,jj,jl) = vnref * ztemp * znut * (zlim_nh4 + (1._wp - zlim_nh4) * zlim_no3) * icetra(ji,jj,jl,jridiac)

                  ! Chl synthesis rate (mg Chl m-3 s-1)
                  zrhoch = zpcmaxidia * lim_PAR(ji,jj,jl) / ( zalphaidia * zbipar + rtrn )
                  chlsyn(ji,jj,jl) = zrhoch * ch2nmax * diaupn(ji,jj,jl)
                  
                  ! Ice algae mortality, off if below threshold min_icedia
                  if_below_min = MAX( 0._wp , SIGN(1._wp, icetra(ji,jj,jl,jridiac) - min_icedia) )
                  ! Linear mortality (mg C m-3 s-1)
                  mortlin_dia(ji,jj,jl)  = if_below_min * r_m1 * exp(t_ia* sst_m(ji,jj)) * icetra(ji,jj,jl,jridiac)
                  ! Quadratic mortality (mg C m-3 s-1)
                  mortquad_dia(ji,jj,jl) = if_below_min * r_m2 * icetra(ji,jj,jl,jridiac)**2._wp

                  ! Remineralization (mmol N m-3 s-1)
                  remin_dia(ji,jj,jl) = f_rm * mortlin_dia(ji,jj,jl) * qnidia(ji,jj,jl) /mmn

                  ! Nitrification, reduced by light
                  nitri(ji,jj,jl) = r_ni / (1.0_wp+zbipar) * icetra(ji,jj,jl,jrinh4)


                  ! Ice algae C biomass dynamics
                  icetra(ji,jj,jl,jridiac) = icetra(ji,jj,jl,jridiac) + rn_Dt * (      &
                              &        - flush_dia(ji,jj,jl)                           & ! loss from flushing 
                              &        - slough_dia(ji,jj,jl)                          & ! loss from sloughing 
                              &        - lamloss_dia(ji,jj,jl)                         & ! loss from lateral melting of ice
                              &        - heatexp_dia(ji,jj,jl)                         & ! loss from heating export
                              &        + bogup_dia(ji,jj,jl)                           & ! Uptake from bottom ice growth
                              &        + lagup_dia(ji,jj,jl)                           & ! Uptake from lateral ice growth
                              &        + phot_dia(ji,jj,jl)                            & ! Photosynthesis
                              &        - etares * diaupn(ji,jj,jl)                     & ! Respiration
                              &        - mortlin_dia(ji,jj,jl)                         & ! linear mortality
                              &        - mortquad_dia(ji,jj,jl)                        & ! quadratic mortality
                              )
                  ! guarantee positive concentration
                  icetra(ji,jj,jl,jridiac) = MAX(0._wp, icetra(ji,jj,jl,jridiac) )
                  ! icetra(ji,jj,jl,jridiac) = MIN(1000._wp, icetra(ji,jj,jl,jridiac) ) 


                  ! Ice algae N biomass dynamics
                  zphyn2c = MIN(qnidiamax(ji,jj,jl), tr(ji,jj,1,jrdn,Kbb)*mmn / ( tr(ji,jj,1,jrdia,Kbb)*mmc +rtrn) ) ! for uptake from ice growth, to limit ice algae N:C to qnidiamx 
                  icetra(ji,jj,jl,jridian) = icetra(ji,jj,jl,jridian) + rn_Dt * (            &
                              &    - flush_dia(ji,jj,jl) * qnidia(ji,jj,jl)                  & ! loss from flushing 
                              &    - slough_dia(ji,jj,jl) * qnidia(ji,jj,jl)                 & ! loss from sloughing 
                              &    - lamloss_dia(ji,jj,jl) * qnidia(ji,jj,jl)                & ! loss from lateral melting of ice 
                              &    - heatexp_dia(ji,jj,jl) * qnidia(ji,jj,jl)                & ! loss from heating export
                              &    + bogup_dia(ji,jj,jl) * zphyn2c                           & ! Uptake from bottom ice growth
                              &    + lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrdia,Kbb)*zphyn2c *mmn - icetra(ji,jj,jl,jridian) * da_lag_cat(ji,jj,jl)   & ! Uptake from lateral ice growth
                              &    + diaupn(ji,jj,jl)                                        & ! NO3+NH4 uptake
                              &    - mortlin_dia(ji,jj,jl) * qnidia(ji,jj,jl)                & ! linear mortality
                              &    - mortquad_dia(ji,jj,jl) * qnidia(ji,jj,jl)               & ! quadratic mortality
                              )
                  ! guarantee positive concentration
                  icetra(ji,jj,jl,jridian) = MAX(0._wp, icetra(ji,jj,jl,jridian) )


                  ! Ice algae Chl biomass dynamics
                  icetra(ji,jj,jl,jridiach) = icetra(ji,jj,jl,jridiach) + rn_Dt * (          &
                              &    - flush_dia(ji,jj,jl) * qchidia(ji,jj,jl)                 & ! loss from flushing
                              &    - slough_dia(ji,jj,jl) * qchidia(ji,jj,jl)                & ! loss from sloughing
                              &    - lamloss_dia(ji,jj,jl) * qchidia(ji,jj,jl)               & ! loss from lateral melting of ice 
                              &    - heatexp_dia(ji,jj,jl) * qchidia(ji,jj,jl)               & ! loss from heating export
                              &    + bogup_dia(ji,jj,jl) * qchidiaref                        & ! Uptake from bottom ice growth 
                              &    + lagup(ji,jj,jl) / z_ia * tr(ji,jj,1,jrdia,Kbb) * qchidiaref - icetra(ji,jj,jl,jridiach) * da_lag_cat(ji,jj,jl)   & ! Uptake from lateral ice growth
                              &    + chlsyn(ji,jj,jl)                                        & ! Chl synthesis
                              &    - mortlin_dia(ji,jj,jl) * qchidia(ji,jj,jl)               & ! linear mortality
                              &    - mortquad_dia(ji,jj,jl) * qchidia(ji,jj,jl)              & ! quadratic mortality
                              )
                  ! guarantee positive concentration
                  icetra(ji,jj,jl,jridiach) = MAX(0._wp, icetra(ji,jj,jl,jridiach) )

                  
                  ! Ice NO3 dynamics
                  icetra(ji,jj,jl,jrino3) = icetra(ji,jj,jl,jrino3) + rn_Dt * (        &
                              &        + bogup_no3(ji,jj,jl)                           & ! Uptake from bottom ice growth
                              &        + lagup_no3(ji,jj,jl)                           & ! Uptake from lateral ice growth
                              &        - flush_no3(ji,jj,jl)                           & ! Loss from flushing 
                              &        - slough_no3(ji,jj,jl)                          & ! Loss from sloughing 
                              &        - lamloss_no3(ji,jj,jl)                         & ! Loss from lateral melting of ice
                              &        - diaupn(ji,jj,jl)/mmn * (1._wp - zlim_nh4) * zlim_no3 / (zlim_nh4 + (1._wp - zlim_nh4) * zlim_no3 +rtrn)   & ! uptake by ice algae
                              &        + nitri(ji,jj,jl)                               & ! nitrification     
                              )
                  ! guarantee positive concentration
                  icetra(ji,jj,jl,jrino3) = MAX(0._wp, icetra(ji,jj,jl,jrino3) )


                  ! Ice NH4 dynamics
                  icetra(ji,jj,jl,jrinh4) = icetra(ji,jj,jl,jrinh4) + rn_Dt * (        &
                              &        + bogup_nh4(ji,jj,jl)                           & ! Uptake from bottom ice growth
                              &        + lagup_nh4(ji,jj,jl)                           & ! Uptake from lateral ice growth
                              &        - flush_nh4(ji,jj,jl)                           & ! Loss from flushing 
                              &        - slough_nh4(ji,jj,jl)                          & ! Loss from sloughing 
                              &        - lamloss_nh4(ji,jj,jl)                         & ! Loss from lateral melting of ice
                              &        - diaupn(ji,jj,jl)/mmn * zlim_nh4 / (zlim_nh4 + (1._wp - zlim_nh4) * zlim_no3 +rtrn)    & ! uptake by ice algae
                              &        - nitri(ji,jj,jl)                               & ! nitrification                
                              &        + remin_dia(ji,jj,jl)                           & ! remineraliztion                
                              )
                  ! guarantee positive concentration
                  icetra(ji,jj,jl,jrinh4) = MAX(0._wp, icetra(ji,jj,jl,jrinh4) )


                  IF (ln_dmsice) THEN
                     zidmspd = icetra(ji,jj,jl,jridmspd) !record DMSPd to have same fluxes for DMSPd and DMS
                     zidms = icetra(ji,jj,jl,jridms)

                     ! nutrient limitation based on monod equation, using h_ni=1 for half-saturation constant
                     zlim_nut_dmspd = (  icetra(ji,jj,jl,jrino3) + icetra(ji,jj,jl,jrinh4)  ) / (  h_ni +  icetra(ji,jj,jl,jrino3) + icetra(ji,jj,jl,jrinh4)  )
                     
                     ! net primary production (Photosynthesis - Respiration)
                     ! NPP (mg C m-3 s-1) = phot_dia (mg C m-3 s-1) - etares (gC gN-1) * diaupn  (mg N m-3 s-1)
                     zinpp =  phot_dia(ji,jj,jl) - etares * diaupn(ji,jj,jl)
                     
                     ! exudation (umolS m-3 s-1) 
                     ! = [f_ei + (1 − f_ei) (1 − lim_nut_dmspd)] * specific growth rate * DMSPp
                     ! =                                         * NPP / ice algal C    * dmsp-to-C * ice algal C
                     ! =                                         * NPP * dmsp-to-C
                     ! q_pi/mmc = umolS mmolC-1 / (g mol-1) = umolS mgC-1
                     dmsp_exud(ji,jj,jl) = (f_ei+(1._wp-f_ei)* (1._wp-zlim_nut_dmspd)) * zinpp * q_pi /mmc 
                     
                     ! lysis (umolS m-3 s-1) 
                     dmsp_lysis(ji,jj,jl) =  q_pi /mmc * mortlin_dia(ji,jj,jl)  * ( 1._wp/(zlim_nut_dmspd+0.1_wp) ) 

                     ! DMS photolysis
                     dms_phot(ji,jj,jl) = k_photoi*icetra(ji,jj,jl,jridms) * zbipar/(zbipar+1._wp)
                     

                     ! Ice DMSP dynamics (umolS/m3)
                     icetra(ji,jj,jl,jridmspd) = icetra(ji,jj,jl,jridmspd) + rn_Dt * (   &
                              &        - flush_dmspd(ji,jj,jl)                           & ! loss from flushing                 (umolS/m3/s)
                              &        - slough_dmspd(ji,jj,jl)                          & ! loss from sloughing                (umolS/m3/s)
                              &        - lamloss_dmspd(ji,jj,jl)                         & ! loss from lateral melting of ice   (umolS/m3/s)
                              &        - k_dmspdi*zidmspd                                & ! loss from bacterial consumption    (1/s*umolS/m3)
                              &        - k_freei*zidmspd                                 & ! loss from free lyase               (1/s*umolS/m3)
                              &        + bogup_dmspd(ji,jj,jl)                           & ! uptake from bottom ice growth      (umolS/m3/s)
                              &        + lagup_dmspd(ji,jj,jl)                           & ! uptake from lateral ice growth     (umolS/m3/s)
                              &        + dmsp_lysis(ji,jj,jl)                            & ! lysis (umolS/m3/s)
                              &        + dmsp_exud(ji,jj,jl)                             & ! exudation (umolS/m3/s)
                              &        )                                                                                       
                     ! guarantee positive concentration
                     icetra(ji,jj,jl,jridmspd) =  MAX(0._wp,icetra(ji,jj,jl,jridmspd))

                     ! Ice DMS dynamics (umolS/m3)
                     icetra(ji,jj,jl,jridms) = icetra(ji,jj,jl,jridms) + rn_Dt * (              &
                              &        - flush_dms(ji,jj,jl)                                    & ! loss from flushing
                              &        - slough_dms(ji,jj,jl)                                   & ! loss from sloughing
                              &        - lamloss_dms(ji,jj,jl)                                  & ! loss from lateral melting of ice
                              &        + f_yieldi*k_dmspdi*zidmspd                              & ! gain from bacterial consumption of dmspd
                              &        + k_freei*zidmspd                                        & ! gain from free lyase
                              &        - k_dmsi*zidms                                           & ! loss from bacterial consumption of dms
                              &        + bogup_dms(ji,jj,jl)                                    & ! uptake from bottom ice growth
                              &        + lagup_dms(ji,jj,jl)                                    & ! uptake from lateral ice growth
                              &        - dms_phot(ji,jj,jl)                                     & ! loss from photlysis 
                              )
                     ! guarantee positive concentration
                     icetra(ji,jj,jl,jridms) =  MAX(0._wp,icetra(ji,jj,jl,jridms))

                  ENDIF ! if DMS


               ! Sea ice C pump
                  IF( sicpump ) THEN
                     ! total water uptake from ice growth: sea ice conc * ( seawater volume uptake from bottom + lateral ice growth per ice area) (m s-1)
                     zsigup_tot = a_i(ji,jj,jl) * ( bogup(ji,jj,jl) + lagup(ji,jj,jl) ) 
                     
                     ! DIC flux from ice growth
                     tr(ji,jj,nmln(ji,jj)-1,jqdic, Krhs) = tr(ji,jj,nmln(ji,jj)-1,jqdic, Krhs) + zsigup_tot * ( tr(ji,jj,nmln(ji,jj)-1,jqdic, Kbb) - icedicref ) * f_dicsw / e3t_0(ji,jj,nmln(ji,jj)-1)
                     ! DIC flux from ice melt
                     tr(ji,jj,1,jqdic, Krhs) = tr(ji,jj,1,jqdic, Krhs) + a_i(ji,jj,jl) * flushrate(ji,jj,jl) * (tr(ji,jj,1,jqdic, Kbb) - icedicref ) * f_dicsw_melt / e3t_0(ji,jj,1)
                     
                     ! TA flux from ice growth
                     tr(ji,jj,nmln(ji,jj)-1,jqtal, Krhs) = tr(ji,jj,nmln(ji,jj)-1,jqtal, Krhs) + zsigup_tot * ( tr(ji,jj,nmln(ji,jj)-1,jqtal, Kbb) - icetalref ) / e3t_0(ji,jj,nmln(ji,jj)-1)
                     ! TA flux from ice melt
                     tr(ji,jj,1,jqtal, Krhs) = tr(ji,jj,1,jqtal, Krhs) + a_i(ji,jj,jl) * flushrate(ji,jj,jl) * (tr(ji,jj,1,jqtal, Kbb) - icetalref ) / e3t_0(ji,jj,1)
                  ENDIF ! if sea ice C pump


               ! Diffusion of N at ice ocean interface: 
                  ! explicit euler scheme can cause numerical problems if large time step 
                  ! instead implicit euler scheme after all other calculations
                  zdtDif = rn_Dt * c_di / c_nu * abs(fric_vel(ji,jj)) / z_ia
                  zscale = a_i(ji,jj,jl) * z_ia / e3t_0(ji,jj,1)
                  ! NO3
                  zno3Old=icetra(ji,jj,jl,jrino3)
                  icetra(ji,jj,jl,jrino3) = ( (1._wp+zdtDif*zscale)*icetra(ji,jj,jl,jrino3) + zdtDif*tr(ji,jj,1,jqno3,Kmm) ) / ( 1._wp + zdtDif*(1._wp+zscale) )
                  tr(ji,jj,1,jqno3,Kmm) = ( zdtDif*zscale*zno3Old + (1._wp+zdtDif)*tr(ji,jj,1,jqno3,Kmm) ) / ( 1._wp + zdtDif*(1._wp+zscale) )
                  ! NH4
                  znh4Old=icetra(ji,jj,jl,jrinh4)
                  icetra(ji,jj,jl,jrinh4) = ( (1._wp+zdtDif*zscale)*icetra(ji,jj,jl,jrinh4) + zdtDif*tr(ji,jj,1,jrnh4,Kmm) ) / ( 1._wp + zdtDif*(1._wp+zscale) )
                  tr(ji,jj,1,jrnh4,Kmm) = ( zdtDif*zscale*znh4Old + (1._wp+zdtDif)*tr(ji,jj,1,jrnh4,Kmm) ) / ( 1._wp + zdtDif*(1._wp+zscale) )

                  ! Diffusion flux, for ouput, as change in ice No3/h4 before and after diffusion flux
                  moldif_no3(ji,jj,jl) = (icetra(ji,jj,jl,jrino3) - zno3Old) /rn_Dt
                  moldif_nh4(ji,jj,jl) = (icetra(ji,jj,jl,jrinh4) - znh4Old) /rn_Dt


               ENDIF ! if ice
               
            ENDDO ! loop jpl ice categories


         ! ocean surface variables affected by exchanges between sea ice and ocean
            ! done after fluxes have been calculated for all ice categories 
            DO jl = 1, jpl ! loop ice cat
         
               IF( a_i(ji,jj,jl) > zsicmin ) THEN ! precence of ice

                  ! scaling factor: conversion of flowrate per sea ice area to flowrate per unit volume
                  ! = sea ice concentration / ocean surface layer height
                  ! = sea ice area / (cell area * ocean surface layer height)
                  zscale = a_i(ji,jj,jl) / e3t_0(ji,jj,1)

               ! Ocean surface large phytoplankton C biomass
                  tr(ji,jj,1,jrdia,Krhs) = tr(ji,jj,1,jrdia,Krhs) + zscale * (    &
                  &        + f_p2 * flush_dia(ji,jj,jl) * z_ia                    & ! flushing of ice algae
                  &        + f_p2 * slough_dia(ji,jj,jl) * z_ia                   & ! sloughing of ice algae
                  &        + f_p2 * lamloss_dia(ji,jj,jl) * z_ia                  & ! ice algae from lateral melting of ice
                  &        + f_p2 * heatexp_dia(ji,jj,jl) * z_ia                  & ! ice algae from heating export
                  &        - bogup_dia(ji,jj,jl) * z_ia                           & ! Uptake from bottom ice growth
                  &        - lagup(ji,jj,jl) * tr(ji,jj,1,jrdia,Kbb)*mmc          & ! Uptake from lateral ice growth
                  ) /mmc ! convert to mmol

               ! Ocean surface large phytoplankton N biomass
                  zphyn2c = MIN(qnidiamax(ji,jj,jl), tr(ji,jj,1,jrdn,Kbb)*mmn / ( tr(ji,jj,1,jrdia,Kbb)*mmc +rtrn) ) ! for uptake from ice growth, to limit ice algae N:C to qnidiamx 
                  tr(ji,jj,1,jrdn,Krhs) = tr(ji,jj,1,jrdn,Krhs) + zscale * (              &
                  &        + f_p2 * flush_dia(ji,jj,jl) * qnidia(ji,jj,jl) * z_ia         & ! flushing of ice algae 
                  &        + f_p2 * slough_dia(ji,jj,jl) * qnidia(ji,jj,jl) * z_ia        & ! sloughing of ice algae 
                  &        + f_p2 * lamloss_dia(ji,jj,jl) * qnidia(ji,jj,jl) * z_ia       & ! ice algae from lateral melting of ice
                  &        + f_p2 * heatexp_dia(ji,jj,jl) * qnidia(ji,jj,jl) * z_ia       & ! ice algae from heating export
                  &        - bogup_dia(ji,jj,jl) * zphyn2c * z_ia                         & ! Uptake from bottom ice growth
                  &        - lagup(ji,jj,jl) * tr(ji,jj,1,jrdia,Kbb)*mmc * zphyn2c        & ! Uptake from lateral ice growth
                  ) /mmn ! convert to mmol


               ! Ocean surface large phytoplankton Chl biomass
                  tr(ji,jj,1,jrdch,Krhs) = tr(ji,jj,1,jrdch,Krhs) + zscale * (         &
                  &        + f_p2 * flush_dia(ji,jj,jl) * qchidia(ji,jj,jl) * z_ia     & ! flushing of ice algae 
                  &        + f_p2 * slough_dia(ji,jj,jl) * qchidia(ji,jj,jl) * z_ia    & ! sloughing of ice algae 
                  &        + f_p2 * lamloss_dia(ji,jj,jl) * qchidia(ji,jj,jl) * z_ia   & ! ice algae from lateral melting of ice
                  &        + f_p2 * heatexp_dia(ji,jj,jl) * qchidia(ji,jj,jl) * z_ia   & ! ice algae from heating export
                  &        - bogup(ji,jj,jl) * tr(ji,jj,1,jrdch,Kbb)                   & ! Uptake from bottom ice growth
                  &        - lagup(ji,jj,jl) * tr(ji,jj,1,jrdch,Kbb)                   & ! Uptake from lateral ice growth
                  )

               ! Ocean surface large POC
                  ztotexp_icediac = (1.0_wp - f_p2) * flush_dia(ji,jj,jl)                  & ! flushing of ice algae
                  &               + (1.0_wp - f_p2) * slough_dia(ji,jj,jl)                 & ! sloughing of ice algae
                  &               + (1.0_wp - f_p2) * lamloss_dia(ji,jj,jl)                & ! ice algae from lateral melting of ice
                  &               + (1.0_wp - f_p2) * heatexp_dia(ji,jj,jl)                & ! ice algae from heating export
                  &               + (1.0_wp - f_rm) * mortlin_dia(ji,jj,jl)                & ! linear mortality 
                  &               + mortquad_dia(ji,jj,jl)                                   ! quadratic mortality 
                  ! C or N excess export
                  ! Detritus are assumed to have fixed C:N redfield ratio (rr_c2n)
                  ! With C export (ztotexp_icediac) there is an implicit N export (= ztotexp_icediac/rr_c2n). 
                  ! If ice algal C:N is different from redfield, there is a mismatch in between the implicit N flux and expected N flux (= ztotexp_icediac*qnidia)
                  ! N excess if qnidia > 1/rr_c2n: difference between expected N export and implicit N export
                  nxsicedia(ji,jj,jl) = MAX(qnidia(ji,jj,jl)-1.0_wp/rr_c2n, 0._wp) * ztotexp_icediac
                  ! C excess otherwise: difference between expected C export and expected N flux converted to a C flux with redfied ration
                  cxsicedia(ji,jj,jl) = (1.0_wp - MIN(qnidia(ji,jj,jl)*rr_c2n, 1.0_wp)) * ztotexp_icediac
                  ! total C export potentially reduced if excess C
                  tr(ji,jj,1,jrgoc,Krhs) = tr(ji,jj,1,jrgoc,Krhs) + zscale * z_ia * MIN(qnidia(ji,jj,jl)*rr_c2n, 1.0_wp) * ztotexp_icediac /mmc

               ! Ocean surface DIC
                  tr(ji,jj,1,jqdic,Krhs) = tr(ji,jj,1,jqdic, Krhs) + zscale * z_ia * cxsicedia(ji,jj,jl)/mmc *1.E-6 ! C excess from export directly remineralized

               ! Ocean surface NO3 
                  tr(ji,jj,1,jqno3,Krhs) = tr(ji,jj,1,jqno3,Krhs) + zscale * (   &
                  &        - bogup_no3(ji,jj,jl) * z_ia                          & ! Uptake from bottom ice growth
                  &        - lagup(ji,jj,jl) * tr(ji,jj,1,jqno3,Kbb)             & ! Uptake from lateral ice growth
                  &        + flush_no3(ji,jj,jl) * z_ia                          & ! Input from flushing 
                  &        + slough_no3(ji,jj,jl) * z_ia                         & ! Input from sloughing 
                  &        + lamloss_no3(ji,jj,jl) * z_ia                        & ! Input from lateral melting of ice
                  )

               ! Ocean surface NH4 
                  tr(ji,jj,1,jrnh4,Krhs) = tr(ji,jj,1,jrnh4,Krhs) + zscale * (   &
                  &        + nxsicedia(ji,jj,jl) * z_ia /mmn                     & ! N excess from export directly remineralized
                  &        - bogup_nh4(ji,jj,jl) * z_ia                          & ! Uptake from bottom ice growth
                  &        - lagup(ji,jj,jl) * tr(ji,jj,1,jrnh4,Kbb)             & ! Uptake from lateral ice growth
                  &        + flush_nh4(ji,jj,jl) * z_ia                          & ! Input from flushing 
                  &        + slough_nh4(ji,jj,jl) * z_ia                         & ! Input from sloughing 
                  &        + lamloss_nh4(ji,jj,jl) * z_ia                        & ! Input from lateral melting of ice
                  )

                  IF (ln_dmsice) THEN
                     ! Ocean surface DMSP
                     tr(ji,jj,1,jrdmspd,Krhs) = tr(ji,jj,1,jrdmspd,Krhs) + zscale * (    &
                     &        + flush_dmspd(ji,jj,jl)   * z_ia                  & ! flushing of DMSPd
                     &        + slough_dmspd(ji,jj,jl)  * z_ia                  & ! sloughing of DMSPd
                     &        + lamloss_dmspd(ji,jj,jl) * z_ia                  & ! DMSPd from lateral melting of ice
                     &        - bogup_dmspd(ji,jj,jl)   * z_ia                  & ! Uptake from bottom ice growth
                     &        - lagup_dmspd(ji,jj,jl)   * z_ia                  & ! Uptake from lateral ice growth
                     )

                     ! Ocean surface DMS
                     tr(ji,jj,1,jrdms,Krhs) = tr(ji,jj,1,jrdms,Krhs) + zscale * (    &
                     &        + flush_dms(ji,jj,jl) * z_ia                    & ! flushing of DMS
                     &        + slough_dms(ji,jj,jl) * z_ia                   & ! sloughing of DMS
                     &        + lamloss_dms(ji,jj,jl) * z_ia                  & ! DMS from lateral melting of ice
                     &        - bogup_dms(ji,jj,jl) * z_ia                    & ! Uptake from bottom ice growth
                     &        - lagup_dms(ji,jj,jl) * z_ia                    & ! Uptake from lateral ice growth
                     )
                  ENDIF

               ENDIF ! if ice
               
            ENDDO ! loop jpl ice categories
               
         ENDDO ! loop jpi
      ENDDO ! loop jpj

      ! record sea ice mean temperature for computation of sea ice temperature tendency for heating export
       DO jj = 1, jpj
         DO ji = 1, jpi
            DO jl = 1, jpl 
               t_i_b(ji,jj,jl) = SUM(t_i(ji,jj,:,jl)) / nlay_i  
            ENDDO
         ENDDO
      ENDDO

      ! Conversion from intensive/equivalent to extensive/global variables for advection by sea ice model
      DO jn = 1,jp_csib
         icetra_gca(:,:,:,jn) = icetra(:,:,:,jn) * a_i(:,:,:)
      ENDDO

      IF( ln_timing )   CALL timing_stop('trc_sms_csib')
      !
   END SUBROUTINE trc_sms_csib



   SUBROUTINE ice_friction_velocity
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_friction_velocity ***
      !!
      !! ** Purpose :   compute friction velocity, for molecular diffusion at sea ice-ocean interface
      !!                code adapted from ICE/icesbc.F90/ice_flx_other
      !!
      !! ** Inputs  :   u_ice, v_ice, ssu_m, ssv_m
      !! ** Outputs :   fric_vel
      !!-----------------------------------------------------------------------
      INTEGER  ::   ji, jj             ! dummy loop indices
      REAL(wp) ::   zu_io, zv_io, zu_iom1, zv_iom1
      !!-----------------------------------------------------------------------
      DO_2D( 0, 0, 0, 0 )
         zu_io   = u_ice(ji  ,jj  ) - ssu_m(ji  ,jj  )
         zu_iom1 = u_ice(ji-1,jj  ) - ssu_m(ji-1,jj  )
         zv_io   = v_ice(ji  ,jj  ) - ssv_m(ji  ,jj  )
         zv_iom1 = v_ice(ji  ,jj-1) - ssv_m(ji  ,jj-1)
         !
         fric_vel(ji,jj) = SQRT( rn_cio * 0.5_wp * ( zu_io*zu_io + zu_iom1*zu_iom1 + zv_io*zv_io + zv_iom1*zv_iom1 ) ) * tmask(ji,jj,1)
      END_2D
   END SUBROUTINE ice_friction_velocity



   SUBROUTINE bottom_ice_PAR
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE bottom_ice_PAR ***
      !!
      !! ** Purpose :   compute bottom ice PAR from qsr_ice 
      !!                   surface scaterring layer (SSL) transmittes only a fraction i0 of light
      !!                   attenuation of light by snow and ice with beer lambert 
      !!                only for sea ice BGC - allows to use different attenuation parameters than used by sea ice model which is focused on thermodynamics
      !!                code adapted from 
      !!                   SBC/sbcblk for surface scatering layer
      !!                   ICE/icethd_zdf_bl99 for attenuation in snow and ice
      !!-----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl         ! dummy loop indices
      REAL(wp) ::   zi0                ! i0 fraction of light transmitted through surface scattering layer 
      REAL(wp) ::   zh0                ! surface scattering layer height
      REAL(wp) ::   zraext_s           ! snow extinction ceof
      REAL(wp) ::   zqtr_ssl           ! light penetrating surface scattering layer
      REAL(wp) ::   zqtr_top_ice       ! light reaching surface of ice
      !!-----------------------------------------------------------------------

       DO jj = 1, jpj
         DO ji = 1, jpi
            DO jl = 1, jpl  ! loop ice cat

               IF( h_i(ji,jj,jl) < 0.1_wp ) THEN ! little or no ice 
                  par_bi_cat(ji,jj,jl) = 0._wp
               
               ELSE ! sea ice is at least 10 cm
                  IF( h_s(ji,jj,jl) > 0._wp ) THEN ! if snow: SSL in snow

                     IF( t_su(ji,jj,jl) < rt0 ) THEN  ! sea ice surface temperature < 0 deg Celsius ->  no surface melting : dry snow
                        zi0 = i0_sdry ! 1._wp (no SSL)
                        zh0 = sslh_sdry  ! 0._wp
                        zraext_s = parext_sdry ! 7._wp 
                     ELSE ! surface melting: wet melting snow
                        zi0 = i0_swet ! 0.45_wp
                        zh0 = sslh_swet ! 0.03_wp
                        zraext_s = parext_swet ! 5._wp  
                     END IF
                     
                     ! light penetrating surface scaterring layer
                     ! qsr_ice : downwelling shortwave * ( 1 - albedo)
                     zqtr_ssl = zi0 * qsr_ice(ji,jj,jl)

                     ! light at top of ice
                     ! effects of melt ponds : weigted sum with melt pond fraction a_ip_frac
                     zqtr_top_ice = (1.0_wp - a_ip_frac(ji,jj,jl)) * zqtr_ssl * EXP( - zraext_s * MAX( 0._wp, h_s(ji,jj,jl) - zh0 ) ) &
                     &              + a_ip_frac(ji,jj,jl) * qsr_ice(ji,jj,jl) * EXP( - zraext_s * h_s(ji,jj,jl)  )  ! under meltpond : no SSL

                     ! PAR in bottom ice layer
                     par_bi_cat(ji,jj,jl) = zqtr_top_ice * EXP( - parext_ice * h_i(ji,jj,jl) ) 
                     
                  ELSE ! no snow: SSL in sea ice 
                     zi0= i0_ice ! 0.26_wp
                     zh0= sslh_ice ! 0.1_wp
                     zqtr_ssl = zi0 * qsr_ice(ji,jj,jl)

                     ! effects of melt ponds : weigted sum with melt pond fraction a_ip_frac
                     par_bi_cat(ji,jj,jl) = (1.0_wp - a_ip_frac(ji,jj,jl)) * zqtr_ssl * EXP( - parext_ice * MAX( 0._wp, h_i(ji,jj,jl)-zh0) ) &
                     &                      + a_ip_frac(ji,jj,jl) * qsr_ice(ji,jj,jl) * EXP( - parext_ice * h_i(ji,jj,jl) ) ! under meltpond : no SSL

                  END IF
               END IF

            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE bottom_ice_PAR


   INTEGER FUNCTION trc_sms_csib_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_csib_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to CSIB
      ! ALLOCATE( tab(...) , STAT=trc_sms_csib_alloc )
      trc_sms_csib_alloc = 0      ! set to zero if no array to be allocated
      
      ALLOCATE( &
      ! variables
         &     icetra(jpi,jpj,jpl,jp_csib) , icetra_gca (jpi,jpj,jpl,jp_csib) , icetragca_2d(jpij,jpl,jp_csib) , &
         &     qnidia(jpi,jpj,jpl), qchidia(jpi,jpj,jpl), qnidiamax(jpi,jpj,jpl),  &
      ! sources and sinks
         &     flushrate (jpi,jpj,jpl) , flush_dia  (jpi,jpj,jpl) , lamloss_dia (jpi,jpj,jpl) , & 
         &     slough_dia(jpi,jpj,jpl) , heatexp_dia(jpi,jpj,jpl) ,                             &
         &     t_i_b     (jpi,jpj,jpl) , dt_i       (jpi,jpj,jpl) ,                             &
         &     bogup     (jpi,jpj,jpl) , bogup_dia  (jpi,jpj,jpl) , lagup       (jpi,jpj,jpl) , &
         &     lagup_dia (jpi,jpj,jpl) , lagup_no3  (jpi,jpj,jpl) , lagup_nh4   (jpi,jpj,jpl) , &
         &     nxsicedia (jpi,jpj,jpl) , cxsicedia  (jpi,jpj,jpl) , par_bi_cat  (jpi,jpj,jpl) , &
         &     phot_dia  (jpi,jpj,jpl) , mortlin_dia(jpi,jpj,jpl) , mortquad_dia(jpi,jpj,jpl) , &
         &     lim_PAR   (jpi,jpj,jpl) , lim_nut    (jpi,jpj,jpl) , lim_ice     (jpi,jpj,jpl) , &
         &     diaupn    (jpi,jpj,jpl) , chlsyn     (jpi,jpj,jpl) ,                             &
         &     remin_dia (jpi,jpj,jpl) , nitri      (jpi,jpj,jpl) ,                             &
         &     flush_no3 (jpi,jpj,jpl) , lamloss_no3(jpi,jpj,jpl) , moldif_no3  (jpi,jpj,jpl) , slough_no3  (jpi,jpj,jpl) , &
         &     flush_nh4 (jpi,jpj,jpl) , lamloss_nh4(jpi,jpj,jpl) , moldif_nh4  (jpi,jpj,jpl) , slough_nh4  (jpi,jpj,jpl) , &
         &     bogup_no3 (jpi,jpj,jpl) , bogup_nh4  (jpi,jpj,jpl) ,                             &
         &     fric_vel  (jpi,jpj)     ,                                                        &
         &     STAT=trc_sms_csib_alloc)
      IF( trc_sms_csib_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_csib_alloc : failed to allocate arrays' )

      IF (ln_dmsice) THEN
         ALLOCATE( &
            & flush_dmspd(jpi,jpj,jpl)  , slough_dmspd(jpi,jpj,jpl), lamloss_dmspd(jpi,jpj,jpl), lagup_dmspd(jpi,jpj,jpl), bogup_dmspd(jpi,jpj,jpl), &
            & flush_dms(jpi,jpj,jpl)    , slough_dms(jpi,jpj,jpl)  , lamloss_dms(jpi,jpj,jpl)  , lagup_dms(jpi,jpj,jpl)  , bogup_dms(jpi,jpj,jpl), &
            & dmsp_exud(jpi,jpj,jpl)    , dmsp_lysis(jpi,jpj,jpl)  , dms_phot(jpi,jpj,jpl)     ,               &
            & STAT=trc_sms_csib_alloc)
         IF( trc_sms_csib_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_csib_alloc : failed to allocate dms arrays' )
      ENDIF

      !
   END FUNCTION trc_sms_csib_alloc

   !!======================================================================
END MODULE trcsms_csib
