MODULE trcsms_dms
   !!======================================================================
   !!                         ***  MODULE trcsms_my_dms  ***
   !! TOP :   Main module of the DMS tracers
   !!======================================================================
   !! History :      ! 2026 (T. Sou, A. Haddon) Ocean DMS 
   !!----------------------------------------------------------------------
   !! trc_sms_dms       : DMS model main routine
   !! trc_sms_dms_alloc : allocate arrays specific to DMS sms
   !!----------------------------------------------------------------------
   USE par_trc                                  ! TOP parameters
   USE par_oce ,  ONLY : jpkm1   
   USE oce_trc                                  ! Ocean variables
   USE trc                                      ! TOP variables
   USE trd_oce
   USE trdtrc                                   ! TESSA TODO:do I need this?
   USE sms_canoe, ONLY : grazing1, grazing2, zmortpn, zmortpd, zprocn, zprocd, xlimdn    ! 
   USE sms_top_canbgc, ONLY : qfact2            ! bgc time-step
   USE sbc_oce , ONLY : sst_m                   ! sea surface temperature (Celsius)
   USE sbc_oce , ONLY : wndm                    ! wind speed module at T-point (=|U10m-Uoce|)  [m/s]
   USE sbc_oce , ONLY : fr_i                    ! ice fraction = 1 - lead fraction      (between 0 to 1)
   USE dom_oce,  ONLY : e3t_0                   !: t- vert. scale factor [m]
   USE trcopt_canbgc, ONLY : par_3bands         ! PAR w/o the mxl averaging and w/ the diurnal cycle if any

   USE ice                                      ! ice variables 
   USE par_dms

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_dms       ! called by trcsms.F90 module
   PUBLIC   trc_sms_dms_alloc ! called by trcini_dms.F90 module

   ! Defined HERE the arrays specific to DMS sms and ALLOCATE them in trc_sms_dms_alloc

   !t ocedms
   REAL(wp), PUBLIC, SAVE :: q_p1, q_p2           ! intracellular DMSPp-to-Carbon ratio for p1,p2 [umol S:mmol C]
   REAL(wp), PUBLIC, SAVE :: f_z1, f_z2           ! sloppy feeding fraction for z1,z2 [-]
   REAL(wp), PUBLIC, SAVE :: f_e1, f_e2           ! exudation fraction for p1,p2 [-]
   REAL(wp), PUBLIC, SAVE :: f_yield              ! DMS yield [-]
   REAL(wp), PUBLIC, SAVE :: k_dmspd, k_dms       ! bacterial DMSPd,DMS consumption rate constant [d-1] (converted to 1/s)
   REAL(wp), PUBLIC, SAVE :: k_free               ! DMSPd free lyase rate constant [d-1] (converted to 1/s)
   REAL(wp), PUBLIC, SAVE :: k_photo              ! photolysis rate constant [d-1] (converted to 1/s)

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zdmsflx                   ! Sea-to-air DMS flux (umolS/m2/s)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: smsdmspd                  ! Uppermost DMSPd SMS (source minus sinks) (umolS/m2/s) 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: smsdms                    ! Uppermost DMS SMS (source minus sinks) (umolS/m2/s)


   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms_dms.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_dms( kt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_dms  ***
      !!
      !! ** Purpose :   main routine of DMS model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      INTEGER  :: ji,jj,jl, jk,jn          ! dummy loop index
      REAL(wp) :: zsch_dms, zk_dms, epsilon15=1.e-15_wp
      REAL(wp) :: zscale               ! scale factor between sea ice skeletal layer and ocean surface layer

      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_dms')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_dms:  OCEAN DMS model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      zdmsflx(:,:)  = 0._wp
      smsdmspd(:,:) = 0._wp
      smsdms(:,:)   = 0._wp

       DO jk = 1, jpkm1
        DO jj = 1, jpj
         DO ji = 1 ,jpi

              ! unit for grazing1, grazing2, zmortpn, zmortpd, zprocn, zprocd  is umolC/L/timestep (same as mmolC/m3/timestep)
              ! unit for qfact2 is seconds in a timestep, xlimndn has no unit
              ! unit for q_p1 and q_p2 is [umol S:mmol C]

              ! Update DMSPd (umolS/m3/s)
              tr(ji,jj,jk,jrdmspd,Krhs) = + q_p1*f_z1*grazing1(ji,jj,jk)/qfact2                           & ! + sloppy feeding      
                                  + q_p2*f_z2*grazing2(ji,jj,jk)/qfact2                                   & ! + sloppy feeding       
                                  + q_p1*zmortpn(ji,jj,jk)/qfact2                                         & ! + cell lysis            
                                  + q_p2*zmortpd(ji,jj,jk)/qfact2*(1./(xlimdn(ji,jj,jk)+0.1))             & ! + cell lysis
                                  + q_p1*zprocn(ji,jj,jk)/qfact2*(f_e1)                                   & ! + exudation 
                                  + q_p2*zprocd(ji,jj,jk)/qfact2*(f_e2+(1.-f_e2)*(1.-xlimdn(ji,jj,jk)))   & ! + exudation
                                  - k_dmspd*tr(ji,jj,jk,jrdmspd,Kbb)                                      & ! - bacterial consumption&           ! TESSA TODO was HH:trn, now Kbb (not Kmm?)
                                  - k_free*tr(ji,jj,jk,jrdmspd,Kbb)                                         ! - free DMSPd-lyase                 ! TESSA TODO was HH:trn, now Kbb

              ! Update DMS
              tr(ji,jj,jk,jrdms,Krhs) = + f_yield*k_dmspd*tr(ji,jj,jk,jrdmspd,Kbb)                                 & ! + bacterial production (yield)     ! TESSA TODO was HH:trn
                                  + k_free*tr(ji,jj,jk,jrdmspd,Kbb)                                                & ! + free DMSPd-lyase                 ! TESSA TODO was HH:trn, now Kbb
                                  - k_dms*tr(ji,jj,jk,jrdms,Kbb)                                                   & ! -bacterial consumption             ! TESSA TODO was HH:trn, now Kbb
                                  - k_photo*par_3bands(ji,jj,jk)/(par_3bands(ji,jj,jk)+1.)*tr(ji,jj,jk,jrdms,Kbb)    !  photolysis  TESSA TODO check par_3bands is correct
                                                                                                                     ! TESSA TODO  was HH:trn (not Kbb), was etot (now par_3bands)
         END DO
        END DO
       END DO

      ! Compute sea-to-air and ice-to-sea fluxes
       DO jj = 1, jpj
        DO ji = 1, jpi
         ! schmidt number calculation
         zsch_dms = 2674.0_wp - 147.12_wp * sst_m(ji,jj)+   3.726_wp *sst_m(ji,jj)**2 - 0.038_wp    * sst_m(ji,jj)**3
         ! 0.01/3600 = [cm/h] to [m/s]; 1-fr_i = lead fraction; Nightingale2000+Luce2011; schmidt number normalization;
         zk_dms =  0.01_wp/3600._wp * ( (1 - fr_i(ji,jj)) **0.4) * (0.222_wp*wndm(ji,jj)**2 + 0.333_wp*wndm(ji,jj)) * SQRT(600./zsch_dms)
         ! Compute DMS flux [umol S m-2 s-1]
         zdmsflx(ji,jj) = tr(ji,jj,1,jrdms,Kbb) * tmask(ji,jj,1) * zk_dms  !TODO was HH:trn
         ! Update DMS due to sea-to-air flux
         tr(ji,jj,1,jrdms,Krhs) = tr(ji,jj,1,jrdms,Krhs) - zdmsflx(ji,jj) / e3t_0(ji,jj,1)  ! TODO HH:tra
        END DO
       END DO

       ! save time derivative (Sink minus source) for dms write
       smsdmspd(:,:) = tr(:,:,1,jrdmspd,Krhs)
       smsdms(:,:) = tr(:,:,1,jrdms,Krhs)

      IF( ln_timing )   CALL timing_stop('trc_sms_dms')
      !
   END SUBROUTINE trc_sms_dms


   INTEGER FUNCTION trc_sms_dms_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_dms_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to DMS
      ! ALLOCATE( tab(...) , STAT=trc_sms_dms_alloc )
      trc_sms_dms_alloc = 0      ! set to zero if no array to be allocated

      ! ocedms
      ALLOCATE(zdmsflx(jpi,jpj), smsdmspd(jpi,jpj), smsdms(jpi,jpj), STAT=trc_sms_dms_alloc)
      IF( trc_sms_dms_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_dms_alloc : failed to allocate arrays' )
      !
   END FUNCTION trc_sms_dms_alloc

   !!======================================================================
END MODULE trcsms_dms
