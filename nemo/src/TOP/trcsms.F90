MODULE trcsms
   !!======================================================================
   !!                         ***  MODULE trcsms  ***
   !! TOP :   Time loop of passive tracers sms
   !!======================================================================
   !! History :   1.0  !  2005-03 (O. Aumont, A. El Moussaoui) F90
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_sms        :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc            !
   USE trc                !

   USE trcsms_canoe       ! CANOE  biogeo-model
   USE trcsms_cmoc        ! CMOC   biogeo-model
   USE par_csib           ! sea ice biogeochemistry model CSIB parameters
   USE trcsms_csib        ! CSIB  tracers
   USE trcsms_dms         ! ocean DMS  tracers
   USE trcsms_pisces      ! PISCES biogeo-model
   USE trcsms_cfc         ! CFC 11 &/or 12
   USE trcsms_c14         ! C14 
   USE trcsms_age         ! AGE
   USE trcsms_my_trc      ! MY_TRC  tracers
   USE prtctl             ! Print control for debbuging
   USE sms_top_canbgc

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms    ! called in trcstp.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms.F90 13286 2020-07-09 15:48:29Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms( kt, Kbb, Kmm , Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms  ***
      !!
      !! ** Purpose :   Managment of the time loop of passive tracers sms 
      !!
      !! ** Method  : -  call the main routine of of each defined tracer model
      !! -------------------------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt        ! ocean time-step index      
      INTEGER, INTENT( in ) ::   Kbb, Kmm, Krhs ! time level indices
      !!
      !!
      INTEGER               ::   jn      ! BGC tracer indexocean time-step index
      INTEGER               ::   jp_tot  ! total number of BGC tracers (shared TOP + activated CanBGC model)
      !!
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms')

      !
      IF( kt == nit000 .AND. lwp) THEN
        WRITE(numout,*)
        WRITE(numout,*) ' trc_sms:  shared BGC processes'
        WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      ENDIF
      !
      ! O Riche Oct 24th 2022
      ! CanBGC BGCMS - adding ln_cmoc/ln_canoe conditional branching
      IF( ln_cmoc .OR. ln_canoe) THEN
        ! total number of shared TOP + CanBGC tracers
        jp_tot = jp_bgc + jp_cmoc                     ! assume CMOC has benen activated
        IF( ln_canoe )  jp_tot = jp_bgc + jp_canoe    ! if assumption above is wrong
        !
        qfact = rDt_trc
        !
        IF( ( ln_top_euler .AND. kt == nittrc000 )  .OR. ( .NOT.ln_top_euler .AND. kt <= nittrc000 + 1 ) ) THEN
          qfactr  = 1. / qfact
          qfact2  = qfact / REAL( qnrdttrc, wp )  ! time split of BGC time step if qnrdttrc is greater than 1.
          qfact2r = 1. / qfact2
          xstepb  = qfact2 / rday    ! time step converted to per day (using in-sec values of time step and day duration)
                                     ! or the fraction of day that is the current time step
          xfactb  = 1.e3 * qfact2r   ! 1 thousand divided by time step for BGC/biology (could be useful?)
          
          IF(lwp) WRITE(numout,*) 
          IF(lwp) WRITE(numout,*) '                    time step    rdt     = ', rdt
          IF(lwp) WRITE(numout,*) '    Passive Tracer  time step    qfact   = ', qfact
          IF(lwp) write(numout,*) '            Biology time step    qfact2  = ', qfact2
          IF(lwp) WRITE(numout,*) '    Passive Tracer  inverse ts   qfactr  = ', qfactr
          IF(lwp) write(numout,*) '            Biology inverse ts   qfact2r = ', qfact2r
          IF(lwp) WRITE(numout,*)
        ENDIF
        ! O Riche Oct 24th 2022 - adding trb/trn swap as appearing in p4zsms.F90 / PISCES BGC
        ! according to comment in p4zsms.F90 this is for restart mode (l_1st_euler) which means
        ! restarts from with Euler forward otherwise leapfrog) and see namelist for OCE component.
        ! ln_top_euler is for TOP, and is like the condition l_1st_euler == 0 for the 1st time step
        ! (but) starting from rest (not from restart).
        IF( (l_1st_euler) .OR. ln_top_euler ) THEN
           DO jn = 1, jp_tot               !   SMS on tracer without Asselin time-filter
              tr(:,:,:,jn,Kbb) = tr(:,:,:,jn,Kmm)
           END DO
        ENDIF
      ENDIF
      ! End of CanBGC BGCMs
      !
      IF( ln_canoe   )   CALL trc_sms_canoe  ( kt, Kbb, Kmm, Krhs )    ! main program of CANOE  
      IF( ln_cmoc    )   CALL trc_sms_cmoc   ( kt, Kbb, Kmm, Krhs )    ! main program of CMOC   
      IF( ln_csib    )   CALL trc_sms_csib   ( kt, Kbb, Kmm, Krhs )    ! main program of CSIB
      IF( ln_dmsoce )    CALL trc_sms_dms    ( kt, Kbb, Kmm, Krhs )    ! main program of DMS
      IF( ln_pisces  )   CALL trc_sms_pisces ( kt, Kbb, Kmm, Krhs )    ! main program of PISCES 
      IF( ll_cfc     )   CALL trc_sms_cfc    ( kt, Kbb, Kmm, Krhs )    ! surface fluxes of CFC
      IF( ln_c14     )   CALL trc_sms_c14    ( kt, Kbb, Kmm, Krhs )    ! surface fluxes of C14
      IF( ln_age     )   CALL trc_sms_age    ( kt, Kbb, Kmm, Krhs )    ! Age tracer
      IF( ln_my_trc  )   CALL trc_sms_my_trc ( kt, Kbb, Kmm, Krhs )    ! MY_TRC  tracers

      IF(sn_cfctl%l_prttrc) THEN                       ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sms ')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
         CALL prt_ctl( tab4d_1=tr(:,:,:,:,Kmm), mask1=tmask, clinfo=ctrcnm )
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('trc_sms')
      !
   END SUBROUTINE trc_sms

#else
   !!======================================================================
   !!  Dummy module :                                     No passive tracer
   !!======================================================================
CONTAINS
   SUBROUTINE trc_sms( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms
#endif 

   !!======================================================================
END MODULE trcsms
