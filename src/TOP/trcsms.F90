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
   USE trcsms_pisces      ! PISCES biogeo-model
   USE trcsms_cfc         ! CFC 11 &/or 12
   USE trcsms_c14         ! C14 
   USE trcsms_age         ! AGE
   USE trcsms_my_trc      ! MY_TRC  tracers
   
   USE prtctl_trc         ! Print control for debbuging

   ! TOP-level processes time integration
   USE par_trc            ! TOP parameters
   USE trd_oce			      ! Ocean trends :   set tracer and momentum trend variables	
                          ! contains logical switches and jp_*** style tracer/var. indices
   USE trdtrc			        ! Dummy module??? O Riche July 5th 2022 not sure what's that doing.

   USE sms_top_canbgc

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms       ! called in trcstp.F90
   PUBLIC   trc_sms_alloc ! called by trcini.F90?

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms  ***
      !!
      !! ** Purpose :   Managment of the time loop of passive tracers sms 
      !!
      !! ** Method  : -  call the main routine of of each defined tracer model
      !! -------------------------------------------------------------------------------------
      !!
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER               ::   jn      ! BGC tracer indexocean time-step index
      INTEGER               ::   jp_tot  ! total number of BGC tracers (shared TOP + activated CanBGC model)
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms')

      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms:  shared BGC processes'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      !
      !
      ! O Riche July 05th 2022
      ! undo this as conditional branching
      ! requires up-to-date kt time step.
      ! perhaps move it higher up
      ! in trcstp.F90. Will it be consistent
      ! with restart state? It won't work
      ! necessarily as rfact is used in other
      ! part of TOP, e.g. rfact2 in TRP/trcnxt.F90
      ! it is possible that var. scope will
      ! locally prevent any conflict with global
      ! variables 
      ! O Riche July 05th 2022
      ! moved from trcsms.F90 to here
      ! to be consistent with its location
      ! in CanESM/CanNEMO code and the migration
      ! of common code to TOP-tier modules
      ! O Riche June 13th 2022
      ! Time step header, select the namelist time step
      ! or double it if Leap-Frog and not Euler has been selected as
      ! the time integration scheme
      !
      qfact = r2dttrc
      ! O Riche July 5th 2022
      ! mitigating the impact of this line for now but might want to keep it
      ! or upgrade it in the final version of the code.
      ! qnrdttrc enables biology components of BGCMs to integrate over extra shorter time steps.
      ! not to confuse with nn_dttrc (lumping physics time steps together drive BGCMs over a longer time step than OCE) and rdttrc
      ! the new time step for BGCM tracers if nn_dtrc/=1
      ! qnrdttrc = 4 ! should be read from namelist_pisces (or _canoe) by trcnam_pisces (or _canoe) or perhaps moved to namelist_top
      !
      ! O Riche Sept 13th 2022
      ! added nrdttrc in namelist_top_* in &namtrc_run section    
      !
      ! O Riche Oct 24th 2022
      ! CanBGC BGCMS - adding ln_cmoc/ln_canoe conditional branching
      IF( ln_cmoc .OR. ln_canoe) THEN
        ! total number of shared TOP + CanBGC tracers
        jp_tot = jp_bgc + jp_cmoc                     ! assume CMOC has benen activated
        IF( ln_canoe )  jp_tot = jp_bgc + jp_canoe    ! if assumption above is wrong
        !
        IF( ( ln_top_euler .AND. kt == nittrc000 )  .OR. ( .NOT.ln_top_euler .AND. kt <= nittrc000 + nn_dttrc ) ) THEN
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
        ! according to comment in p4zsms.F90 this is for restart mode (neuler == 0 which means
        ! restarts from with Euler forward otherwise leapfrog) and see namelist for OCE component.
        ! ln_top_euler is for TOP, and is like the condition neuler == 0 for the 1st time step
        ! (but) starting from rest (not from restart).
        IF( ( neuler == 0 .AND. kt == nittrc000 ) .OR. ln_top_euler ) THEN
           DO jn = 1, jp_tot               !   SMS on tracer without Asselin time-filter
              trb(:,:,:,jn) = trn(:,:,:,jn)
           END DO
        ENDIF
      ENDIF
      ! End of CanBGC BGCMs
      !
      IF( ln_canoe   )   CALL trc_sms_canoe  ( kt )    ! main program of CANOE  
      IF( ln_cmoc    )   CALL trc_sms_cmoc   ( kt )    ! main program of CMOC   
      IF( ln_pisces  )   CALL trc_sms_pisces ( kt )    ! main program of PISCES 
      IF( ll_cfc     )   CALL trc_sms_cfc    ( kt )    ! surface fluxes of CFC
      IF( ln_c14     )   CALL trc_sms_c14    ( kt )    ! surface fluxes of C14
      IF( ln_age     )   CALL trc_sms_age    ( kt )    ! Age tracer
      IF( ln_my_trc  )   CALL trc_sms_my_trc ( kt )    ! MY_TRC  tracers

      ! O Riche Oct 25th 2022
      ! test value of jp_tot to see if jp_age is involved
      IF( lwp .AND. kt == nittrc000 ) THEN
        WRITE(numout,*) 'trc_sms: jp_age and jp_tot check'
        WRITE(numout,*) 'jp_age = ', jp_age
        WRITE(numout,*) 'jp_tot = ', jp_tot
      ENDIF
  
      IF(ln_ctl) THEN      ! print mean trends (used for debugging)
        WRITE(charout, FMT="('sms ')")
        CALL prt_ctl_trc_info( charout )
        CALL prt_ctl_trc( tab4d=trn, mask=tmask, clinfo=ctrcnm )
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('trc_sms')
      !
   END SUBROUTINE trc_sms


   INTEGER FUNCTION trc_sms_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to TOP
      ! ALLOCATE( tab(...) , STAT=trc_sms_alloc )
      trc_sms_alloc = 0      ! set to zero if no array to be allocated
      !
      IF( trc_sms_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_alloc : failed to allocate arrays' )
      !
   END FUNCTION trc_sms_alloc

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
