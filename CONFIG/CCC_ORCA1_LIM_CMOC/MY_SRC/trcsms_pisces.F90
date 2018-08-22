MODULE trcsms_pisces
   !!======================================================================
   !!                         ***  MODULE trcsms_pisces  ***
   !! TOP :   PISCES Source Minus Sink manager
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!           CMOC1  !  2013-15 (O. Riche) code editing for consistency with CMOC1 adaptation
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   trcsms_pisces        :  Time loop of passive tracers sms
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zbio          !  Biological model
   USE p4zche          !  Chemical model
   USE p4zdcy          !  Radioactive decay
   USE p4zflx          !  Gas exchange
   USE p4zint          !  Interpolation and computation of accessory fields 
   USE p4zsed          !  Sedimentation
   USE trdmod_oce      !  Ocean trends variables
   USE trdmod_trc      !  TOP trends variables
   USE sedmodel        !  Sediment model
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_pisces    ! called in trcsms.F90

   LOGICAL ::  ln_check_mass = .false.       !: Flag to check mass conservation 

   INTEGER ::  numno3  !: logical unit for NO3 budget
   INTEGER ::  numalk  !: logical unit for talk budget
   INTEGER ::  numsil  !: logical unit for Si budget

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_pisces.F90 3320 2012-03-05 16:37:52Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_sms_pisces( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sms_pisces  ***
      !!
      !! ** Purpose :   Managment of the call to Biological sources and sinks 
      !!              routines of PISCES bio-model
      !!
      !! ** Method  : - at each new day ...
      !!              - several calls of bio and sed ???
      !!              - ...
      !!---------------------------------------------------------------------
      !
       INTEGER, INTENT( in ) ::   kt               ! ocean time-step index      
       !!
       INTEGER ::   jnt, jn, jl
       CHARACTER (len=25) :: charout
       REAL(wp), POINTER, DIMENSION(:,:,:,:)  :: ztrdpis
       !!---------------------------------------------------------------------
       !
       IF( nn_timing == 1 )  CALL timing_start('trc_sms_pisces')
       !
       IF( ln_pisdmp .AND. MOD( kt - nn_dttrc, nn_pisdmp ) == 0 )   CALL trc_sms_pisces_dmp( kt )  ! Relaxation of some tracers
       IF( l_trdtrc )  THEN
          CALL wrk_alloc( jpi, jpj, jpk, jp_pisces, ztrdpis )
          DO jn = 1, jp_pisces
             jl = jn + jp_pcs0 - 1
             ztrdpis(:,:,:,jn) = trn(:,:,:,jl)
          ENDDO
       ENDIF
 
       IF( ndayflxtr /= nday_year ) THEN      ! New days
          !
          ndayflxtr = nday_year
 
          IF(lwp) write(numout,*)
          IF(lwp) write(numout,*) ' New chemical constants and various rates for biogeochemistry at new day : ', nday_year
          IF(lwp) write(numout,*) '~~~~~~'
 
         CALL p4z_che              ! computation of chemical constants
         !
      ENDIF

      DO jnt = 1, nrdttrc          ! Potential time splitting if requested
         !
         CALL p4z_bio (kt, jnt)    ! Compute soft tissue production (POC)
         CALL p4z_sed (kt, jnt)    ! compute soft tissue remineralisation
         !
         DO jn = jp_pcs0, jp_pcs1
            trb(:,:,:,jn) = trn(:,:,:,jn)
         ENDDO
         !
      END DO

      IF( l_trdtrc )  THEN
         DO jn = 1, jp_pisces
            jl = jn + jp_pcs0 - 1
            ztrdpis(:,:,:,jn) = ( ztrdpis(:,:,:,jn) - trn(:,:,:,jl) ) * rfact2r
         ENDDO
      ENDIF

      CALL p4z_flx( kt )             ! Compute surface fluxes
      CALL p4z_dcy( kt )

      DO jn = jp_pcs0, jp_pcs1
        CALL lbc_lnk( trn(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk( trb(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk( tra(:,:,:,jn), 'T', 1. )
      END DO

      IF( l_trdtrc ) THEN
         DO jn = 1, jp_pisces
            jl = jn + jp_pcs0 - 1
             ztrdpis(:,:,:,jn) = ztrdpis(:,:,:,jn) + tra(:,:,:,jl)
             CALL trd_mod_trc( ztrdpis(:,:,:,jn), jn, jptra_trd_sms, kt )   ! save trends
          END DO
          CALL wrk_dealloc( jpi, jpj, jpk, jp_pisces, ztrdpis ) 
      END IF

      IF( lk_sed ) THEN 
         !
         CALL sed_model( kt )     !  Main program of Sediment model
         !
         DO jn = jp_pcs0, jp_pcs1
           CALL lbc_lnk( trn(:,:,:,jn), 'T', 1. )
         END DO
         !
      ENDIF
      
      ! Update average of sea surface salinity
      CALL update_salt_avg_2d(trn(:,:,1,jp_sal), rdt, salt_avg, salt_dtsum)  
      ! If this is the last timestep of the year, reset all the annual
      ! salt-related variables
      IF ( (salt_dtsum + rdt) == nyear_len(1)*86400.) THEN
         sss_glob_avg = glob_avg_area_wt( salt_avg ) 
         salt_avg(:,:) = 0.
         salt_dtsum = 0.
      ENDIF
      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_pisces')
      !
   END SUBROUTINE trc_sms_pisces

   SUBROUTINE trc_sms_pisces_dmp( kt )
      !!----------------------------------------------------------------------
      !!                    ***  trc_sms_pisces_dmp  ***
      !!
      !! ** purpose  : Relaxation of some tracers
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT( in )  ::     kt ! time step
      !
      REAL(wp) ::  alkmean = 2426.     ! mean value of alkalinity ( Glodap ; for Goyet 2391. )
      REAL(wp) ::  po4mean = 2.165     ! mean value of phosphates
      REAL(wp) ::  no3mean = 30.90     ! mean value of nitrate
      REAL(wp) ::  silmean = 91.51     ! mean value of silicate
      !
      REAL(wp) :: zarea, zalksum, zpo4sum, zno3sum, zsilsum
      !!---------------------------------------------------------------------


      IF(lwp)  WRITE(numout,*)
      IF(lwp)  WRITE(numout,*) ' trc_sms_pisces_dmp : Relaxation of nutrients at time-step kt = ', kt
      IF(lwp)  WRITE(numout,*)

      IF( cp_cfg == "orca" .AND. .NOT. lk_c1d ) THEN      ! ORCA configuration (not 1D) !
         !                                                    ! --------------------------- !
         ! set total alkalinity, phosphate, nitrate & silicate
         zarea          = 1._wp / glob_sum( cvol(:,:,:) ) * 1e6              

         zalksum = glob_sum( trn(:,:,:,jptal) * cvol(:,:,:)  ) * zarea
         zno3sum = glob_sum( trn(:,:,:,jpno3) * cvol(:,:,:)  ) * zarea / 7.6
 
         IF(lwp) WRITE(numout,*) '       TALK mean : ', zalksum
         trn(:,:,:,jptal) = trn(:,:,:,jptal) * alkmean / zalksum

         IF(lwp) WRITE(numout,*) '       NO3  mean : ', zno3sum
         trn(:,:,:,jpno3) = trn(:,:,:,jpno3) * no3mean / zno3sum
         !
      ENDIF

   END SUBROUTINE trc_sms_pisces_dmp

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE trc_sms_pisces( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_pisces: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_pisces
#endif 

   !!======================================================================
END MODULE trcsms_pisces 
