MODULE cmocrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP/CANBGC :   CMOC Compute remineralization/scavenging of organic compounds
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!          CMOC 1  !  2013-2015(O. Riche) remineralization, PIC export and nitrogen fixation, based on Zahariev et al 2008
   !!          CMOC 1  !  2016-02  (N. Swart) Bugfixes and moves calcite flux to p4zsink; DNF to p4zsed.
   !!          CMOC 2  !  2022-23  (O. Riche) NEMO4 remineralization   
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_cmoc'?                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !!   cmoc_rem       :  Compute remineralization/scavenging of organic compounds
   !!   cmoc_rem_init  :  Initialisation of parameters for remineralisation
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 

   USE sms_top_canbgc     !  TOP Source Minus Sink variables
   USE sms_cmoc           !  CMOC specific parameters declaration

   USE trc_closea_canbgc  !  tmask_bgc_closea

   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   ! timing modules
   USE in_out_manager  ! nn_timing integer
   USE timing          ! *_timing subroutines

   IMPLICIT NONE
   PRIVATE

   PUBLIC   cmoc_rem          ! called in trcsms_cmoc.F90
   PUBLIC   cmoc_rem_denit    ! called in trcsms_cmoc.F90
   PUBLIC   cmoc_rem_init     ! called in trcini_cmoc.F90    

   !!* Substitution
!#  include "top_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmort.F90 3295 2012-01-30 15:49:07Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS


  SUBROUTINE cmoc_rem( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Calls the different subroutine to initialize and compute
      !!                the remineralization term
      !!
      !! ** Method  : - forward time integration (Euler or Leapfrog)
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk                                      ! loop indices 
      REAL(wp) :: zcompaph , ztortp , zrespp , zmortp , zfactch
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('cmoc_rem')
      !
      !
      IF( lwp ) THEN
        WRITE(numout,*)
        WRITE(numout,*), 'cmoc_rem: compute organic remineralization'
        WRITE(numout,*), '          and prime denitrication rate    '
        WRITE(numout,*), '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*)
        CALL FLUSH(numout)
      END IF      
      ! Initialization of CMOC arrays
       redet   (:,:,:) = 0._wp
       redettot(:,:)   = 0._wp
      !
      ! Remineralisation rate of detritus
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               redet (ji,jj,jk) = reref_cmoc * xstepb &
               &                 * exp ( -ed_cmoc * 1e3_wp / 8.31_wp *        &
               &                ( 1._wp / ( tsn(ji,jj,jk,jp_tem) + 273.15_wp  &
               &                 + rtrn ) - 1._wp / ( tvm_cmoc + 273.15_wp )  &
               &                )      ) * trn(ji,jj,jk,jqpoc) * tmask_bgc_closea(ji,jj,jk)
            END DO
         END DO 
      END DO         
      !
      !     --------------------------------------------------------------------
      !     Update the arrays TRA which contain the biological sources and sinks
      !     --------------------------------------------------------------------
      DO jk = 1, jpkm1
         tra(:,:,jk,jqpoc) = tra(:,:,jk,jqpoc) - redet(:,:,jk) 
         tra(:,:,jk,jqno3) = tra(:,:,jk,jqno3) + redet(:,:,jk) 
         tra(:,:,jk,jqoxy) = tra(:,:,jk,jqoxy) - redet(:,:,jk) 
         tra(:,:,jk,jqdic) = tra(:,:,jk,jqdic) + redet(:,:,jk) 
         tra(:,:,jk,jqtal) = tra(:,:,jk,jqtal) - redet(:,:,jk) * ncrr_cmoc
         ! tra(:,:,jk,jqdnt) = tra(:,:,jk,jqdnt) + redet(:,:,jk) 
      END DO

      ! print mean trends (used for debugging)
      IF(ln_ctl)   THEN
       WRITE(charout, FMT="('rem')")
       CALL prt_ctl_trc_info(charout)
       CALL prt_ctl_trc(tab4d=tra, mask=tmask_bgc_closea, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )  CALL timing_stop('cmoc_rem')
      !  
  END SUBROUTINE cmoc_rem

  
  SUBROUTINE cmoc_rem_denit( redet0, redettot0 )
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in    ) ::    redet0
      REAL(wp), DIMENSION(jpi,jpj    ), INTENT(   out ) :: redettot0
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk         
      !!!
      ! O Riche Oct 27th 2022
      ! This is a block for denitrification, using the rem. rate as a proxy
      ! The code also requires jk_eud_cmoc the z-level for the bottom of
      ! euphotic layer
      !!!
      ! Integration of remineralization below the euphotic zone (used for dentrification scaling)
      DO jk = jk_eud_cmoc+1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
                redettot0(ji,jj) = redettot0(ji,jj) + redet0(ji,jj,jk)        &
                &                                   * e3t_n(ji,jj,jk)         &
                &                                   * tmask_bgc_closea(ji,jj,jk)
            END DO
          END DO 
      END DO  
      !!!
  
  END SUBROUTINE cmoc_rem_denit

  
  SUBROUTINE cmoc_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the namcmocpoc namelist and check the parameters
      !!                called at the first timestep
      !!
      !! ** input   :   Namelist namcmocpoc
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios       ! Local integer
      ! <CMOC code OR 10/15/2015> CMOC namelist
      NAMELIST/namcmocpoc/ ed_cmoc, reref_cmoc
      NAMELIST/namcmocdeu/ jk_eud_cmoc, nk_bal_cmoc

      REWIND( numnatp_refb )              ! Namelist namcmocpoc in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcmocpoc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcmocpoc in reference namelist_cmoc' )
      REWIND( numnatp_cfgb )              ! Namelist namcmocpoc in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcmocpoc, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcmocpoc in configuration namelist_cmoc' )

      IF(lwm) WRITE( numonpb, namcmocpoc )

      REWIND( numnatp_refb )              ! Namelist namcmocdeu in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcmocdeu, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcmocdeu in reference namelist_cmoc' )
      REWIND( numnatp_cfgb )              ! Namelist namcmocdeu in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcmocdeu, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcmocdeu in configuration namelist_cmoc' )

      IF(lwm) WRITE( numonpb, namcmocdeu )


      ! control print
      IF(lwp) THEN
         WRITE(numout,*) ' Namelist parameters for remineralization, namcmocpoc'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Remineralisation rate of POC              reref_cmoc=', reref_cmoc
         WRITE(numout,*) '    Activation energy for remineralization    ed_cmoc   =', ed_cmoc   
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for remineralization, namcmocdeu'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Scale of euphotic zone ~100 m with jk = 25              jk_eud_cmoc=', jk_eud_cmoc
         WRITE(numout,*) '    Open ocean criterion:  '
         WRITE(numout,*) '    Number of vertical layers required to calculate balance'
         WRITE(numout,*) '    for the calcite fluxes and DNF/denitrification ~ 150 m '
         WRITE(numout,*) '    with nk_bal_cmoc = 28    nk_bal_cmoc   =', nk_bal_cmoc   
         WRITE(numout,*) ' '
      ENDIF
      ! 
      ! Allocate arrays
      ALLOCATE(    redet( jpi, jpj, jpk ) )   ! Remineralization rate
      ALLOCATE( redettot( jpi, jpj ) )        ! Denitrification rate
      !
  END SUBROUTINE cmoc_rem_init

END MODULE cmocrem