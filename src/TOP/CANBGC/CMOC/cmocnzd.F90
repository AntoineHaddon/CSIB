MODULE cmocnzd
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP/CANBGC :   CMOC Compute remineralization/scavenging of organic compounds
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!          CMOC 1  !  2013-2015(O. Riche) remineralization, PIC export and nitrogen fixation, based on Zahariev et al 2008
   !!          CMOC 1  !  2016-02  (N. Swart) Bugfixes and moves calcite flux to p4zsink; DNF to p4zsed.
   !!          CMOC 2  !  2022-23  (O. Riche) NEMO4 remineralization (N), zooplankton (Z) grazing, and mortality (D)   
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_cmoc'?                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !!   cmoc_rem        :  Compute remineralization/scavenging of organic compounds
   !!   cmoc_rem_init   :  Initialisation of parameters for remineralisation
   !!   cmoc_zoo        :  Compute the sources/sinks for microzooplankton
   !!   cmoc_zoo_init   :  Initialize and read the appropriate namelist
   !!   cmoc_mort       :  Compute the phytoplankton mortality terms
   !!   cmoc_mort_init  :  Initialization of the mortality parameters
   !!----------------------------------------------------------------------
   
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

   ! remineralisation to nitrate
   PUBLIC   cmoc_rem          ! called in trcsms_cmoc.F90
   PUBLIC   cmoc_rem_denit    ! called in trcsms_cmoc.F90
   PUBLIC   cmoc_rem_init     ! called in trcini_cmoc.F90    
   ! zooplankton grazing
   PUBLIC cmoc_zoo
   PUBLIC cmoc_zoo_init
   ! phytoplankton mortality
   PUBLIC   cmoc_mort          ! called in trcsms_cmoc.F90
   PUBLIC   cmoc_mort_init     ! called in trcini_cmoc.F90    
   !
   !!* Substitution
!#  include "top_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmort.F90 3295 2012-01-30 15:49:07Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
   !
   !!----------------------------------------------------------------------
   !! Remineralization into nitrate
   !!----------------------------------------------------------------------
   !
   SUBROUTINE cmoc_rem( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Calls the different subroutine to initialize and compute
      !!                the remineralization term
      !!
      !! ** Method  : - forward time integration (Euler or Leapfrog)
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk                                      ! loop indices 
      REAL(wp) :: zcompaph , ztortp , zrespp , zmortp , zfactch
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('cmoc_rem')
      !
      IF( lwp ) THEN
        WRITE(numout,*)
        WRITE(numout,*), 'cmoc_rem: compute organic remineralization'
        WRITE(numout,*), '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*)
        CALL FLUSH(numout)
      END IF
      !
      ! Initialization
      redet(:,:,:) = 0._wp
      !
      ! Remineralisation rate of detritus
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               redet(ji,jj,jk) = reref_cmoc * xstepb &
               &                 * exp ( -ed_cmoc * 1e3_wp / 8.31_wp *         &
               &                 ( 1._wp / ( tsn(ji,jj,jk,jp_tem) + 273.15_wp  &
               &                 + rtrn ) - 1._wp / ( tvm_cmoc + 273.15_wp )   &
               &                 )      ) * trn(ji,jj,jk,jqpoc) * tmask_bgc_closea(ji,jj,jk)
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

  
  SUBROUTINE cmoc_rem_denit
      !
      INTEGER  :: ji, jj, jk         
      !!!
      ! O Riche Oct 27th 2022
      ! This is a block for denitrification, using the rem. rate as a proxy
      ! The code also requires jk_eud_cmoc the z-level for the bottom of
      ! euphotic layer
      !
      IF( lwp ) THEN
        WRITE(numout,*)
        WRITE(numout,*), 'cmoc_rem_denit: compute  denitrication rate'
        WRITE(numout,*), '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*)
        CALL FLUSH(numout)
      END IF      
      !
      ! Initialization
      redettot(:,:) = 0._wp
      !
      ! Integration of remineralization below the euphotic zone (used for dentrification scaling)
      DO jk = jk_eud_cmoc+1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
                redettot(ji,jj) = redettot(ji,jj) + redet(ji,jj,jk)          &
                &                                 * e3t_n(ji,jj,jk)          &
                &                                 * tmask_bgc_closea(ji,jj,jk)
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
   !
   !!----------------------------------------------------------------------
   !! Zooplankton grazing
   !!----------------------------------------------------------------------
   !
  SUBROUTINE cmoc_zoo( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - forward time integration (Euler or Leapfrog)
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph
      REAL(wp) :: zgrapoc
      REAL(wp) :: zgrazpcmoc
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('cmoc_zoo')
      !
      IF( lwp ) THEN
        WRITE(numout,*)
        WRITE(numout,*), 'cmoc_zoo: compute zooplankton grazing'
        WRITE(numout,*), '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*)
        CALL FLUSH(numout)
      END IF
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! Conserve the PISCES code principle of a minimum phytoplankton biomass
               zcompaph  = MAX( ( trn(ji,jj,jk,jqphy) - xthreshphy ), 0.e0 )
               !
               ! Convert kp_cmoc from uM N (mmol N m^-3) to mol C L^-1 with 1e-6_wp * cnrr_cmoc
               ! lambda formula in Zahariev et al 2008
               ! grazing tendency
               zgrazpcmoc = xstepb  * rm_cmoc * zcompaph  * trn(ji,jj,jk,jqphy)  /       &
               &            ( kp_cmoc * 1e-6_wp * cnrr_cmoc * kp_cmoc * 1e-6_wp *        &
               &             cnrr_cmoc + trn(ji,jj,jk,jqphy)                             &
               &            * trn(ji,jj,jk,jqphy) + rtrn ) * trn(ji,jj,jk,jqzoo)
               ! POC tendency due to detritus fraction of grazed phytoplankton
               zgrapoc   = ( 1._wp - ga_cmoc ) * zgrazpcmoc

               ! zooplankton tendencies
               tra(ji,jj,jk,jqzoo) = tra(ji,jj,jk,jqzoo) & 
               ! grazing
               &                     +  ga_cmoc * zgrazpcmoc &
               ! linear mortality and loss to POC
               &                     - ( mzn_cmoc + mzd_cmoc ) * xstepb * trn(ji,jj,jk,jqzoo) &
               ! quadratic mortality ( convert mz2_cmoc from (molN m^-3)^-1 to (molC L^-1)^-1 )
               &                     - mz2_cmoc * ncrr_cmoc * 1e3_wp                          &
               &                      * xstepb * trn(ji,jj,jk,jqzoo) * trn(ji,jj,jk,jqzoo)

               ! contribution to phytoplankton and POC 
               tra(ji,jj,jk,jqphy) = tra(ji,jj,jk,jqphy) - zgrazpcmoc
               tra(ji,jj,jk,jqnch) = tra(ji,jj,jk,jqnch) - zgrazpcmoc * trn(ji,jj,jk,jqnch)/(trn(ji,jj,jk,jqphy)+rtrn)
               tra(ji,jj,jk,jqpoc) = tra(ji,jj,jk,jqpoc) + zgrapoc

               ! mortality contribution to nutrients, carbon and oxygen cycle
               tra(ji,jj,jk,jqno3) = tra(ji,jj,jk,jqno3) + mzn_cmoc * xstepb * trn(ji,jj,jk,jqzoo)
               tra(ji,jj,jk,jqoxy) = tra(ji,jj,jk,jqoxy) - mzn_cmoc * xstepb * trn(ji,jj,jk,jqzoo)
               tra(ji,jj,jk,jqdic) = tra(ji,jj,jk,jqdic) + mzn_cmoc * xstepb * trn(ji,jj,jk,jqzoo)
               tra(ji,jj,jk,jqtal) = tra(ji,jj,jk,jqtal) - mzn_cmoc * xstepb * trn(ji,jj,jk,jqzoo) * ncrr_cmoc               
               tra(ji,jj,jk,jqpoc) = tra(ji,jj,jk,jqpoc) + mzd_cmoc * xstepb * trn(ji,jj,jk,jqzoo) &
               &                    + ncrr_cmoc * 1e3_wp * mz2_cmoc * xstepb * trn(ji,jj,jk,jqzoo) * trn(ji,jj,jk,jqzoo)
               ! O Riche Oct 28th 2022 ! This is not yet implemented.
               ! tra(ji,jj,jk,jqdnt) = tra(ji,jj,jk,jqdnt) + mzn_cmoc * xstepb * trn(ji,jj,jk,jqzoo)
               !
               
            END DO
         END DO
      END DO
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('zoo')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask_bgc_closea, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )  CALL timing_stop('cmoc_zoo')
      !
  END SUBROUTINE cmoc_zoo
  
  
  SUBROUTINE cmoc_zoo_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE cmoc_zoo_init  ***
      !!
      !! ** Purpose :   Initialization of zooplankton parameters
      !!
      !! ** Method  :   Read the namcmoczoo namelist and check the parameters
      !!                called at the first timestep
      !!
      !! ** input   :   Namelist namcmoczoo
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios       ! Local integer
      ! <CMOC code OR 10/20/2015> CMOC namelist
      ! O Riche October 28th 2022 ! move xthreshphy from namelist_pisces* to namelist_cmoc_ref
      NAMELIST/namcmoczoo/ rm_cmoc, kp_cmoc, ga_cmoc, mzn_cmoc, mzd_cmoc, mz2_cmoc, xthreshphy      
      ! <CMOC code OR 10/20/2015> CMOC namelist end 

      REWIND( numnatp_refb )              ! Namelist namcmoczoo in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcmoczoo, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcmoczoo in reference namelist_cmoc' )
      REWIND( numnatp_cfgb )              ! Namelist namcmoczoo in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcmoczoo, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcmoczoo in configuration namelist_cmoc' )

      IF(lwm) WRITE( numonpb, namcmoczoo )     

      ! control print
      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for microzooplankton, namcmoczoo'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Maximum grazing rate                            rm_cmoc     =', rm_cmoc
         WRITE(numout,*) '    Grazing half-saturation constant                kp_cmoc     =', kp_cmoc
         WRITE(numout,*) '    Grazing efficiency                              ga_cmoc     =', ga_cmoc
         WRITE(numout,*) '    Loss to nitrogen                                mzn_cmoc    =', mzn_cmoc
         WRITE(numout,*) '    Loss to detritus                                mzd_cmoc    =', mzd_cmoc
         WRITE(numout,*) '    Quadratic mortality                             mz2_cmoc    =', mz2_cmoc
         WRITE(numout,*) '    Phyto biomass threshold on zooplankton feeding  xthreshphy  =', xthreshphy
      ENDIF
      
  END SUBROUTINE cmoc_zoo_init
   !
   !!----------------------------------------------------------------------
   !! Phytoplankton mortality
   !!----------------------------------------------------------------------
   !

  SUBROUTINE cmoc_mort( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to initialize and compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - forward time integration (Euler or Leapfrog)
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk                                      ! loop indices
      REAL(wp) :: zcompaph , ztortp , zrespp , zmortp , zfactch   ! working variables
      CHARACTER (len=25) :: charout
      !
      IF( ln_timing )  CALL timing_start('cmoc_mort')
      !
      !
      IF( lwp ) THEN
        WRITE(numout,*)
        WRITE(numout,*), 'cmoc_mort: compute phytoplankton mortality terms'
        WRITE(numout,*), '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*)
        CALL FLUSH(numout)
      END IF      ! <CMOC code OR 10/19/2015> Note: replacing zstep (and removing facvol/key_grad instances) by xstep time step in days
      ! O Riche Oct 27th 2022, xstep is xstepb in TOP/CANBGC
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Note: conserve the PISCES practice of minimum phyto biomass (see zcompaph below)
               zcompaph = MAX( ( trn(ji,jj,jk,jqphy) - 1e-8 ), 0.e0 )
               
               ! Quadratic mortality
               ! <CMOC code OR 10/19/2015> 1e.3_wp convert (umol? OR Nov 2022) L^-1 to (mol ? OR Nov 2022) m^-3 ; use xstepb the global constant to convert to d^-1
               zrespp = mpd2_cmoc * ncrr_cmoc * 1.e3_wp * xstepb * zcompaph * trn(ji,jj,jk,jqphy)

               !  Linear mortality
               ztortp = mpd_cmoc * xstepb * zcompaph

               !  Total mortality
               zmortp = zrespp + ztortp

               !   Update the arrays TRA which contains the biological sources and sinks
               !   Calculate the chlorophyll to phytoplankton ratio
               zfactch = trn(ji,jj,jk,jqnch)/(trn(ji,jj,jk,jqphy)+rtrn)

               tra(ji,jj,jk,jqphy) = tra(ji,jj,jk,jqphy) - zmortp
               tra(ji,jj,jk,jqnch) = tra(ji,jj,jk,jqnch) - zmortp * zfactch
               tra(ji,jj,jk,jqpoc) = tra(ji,jj,jk,jqpoc) + zmortp      
               
            END DO
         END DO
      END DO      
      ! print mean trends (used for debugging)
      IF(ln_ctl)   THEN
       WRITE(charout, FMT="('mort')")
       CALL prt_ctl_trc_info(charout)
       CALL prt_ctl_trc(tab4d=tra, mask=tmask_bgc_closea, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )  CALL timing_stop('cmoc_mort')
      !
  END SUBROUTINE cmoc_mort


  SUBROUTINE cmoc_mort_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the namcmocmor namelist and check the parameters
      !!                called at the first timestep
      !!
      !! ** input   :   Namelist namcmocmor
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios       ! Local integer
      ! <CMOC code OR 10/19/2015> CMOC namelist
      NAMELIST/namcmocmor/ mpd_cmoc, mpd2_cmoc
      ! <CMOC code OR 10/19/2015> CMOC namelist end 
      !!----------------------------------------------------------------------

      REWIND( numnatp_refb )              ! Namelist namcmocmor in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcmocmor, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcmocmor in reference namelist_cmoc' )
      REWIND( numnatp_cfgb )              ! Namelist namcmocmor in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcmocmor, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcmocmor in configuration namelist_cmoc' )

      IF(lwm) WRITE( numonpb, namcmocmor )


      ! control print
      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton mortality, namcmocmor'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Phytoplankton mortality to detritus       mpd_cmoc  =', mpd_cmoc
         WRITE(numout,*) '    Phytoplankton quadratic mortality         mpd2_cmoc =', mpd2_cmoc
      ENDIF

  END SUBROUTINE cmoc_mort_init


END MODULE cmocnzd