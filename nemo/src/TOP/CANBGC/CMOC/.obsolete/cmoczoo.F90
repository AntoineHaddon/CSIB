MODULE cmoczoo
   !!======================================================================
   !!                         ***  MODULE p4zmicro  ***
   !! TOP :   PISCES Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!           CMOC1  !  2013-15  (O. Riche) Zooplankton: phyto grazing and zoop. mortality
   !!          CMOC 2  !  2022-23  (O. Riche) NEMO4 zooplankton grazing/mortality     
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_cmoc'?                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !!   cmoc_zoo       :   Compute the sources/sinks for microzooplankton
   !!   cmoc_zoo_init  :   Initialize and read the appropriate namelist
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

   PUBLIC cmoc_zoo
   PUBLIC cmoc_zoo_init

   !!* Substitution
!#  include "top_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmort.F90 3295 2012-01-30 15:49:07Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS


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

END MODULE cmoczoo