MODULE cmocmort
   !!======================================================================
   !!                         ***  MODULE p4zmort  ***
   !! TOP/CANBGC :   CMOC Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!          CMOC 1  !  2013-2015(O. Riche) phytoplankton mortality, based on Zahariev et al 2008
   !!          CMOC 2  !  2022-23  (O. Riche) NEMO4 phytoplankton mortality   
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_cmoc'?                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !!   cmoc_mort       :   Compute the phytoplankton mortality terms
   !!   cmoc_mort_init  :   Initialization of the mortality parameters   !!----------------------------------------------------------------------
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

   PUBLIC   cmoc_mort          ! called in trcsms_cmoc.F90
   PUBLIC   cmoc_mort_init     ! called in trcini_cmoc.F90    

   !!* Substitution
!#  include "top_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmort.F90 3295 2012-01-30 15:49:07Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

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



END MODULE cmocmort