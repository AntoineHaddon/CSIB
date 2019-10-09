MODULE p4zmort
   !!======================================================================
   !!                         ***  MODULE p4zmort  ***
   !! TOP :   PISCES Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!          CMOC 1  !  2013-2015(O. Riche) phytoplankton mortality, based on Zahariev et al 2008
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_mort       :   Compute the mortality terms for phytoplankton
   !!   p4z_mort_init  :   Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_mort    
   PUBLIC   p4z_mort_init    

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmort.F90 3295 2012-01-30 15:49:07Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_mort( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to initialize and compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      !!---------------------------------------------------------------------

      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph , ztortp , zrespp , zmortp , zfactch
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_mort')
      !
      
      ! <CMOC code OR 10/19/2015> Note: replacing zstep (and removing facvol/key_grad instances) by xstep time step in days
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Note: conserve the PISCES practice of minimum phyto biomass (see zcompaph below)
               zcompaph = MAX( ( trn(ji,jj,jk,jpphy) - 1e-8 ), 0.e0 )
               
               ! Quadratic mortality
               ! <CMOC code OR 10/19/2015> 1e+3_wp convert L^-1 to m^-3 (); use xstep the global constant to convert to d^-1
               zrespp = mpd2_cmoc * ncrr_cmoc * 1.0e3_wp * xstep * zcompaph * trn(ji,jj,jk,jpphy)

               !  Linear mortality
               ztortp = mpd_cmoc * xstep * zcompaph

               !  Total mortality
               zmortp = zrespp + ztortp

               !   Update the arrays TRA which contains the biological sources and sinks
               !   Calculate the chlorophyll to phytoplankton ratio
               zfactch = trn(ji,jj,jk,jpnch)/(trn(ji,jj,jk,jpphy)+rtrn)

               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zmortp
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zmortp * zfactch
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortp
               
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
      IF( nn_timing == 1 )  CALL timing_stop('p4z_mort')
      
   END SUBROUTINE p4z_mort


   SUBROUTINE p4z_mort_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the nampismort namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampismort
      !!
      !!----------------------------------------------------------------------

      ! <CMOC code OR 10/19/2015> CMOC namelist
      NAMELIST/namcmocmor/ mpd_cmoc, mpd2_cmoc
      ! <CMOC code OR 10/19/2015> CMOC namelist end 
      !!----------------------------------------------------------------------

      REWIND( numcmoc )
      READ  ( numcmoc, namcmocmor )

      ! control print
      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton mortality, namcmocmor'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Phytoplankton mortality to detritus       mpd_cmoc  =', mpd_cmoc
         WRITE(numout,*) '    Phytoplankton quadratic mortality         mpd2_cmoc =', mpd2_cmoc
      ENDIF

   END SUBROUTINE p4z_mort_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_mort                    ! Empty routine
   END SUBROUTINE p4z_mort
#endif 

   !!======================================================================
END MODULE  p4zmort
