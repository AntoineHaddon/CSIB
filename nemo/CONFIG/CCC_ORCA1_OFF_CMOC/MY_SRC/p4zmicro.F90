MODULE p4zmicro
   !!======================================================================
   !!                         ***  MODULE p4zmicro  ***
   !! TOP :   PISCES Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!           CMOC1  !  2013-15  (O. Riche) Zooplankton: phyto grazing and zoop. mortality
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_micro       :   Compute the sources/sinks for microzooplankton
   !!   p4z_micro_init  :   Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zprod         !  production
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_micro         ! called in p4zbio.F90
   PUBLIC   p4z_micro_init    ! called in trcsms_pisces.F90
   REAL(wp), PUBLIC ::  xthreshphy = 2E-7_wp    !: nanophyto threshold for microzooplankton 

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmicro.F90 3295 2012-01-30 15:49:07Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_micro( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph
      REAL(wp) :: zgrapoc
      REAL(wp) :: zgrazpcmoc
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_micro')
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! Conserve the PISCES code principle of a minimum phytoplankton biomass
               zcompaph  = MAX( ( trn(ji,jj,jk,jpphy) - xthreshphy ), 0.e0 )
               !
               ! Convert kp_cmoc from uM N (mmol N m^-3) to mol C L^-1 with 1e-6_wp * cnrr_cmoc
               ! lambda formula in Zahariev et al 2008
               ! grazing tendency
               zgrazpcmoc = xstep   * rm_cmoc * zcompaph  * trn(ji,jj,jk,jpphy)  /       &
               &            ( kp_cmoc * 1e-6_wp * cnrr_cmoc * kp_cmoc * 1e-6_wp *        &
               &             cnrr_cmoc + trn(ji,jj,jk,jpphy)                             &
               &            * trn(ji,jj,jk,jpphy) + rtrn ) * trn(ji,jj,jk,jpzoo)
               ! POC tendency due to detritus fraction of grazed phytoplankton
               zgrapoc   = ( 1._wp - ga_cmoc ) * zgrazpcmoc

               ! zooplankton tendencies
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) & 
               ! grazing
               &                     +  ga_cmoc * zgrazpcmoc &
               ! linear mortality and loss to POC
               &                     - ( mzn_cmoc + mzd_cmoc ) * xstep * trn(ji,jj,jk,jpzoo) &
               ! quadratic mortality ( convert mz2_cmoc from (molN m^-3)^-1 to (molC L^-1)^-1 )
               &                     - mz2_cmoc * ncrr_cmoc * 1e3_wp                         &
               &                      * xstep * trn(ji,jj,jk,jpzoo) * trn(ji,jj,jk,jpzoo)

               ! contribution to phytoplankton and POC 
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgrazpcmoc
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgrazpcmoc * trn(ji,jj,jk,jpnch)/(trn(ji,jj,jk,jpphy)+rtrn)
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zgrapoc

               ! mortality contribution to nutrients, carbon and oxygen cycle
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + mzn_cmoc * xstep * trn(ji,jj,jk,jpzoo)
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - mzn_cmoc * xstep * trn(ji,jj,jk,jpzoo)
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + mzn_cmoc * xstep * trn(ji,jj,jk,jpzoo)
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - mzn_cmoc * xstep * trn(ji,jj,jk,jpzoo) * ncrr_cmoc               
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + mzd_cmoc * xstep * trn(ji,jj,jk,jpzoo) &
               &                    + ncrr_cmoc * 1e3_wp * mz2_cmoc * xstep * trn(ji,jj,jk,jpzoo) * trn(ji,jj,jk,jpzoo)
               tra(ji,jj,jk,jpdnt) = tra(ji,jj,jk,jpdnt) + mzn_cmoc * xstep * trn(ji,jj,jk,jpzoo)
               !

            END DO
         END DO
      END DO
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_micro')
      !
   END SUBROUTINE p4z_micro


   SUBROUTINE p4z_micro_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the nampiszoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiszoo
      !!
      !!----------------------------------------------------------------------

      ! <CMOC code OR 10/20/2015> CMOC namelist
      NAMELIST/nampiszoo/ xthreshphy
      NAMELIST/namcmoczoo/ rm_cmoc, kp_cmoc, ga_cmoc, mzn_cmoc, mzd_cmoc, mz2_cmoc
      ! <CMOC code OR 10/20/2015> CMOC namelist end 

      REWIND( numnatp )
      READ  ( numnatp, nampiszoo )
      REWIND( numcmoc )
      READ  ( numcmoc, namcmoczoo )

      ! control print
      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for microzooplankton, nampiszoo'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    nanophyto feeding threshold for microzoo        xthreshphy  =', xthreshphy
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for microzooplankton, namcmoczoo'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Maximum grazing rate                            rm_cmoc     =', rm_cmoc
         WRITE(numout,*) '    Grazing half-saturation constant                kp_cmoc     =', kp_cmoc
         WRITE(numout,*) '    Grazing efficiency                              ga_cmoc     =', ga_cmoc
         WRITE(numout,*) '    Loss to nitrogen                                mzn_cmoc    =', mzn_cmoc
         WRITE(numout,*) '    Loss to detritus                                mzd_cmoc    =', mzd_cmoc
         WRITE(numout,*) '    Quadratic mortality                             mz2_cmoc    =', mz2_cmoc
      ENDIF

   END SUBROUTINE p4z_micro_init


#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_micro                    ! Empty routine
   END SUBROUTINE p4z_micro
#endif 

   !!======================================================================
END MODULE  p4zmicro
