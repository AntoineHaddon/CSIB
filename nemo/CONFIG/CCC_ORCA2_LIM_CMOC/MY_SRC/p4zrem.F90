MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   PISCES Compute remineralization/scavenging of organic compounds
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!          CMOC 1  !  2013-2015(O. Riche) remineralization, PIC export and nitrogen fixation, based on Zahariev et al 2008
   !!          CMOC 1  !  2016-02  (N. Swart) Bugfixes and moves calcite flux to p4zsink; DNF to p4zsed.
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_top'       and                                      TOP models
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/scavenging of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zopt          !  optical model
   USE prtctl_trc      !  print control for debugging
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
   USE iom             !
   USE fldread         !
   USE p4zsink         !
   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p4z_rem_init    ! called in trcsms_pisces.F90

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zrem.F90 3558 2012-11-14 19:15:05Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_rem( kt, jnt ) ! <CMOC code OR 10/15/2015> jnt is to be i
                                 ! known to prevent iom_put('Nfix') to cause 
                                 ! an error (see below) 
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization of detritus
      !!
      !! ** Method  : Temperature dependent remineralization based on
      !!              Zahariev et al. (2008)
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, jnt ! <CMOC code OR 10/15/2015> add the time 
                                       ! split index jnt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_rem')
      !

      ! Initialization of CMOC arrays
       redet   (:,:,:) = 0._wp
       redettot(:,:)   = 0._wp

      ! Remineralisation rate of detritus
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               redet (ji,jj,jk) = reref_cmoc * xstep &
               &                 * exp ( -ed_cmoc * 1e3_wp / 8.31_wp *        &
               &                ( 1._wp / ( tsn(ji,jj,jk,jp_tem) + 273.15_wp  &
               &                 + rtrn ) - 1._wp / ( tvm_cmoc + 273.15_wp )  &
               &                )      ) * tmask(ji,jj,jk)
            END DO
          END DO 
      END DO      

      ! Integration of remineralization below the euphotic zone (used for dentrification scaling)
      DO jk = jk_eud_cmoc+1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
                redettot(ji,jj) = redettot(ji,jj) + redet(ji,jj,jk)        &
                &                                   *    trn(ji,jj,jk,jppoc)  &
                &                                   * fse3t(ji,jj,jk)         &
                &                                   * tmask(ji,jj,jk)
            END DO
          END DO 
      END DO      

      !     --------------------------------------------------------------------
      !     Update the arrays TRA which contain the biological sources and sinks
      !     --------------------------------------------------------------------
      DO jk = 1, jpkm1
         tra(:,:,jk,jppoc) = tra(:,:,jk,jppoc) - redet(:,:,jk) * trn(:,:,jk,jppoc) 
         tra(:,:,jk,jpno3) = tra(:,:,jk,jpno3) + redet(:,:,jk) * trn(:,:,jk,jppoc) 
         tra(:,:,jk,jpoxy) = tra(:,:,jk,jpoxy) - redet(:,:,jk) * trn(:,:,jk,jppoc)
         tra(:,:,jk,jpdic) = tra(:,:,jk,jpdic) + redet(:,:,jk) * trn(:,:,jk,jppoc) 
         tra(:,:,jk,jptal) = tra(:,:,jk,jptal) - redet(:,:,jk) * trn(:,:,jk,jppoc) * ncrr_cmoc
      END DO


      ! print mean trends (used for debugging)
      IF(ln_ctl)   THEN
         WRITE(charout, FMT="('rem6')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_rem')
      !
   END SUBROUTINE p4z_rem


   SUBROUTINE p4z_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------

      ! <CMOC code OR 10/15/2015> CMOC namelist
      NAMELIST/namcmocpoc/ ed_cmoc, reref_cmoc
      NAMELIST/namcmocdeu/ jk_eud_cmoc, nk_bal_cmoc

      REWIND( numcmoc )            
      READ  ( numcmoc, namcmocpoc )

      ! <CMOC code OR 10/15/2015> CMOC namelist end 
      !!----------------------------------------------------------------------
      
      ! control print
      IF(lwp) THEN
         WRITE(numout,*) ' Namelist parameters for remineralization, namcmocpoc'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Remineralisation rate of POC              reref_cmoc=', reref_cmoc
         WRITE(numout,*) '    Activation energy for remineralization    ed_cmoc   =', ed_cmoc   
         WRITE(numout,*) ' '
      ENDIF
      !
   END SUBROUTINE p4z_rem_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_rem                    ! Empty routine
   END SUBROUTINE p4z_rem
#endif 

   !!======================================================================
END MODULE p4zrem
