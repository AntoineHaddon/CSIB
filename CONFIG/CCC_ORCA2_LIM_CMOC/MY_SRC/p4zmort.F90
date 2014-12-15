MODULE p4zmort
   !!======================================================================
   !!                         ***  MODULE p4zmort  ***
   !! TOP :   PISCES Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
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
   ! <CMOC OR 06/13/2014> Code trimming !  USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_mort    
   PUBLIC   p4z_mort_init    

   !! * Shared module variables
! <CMOC OR 06/13/2014> Code trimming !    REAL(wp), PUBLIC :: wchl   = 0.001_wp  !:
! <CMOC OR 06/13/2014> Code trimming !    REAL(wp), PUBLIC :: wchld  = 0.02_wp   !:
! <CMOC OR 06/13/2014> Code trimming !    REAL(wp), PUBLIC :: mprat  = 0.01_wp   !:
! <CMOC OR 06/13/2014> Code trimming !    REAL(wp), PUBLIC :: mprat2 = 0.01_wp   !:
! <CMOC OR 06/13/2014> Code trimming !    REAL(wp), PUBLIC :: mpratm = 0.01_wp   !:


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

! <CMOC OR 06/13/2014> Code trimming !       CALL p4z_nano            ! nanophytoplankton

! <CMOC OR 06/13/2014> Code trimming !       CALL p4z_diat            ! diatoms

! <CMOC OR 06/13/2014> Code trimming !    END SUBROUTINE p4z_mort


! <CMOC OR 06/13/2014> Code trimming !    SUBROUTINE p4z_nano
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph
      ! <CMOC OR 06/13/2014> Code trimming !  REAL(wp) :: zfactfe, zfactch, zprcaca, zfracal
      REAL(wp) :: ztortp , zrespp , zmortp , zstep, zfactch
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_mort')! <CMOC OR 06/13/2014> Code trimming !  'p4z_nano')
      !
      ! <CMOC OR 06/13/2014> Code trimming !  prodcal(:,:,:) = 0.  !: calcite production variable set to zero
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompaph = MAX( ( trn(ji,jj,jk,jpphy) - 1e-8 ), 0.e0 )
               zstep    = xstep ! xstep convert from d-1 to s-1 and then integrate over the time step (standard time step is 5760 s)
# if defined key_degrad
               zstep    = zstep * facvol(ji,jj,jk)
# endif
               !     Squared mortality of Phyto similar to a sedimentation term during
               !     blooms (Doney et al. 1996)
               zrespp = mpd2_cmoc * ncrr_cmoc * 1.0e3_wp * zstep * zcompaph * trn(ji,jj,jk,jpphy) ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones ! <CMOC OR 12/05/2013> fixed volume in 0.1_wp parameter (from m^-3 to L^-1) ! <CMOC OR 09/26/2013> wchl * 1.e6 * zstep * xdiss(ji,jj,jk) * zcompaph * trn(ji,jj,jk,jpphy) 

               !     Phytoplankton mortality. This mortality loss is slightly
               !     increased when nutrients are limiting phytoplankton growth
               !     as observed for instance in case of iron limitation.
               ztortp = mpd_cmoc * zstep * zcompaph ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones ! <CMOC OR 09/26/2013> mprat * xstep * trn(ji,jj,jk,jpphy) / ( xkmort + trn(ji,jj,jk,jpphy) ) * zcompaph 

               zmortp = zrespp + ztortp

               !   Update the arrays TRA which contains the biological sources and sinks

! <CMOC OR 06/13/2014> Code trimming !                 zfactfe = trn(ji,jj,jk,jpnfe)/(trn(ji,jj,jk,jpphy)+rtrn)
                 zfactch = trn(ji,jj,jk,jpnch)/(trn(ji,jj,jk,jpphy)+rtrn)

               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zmortp
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zmortp * zfactch
! <CMOC OR 06/13/2014> Code trimming !                 tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zmortp * zfactfe
! <CMOC OR 06/13/2014> Code trimming !                 zprcaca = xfracal(ji,jj,jk) * zmortp
               !
! <CMOC OR 06/13/2014> Code trimming !                 prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
! <CMOC OR 06/13/2014> Code trimming !                 zfracal = 0.5 * xfracal(ji,jj,jk)
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) ! <CMOC OR 02/19/2014> - zprcaca       ! <CMOC OR 01/07/2014> PISCES calcification
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) ! <CMOC OR 02/19/2014> - 2. * zprcaca  ! <CMOC OR 01/07/2014> PISCES calcification
! <CMOC OR 06/13/2014> Code trimming !                 tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
! <CMOC OR 06/13/2014> Code trimming !  #if defined key_kriest
! <CMOC OR 06/13/2014> Code trimming !                 tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortp
! <CMOC OR 06/13/2014> Code trimming !                 tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + ztortp * xkr_dnano + zrespp * xkr_ddiat
! <CMOC OR 06/13/2014> Code trimming !                 tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zmortp * zfactfe
! <CMOC OR 06/13/2014> Code trimming !  #else
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) +  zfracal * zmortp 
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortp! + ( 1. - zfracal ) * zmortp <CMOC OR 10/04/2013> in CMOC all decaying phytoplankton goes to detritus
! <CMOC OR 06/13/2014> Code trimming !                 tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ( 1. - zfracal ) * zmortp * zfactfe
! <CMOC OR 06/13/2014> Code trimming !                 tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zfracal * zmortp * zfactfe
! <CMOC OR 06/13/2014> Code trimming !  #endif
            END DO
         END DO
      END DO
      !
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('mort')") ! <CMOC OR 06/13/2014> Code trimming !  nano')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_mort') ! <CMOC OR 06/13/2014> Code trimming !  p4z_nano')
      !
! <CMOC OR 06/13/2014> Code trimming !    END SUBROUTINE p4z_nano
   END SUBROUTINE p4z_mort

! <CMOC OR 06/13/2014> Code trimming !    SUBROUTINE p4z_diat
!      !!---------------------------------------------------------------------
!      !!                     ***  ROUTINE p4z_diat  ***
!      !!
!      !! ** Purpose :   Compute the mortality terms for diatoms
!      !!
!      !! ** Method  : - ???
!      !!---------------------------------------------------------------------
!      INTEGER  ::  ji, jj, jk
!      REAL(wp) ::  zfactfe,zfactsi,zfactch, zcompadi
!      REAL(wp) ::  zrespp2, ztortp2, zmortp2, zstep
!      CHARACTER (len=25) :: charout
!      !!---------------------------------------------------------------------
!      !
!      IF( nn_timing == 1 )  CALL timing_start('p4z_diat')
!      !
!
!      !    Aggregation term for diatoms is increased in case of nutrient
!      !    stress as observed in reality. The stressed cells become more
!      !    sticky and coagulate to sink quickly out of the euphotic zone
!      !     ------------------------------------------------------------
!
!      DO jk = 1, jpkm1
!         DO jj = 1, jpj
!            DO ji = 1, jpi
!
!               zcompadi = MAX( ( trn(ji,jj,jk,jpdia) - 1e-8), 0. )
!
!               !    Aggregation term for diatoms is increased in case of nutrient
!               !    stress as observed in reality. The stressed cells become more
!               !    sticky and coagulate to sink quickly out of the euphotic zone
!               !     ------------------------------------------------------------
!               zstep   = xstep
!# if defined key_degrad
!               zstep = zstep * facvol(ji,jj,jk)
!# endif
!               !  Phytoplankton respiration 
!               !     ------------------------
!               zrespp2  = 1.e6 * zstep * (  wchl + wchld * ( 1.- xlimdia(ji,jj,jk) )  )    &
!                  &       * xdiss(ji,jj,jk) * zcompadi * trn(ji,jj,jk,jpdia)
!
!               !     Phytoplankton mortality. 
!               !     ------------------------
!               ztortp2  = mprat2 * zstep * trn(ji,jj,jk,jpdia)  / ( xkmort + trn(ji,jj,jk,jpdia) ) * zcompadi 
!
!               zmortp2 = zrespp2 + ztortp2
!
!               !   Update the arrays tra which contains the biological sources and sinks
!               !   ---------------------------------------------------------------------
!               zfactch = trn(ji,jj,jk,jpdch) / ( trn(ji,jj,jk,jpdia) + rtrn )
!               zfactfe = trn(ji,jj,jk,jpdfe) / ( trn(ji,jj,jk,jpdia) + rtrn )
!               zfactsi = trn(ji,jj,jk,jpdsi) / ( trn(ji,jj,jk,jpdia) + rtrn )
!
!               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zmortp2 
!               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zmortp2 * zfactch
!               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zmortp2 * zfactfe
!               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zmortp2 * zfactsi
!               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zmortp2 * zfactsi
!#if defined key_kriest
!               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) !+ zmortp2 <CMOC OR 11/13/2013> 
!               tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + ztortp2 * xkr_ddiat + zrespp2 * xkr_daggr
!               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zmortp2 * zfactfe
!#else
!               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zrespp2 + 0.5 * ztortp2
!               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) !+ 0.5 * ztortp2 <CMOC OR 11/13/2013>
!               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + 0.5 * ztortp2 * zfactfe
!               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + ( zrespp2 + 0.5 * ztortp2 ) * zfactfe
!#endif
!            END DO
!         END DO
!      END DO
!      !
!      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
!         WRITE(charout, FMT="('diat')")
!         CALL prt_ctl_trc_info(charout)
!         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
!      ENDIF
!      !
!      IF( nn_timing == 1 )  CALL timing_stop('p4z_diat')
!      !
! <CMOC OR 06/13/2014> Code trimming !     END SUBROUTINE p4z_diat

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

! <CMOC OR 06/13/2014> Code trimming !       NAMELIST/nampismort/ wchl, wchld, mprat, mprat2, mpratm

      ! <CMOC OR 03/08/2014> CMOC namelist
      NAMELIST/namcmocmor/ mpd_cmoc, mpd2_cmoc
      ! <CMOC OR 03/08/2014> CMOC namelist end 
      !!----------------------------------------------------------------------

! <CMOC OR 06/13/2014> Code trimming !       REWIND( numnatp )                     ! read numnatp
! <CMOC OR 06/13/2014> Code trimming !       READ  ( numnatp, nampismort )

      REWIND( numcmoc )                    ! <CMOC OR 03/10/2014> ! read numcmoc, cmocphy
      READ  ( numcmoc, namcmocmor )

      
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
! <CMOC OR 06/13/2014> Code trimming !          WRITE(numout,*) ' Namelist parameters for phytoplankton mortality, nampismort'
! <CMOC OR 06/13/2014> Code trimming !          WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
! <CMOC OR 06/13/2014> Code trimming !          WRITE(numout,*) '    quadratic mortality of phytoplankton      wchl      =', wchl
! <CMOC OR 06/13/2014> Code trimming !          WRITE(numout,*) '    maximum quadratic mortality of diatoms    wchld     =', wchld
! <CMOC OR 06/13/2014> Code trimming !          WRITE(numout,*) '    phytoplankton mortality rate              mprat     =', mprat
! <CMOC OR 06/13/2014> Code trimming ! WRITE(numout,*) '    Diatoms mortality rate                    mprat2    =', mprat2
! <CMOC OR 06/13/2014> Code trimming !          WRITE(numout,*) '    Phytoplankton minimum mortality rate      mpratm    =', mpratm
! <CMOC OR 06/13/2014> Code trimming !          WRITE(numout,*) ' '
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
