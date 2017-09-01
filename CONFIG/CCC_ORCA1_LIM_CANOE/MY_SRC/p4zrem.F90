MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   PISCES Compute remineralization/scavenging of organic compounds
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_top'       and                                      TOP models
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/scavenging of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!   p4z_rem_alloc :  Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zopt          !  optical model
   USE p4zche          !  chemical model
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zmeso         !  Sources and sinks of mesozooplankton
   USE p4zint          !  interpolation and computation of various fields
   USE prtctl_trc      !  print control for debugging
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p4z_rem_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_rem_alloc

   !! * Shared module variables
   REAL(wp), PUBLIC ::  xremik    = 0.25_wp    !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::  xremip    = 0.025_wp   !: remineralisation rate of DOC
   REAL(wp), PUBLIC ::  nitrif    = 0.05_wp    !: NH4 nitrification rate 
   REAL(wp), PUBLIC ::  xlam1     = 0.0001_wp  !: scavenging rate of iron (low concentrations)
   REAL(wp), PUBLIC ::  xlam2     = 2.5_wp     !: scavenging rate of iron (high concentrations)
   REAL(wp), PUBLIC ::  ligand    = 6.0E+2_wp  !: ligand concentration
   REAL(wp), PUBLIC ::  pocfctr   = 0.65574_wp !: multiplier for POC-dependent scavenging
   REAL(wp), PUBLIC ::  o2thresh  = 6._wp      !: O2 threshold for denitrification
   REAL(wp), PUBLIC ::  nh4frx    = 02.5_wp    !: annamox fraction of denitrification
   REAL(wp), PUBLIC ::  oxymin    = 1._wp      !: half saturation constant for anoxia 
   REAL(wp), PUBLIC ::  nyld      = 0.8_wp     !: denitrification stoichiometric coefficient

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zrem.F90 3558 2012-11-14 19:15:05Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_rem( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zremip, zremik, Tf
      REAL(wp) ::   zkeq, zfeequi
      REAL(wp) ::   zsatur, zsatur2, zdep, zfactdep
      REAL(wp) ::   zorem, zorem2, zofer, zofer2
      REAL(wp) ::   zscave, zscavex, fexs, zcoag
      REAL(wp) ::   zlamfac, zonitr, zstep
      REAL(wp) ::   zrfact2
      CHARACTER (len=25) :: charout
      !REAL(wp), POINTER, DIMENSION(:,:,:) :: zolimi, zolimi2, zwork
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_rem')
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! nitrification rate depends on O2 concentration
               nitrfac(ji,jj,jk) = MAX(  0.e0, 0.4 * ( 6.  - trn(ji,jj,jk,jpoxy) )    &
                  &                                / ( oxymin + trn(ji,jj,jk,jpoxy) )  )
               nitrfac(ji,jj,jk) = MIN( 1., nitrfac(ji,jj,jk) )
            END DO
         END DO
      END DO

      nh4ox(:,:,:)=0.
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
               !    NH4 nitrification to NO3. Ceased for oxygen concentrations
               !    below 2 umol/L. Inhibited at strong light 
               !    ----------------------------------------------------------
               zonitr  =nitrif * zstep * trn(ji,jj,jk,jpnh4) / ( 1.+ emoy(ji,jj,jk) ) * ( 1.- nitrfac(ji,jj,jk) ) 
               !denitnh4(ji,jj,jk) = nitrif * zstep * trn(ji,jj,jk,jpnh4) * nitrfac(ji,jj,jk) 
               !   Update of the tracers trends
               !   ----------------------------
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zonitr
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + zonitr
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - 2. * zonitr
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2.e-6 * zonitr
               nh4ox(ji,jj,jk) = zonitr
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem2')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      denitr(:,:,:)=0.
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               Tf = tgfuncr(ji,jj,jk)
               zorem  = xremik * xstep * Tf * trn(ji,jj,jk,jppoc)
               zofer  = zorem * rr_fe2c
               zorem2 = xremik * xstep * Tf * trn(ji,jj,jk,jpgoc)
               zofer2 = zorem2 * rr_fe2c

! denitrification is assumed to remove NO3 as a fraction of remineralization increasing linearly from 0 to 1 with declining [O2] for [O2]<10 uM
! NO3 fraction is then divided 0.75/0.25 between NO3 and NH4 (50% classical denitrification and 50% annamox)
               zonitr=1.-MIN(trn(ji,jj,jk,jpoxy),o2thresh)/o2thresh
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + (zorem + zorem2)*rr_n2c - (zorem + zorem2)*nyld*zonitr*0.5*nh4frx
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - (zorem + zorem2)*nyld*zonitr*(1.-0.5*nh4frx)
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - (zorem + zorem2)*(1.-zonitr)
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zofer + zofer2
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + (zorem + zorem2)*1.e-6
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zorem
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) - zorem2
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 1.e-6 * (zorem + zorem2)*rr_n2c                         ! 1 mol of alkalinity per mol of N
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 1.e-6 * (zorem + zorem2)*nyld*zonitr*(1.-nh4frx)        ! +1 mol if denitrification, 0 if annamox
               denitr(ji,jj,jk) = (zorem + zorem2)*zonitr*nyld

            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem3')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem4')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji = 1, jpi
               zstep   = xstep
! irreversible scavenging as in Christian et al 2002
               zcoag = MIN((trn(ji,jj,jk,jppoc)+trn(ji,jj,jk,jpgoc))*pocfctr,1.)
               fexs = MAX(trn(ji,jj,jk,jpfer)-ligand,0.)
               zscave = xlam1 * xstep * (trn(ji,jj,jk,jpfer)-fexs) * zcoag
               zscavex = xlam2 * xstep * fexs
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - (zscave+zscavex)
            END DO
         END DO
      END DO
      !

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem5')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      !     Update the arrays TRA which contain the biological sources and sinks
      !     --------------------------------------------------------------------

      IF( ln_diatrc ) THEN
         zrfact2 = 1.e-3 * rfact2r  ! conversion from umol/L/timestep into mol/m3/s
         denitr(:,:,:) = denitr(:,:,:) * zrfact2
         nh4ox(:,:,:) = nh4ox(:,:,:) * zrfact2
         IF( jnt == nrdttrc ) THEN
          CALL iom_put( "Denitr"   , denitr(:,:,:) * tmask(:,:,:) )  ! rate of denitrification
          CALL iom_put( "Nitrif"   , nh4ox(:,:,:) * tmask(:,:,:) )  ! rate of nitrification
         ENDIF
      ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
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
      NAMELIST/nampisrem/ xremik, xremip, nitrif, xlam1, xlam2, ligand, pocfctr, o2thresh, &
                        & nh4frx, oxymin

      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampisrem )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for remineralization, nampisrem'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    remineralisation rate of POC(1)           xremik    =', xremik
         WRITE(numout,*) '    remineralisation rate of POC(2)           xremip    =', xremip
         WRITE(numout,*) '    Nitrification maximum rate                nitrif    =', nitrif
         WRITE(numout,*) '    Fe scavenging rate (low concentrations)   xlam1     =', xlam1
         WRITE(numout,*) '    Fe scavenging rate (high concentrations)  xlam2     =', xlam2
         WRITE(numout,*) '    Fe-binding ligand concentration           ligand    =', ligand
         WRITE(numout,*) '    POC-dependence of Fe scavenging           pocfctr   =', pocfctr
         WRITE(numout,*) '    O2 threshold for denitrification          o2thresh  =', o2thresh
         WRITE(numout,*) '    Annamox fraction of denitrification       nh4frx    =', nh4frx
         WRITE(numout,*) '    O2 dependence of nitrification            oxymin    =', oxymin
      ENDIF
      !
      nitrfac (:,:,:) = 0._wp
      denitr  (:,:,:) = 0._wp
      nh4ox   (:,:,:) = 0._wp
      !denitnh4(:,:,:) = 0._wp
      !
   END SUBROUTINE p4z_rem_init


   INTEGER FUNCTION p4z_rem_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem_alloc  ***
      !!----------------------------------------------------------------------
!      ALLOCATE( denitr(jpi,jpj,jpk), STAT=p4z_rem_alloc )
      !
      IF( p4z_rem_alloc /= 0 )   CALL ctl_warn('p4z_rem_alloc: failed to allocate arrays')
      !
   END FUNCTION p4z_rem_alloc

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
