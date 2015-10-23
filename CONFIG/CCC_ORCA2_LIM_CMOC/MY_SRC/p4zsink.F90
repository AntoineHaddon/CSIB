MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  PISCES  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!           CMOC1  !  2013-15  (O. Riche) POC sinking
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_init  :  Unitialisation of sinking speed parameters
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_sink_alloc

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking  !: POC sinking fluxes

   INTEGER  :: iksed  = 10

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zsink.F90 3685 2012-11-27 15:39:02Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_sink ( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, jnt
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zfact, zwsmax, zmax, zstep
      REAL(wp) ::   zrfact2
      INTEGER  ::   ik1
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink')
      !
      !    Sinking speeds of detritus is increased with depth as shown
      !    by data and from the coagulation theory
      !    -----------------------------------------------------------

      ! limit the values of the sinking speeds to avoid numerical instabilities
      wsbio3(:,:,:) = wsbio
      !

      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zwsmax = 0.8 * fse3t(ji,jj,jk) / xstep
               wsbio3(ji,jj,jk) = MIN( wsbio3(ji,jj,jk), zwsmax )
            END DO
         END DO
      END DO

      !  Initialization to zero all the sinking arrays 
      !   -----------------------------------------

      sinking (:,:,:) = 0.e0

      !   Compute the sedimentation term using p4zsink2 for all the sinking particles
      !   -----------------------------------------------------

      CALL p4z_sink2( wsbio3, sinking , jppoc )

      IF( ln_diatrc ) THEN
         zrfact2 = 1.e3 * rfact2r
         ik1  = iksed + 1
         IF( lk_iomput ) THEN
           IF( jnt == nrdttrc ) THEN
              CALL iom_put( "EPC100"  ,   sinking(:,:,ik1)                       * zrfact2 * tmask(:,:,1) )
              CALL iom_put( "EPCALC100",  wsbio3(:,:,11) / rday * trn(:,:,11,jppoc) * 1e3_wp * xrcico(:,:) * tmask(:,:,1) ) ! <CMOC code OR 10/22/2015> PIC diagnostics
              
           ENDIF
         ELSE
           trc2d(:,:,jp_pcs0_2d + 4) = sinking (:,:,ik1) * zrfact2 * tmask(:,:,1)

         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink

   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!----------------------------------------------------------------------
 
      ! <CMOC code OR 10/22/2015> CMOC namelist
      NAMELIST/namcmocrr/  cnrr_cmoc, ncrr_cmoc
      !
      !!----------------------------------------------------------------------

      REWIND( numcmoc )
      READ  ( numcmoc, namcmocrr  )
      ! <CMOC code OR 10/22/2015> CMOC namelist end
      
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters, Redfield ratios, namcmocrr'    
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    C:N ratio                                 cnrr_cmoc =', cnrr_cmoc
         WRITE(numout,*) '    N:C ratio                                 ncrr_cmoc =', ncrr_cmoc
      ENDIF

   END SUBROUTINE p4z_sink_init


   SUBROUTINE p4z_sink2( pwsink, psinkflx, jp_tra )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink2  ***
      !!
      !! ** Purpose :   Compute the sedimentation terms for the various sinking
      !!     particles. The scheme used to compute the trends is based
      !!     on MUSCL.
      !!
      !! ** Method  : - this ROUTINE compute not exactly the advection but the
      !!      transport term, i.e.  div(u*tra).
      !!---------------------------------------------------------------------
      !
      INTEGER , INTENT(in   )                         ::   jp_tra    ! tracer index index
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj,jpk) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   psinkflx  ! sinking fluxe
      !!
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zigma,zew,zign, zflx, zstep
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztraz, zakz, zwsink2, ztrb 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink2')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, jpk, ztraz, zakz, zwsink2, ztrb )

      zstep = rfact2 / 2.

      ztraz(:,:,:) = 0.e0
      zakz (:,:,:) = 0.e0
      ztrb (:,:,:) = trn(:,:,:,jp_tra)

      DO jk = 1, jpkm1
         zwsink2(:,:,jk+1) = -pwsink(:,:,jk) / rday * tmask(:,:,jk+1) 
      END DO
      zwsink2(:,:,1) = 0.e0


      ! Vertical advective flux
      DO jn = 1, 2
         !  first guess of the slopes interior values
         DO jk = 2, jpkm1
            ztraz(:,:,jk) = ( trn(:,:,jk-1,jp_tra) - trn(:,:,jk,jp_tra) ) * tmask(:,:,jk)
         END DO
         ztraz(:,:,1  ) = 0.0
         ztraz(:,:,jpk) = 0.0

         ! slopes
         DO jk = 2, jpkm1
            DO jj = 1,jpj
               DO ji = 1, jpi
                  zign = 0.25 + SIGN( 0.25, ztraz(ji,jj,jk) * ztraz(ji,jj,jk+1) )
                  zakz(ji,jj,jk) = ( ztraz(ji,jj,jk) + ztraz(ji,jj,jk+1) ) * zign
               END DO
            END DO
         END DO
         
         ! Slopes limitation
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zakz(ji,jj,jk) = SIGN( 1., zakz(ji,jj,jk) ) *        &
                     &             MIN( ABS( zakz(ji,jj,jk) ), 2. * ABS(ztraz(ji,jj,jk+1)), 2. * ABS(ztraz(ji,jj,jk) ) )
               END DO
            END DO
         END DO
         
         ! vertical advective flux
         DO jk = 1, jpkm1
            DO jj = 1, jpj      
               DO ji = 1, jpi    
                  zigma = zwsink2(ji,jj,jk+1) * zstep / fse3w(ji,jj,jk+1)
                  zew   = zwsink2(ji,jj,jk+1)
                  psinkflx(ji,jj,jk+1) = -zew * ( trn(ji,jj,jk,jp_tra) - 0.5 * ( 1 + zigma ) * zakz(ji,jj,jk) ) * zstep
               END DO
            END DO
         END DO
         !
         ! Boundary conditions
         psinkflx(:,:,1  ) = 0.e0
         psinkflx(:,:,jpk) = 0.e0
         
         DO jk=1,jpkm1
            DO jj = 1,jpj
               DO ji = 1, jpi
                  zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
                  trn(ji,jj,jk,jp_tra) = trn(ji,jj,jk,jp_tra) + zflx
               END DO
            END DO
         END DO

      ENDDO

      DO jk=1,jpkm1
         DO jj = 1,jpj
            DO ji = 1, jpi
               zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
               ztrb(ji,jj,jk) = ztrb(ji,jj,jk) + 2. * zflx
            END DO
         END DO
      END DO

      trn     (:,:,:,jp_tra) = ztrb(:,:,:)
      psinkflx(:,:,:)        = 2. * psinkflx(:,:,:)
      !
      CALL wrk_dealloc( jpi, jpj, jpk, ztraz, zakz, zwsink2, ztrb )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sink2')
      !
   END SUBROUTINE p4z_sink2


   INTEGER FUNCTION p4z_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink_alloc  ***
      !!----------------------------------------------------------------------
      !
      ALLOCATE( wsbio3 (jpi,jpj,jpk) ,        &
         &      sinking(jpi,jpj,jpk) ,        &
         &                                    STAT=p4z_sink_alloc )
         !
      IF( p4z_sink_alloc /= 0 ) CALL ctl_warn('p4z_sink_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_sink_alloc
   
#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sink                    ! Empty routine
   END SUBROUTINE p4z_sink
#endif 

   !!======================================================================
END MODULE  p4zsink
