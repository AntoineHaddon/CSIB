MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  PISCES  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!           CMOC1  !  2013-15  (O. Riche) POC sinking, Oct 28th 2015 simplified sink scheme according to Jim's own modifications for CMOC2
   !!           CMOC1  !  2016-02  (N. Swart) Adds calcite sinking.
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
      ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
      !REAL(wp) ::   zfact, zwsmax, zmax, zstep
      REAL(wp) ::   zwsmax, zmax
      REAL(wp) ::   zrfact2
      REAL(wp), POINTER, DIMENSION(:,:  ) ::   zfpon         ! Calcite export flx at the bottom of the euphotic zone
      REAL(wp), POINTER, DIMENSION(:,:  ) ::   zcalbotflx    ! Calcite flux to sediments
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zcalflxexp    ! exponential decay of calcite flux with depth
      REAL(wp) ::   zcaldiv                                  ! divergence of the calcite flux
      REAL(wp) ::   globvol, globtal                         ! TAL conservation diagnostics
      REAL(wp) ::   ztaleuz, ztalapz, ztalflxsum,ztalapb     ! TAL conservation diagnostics
      REAL(wp), POINTER, DIMENSION(:) :: zdepw ! computation of depths between t-grid cells.
      REAL(wp) ::   r_dci_cmoc                               ! inverse length of calcite dissolution.
      REAL(wp) ::   zdeup, zideup                            ! Euphotic zone depth, inverse depth
      INTEGER  ::   jk_eud_cmoc_p1                           ! level below the euphotic zone
      INTEGER  ::   ikt_p1, ikt                              ! bottom index / plus 1

      
      INTEGER  ::   ik1
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink')

      CALL wrk_alloc( jpi, jpj, zfpon, zcalbotflx)
      CALL wrk_alloc( jpi, jpj, jpk, zcalflxexp )
      CALL wrk_alloc( jpk, zdepw)

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

      sinking (:,:,:) = 0.0e0_wp
      zfpon   (:,:  ) = 0.0_wp


      !   Compute the sedimentation term using p4zsink2 for POC
      !   -----------------------------------------------------

      CALL p4z_sink2( wsbio3, sinking , jppoc )

      !     Calcite sinking flux
      !     --------------------------------------------------------------------
      ! Define an open ocean mask, based on where mbkt > nk_bal_cmoc == 15 (or as define in namelist)
      oomask(:,:) = 0.0_wp
      !WHERE ( mbkt(:,:) >= nk_bal_cmoc ) oomask = 1._wp
      DO jj = 1, jpj
         DO ji = 1,jpi
            ikt = mbkt(ji,jj)
            IF ( ikt > nk_bal_cmoc ) THEN
               oomask(ji, jj) = 1.0_wp
            ENDIF
         ENDDO                                                               
      ENDDO

      ! The euphotic zone depth and inverse depth, used for computing averages.
      jk_eud_cmoc_p1 = jk_eud_cmoc+1  ! t-grid index below the mixed layer 

      DO jj = 1, jpj
         DO ji = 1,jpi
            !  Rain ratio at level jk_eud_cmoc - bottom of the euphotic zone:
            !  Temperature is however taken from the 1st layer (confirmed with RJC, 16/02/2016)
            xrcico(ji,jj) = rmcico_cmoc * exp(aci_cmoc * ( tsn(ji,jj,1,jp_tem)  - trcico_cmoc ) )                  &
            &                         / (1.0_wp + exp(aci_cmoc *( tsn(ji,jj,1,jp_tem) - trcico_cmoc ) + rtrn ) )

           ! PIC export at the bottom of the euphotic zone based on Zahariev et al 2008 p.59
           ! Time stepping is included with xstep, so units are in mol/m2/step
           zfpon(ji,jj) = xrcico(ji,jj) * wsbio3(ji,jj,jk_eud_cmoc) * xstep                                & 
           &                        * trn(ji,jj,jk_eud_cmoc,jppoc)                                         &
           &                        *  tmask_bgc_closea(ji,jj,jk_eud_cmoc) * oomask(ji,jj)

         ENDDO
      ENDDO
      ! Exponential decay of calcite flux with depth, at w-points.
      r_dci_cmoc = 1.0_wp / dci_cmoc
      zcalflxexp(:,:,:) = 0.0_wp
      DO jk = jk_eud_cmoc_p1, jpk                                                         
         DO jj = 1, jpj
            DO ji = 1,jpi
               zcalflxexp(ji,jj,jk) = zfpon(ji,jj) * exp(-1.0_wp*(gdepw(ji,jj,jk)-gdepw(ji,jj,jk_eud_cmoc_p1)) * r_dci_cmoc)      &
               &                                   * tmask_bgc_closea(ji,jj,jk-1)                ! Mask at jk-1 ensures bottom flux
            ENDDO                                                                     ! is included.
         ENDDO
      ENDDO
                
      ! Bottom PIC flux into sediments, which is removed from the deepest layer and
      ! added back to the surface (below).
      DO jj = 1, jpj
         DO ji = 1,jpi
         ikt_p1 = mbkt(ji,jj)+1
         zcalbotflx(ji,jj) = zcalflxexp(ji,jj,ikt_p1)

         ! Set the bottom boundary condition on the calcite flux to zero (no flux to sediment),
         ! and deal with the sediment flux separately from the sinking in the sections below.
           zcalflxexp(ji,jj,ikt_p1) = 0.0_wp
         ENDDO
      ENDDO

      ! Over the levels of the euphotic zone, remove the euphotic-zone averaged
      ! PIC flux (mol/m3) from each level. 
      !ztaleuz = 0._wp
      DO jk =1, jk_eud_cmoc
         DO jj = 1, jpj
            DO ji = 1,jpi
                zdeup = gdepw(ji,jj,jk_eud_cmoc_p1) ! w-grid depth at jk_eud_cmoc_p1 defines bottom
                                      ! boundary of the mixed layer.                           
                zideup = 1.0_wp / zdeup    

               trn(ji,jj,jk,jpdic) = trn(ji,jj,jk,jpdic) -                                   &
               &                              zfpon(ji,jj) * zideup 
               trn(ji,jj,jk,jpdnt) = trn(ji,jj,jk,jpdnt) -                                   &
               &                              zfpon(ji,jj) * zideup 

               trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal) -                                   &
               &                      2.0_wp * zfpon(ji,jj) * zideup 
               !   ztaleuz = ztaleuz + SUM(2._wp * zfpon(:,:) * zideup * fse3t(:,:,jk))
            ENDDO                                                                     ! is included.
         ENDDO
      END DO

      ! Below the euphotic zone:
      ! Compute the divergence of the calcite flux and distribute it over the t-cell. 
      ! No sinking flux through the bottom here, given the masking above.
      !ztalapz = 0._wp
      DO jk = jk_eud_cmoc_p1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zcaldiv =  ( zcalflxexp(ji,jj,jk) - zcalflxexp(ji,jj,jk+1) ) / fse3t(ji,jj,jk) * tmask_bgc_closea(ji,jj,jk)

               trn(ji,jj,jk,jpdic) = trn(ji,jj,jk,jpdic) +         zcaldiv 
               trn(ji,jj,jk,jpdnt) = trn(ji,jj,jk,jpdnt) +         zcaldiv 
               trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal) + 2.0_wp * zcaldiv                      
      !         ztalapz = ztalapz + 2._wp * zcaldiv * fse3t(ji,jj,jk) 
            ENDDO
         ENDDO
      ENDDO

      !ztalflxsum = 0._wp
      !ztalflxsum = SUM(2._wp * zfpon(:,:))
      !ztalapb = ztalapz + SUM(2._wp*zcalbotflx(:,:))
      !WRITE(numout,*) 'talflxsum', ztalflxsum
      !WRITE(numout,*) 'taleuzsum', ztaleuz
      !WRITE(numout,*) 'talapzsum', ztalapz
      !WRITE(numout,*) 'talapbsum', ztalapb
      !WRITE(numout,*) 'sumspace' 

    
      ! Do the bottom sedimentation of calcite. The sedimenting flux is added back
      ! to the surface layer (psuedo "river flux") for conservation.
      DO jj = 1, jpj
         DO ji = 1,jpi
            ikt = mbkt(ji,jj)
            trn(ji,jj,ikt,jpdic) = trn(ji,jj,ikt,jpdic) - zcalbotflx(ji,jj) / fse3t(ji,jj, ikt)
            trn(ji,jj,1,jpdic) = trn(ji,jj,1,jpdic)  + zcalbotflx(ji,jj) / fse3t(ji,jj, 1) 
            trn(ji,jj,ikt,jpdnt) = trn(ji,jj,ikt,jpdnt) - zcalbotflx(ji,jj) / fse3t(ji,jj, ikt)
            trn(ji,jj,1,jpdnt) = trn(ji,jj,1,jpdnt)  + zcalbotflx(ji,jj) / fse3t(ji,jj, 1) 
            trn(ji,jj,ikt,jptal) = trn(ji,jj,ikt,jptal) - 2.0_wp * zcalbotflx(ji,jj) / fse3t(ji,jj,ikt)
            trn(ji,jj,1,jptal) = trn(ji,jj,1,jptal)  + 2.0_wp * zcalbotflx(ji,jj) / fse3t(ji,jj, 1) 
         ENDDO
      ENDDO

      !globvol = glob_sum( cvol(:,:,:) )
      !globtal = glob_sum( trn(:,:,:,jptal) * cvol(:,:,:) ) / globvol
      !WRITE(numout,*) 'TAL integral : ', globtal*1000._wp
      !     --------------------------------------------------------------------

      IF( ln_diatrc ) THEN
         zrfact2 = 1.e3 * rfact2r
         ik1  = iksed + 1
         IF( lk_iomput ) THEN
           IF( jnt == nrdttrc ) THEN
              CALL iom_put( "oomask", oomask(:,:))
              CALL iom_put( "EPC100", sinking(:,:,ik1) * zrfact2 * tmask_bgc_closea(:,:,1) )
              CALL iom_put( "EPCALC100",    zfpon(:,:) * zrfact2 * tmask_bgc_closea(:,:,1) ) !
              ! <CMOC code OR 12/11/2015> denitrification ! CALL iom_put( "BUPOC"  , wsbio3(:,:,11) /rday * zbpoc(:,:) * 1e+3_wp  )  ! POC burial flux
              ! <CMOC code OR 12/11/2015> denitrification ! CALL iom_put( "BUCALC" , zfpon(:,:) * 1e+3_wp * rfact2r * zbpon(:,:)  )  ! <CMOC code OR 12/11/2015> *rfact2r replaces /rfact2 ! PIC burial flux
           ENDIF
         ELSE
           trc2d(:,:,jp_pcs0_2d + 4) = sinking (:,:,ik1) * zrfact2 * tmask_bgc_closea(:,:,1)

         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask_bgc_closea, clinfo=ctrcnm)
      ENDIF
      !

      CALL wrk_dealloc( jpi, jpj, zfpon, zcalbotflx)
      CALL wrk_dealloc( jpi, jpj, jpk, zcalflxexp )
      CALL wrk_dealloc( jpk, zdepw)
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
      NAMELIST/namcmoccal/ rmcico_cmoc, trcico_cmoc, aci_cmoc, dci_cmoc
      NAMELIST/namcmocdeu/ jk_eud_cmoc, nk_bal_cmoc

      !
      !!----------------------------------------------------------------------

      REWIND( numcmoc )
      READ  ( numcmoc, namcmocrr  )
      REWIND( numcmoc )
      READ  ( numcmoc, namcmoccal )
      REWIND( numcmoc )
      READ  ( numcmoc, namcmocdeu )

      ! <CMOC code OR 10/22/2015> CMOC namelist end
      
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters, Redfield ratios, namcmocrr'    
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    C:N ratio                                 cnrr_cmoc =', cnrr_cmoc
         WRITE(numout,*) '    N:C ratio                                 ncrr_cmoc =', ncrr_cmoc
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for calcite export  , namcmoccal'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Maximum rain ratio                      rmcico_cmoc =', rmcico_cmoc
         WRITE(numout,*) '    Rain ratio half-point temperature       trcico_cmoc =', trcico_cmoc
         WRITE(numout,*) '    Rain ratio scaling factor                  aci_cmoc =',    aci_cmoc
         WRITE(numout,*) '    CaCO3 redissolution depth scale            dci_cmoc =',    dci_cmoc
         WRITE(numout,*) ' '
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters, Euphotic zone depth , namcmocdeu'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Index of the bottom of the euphotic zone jk_eud_cmoc =', jk_eud_cmoc
         WRITE(numout,*) '    Open ocean, min. number of vert. layers  nk_bal_cmoc =', nk_bal_cmoc
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
      ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
      !REAL(wp) ::   zigma,zew,zign, zflx, zstep
      REAL(wp) ::   zew, zflx
      ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
      !REAL(wp), POINTER, DIMENSION(:,:,:) :: ztraz, zakz, zwsink2, ztrb
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zwsink2, ztrb 
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sink2')
      !
      ! Allocate temporary workspace
      ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
      !CALL wrk_alloc( jpi, jpj, jpk, ztraz, zakz, zwsink2, ztrb )
      CALL wrk_alloc( jpi, jpj, jpk, zwsink2, ztrb )
      
      ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
      !zstep = rfact2 / 2.

      ! <CMOC code OR Oct 28th 2015> Simple sinking scheme      
!       ztraz(:,:,:) = 0.e0
!       zakz (:,:,:) = 0.e0
      ztrb (:,:,:) = trn(:,:,:,jp_tra)

      DO jk = 1, jpkm1
         zwsink2(:,:,jk+1) = -pwsink(:,:,jk) / rday * tmask_bgc_closea(:,:,jk+1) 
      END DO
      zwsink2(:,:,1) = 0.e0


      ! Vertical advective flux
      ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
      !DO jn = 1, 2
         !  first guess of the slopes interior values
!          DO jk = 2, jpkm1
!             ztraz(:,:,jk) = ( trn(:,:,jk-1,jp_tra) - trn(:,:,jk,jp_tra) ) * tmask_bgc_closea(:,:,jk)
!          END DO
!          ztraz(:,:,1  ) = 0.0
!          ztraz(:,:,jpk) = 0.0
! 
!          ! slopes
!          DO jk = 2, jpkm1
!             DO jj = 1,jpj
!                DO ji = 1, jpi
!                   zign = 0.25 + SIGN( 0.25, ztraz(ji,jj,jk) * ztraz(ji,jj,jk+1) )
!                   zakz(ji,jj,jk) = ( ztraz(ji,jj,jk) + ztraz(ji,jj,jk+1) ) * zign
!                END DO
!             END DO
!          END DO
!          
!          ! Slopes limitation
!          DO jk = 2, jpkm1
!             DO jj = 1, jpj
!                DO ji = 1, jpi
!                   zakz(ji,jj,jk) = SIGN( 1., zakz(ji,jj,jk) ) *        &
!                      &             MIN( ABS( zakz(ji,jj,jk) ), 2. * ABS(ztraz(ji,jj,jk+1)), 2. * ABS(ztraz(ji,jj,jk) ) )
!                END DO
!             END DO
!          END DO
         
         ! vertical advective flux
      DO jk = 1, jpkm1
	DO jj = 1, jpj      
	    DO ji = 1, jpi    
	      ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
	      !zigma = zwsink2(ji,jj,jk+1) * zstep / fse3w(ji,jj,jk+1)
	      zew   = zwsink2(ji,jj,jk+1)
	      ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
	      !psinkflx(ji,jj,jk+1) = -zew * ( trn(ji,jj,jk,jp_tra) - 0.5 * ( 1 + zigma ) * zakz(ji,jj,jk) ) * zstep
	      psinkflx(ji,jj,jk+1) = -zew * trn(ji,jj,jk,jp_tra) * rfact2
	    END DO
	END DO
      END DO
      !
      ! Boundary conditions
      psinkflx(:,:,1  ) = 0.e0
      psinkflx(:,:,jpk) = 0.e0
         
         ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
!          DO jk=1,jpkm1
!             DO jj = 1,jpj
!                DO ji = 1, jpi
!                   zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
!                   trn(ji,jj,jk,jp_tra) = trn(ji,jj,jk,jp_tra) + zflx
!                END DO
!             END DO
!          END DO

      !ENDDO

      DO jk=1,jpkm1
         DO jj = 1,jpj
            DO ji = 1, jpi
               zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / fse3t(ji,jj,jk)
	       ! <CMOC code OR Oct 28th 2015> Simple sinking scheme
               !ztrb(ji,jj,jk) = ztrb(ji,jj,jk) + 2. * zflx
               ztrb(ji,jj,jk) = ztrb(ji,jj,jk) + zflx
            END DO
         END DO
      END DO

      trn     (:,:,:,jp_tra) = ztrb(:,:,:)
      ! <CMOC code OR Oct 28th 2015> Simple sinking scheme      
      ! psinkflx(:,:,:)        = 2. * psinkflx(:,:,:)
      !
      !CALL wrk_dealloc( jpi, jpj, jpk, ztraz, zakz, zwsink2, ztrb )
      CALL wrk_dealloc( jpi, jpj, jpk, zwsink2, ztrb )      
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
