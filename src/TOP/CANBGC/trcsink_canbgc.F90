MODULE trcsink_canbgc
   !!======================================================================
   !!                         ***  MODULE cmocsink  ***
   !! TOP :  CANBGC  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!           CMOC1  !  2013-15  (O. Riche) POC sinking, Oct 28th 2015 simplified sink scheme according to Jim's own modifications for CMOC2
   !!           CMOC1  !  2016-02  (N. Swart) Adds calcite sinking.
   !!           CMOC2  !  2022-23  NEMO4 implementation
   !!----------------------------------------------------------------------
   !!   cmoc_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   cmoc_sink_init  :  Unitialisation of sinking speed parameters
   !!   cmoc_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   USE in_out_manager  ! I/O manager
   USE dom_oce         ! ocean space and time domain 
   USE timing          ! Timing
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   USE sms_top_canbgc
   USE sms_cmoc
   USE sms_canoe
   USE trc_closea_canbgc ! tmask_bgc_closea

   IMPLICIT NONE
   PRIVATE

   PUBLIC cmoc_sink         ! called in trcsms_cmoc.F90
   PUBLIC cmoc_sink_init    ! called in trcsms_cmoc.F90
   PUBLIC cmoc_sink_alloc

   PUBLIC canoe_sink         ! called in trcsms_canoe.F90
   PUBLIC canoe_sink_init    ! called in trcsms_canoe.F90
   PUBLIC canoe_sink_alloc

   ! Common CanBGC arrays
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking  !: POC sinking fluxes
   ! CanOE specific arrays
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio4   !: GOC sinking speed
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wscal    !: Calcite sinking speed
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking2 !: POC sinking fluxes 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkcal  !: CaCO3 sinking flux

   INTEGER  :: iksed  = 10

   !!* Substitution
!#  include "top_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: cmocsink.F90 3685 2012-11-27 15:39:02Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE cmoc_sink_init  ***
      !!----------------------------------------------------------------------
 
  SUBROUTINE cmoc_sink( kt , jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - NEED TO BE DETAILED FOR ONCE
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt, jnt
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zwsmax, zmax
      REAL(wp) ::   zrfact2
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) ::   zfpon         ! Calcite export flx at the bottom of the euphotic zone
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) ::   zcalbotflx    ! Calcite flux to sediments
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zcalflxexp    ! exponential decay of calcite flux with depth
      REAL(wp) ::   zcaldiv                                      ! divergence of the calcite flux
      !
      REAL(wp) ::   globvol, globtal                             ! TAL conservation diagnostics
      REAL(wp) ::   ztaleuz, ztalapz, ztalflxsum,ztalapb         ! TAL conservation diagnostics
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:) :: zdepw               ! computation of depths between t-grid cells.
      REAL(wp) ::   r_dci_cmoc                                   ! inverse length of calcite dissolution.
      REAL(wp) ::   zdeup, zideup                                ! Euphotic zone depth, inverse depth
      INTEGER  ::   jk_eud_cmoc_p1                               ! level below the euphotic zone
      INTEGER  ::   ikt_p1, ikt                                  ! bottom index / plus 1
      !
      INTEGER  ::   ik1
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('cmoc_sink')
      !
      ALLOCATE( zfpon( jpi, jpj ), zcalbotflx( jpi, jpj ))
      ALLOCATE( zcalflxexp ( jpi, jpj, jpk ) )
      ALLOCATE( zdepw( jpk ) )
      !
      !    Sinking speeds of detritus is increased with depth as shown
      !    by data and from the coagulation theory
      !    -----------------------------------------------------------
      ! limit the values of the sinking speeds to avoid numerical instabilities
      wsbio3(:,:,:) = ws_cmoc  ! = wsbio   ! OR Nov 18th 2022 ! Was wsbio declared?
      !
      DO jk = 1,jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zwsmax = 0.8_wp * e3t_n(ji,jj,jk) / xstepb
               wsbio3(ji,jj,jk) = MIN( wsbio3(ji,jj,jk), zwsmax )
            END DO
         END DO
      END DO
      !
      !  Initialization to zero all the sinking arrays 
      !   -----------------------------------------
      !
      sinking (:,:,:) = 0._wp
      zfpon   (:,:  ) = 0._wp
      !
      !   Compute the sedimentation term using cmocsink2 for POC
      !   -----------------------------------------------------
      !
      CALL cmoc_sink2( wsbio3, sinking , jppoc )
      !
      !     Calcite sinking flux
      !     --------------------------------------------------------------------
      ! Define an open ocean mask, based on where mbkt > nk_bal_cmoc == 15 (or as define in namelist)
      oomask(:,:) = 0.0_wp
      !
      DO jj = 1, jpj
         DO ji = 1,jpi
            ikt = mbkt(ji,jj)
            IF ( ikt > nk_bal_cmoc ) THEN
               oomask(ji, jj) = 1.0_wp
            ENDIF
         ENDDO                                                               
      ENDDO
      !
      ! The euphotic zone depth and inverse depth, used for computing averages.
      jk_eud_cmoc_p1 = jk_eud_cmoc+1  ! t-grid index below the mixed layer 
      !
      DO jj = 1, jpj
         DO ji = 1,jpi
            !  Rain ratio at level jk_eud_cmoc - bottom of the euphotic zone:
            !  Temperature is however taken from the 1st layer (confirmed with RJC, 16/02/2016)
            xrcico(ji,jj) = rmcico_cmoc * exp(aci_cmoc * ( tsn(ji,jj,1,jp_tem)  - trcico_cmoc ) )                  &
            &                         / (1.0_wp + exp(aci_cmoc *( tsn(ji,jj,1,jp_tem) - trcico_cmoc ) + rtrn ) )
            !
            ! PIC export at the bottom of the euphotic zone based on Zahariev et al 2008 p.59
            ! Time stepping is included with xstepb, so units are in mol/m2/step
           zfpon(ji,jj) = xrcico(ji,jj) * wsbio3(ji,jj,jk_eud_cmoc) * xstepb                                & 
           &                        * trn(ji,jj,jk_eud_cmoc,jppoc)                                          &
           &                        * tmask_bgc_closea(ji,jj,jk_eud_cmoc) * oomask(ji,jj)
           !
         ENDDO
      ENDDO
      ! Exponential decay of calcite flux with depth, at w-points.
      !
      r_dci_cmoc = 1.0_wp / dci_cmoc
      zcalflxexp(:,:,:) = 0.0_wp
      DO jk = jk_eud_cmoc_p1, jpk                                                         
         DO jj = 1, jpj
            DO ji = 1,jpi
               zcalflxexp(ji,jj,jk) = zfpon(ji,jj) * exp(-1.0_wp*(gdepw_n(ji,jj,jk)-gdepw_n(ji,jj,jk_eud_cmoc_p1)) * r_dci_cmoc)      &
               &                                   * tmask_bgc_closea(ji,jj,jk-1)                ! Mask at jk-1 ensures bottom flux
            ENDDO                                                                                ! is included.
         ENDDO
      ENDDO
      !          
      ! Bottom PIC flux into sediments, which is removed from the deepest layer and
      ! added back to the surface (below).
      DO jj = 1, jpj
         DO ji = 1,jpi
           ikt_p1 = mbkt(ji,jj)+1
           zcalbotflx(ji,jj) = zcalflxexp(ji,jj,ikt_p1)
           !
           ! Set the bottom boundary condition on the calcite flux to zero (no flux to sediment),
           ! and deal with the sediment flux separately from the sinking in the sections below.
           zcalflxexp(ji,jj,ikt_p1) = 0.0_wp
         ENDDO
      ENDDO
      !
      ! Over the levels of the euphotic zone, remove the euphotic-zone averaged
      ! PIC flux (mol/m3) from each level. 
      !ztaleuz = 0._wp
      !
      DO jk =1, jk_eud_cmoc
         DO jj = 1, jpj
            DO ji = 1,jpi
               zdeup = gdepw_n(ji,jj,jk_eud_cmoc_p1) ! w-grid depth at jk_eud_cmoc_p1 defines bottom
                                                     ! boundary of the mixed layer.                           
               zideup = 1.0_wp / zdeup    
               !
               trn(ji,jj,jk,jpdic) = trn(ji,jj,jk,jpdic) -                                   &
               &                              zfpon(ji,jj) * zideup 
               ! trn(ji,jj,jk,jpdnt) = trn(ji,jj,jk,jpdnt) -                                   &
               ! &                              zfpon(ji,jj) * zideup 
               !
               trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal) -                                   &
               &                      2.0_wp * zfpon(ji,jj) * zideup 
               !
            ENDDO
         ENDDO
      END DO
      !
      ! Below the euphotic zone:
      ! Compute the divergence of the calcite flux and distribute it over the t-cell. 
      !
      DO jk = jk_eud_cmoc_p1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zcaldiv =  ( zcalflxexp(ji,jj,jk) - zcalflxexp(ji,jj,jk+1) ) / e3t_n(ji,jj,jk) * tmask_bgc_closea(ji,jj,jk)
               !
               trn(ji,jj,jk,jpdic) = trn(ji,jj,jk,jpdic) +          zcaldiv 
               ! trn(ji,jj,jk,jpdnt) = trn(ji,jj,jk,jpdnt) +          zcaldiv 
               trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal) + 2.0_wp * zcaldiv                      
               !
            ENDDO
         ENDDO
      ENDDO
      !
      ! Do the bottom sedimentation of calcite. The sedimenting flux is added back
      ! to the surface layer (psuedo "river flux") for conservation.
      !
      DO jj = 1, jpj
         DO ji = 1,jpi
            ikt = mbkt(ji,jj)
            trn(ji,jj,ikt,jpdic) = trn(ji,jj,ikt,jpdic) - zcalbotflx(ji,jj)          / e3t_n(ji,jj, ikt)
            trn(ji,jj,1,jpdic)   = trn(ji,jj,1,jpdic)   + zcalbotflx(ji,jj)          / e3t_n(ji,jj, 1) 
            ! trn(ji,jj,ikt,jpdnt) = trn(ji,jj,ikt,jpdnt) - zcalbotflx(ji,jj)          / e3t_n(ji,jj, ikt)
            ! trn(ji,jj,1,jpdnt)   = trn(ji,jj,1,jpdnt)   + zcalbotflx(ji,jj)          / e3t_n(ji,jj, 1) 
            trn(ji,jj,ikt,jptal) = trn(ji,jj,ikt,jptal) - 2.0_wp * zcalbotflx(ji,jj) / e3t_n(ji,jj,ikt)
            trn(ji,jj,1,jptal)   = trn(ji,jj,1,jptal)   + 2.0_wp * zcalbotflx(ji,jj) / e3t_n(ji,jj, 1) 
         ENDDO
      ENDDO
      !
      !
      ! IF( ln_diatrc ) THEN
         ! zrfact2 = 1.e3 * rfact2r
         ! ik1  = iksed + 1
         ! IF( lk_iomput ) THEN
           ! IF( jnt == qnrdttrc ) THEN
              ! CALL iom_put( "oomask", oomask(:,:))
              ! CALL iom_put( "EPC100", sinking(:,:,ik1) * zrfact2 * tmask_bgc_closea(:,:,1) )
              ! CALL iom_put( "EPCALC100",    zfpon(:,:) * zrfact2 * tmask_bgc_closea(:,:,1) ) !
           ! ENDIF
         ! ELSE
           ! trc2d(:,:,jp_pcs0_2d + 4) = sinking (:,:,ik1) * zrfact2 * tmask_bgc_closea(:,:,1)
         ! !
         ! ENDIF
      ! ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('cmocsink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask_bgc_closea, clinfo=ctrcnm)
      ENDIF
      !
      !
      DEALLOCATE( zfpon, zcalbotflx, zcalflxexp, zdepw)
      !
      IF( ln_timing )  CALL timing_stop('cmoc_sink')
      ! 
  END SUBROUTINE cmoc_sink


  SUBROUTINE cmoc_sink2( pwsink, psinkflx, jp_tra )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_sink2  ***
      !!
      !! ** Purpose :   Compute the sedimentation terms for the various sinking
      !!                particles. The scheme used to compute the trends is based
      !!                on MUSCL.
      !!
      !! ** Method  : - this ROUTINE computes not exactly the advection but the
      !!                transport term, i.e.  div(u*tra).
      !!---------------------------------------------------------------------
      !
      INTEGER , INTENT(in   )                         ::   jp_tra    ! tracer index index
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj,jpk) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   psinkflx  ! sinking fluxe
      !!
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zew, zflx
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwsink2, ztrb 
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('cmoc_sink2')
      !
      ! Allocate temporary workspace
      ALLOCATE( zwsink2(jpi, jpj, jpk), ztrb(jpi, jpj, jpk) )
      !
      ztrb (:,:,:) = trn(:,:,:,jp_tra)
      !
      DO jk = 1, jpkm1
         zwsink2(:,:,jk+1) = -pwsink(:,:,jk) / rday * tmask_bgc_closea(:,:,jk+1) 
      END DO
      zwsink2(:,:,1) = 0._wp
      !
      ! vertical advective flux
      DO jk = 1, jpkm1
        DO jj = 1, jpj      
            DO ji = 1, jpi    
              zew   = zwsink2(ji,jj,jk+1)
              psinkflx(ji,jj,jk+1) = -zew * trn(ji,jj,jk,jp_tra) * qfact2
            END DO
        END DO
      END DO
      !
      ! Boundary conditions
      psinkflx(:,:,1  ) = 0._wp
      psinkflx(:,:,jpk) = 0._wp
      !
      DO jk=1,jpkm1
         DO jj = 1,jpj
            DO ji = 1, jpi
               zflx = ( psinkflx(ji,jj,jk) - psinkflx(ji,jj,jk+1) ) / e3t_n(ji,jj,jk)
               ztrb(ji,jj,jk) = ztrb(ji,jj,jk) + zflx
            END DO
         END DO
      END DO
      !
      trn(:,:,:,jp_tra) = ztrb(:,:,:)
      !
      CALL DEALLOCATE( zwsink2, ztrb )      
      !
      IF( ln_timing )  CALL timing_stop('cmoc_sink2')
      !
  END SUBROUTINE cmoc_sink2


  SUBROUTINE cmoc_sink_init
      !
      INTEGER :: ji, jj, ikt, ios     !: working integers for loops and I/O
      !
      NAMELIST/namcmoccal/ rmcico_cmoc, trcico_cmoc, aci_cmoc, dci_cmoc
      !
      REWIND( numnatp_refb )              ! Namelist namcmoccal in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcmoccal, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcmoccal in reference namelist_cmoc' )
      REWIND( numnatp_cfgb )              ! Namelist namcmoccal in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcmoccal, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcmoccal in configuration namelist_cmoc' )

      IF(lwm) WRITE( numonpb, namcmoccal )

      !
      IF( lwp ) THEN
        WRITE(numout,*) ' Namelist parameters for calcite export  , namcmoccal'
        WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*) '    Maximum rain ratio                      rmcico_cmoc =', rmcico_cmoc
        WRITE(numout,*) '    Rain ratio half-point temperature       trcico_cmoc =', trcico_cmoc
        WRITE(numout,*) '    Rain ratio scaling factor                  aci_cmoc =',    aci_cmoc
        WRITE(numout,*) '    CaCO3 redissolution depth scale            dci_cmoc =',    dci_cmoc
        WRITE(numout,*) ' '      
      END IF
      !
      ! Define an open ocean mask, based on where mbkt > nk_bal_cmoc == 28 (or as define in namelist_cmoc_ref)
      ! Move to here instead of cmoc_sink2 in the original CanESM5/CMOC code, since this only needs to be defined
      ! once
      !
      ALLOCATE( oomask(jpi, jpj) )
      oomask(:,:) = 0.0_wp
      !
      DO jj = 1, jpj
         DO ji = 1,jpi
            ikt = mbkt(ji,jj)
            IF ( ikt > nk_bal_cmoc ) THEN
               oomask(ji, jj) = 1.0_wp
            ENDIF
         ENDDO                                                               
      ENDDO
      !
      !
  END SUBROUTINE cmoc_sink_init
  
  
  INTEGER FUNCTION cmoc_sink_alloc()
    !!----------------------------------------------------------------------
    !!                     ***  ROUTINE cmoc_sink_alloc  ***
    !!----------------------------------------------------------------------
    !
    ALLOCATE( wsbio3 (jpi,jpj,jpk) ,        &
       &      sinking(jpi,jpj,jpk) ,        &
       &                                    STAT=cmoc_sink_alloc )
       !
    IF( cmoc_sink_alloc /= 0 ) CALL ctl_warn('cmoc_sink_alloc : failed to allocate arrays.')
    !
  END FUNCTION cmoc_sink_alloc


  SUBROUTINE canoe_sink ( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE canoe_sink  ***
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
      IF( ln_timing )  CALL timing_start('canoe_sink')
      !
      !  Initialize to zero all the sinking arrays 
      !   -----------------------------------------
      !
      sinking (:,:,:) = 0.e0
      sinking2(:,:,:) = 0.e0
      sinkcal (:,:,:) = 0.e0
      !
      !   Compute the sedimentation term using canoesink2 for all the sinking particles
      !   -----------------------------------------------------
      !
      CALL canoe_sink2( wsbio3, sinking , jppoc )
      CALL canoe_sink2( wsbio4, sinking2, jpgoc )
      CALL canoe_sink2( wscal , sinkcal , jpcal )
      !
      ! IF( ln_diatrc ) THEN
         ! zrfact2 = 1.e-3 * qfact2r
         ! ik1  = iksed + 1
         ! IF( lk_iomput ) THEN
           ! IF( jnt == nrdttrc ) THEN
              ! CALL iom_put( "EPC100"  , ( sinking(:,:,ik1) + sinking2(:,:,ik1) ) * zrfact2 * tmask_bgc_closea(:,:,1) ) ! Export of carbon at 100m
              ! CALL iom_put( "EPCALC100",  sinkcal(:,:,ik1)                       * zrfact2 * tmask_bgc_closea(:,:,1) ) ! Export of calcite  at 100m
           ! ENDIF
         ! ELSE
           ! trc2d(:,:,jp_pcs0_2d + 4) = sinking (:,:,ik1) * zrfact2 * tmask_bgc_closea(:,:,1)
           ! trc2d(:,:,jp_pcs0_2d + 5) = sinking2(:,:,ik1) * zrfact2 * tmask_bgc_closea(:,:,1)
           ! trc2d(:,:,jp_pcs0_2d + 9) = sinkcal (:,:,ik1) * zrfact2 * tmask_bgc_closea(:,:,1)
         ! ENDIF
      ! ENDIF
      ! !
      ! !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('canoe_sink')
      !
  END SUBROUTINE canoe_sink

  
  SUBROUTINE canoe_sink2
   
   
  END SUBROUTINE canoe_sink2


  SUBROUTINE canoe_sink_init
   
   
  END SUBROUTINE canoe_sink_init
   
   
  INTEGER FUNCTION canoe_sink_alloc()
    !!----------------------------------------------------------------------
    !!                     ***  ROUTINE cmoc_sink_alloc  ***
    !!----------------------------------------------------------------------
    !
    ALLOCATE( wsbio3 (jpi,jpj,jpk),  wsbio4  (jpi,jpj,jpk),        &
       &      sinking(jpi,jpj,jpk),  sinking2(jpi,jpj,jpk),        &
       &      sinkcal(jpi,jpj,jpk),          STAT=canoe_sink_alloc )
       !
    IF( canoe_sink_alloc /= 0 ) CALL ctl_warn('canoe_sink_alloc : failed to allocate arrays.')
    !
  END FUNCTION canoe_sink_alloc
   


END MODULE trcsink_canbgc