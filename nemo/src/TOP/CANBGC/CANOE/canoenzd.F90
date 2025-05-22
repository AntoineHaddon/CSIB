MODULE canoenzd
   
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 

   USE sms_top_canbgc     ! TOP Source Minus Sink variables
   USE sms_canoe          ! CanOE specific parameters declaration
   USE trc_closea_canbgc  !  tmask_bgc_closea
   USE trcopt_canbgc      ! PAR attenuation
   USE canoetemp          ! CanOE temperature dependencies module

   USE prtctl          !  print control for debugging
   USE lib_mpp         !  ctl_stop on failed mem allocate check
   USE lib_fortran     !  access glob_sum function
   USE iom             !  I/O manager

   ! timing modules
   USE in_out_manager  ! nn_timing integer
   USE timing          ! *_timing subroutines

   IMPLICIT NONE
   PRIVATE

   !!----------------------------------------------------------------------
   PUBLIC canoe_meso
   PUBLIC canoe_mzoo
   PUBLIC canoe_mort1
   PUBLIC canoe_mort2
   PUBLIC canoe_rem
   PUBLIC canoe_nzd_init
   PUBLIC canoe_nzd_alloc
 
   REAL(wp), PUBLIC :: mpratps = 5.E-2_wp   !: small phytoplankton mortality rate
   REAL(wp), PUBLIC :: mpratpl = 5.E-2_wp   !: large phytoplankton mortality rate
   REAL(wp), PUBLIC :: mpratzs = 5.E-2_wp   !: microzooplankton mortality rate
   REAL(wp), PUBLIC :: mpratzl = 2.E-1_wp   !: mesozooplankton mortality rate
   REAL(wp), PUBLIC :: mpquaps = 6.E-2_wp   !: quadratic mortality of small phytoplankton
   REAL(wp), PUBLIC :: mpquapl = 6.E-2_wp   !: quadratic mortality of large phytoplankton
   REAL(wp), PUBLIC :: mpquazs = 6.E-2_wp   !: quadratic mortality of microzooplankton
   REAL(wp), PUBLIC :: mpquazl = 6.E-2_wp   !: quadratic mortality of mesozooplankton
!   REAL(wp), PUBLIC :: mprat   = 5.E-2_wp   !: phytoplankton mortality rate 
!   REAL(wp), PUBLIC :: mprat2  = 2.E-1_wp   !: Diatoms mortality rate
!   REAL(wp), PUBLIC :: mpratm  = 5.E-2_wp   !: Phytoplankton minimum mortality rate
!   REAL(wp), PUBLIC :: mpqua   = 1.E-09_wp  !: quadratic mortality of phytoplankton
!   REAL(wp), PUBLIC :: mpquad  = 2.E-08_wp  !: maximum quadratic mortality of diatoms
   REAL(wp), PUBLIC :: chldegr = 2.E-2_wp   !: Chlorophyll photooxidation rate
   REAL(wp), PUBLIC :: picfrx  = 1.E-1_wp   !: CaCO3 fraction of mortality (0.1 implies 1 mol caCO3 for each 10 mol POC)
   REAL(wp), PUBLIC :: xminp   = 0.01       !: minimum phytoplankton concentration for linear mortality
   REAL(wp), PUBLIC :: part    = 0.5_wp     !: part of calcite not dissolved in microzoo guts (not used)
   REAL(wp), PUBLIC :: gmax1   = 1.7_wp     !: maximum grazing rate
   REAL(wp), PUBLIC :: aps     = 0.075_wp   !: small zooplankton functional response parameter
   REAL(wp), PUBLIC :: zsr1    = 0.3_wp     !: specific respiration rate
   REAL(wp), PUBLIC :: lambda1 = 0.8_wp     !: assimilation efficiency
   REAL(wp), PUBLIC :: part2   = 0.5_wp     !: part of calcite not dissolved in mesozoo guts (not used)
   REAL(wp), PUBLIC :: gmax2   = 0.85_wp    !: maximum grazing rate rate
   REAL(wp), PUBLIC :: apl     = 0.075_wp   !: large zooplankton functional response parameter
   REAL(wp), PUBLIC :: zsr2    = 0.3_wp     !: specific respiration rate
   REAL(wp), PUBLIC :: lambda2 = 0.8_wp     !: assimilation efficiency
   REAL(wp), PUBLIC :: xremik = 0.25_wp     !: remineralisation rate of POC 
   REAL(wp), PUBLIC :: xremip = 0.025_wp    !: remineralisation rate of DOC (not used)
   REAL(wp), PUBLIC :: nitrif = 0.05_wp     !: NH4 nitrification rate 
   REAL(wp), PUBLIC :: xlam1  = 0.0001_wp   !: scavenging rate of iron (low concentrations)
   REAL(wp), PUBLIC :: xlam2  = 2.5_wp      !: scavenging rate of iron (high concentrations)
   REAL(wp), PUBLIC :: ligand = 6.0E+2_wp   !: ligand concentration
   REAL(wp), PUBLIC :: pocfctr = 0.65574_wp !: multiplier for POC-dependent Fe scavenging
   REAL(wp), PUBLIC :: o2thresh = 6._wp     !: O2 threshold for denitrification
   REAL(wp), PUBLIC :: nh4frx = 0.25_wp     !: anammox fraction of denitrification
   REAL(wp), PUBLIC :: oxymin = 1._wp       !: half saturation constant for O2 inhibition of nitrification
   REAL(wp), PUBLIC :: nyld   = 0.8_wp      !: denitrification stoichiometric coefficient
   REAL(wp), PUBLIC :: kdca   = 0.0074_wp   !: dissolution rate of CaCO3
   REAL(wp), PUBLIC :: nca    = 1._wp       !: order of dissolution reaction (not used)

   !REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr

#  include "vectopt_loop_substitute.h90"

CONTAINS
   !
   SUBROUTINE canoe_meso( kt, jnt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE canoe_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      INTEGER  :: ji, jj, jk
      REAL(wp) :: ztn,Tf,lpc,lpn,lpf,chl,grazt,grazp,grazz,R,szc,itfc
      REAL(wp) :: c2n,n2c,c2fe,fe2c,n2fe,fe2n
      REAL(wp) :: cxs,nxs1,fexs1,nxs2,fexs2
      REAL(wp) :: csw1,csw2
      CHARACTER (len=25) :: charout
      REAL(wp) :: zrfact2
      !!---------------------------------------------------------------------
      !
      !IF( nn_timing == 1 )  CALL timing_start('canoe_meso')
      !
      !grazing2(:,:,:) = 0.  !: grazing set to zero
      !grazing3(:,:,:) = 0.  !: grazing set to zero

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztn = ts(ji,jj,jk,jp_tem, Kmm)
               Tf = tgfuncz20(ji,jj,jk)

               lpc = MAX(tr(ji,jj,jk,jrdia,Kbb),0.)
               lpn = MAX(tr(ji,jj,jk,jrdn,Kbb),0.)
               lpf = MAX(tr(ji,jj,jk,jrdfe,Kbb),0.)
               chl = MAX(tr(ji,jj,jk,jrdch,Kbb),0.)
               szc = MAX(tr(ji,jj,jk,jrzoo,Kbb),0.)
               itfc=1./(lpc+szc+rtrn)

               c2n=lpc/(lpn+rtrn)
               n2c=lpn/(lpc+rtrn)
               c2fe=lpc/(lpf+rtrn)
               fe2c=lpf/(lpc+rtrn)
               n2fe=lpn/(lpf+rtrn)
               fe2n=lpf/(lpn+rtrn)

! assume grazing hyperbola is determined by total food concentration and the two food types are consumed in proportion to their concentrations (in C units)
               grazt=gmax2*(1.-EXP(-apl*(lpc+szc)))*tr(ji,jj,jk,jrmes,Kbb)*xstepb
               grazz=grazt*szc*itfc
               grazp=grazt*lpc*itfc
! reduce phytoplankton fraction to what can support grazer biomass production based on the least abundant element: the MIN(...) term should be 1 if N and Fe are in excess of the grazer ratio
               grazp=grazp*MIN(n2c*rr_c2n,fe2c*rr_c2fe,1.)
! calculate "excess" relative to grazer RR
               cxs=grazp*MAX(c2n*rr_n2c-1.,c2fe*rr_fe2c-1.,0.)
               nxs1=grazp*(n2c-rr_n2c)
               nxs1=MAX(nxs1,0.)
               nxs2=grazp*rr_n2c*(n2fe*rr_fe2n-1.)
               nxs2=MAX(nxs2,0.)
               fexs1=grazp*(fe2c-rr_fe2c)
               fexs1=MAX(fexs1,0.)
               fexs2=grazp*rr_fe2c*(fe2n*rr_n2fe-1.)
               fexs2=MAX(fexs2,0.)
               csw1=MAX(cxs,0.)
               csw1=csw1/(csw1+rtrn)   !!! csw1 is 1 when cxs>0 and 0 otherwise
               csw2=1.-csw1            !!! csw2 is 0 when cxs>0 and 1 otherwise
               !!! apply csw1 switch on nxs2 and fexs2 terms
               nxs1 = csw2*nxs1
               fexs1= csw2*fexs1
               !!! apply csw2 switch on nxs1 and fexs1 terms
               nxs2 = csw1*nxs2
               fexs2= csw1*fexs2
! calculate zooplankton respiration (in carbon units)
               R = MAX(zsr2*Tf*tr(ji,jj,jk,jrmes,Kbb)*xstepb-cxs,0.)

               !   Update the arrays TRA which contain the biological sources and sinks
               tr(ji,jj,jk,jqdic, Krhs) = tr(ji,jj,jk,jqdic, Krhs) + (R + cxs)*1.E-6
               tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) + (R*rr_n2c + nxs1 + nxs2)*1.E-6
               tr(ji,jj,jk,jqoxy, Krhs) = tr(ji,jj,jk,jqoxy, Krhs) - R - cxs
               tr(ji,jj,jk,jrnh4, Krhs) = tr(ji,jj,jk,jrnh4, Krhs) + R*rr_n2c + nxs1 + nxs2
               tr(ji,jj,jk,jrfer, Krhs) = tr(ji,jj,jk,jrfer, Krhs) + R*rr_fe2c + fexs1 + fexs2
               tr(ji,jj,jk,jrzoo, Krhs) = tr(ji,jj,jk,jrzoo, Krhs) - grazz
               tr(ji,jj,jk,jrmes, Krhs) = tr(ji,jj,jk,jrmes, Krhs) + lambda2*(grazz+grazp) - R
               tr(ji,jj,jk,jrdia, Krhs) = tr(ji,jj,jk,jrdia, Krhs) - grazp - cxs
               tr(ji,jj,jk,jrdn, Krhs) = tr(ji,jj,jk,jrdn, Krhs) - grazp*rr_n2c - (nxs1+nxs2)
               tr(ji,jj,jk,jrdfe, Krhs) = tr(ji,jj,jk,jrdfe, Krhs) - grazp*rr_fe2c - (fexs1+fexs2)
               tr(ji,jj,jk,jrdch, Krhs) = tr(ji,jj,jk,jrdch, Krhs) - (grazp+cxs)*chl/(lpc+rtrn)
               tr(ji,jj,jk,jrgoc, Krhs) = tr(ji,jj,jk,jrgoc, Krhs) + (1.-lambda2)*(grazz+grazp)

! I am excluding cxs from grazing diagnostics (represent zooplankton gains rather than phytoplankton losses)
               !grazing2(ji,jj,jk) = grazp
               !grazing3(ji,jj,jk) = grazz

            END DO
         END DO
      END DO
      !
!      IF( lk_iomput ) THEN
!         zrfact2 = 1.e-3 * qfact2r
!         grazing(:,:,:) = grazing(:,:,:) * zrfact2 * tmask_bgc_closea(:,:,:)   ! Total grazing of phyto by zoo
!         IF( jnt == qnrdttrc ) THEN
!            CALL iom_put( "GRAZ2" , grazing2 * zrfact2 * tmask_bgc_closea(:,:,:) )  ! Grazing of large phytoplankton
!            CALL iom_put( "GRAZ3" , grazing3 * zrfact2 * tmask_bgc_closea(:,:,:) )  ! Grazing of microzooplankton
!         ENDIF
!      ENDIF
!      !
      IF( sn_cfctl%l_prttrc )   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('meso')")
         CALL prt_ctl_info(charout, cdcomp = 'top')
         CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      !IF( nn_timing == 1 )  CALL timing_stop('canoe_meso')
      !
   END SUBROUTINE canoe_meso

   SUBROUTINE canoe_mzoo( kt, jnt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE canoe_mzoo  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      INTEGER  :: ji, jj, jk
      REAL(wp) :: ztn,Tf,spc,spn,spf,chl,grazp,R
      REAL(wp) :: c2n,n2c,c2fe,fe2c,n2fe,fe2n
      REAL(wp) :: cxs,nxs1,fexs1,nxs2,fexs2
      REAL(wp) :: csw1,csw2
      REAL(wp) :: zrfact2
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      !IF( nn_timing == 1 )  CALL timing_start('canoe_mzoo')
      !
      !grazing1(:,:,:) = 0.  !: grazing set to zero

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ztn = ts(ji,jj,jk,jp_tem, Kmm)
               Tf = tgfuncz0(ji,jj,jk) 

               spc = MAX(tr(ji,jj,jk,jrphy,Kbb),0.)
               spn = MAX(tr(ji,jj,jk,jrnn,Kbb),0.)
               spf = MAX(tr(ji,jj,jk,jrnfe,Kbb),0.)
               chl = MAX(tr(ji,jj,jk,jrnch,Kbb),0.)

               c2n=spc/(spn+rtrn)
               n2c=spn/(spc+rtrn)
               c2fe=spc/(spf+rtrn)
               fe2c=spf/(spc+rtrn)
               n2fe=spn/(spf+rtrn)
               fe2n=spf/(spn+rtrn)

! Micrograzer functional response is determined by phytoplankton C
               grazp=gmax1*(1.-EXP(-aps*spc))*tr(ji,jj,jk,jrzoo,Kbb)*xstepb
! reduce phytoplankton consumption to what can support grazer biomass production based on the least abundant element: the MIN(...) term should be 1 if N and Fe are in excess of the grazer ratio
               grazp=grazp*MIN(n2c*rr_c2n,fe2c*rr_c2fe,1.)
! calculate "excess" relative to grazer RR
               cxs=grazp*MAX(c2n*rr_n2c-1.,c2fe*rr_fe2c-1.,0.)
               nxs1=grazp*(n2c-rr_n2c)
               nxs1=MAX(nxs1,0.)
               nxs2=grazp*rr_n2c*(n2fe*rr_fe2n-1.)
               nxs2=MAX(nxs2,0.)
               fexs1=grazp*(fe2c-rr_fe2c)
               fexs1=MAX(fexs1,0.)
               fexs2=grazp*rr_fe2c*(fe2n*rr_n2fe-1.)
               fexs2=MAX(fexs2,0.)
               csw1=MAX(cxs,0.)
               csw1=csw1/(csw1+rtrn)   !!! csw1 is 1 when cxs>0 and 0 otherwise
               csw2=1.-csw1            !!! csw2 is 0 when cxs>0 and 1 otherwise
               !!! apply csw1 switch on nxs2 and fexs2 terms
               nxs1 = csw2*nxs1
               fexs1= csw2*fexs1
               !!! apply csw2 switch on nxs1 and fexs1 terms
               nxs2 = csw1*nxs2
               fexs2= csw1*fexs2
! calculate zooplankton respiration (in carbon units)
               R = MAX(zsr1*Tf*tr(ji,jj,jk,jrzoo,Kbb)*xstepb-cxs,0.)

               ! Grazing by microzooplankton
               !grazing1(ji,jj,jk) = grazp

               !  Update of the TRA arrays
               !  ------------------------
               !zgrarsig  = zgrarem * sigma1
               tr(ji,jj,jk,jqdic, Krhs) = tr(ji,jj,jk,jqdic, Krhs) + (R + cxs)*1.E-6
               tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) + (R*rr_n2c + nxs1 + nxs2)*1.E-6
               tr(ji,jj,jk,jqoxy, Krhs) = tr(ji,jj,jk,jqoxy, Krhs) - R - cxs
               tr(ji,jj,jk,jrfer, Krhs) = tr(ji,jj,jk,jrfer, Krhs) + R*rr_fe2c + fexs1 + fexs2
               tr(ji,jj,jk,jrnh4, Krhs) = tr(ji,jj,jk,jrnh4, Krhs) + R*rr_n2c + nxs1 + nxs2
               !   Update the arrays TRA which contain the biological sources and sinks
               !   --------------------------------------------------------------------
               tr(ji,jj,jk,jrzoo, Krhs) = tr(ji,jj,jk,jrzoo, Krhs) + lambda1*grazp - R
               tr(ji,jj,jk,jrphy, Krhs) = tr(ji,jj,jk,jrphy, Krhs) - grazp - cxs
               tr(ji,jj,jk,jrnn, Krhs) = tr(ji,jj,jk,jrnn, Krhs) - grazp*rr_n2c - (nxs1+nxs2)
               tr(ji,jj,jk,jrnfe, Krhs) = tr(ji,jj,jk,jrnfe, Krhs) - grazp*rr_fe2c - (fexs1+fexs2)
               tr(ji,jj,jk,jrnch, Krhs) = tr(ji,jj,jk,jrnch, Krhs) - (grazp+cxs)*chl/(spc+rtrn)
               tr(ji,jj,jk,jrpoc, Krhs) = tr(ji,jj,jk,jrpoc, Krhs) + (1.-lambda1)*grazp
            END DO
         END DO
      END DO
      !
!      IF( lk_iomput ) THEN
!         zrfact2 = 1.e-3 * rfact2r  ! conversion from umol/L/timestep into mol/m3/s
!         IF( jnt == qnrdttrc ) THEN
!          CALL iom_put( "GRAZ1"   , grazing1(:,:,:) * zrfact2 * tmask_bgc_closea(:,:,:) )  ! microzooplankton grazing on nanophytoplankton
!         ENDIF
!      ENDIF

      IF( sn_cfctl%l_prttrc )   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_info(charout, cdcomp = 'top')
         CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      !IF( nn_timing == 1 )  CALL timing_stop('canoe_mzoo')
      !
   END SUBROUTINE canoe_mzoo
      !
   SUBROUTINE canoe_mort1( kt, jnt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE canoe_mort1  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zmortp,zmortz
      REAL(wp) :: spc,spn,spf,szc,chl
      REAL(wp) :: c2n,n2c,c2fe,fe2c,n2fe,fe2n,thetac
      REAL(wp) :: cxs,nxs1,nxs2,fexs1,fexs2
      REAL(wp) :: csw1,csw2
      REAL(wp) :: zrfact2
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zmortpn
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodcal
      CHARACTER (len=25) :: charout

      ALLOCATE( zmortpn(  jpi, jpj, jpk ), prodcal(  jpi, jpj, jpk ) )
      zmortpn(:,:,:) = 0._wp

      !!---------------------------------------------------------------------
      !
      !IF( nn_timing == 1 )  CALL timing_start('canoe_mort1')
      !

      prodcal(:,:,:) = 0.  !: calcite production variable set to zero

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               spc = MAX(tr(ji,jj,jk,jrphy,Kbb),0.)
               spn = MAX(tr(ji,jj,jk,jrnn,Kbb),0.)
               spf = MAX(tr(ji,jj,jk,jrnfe,Kbb),0.)
               szc = MAX(tr(ji,jj,jk,jrzoo,Kbb),0.)
               chl = MAX(tr(ji,jj,jk,jrnch,Kbb),0.)

               c2n=spc/(spn+rtrn)
               n2c=spn/(spc+rtrn)
               c2fe=spc/(spf+rtrn)
               fe2c=spf/(spc+rtrn)
               n2fe=spn/(spf+rtrn)
               fe2n=spf/(spn+rtrn)
               thetac=chl/(spc+rtrn)

! simplified CMOC type mortality: sum of linear and quadratic terms
               zmortp = mpratps * xstepb * spc + mpquaps * xstepb * spc * spc
               if (spc.le.xminp) zmortp = mpquaps * xstepb * spc * spc           ! no linear mortality below biomass threshold xminp
               zmortz = mpratzs * xstepb * szc + mpquazs * xstepb * szc * szc
               if (szc.le.xminp) zmortz = mpquazs * xstepb * szc * szc
! reduce mortality to what can support detritus production based on the least abundant element: the MIN(...) term should be 1 if N and Fe are in excess of the detritus ratio
               zmortp=zmortp*MIN(n2c*rr_c2n,fe2c*rr_c2fe,1.)
               zmortpn(ji,jj,jk) = zmortp
! calculate "excess" relative to grazer RR
               cxs=zmortp*MAX(c2n*rr_n2c-1.,c2fe*rr_fe2c-1.,0.)
               nxs1=zmortp*(n2c-rr_n2c)
               nxs1=MAX(nxs1,0.)
               nxs2=zmortp*rr_n2c*(n2fe*rr_fe2n-1.)
               nxs2=MAX(nxs2,0.)
               fexs1=zmortp*(fe2c-rr_fe2c)
               fexs1=MAX(fexs1,0.)
               fexs2=zmortp*rr_fe2c*(fe2n*rr_n2fe-1.)
               fexs2=MAX(fexs2,0.)

               csw1=MAX(cxs,0.)
               csw1=csw1/(csw1+rtrn)   !!! csw1 is 1 when cxs>0 and 0 otherwise
               csw2=1.-csw1            !!! csw2 is 0 when cxs>0 and 1 otherwise
               !!! apply csw1 switch on nxs2 and fexs2 terms
               nxs1 = csw2*nxs1
               fexs1= csw2*fexs1
               !!! apply csw2 switch on nxs1 and fexs1 terms
               nxs2 = csw1*nxs2
               fexs2= csw1*fexs2

               !   Update the arrays TRA which contains the biological sources and sinks

               tr(ji,jj,jk,jqdic, Krhs) = tr(ji,jj,jk,jqdic, Krhs) + cxs*1.E-6
               tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) + (nxs1 + nxs2)*1.E-6
               tr(ji,jj,jk,jqoxy, Krhs) = tr(ji,jj,jk,jqoxy, Krhs) - cxs
               tr(ji,jj,jk,jrphy, Krhs) = tr(ji,jj,jk,jrphy, Krhs) - zmortp - cxs
               tr(ji,jj,jk,jrnn, Krhs) = tr(ji,jj,jk,jrnn, Krhs) - zmortp*rr_n2c - (nxs1+nxs2)
               tr(ji,jj,jk,jrnfe, Krhs) = tr(ji,jj,jk,jrnfe, Krhs) - zmortp*rr_fe2c - (fexs1+fexs2)
               tr(ji,jj,jk,jrnch, Krhs) = tr(ji,jj,jk,jrnch, Krhs) - (zmortp+cxs) * thetac - chl*chldegr*xstepb
               tr(ji,jj,jk,jrnh4, Krhs) = tr(ji,jj,jk,jrnh4, Krhs) + nxs1 + nxs2
               tr(ji,jj,jk,jrfer, Krhs) = tr(ji,jj,jk,jrfer, Krhs) + fexs1 + fexs2
               tr(ji,jj,jk,jrzoo, Krhs) = tr(ji,jj,jk,jrzoo, Krhs) - zmortz
               tr(ji,jj,jk,jrpoc, Krhs) = tr(ji,jj,jk,jrpoc, Krhs) + zmortp + zmortz
! Calcification
               tr(ji,jj,jk,jqdic, Krhs) = tr(ji,jj,jk,jqdic, Krhs) - picfrx*(zmortp + zmortz)*1.E-6
               tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) - 2.*picfrx*(zmortp + zmortz)*1.E-6
               tr(ji,jj,jk,jrcal, Krhs) = tr(ji,jj,jk,jrcal, Krhs) + picfrx*(zmortp + zmortz)
               prodcal(ji,jj,jk) = picfrx*(zmortp + zmortz)         ! diagnostic array should be in mmol/m^-3/s but conversion is in p4zmeso for now
            END DO
         END DO
      END DO
      !
      IF( lk_iomput ) THEN
         zrfact2 = 1.e-3 * qfact2r
         prodcal(:,:,:) = prodcal(:,:,:) * zrfact2 * tmask_bgc_closea(:,:,:)   ! Calcite production
         IF( jnt == qnrdttrc ) THEN
            CALL iom_put( "PCAL" , prodcal  )  ! Calcite production
         ENDIF
      ENDIF
!      !
      DEALLOCATE( zmortpn, prodcal )
      !
      IF( sn_cfctl%l_prttrc )   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nano')")
         CALL prt_ctl_info(charout, cdcomp = 'top')
         CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      !IF( nn_timing == 1 )  CALL timing_stop('canoe_mort1')
      !
   END SUBROUTINE canoe_mort1

   SUBROUTINE canoe_mort2( kt, jnt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE canoe_mort2  ***
      !!
      !! ** Purpose :   Compute the mortality terms for diatoms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      INTEGER  ::  ji, jj, jk
      REAL(wp) :: zmortp,zmortz
      REAL(wp) :: spc,spn,spf,szc,chl
      REAL(wp) :: c2n,n2c,c2fe,fe2c,n2fe,fe2n,thetac
      REAL(wp) :: cxs,nxs1,nxs2,fexs1,fexs2
      REAL(wp) :: csw1,csw2
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zmortpd
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      ALLOCATE( zmortpd(  jpi, jpj, jpk ) )
      zmortpd(:,:,:) = 0._wp
      !
      !IF( nn_timing == 1 )  CALL timing_start('canoe_mort2')
      !

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               !     Phytoplankton mortality. 
               !     ------------------------
               spc = MAX(tr(ji,jj,jk,jrdia,Kbb),0.)
               spn = MAX(tr(ji,jj,jk,jrdn,Kbb),0.)
               spf = MAX(tr(ji,jj,jk,jrdfe,Kbb),0.)
               szc = MAX(tr(ji,jj,jk,jrmes,Kbb),0.)
               chl = MAX(tr(ji,jj,jk,jrdch,Kbb),0.)

               c2n=spc/(spn+rtrn)
               n2c=spn/(spc+rtrn)
               c2fe=spc/(spf+rtrn)
               fe2c=spf/(spc+rtrn)
               n2fe=spn/(spf+rtrn)
               fe2n=spf/(spn+rtrn)
               thetac=chl/(spc+rtrn)

               zmortp = mpratpl * xstepb * spc + mpquapl * xstepb * spc * spc
               if (spc.le.xminp) zmortp = mpquapl * xstepb * spc * spc           ! no linear mortality below biomass threshold xminp
               zmortz = mpratzl * xstepb * szc + mpquazl * xstepb * szc * szc
               if (szc.le.xminp) zmortz = mpquazl * xstepb * szc * szc
! reduce mortality to what can support detritus production based on the least abundant element: the MIN(...) term should be 1 if N and Fe are in excess of the detritus ratio
               zmortp=zmortp*MIN(n2c*rr_c2n,fe2c*rr_c2fe,1.)
               zmortpd(ji,jj,jk) = zmortp
! calculate "excess" relative to grazer RR
               cxs=zmortp*MAX(c2n*rr_n2c-1.,c2fe*rr_fe2c-1.,0.)
               nxs1=zmortp*(n2c-rr_n2c)
               nxs1=MAX(nxs1,0.)
               nxs2=zmortp*rr_n2c*(n2fe*rr_fe2n-1.)
               nxs2=MAX(nxs2,0.)
               fexs1=zmortp*(fe2c-rr_fe2c)
               fexs1=MAX(fexs1,0.)
               fexs2=zmortp*rr_fe2c*(fe2n*rr_n2fe-1.)
               fexs2=MAX(fexs2,0.)

               csw1=MAX(cxs,0.)
               csw1=csw1/(csw1+rtrn)   !!! csw1 is 1 when cxs>0 and 0 otherwise
               csw2=1.-csw1            !!! csw2 is 0 when cxs>0 and 1 otherwise
               !!! apply csw1 switch on nxs2 and fexs2 terms
               nxs1 = csw2*nxs1
               fexs1= csw2*fexs1
               !!! apply csw2 switch on nxs1 and fexs1 terms
               nxs2 = csw1*nxs2
               fexs2= csw1*fexs2

               !   Update the arrays tra which contains the biological sources and sinks
               !   ---------------------------------------------------------------------
               tr(ji,jj,jk,jqdic, Krhs) = tr(ji,jj,jk,jqdic, Krhs) + cxs*1.E-6
               tr(ji,jj,jk,jqoxy, Krhs) = tr(ji,jj,jk,jqoxy, Krhs) - cxs
               tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) + (nxs1 + nxs2)*1.E-6
               tr(ji,jj,jk,jrdia, Krhs) = tr(ji,jj,jk,jrdia, Krhs) - zmortp - cxs 
               tr(ji,jj,jk,jrdn, Krhs) = tr(ji,jj,jk,jrdn, Krhs) - zmortp*rr_n2c - (nxs1+nxs2)
               tr(ji,jj,jk,jrdfe, Krhs) = tr(ji,jj,jk,jrdfe, Krhs) - zmortp*rr_fe2c - (fexs1+fexs2)
               tr(ji,jj,jk,jrdch, Krhs) = tr(ji,jj,jk,jrdch, Krhs) - (zmortp+cxs) * thetac - chl*chldegr*xstepb
               tr(ji,jj,jk,jrnh4, Krhs) = tr(ji,jj,jk,jrnh4, Krhs) + nxs1 + nxs2
               tr(ji,jj,jk,jrfer, Krhs) = tr(ji,jj,jk,jrfer, Krhs) + fexs1 + fexs2
               tr(ji,jj,jk,jrmes, Krhs) = tr(ji,jj,jk,jrmes, Krhs) - zmortz 
               tr(ji,jj,jk,jrgoc, Krhs) = tr(ji,jj,jk,jrgoc, Krhs) + zmortp + zmortz
            END DO
         END DO
      END DO
      !
      DEALLOCATE( zmortpd )
      !
      IF( sn_cfctl%l_prttrc )   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diat')")
         CALL prt_ctl_info(charout, cdcomp = 'top')
         CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      !IF( nn_timing == 1 )  CALL timing_stop('canoe_mort2')
      !
   END SUBROUTINE canoe_mort2

   SUBROUTINE canoe_rem( kt, jnt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------

      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zremip, zremik, Tf
      REAL(wp) ::   zkeq, zfeequi
      REAL(wp) ::   zsatur, zsatur2, zdep, zfactdep
      REAL(wp) ::   zorem, zorem2, zofer, zofer2
      REAL(wp) ::   zscave, zscavex, fexs, zcoag
      REAL(wp) ::   zlamfac, zonitr, zstep, znitro2dep
      REAL(wp) ::   zrfact2
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nh4ox
      CHARACTER (len=25) :: charout

      ALLOCATE( nh4ox(  jpi, jpj, jpk ) )
      !REAL(wp), POINTER, DIMENSION(:,:,:) :: zolimi, zolimi2, zwork
      !!---------------------------------------------------------------------
      !
      !IF( nn_timing == 1 )  CALL timing_start('canoe_rem')
      !

      nh4ox(:,:,:)=0.
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep = xstepb
               !    NH4 nitrification to NO3. Ceased for oxygen concentrations
               !    below 2 umol/L. Inhibited at strong light 
               !    ----------------------------------------------------------
               ! nitrification rate depends on O2 concentration
               znitro2dep = MAX(0., 0.4*(6.-tr(ji,jj,jk,jqoxy,Kbb))/(oxymin + tr(ji,jj,jk,jqoxy,Kbb)))
               znitro2dep = MIN(1., znitro2dep )
               zonitr = nitrif * zstep * tr(ji,jj,jk,jrnh4,Kbb) / (1.+ par_3bands(ji,jj,jk)) * (1.- znitro2dep) 
               !denitnh4(ji,jj,jk) = nitrif * zstep * tr(ji,jj,jk,jpnh4,Kbb) * nitrfac(ji,jj,jk) 
               !   Update of the tracers trends
               !   ----------------------------
               tr(ji,jj,jk,jrnh4, Krhs) = tr(ji,jj,jk,jrnh4, Krhs) - zonitr
               tr(ji,jj,jk,jqno3, Krhs) = tr(ji,jj,jk,jqno3, Krhs) + zonitr
               tr(ji,jj,jk,jqoxy, Krhs) = tr(ji,jj,jk,jqoxy, Krhs) - 2. * zonitr
               tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) - 2.e-6 * zonitr
               nh4ox(ji,jj,jk) = zonitr
            END DO
         END DO
      END DO

      IF( sn_cfctl%l_prttrc )   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_info(charout, cdcomp = 'top')
         CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF

      IF( sn_cfctl%l_prttrc )   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem2')")
         CALL prt_ctl_info(charout, cdcomp = 'top')
         CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF

      denitr(:,:,:)=0.
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               Tf = tgfuncr0(ji,jj,jk)
               zorem  = xremik * xstepb * Tf * tr(ji,jj,jk,jrpoc,Kbb)
               zofer  = zorem * rr_fe2c
               zorem2 = xremik * xstepb * Tf * tr(ji,jj,jk,jrgoc,Kbb)
               zofer2 = zorem2 * rr_fe2c

! denitrification is assumed to remove NO3 as a fraction of remineralization increasing linearly from 0 to 1 with declining [O2] for [O2]<6 uM
! NO3 fraction is then divided between NO3 and NH4 according to the parameter nh4frx (for anammox 50% of N comes from NO3 and 50% from NH4)
               zonitr=1.-MIN(tr(ji,jj,jk,jqoxy,Kbb),o2thresh)/o2thresh
               tr(ji,jj,jk,jrnh4, Krhs) = tr(ji,jj,jk,jrnh4, Krhs) + (zorem + zorem2)*rr_n2c - (zorem + zorem2)*nyld*zonitr*0.5*nh4frx
               tr(ji,jj,jk,jqno3, Krhs) = tr(ji,jj,jk,jqno3, Krhs) - (zorem + zorem2)*nyld*zonitr*(1.-0.5*nh4frx)
               tr(ji,jj,jk,jqoxy, Krhs) = tr(ji,jj,jk,jqoxy, Krhs) - (zorem + zorem2)*(1.-zonitr)
               tr(ji,jj,jk,jrfer, Krhs) = tr(ji,jj,jk,jrfer, Krhs) + zofer + zofer2
               tr(ji,jj,jk,jqdic, Krhs) = tr(ji,jj,jk,jqdic, Krhs) + (zorem + zorem2)*1.e-6
               tr(ji,jj,jk,jrpoc, Krhs) = tr(ji,jj,jk,jrpoc, Krhs) - zorem
               tr(ji,jj,jk,jrgoc, Krhs) = tr(ji,jj,jk,jrgoc, Krhs) - zorem2
               tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) + 1.e-6 * (zorem + zorem2)*rr_n2c                         ! 1 mol of alkalinity per mol of N
               tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) + 1.e-6 * (zorem + zorem2)*nyld*zonitr*(1.-nh4frx)        ! +1 mol if denitrification, 0 if anammox
               denitr(ji,jj,jk) = (zorem + zorem2)*zonitr*nyld

! CaCO3 dissolution
               zorem2 = kdca * xstepb * tr(ji,jj,jk,jrcal,Kbb)
               tr(ji,jj,jk,jrcal, Krhs) = tr(ji,jj,jk,jrcal, Krhs) - zorem2 
               tr(ji,jj,jk,jqdic, Krhs) = tr(ji,jj,jk,jqdic, Krhs) + zorem2 * 1.e-6
               tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) + zorem2 * 2.e-6

            END DO
         END DO
      END DO

      IF( sn_cfctl%l_prttrc )   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem3')")
         CALL prt_ctl_info(charout, cdcomp = 'top')
         CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF

      IF( sn_cfctl%l_prttrc )   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem4')")
         CALL prt_ctl_info(charout, cdcomp = 'top')
         CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji = 1, jpi
               zstep = xstepb
! irreversible scavenging as in Christian et al 2002
               zcoag = MIN((tr(ji,jj,jk,jrpoc,Kbb)+tr(ji,jj,jk,jrgoc,Kbb))*pocfctr,1.)
               fexs = MAX(tr(ji,jj,jk,jrfer,Kbb)-ligand,0.)
               zscave = xlam1 * xstepb * (tr(ji,jj,jk,jrfer,Kbb)-fexs) * zcoag
               zscavex = xlam2 * xstepb * fexs
               tr(ji,jj,jk,jrfer, Krhs) = tr(ji,jj,jk,jrfer, Krhs) - (zscave+zscavex)
            END DO
         END DO
      END DO
      !

      IF( sn_cfctl%l_prttrc )   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem5')")
         CALL prt_ctl_info(charout, cdcomp = 'top')
         CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF

      !     Update the arrays TRA which contain the biological sources and sinks
      !     --------------------------------------------------------------------

      IF( lk_iomput ) THEN
         zrfact2 = 1.e-3 * qfact2r  ! conversion from umol/L/timestep into mol/m3/s
         IF( jnt == qnrdttrc ) THEN
           CALL iom_put( "Denitr"   , denitr(:,:,:) * zrfact2 * tmask_bgc_closea(:,:,:) )  ! rate of denitrification
           CALL iom_put( "Nitrif"   , nh4ox(:,:,:) * zrfact2 * tmask_bgc_closea(:,:,:) )  ! rate of nitrification
         ENDIF
      ENDIF

      DEALLOCATE( nh4ox )

      IF( sn_cfctl%l_prttrc )   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem6')")
         CALL prt_ctl_info(charout, cdcomp = 'top')
         CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      !IF( nn_timing == 1 )  CALL timing_stop('canoe_rem')
      !
   END SUBROUTINE canoe_rem

   SUBROUTINE canoe_nzd_init

      !!----------------------------------------------------------------------
      !!
      !! ** Purpose :   Initialization of parameters
      !!
      !! ** Method  :   Read the namelists and check the parameters
      !!                called at the first timestep
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios, ierr ! Local integers
      NAMELIST/namcanmort/ mpratps, mpratzs, mpratpl, mpratzl, mpquaps, mpquazs, mpquapl, mpquazl, chldegr, picfrx, xminp
      NAMELIST/namcanzoo/ part, gmax1, aps, zsr1, lambda1
      NAMELIST/namcanmes/ part2, gmax2, apl, zsr2, lambda2
      NAMELIST/namcanrem/ xremik, xremip, nitrif, xlam1, xlam2, ligand, pocfctr, o2thresh, &
                        & nh4frx, oxymin
      NAMELIST/namcancal/ kdca, nca

      !!----------------------------------------------------------------------

      REWIND( numnatp_refb )              ! Namelist namcanmort in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcanmort, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcanmort in reference namelist_canoe' )
      REWIND( numnatp_cfgb )              ! Namelist namcanmort in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcanmort, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcanmort in configuration namelist_canoe' )
      IF(lwp) WRITE( numonpb, namcanmort )

      REWIND( numnatp_refb )              ! Namelist namcanzoo in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcanzoo, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcanzoo in reference namelist_canoe' )
      REWIND( numnatp_cfgb )              ! Namelist namcanzoo in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcanzoo, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcanzoo in configuration namelist_canoe' )
      IF(lwp) WRITE( numonpb, namcanzoo )

      REWIND( numnatp_refb )              ! Namelist namcanmes in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcanmes, IOSTAT = ios, ERR = 905)
905   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcanmes in reference namelist_canoe' )
      REWIND( numnatp_cfgb )              ! Namelist namcanmes in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcanmes, IOSTAT = ios, ERR = 906 )
906   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcanmes in configuration namelist_canoe' )
      IF(lwp) WRITE( numonpb, namcanmes )

      REWIND( numnatp_refb )              ! Namelist namcanrem in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcanrem, IOSTAT = ios, ERR = 907)
907   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcanrem in reference namelist_canoe' )
      REWIND( numnatp_cfgb )              ! Namelist namcanrem in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcanrem, IOSTAT = ios, ERR = 908 )
908   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcanrem in configuration namelist_canoe' )
      IF(lwp) WRITE( numonpb, namcanrem )

      REWIND( numnatp_refb )              ! Namelist namcancal in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcancal, IOSTAT = ios, ERR = 909)
909   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcancal in reference namelist_canoe' )
      REWIND( numnatp_cfgb )              ! Namelist namcancal in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcancal, IOSTAT = ios, ERR = 910 )
910   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcancal in configuration namelist_canoe' )
      IF(lwp) WRITE( numonpb, namcancal )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton mortality: '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    small phytoplankton mortality rate          mpratps   =', mpratps
         WRITE(numout,*) '    microzooplankton mortality rate             mpratzs   =', mpratzs
         WRITE(numout,*) '    large phytoplankton mortality rate          mpratpl   =', mpratpl
         WRITE(numout,*) '    mesozooplankton mortality rate              mpratzl   =', mpratzl
         WRITE(numout,*) '    quadratic mortality of small phytoplankton  mpquaps   =', mpquaps
         WRITE(numout,*) '    quadratic mortality of microzooplankton     mpquazs   =', mpquazs
         WRITE(numout,*) '    quadratic mortality of large phytoplankton  mpquapl   =', mpquapl
         WRITE(numout,*) '    quadratic mortality of mesozooplankton      mpquazl   =', mpquazl
         WRITE(numout,*) '    Chlorophyll photooxidation rate             chldegr   =', chldegr
         WRITE(numout,*) '    CaCO3 production rate                       picfrx    =', picfrx
         WRITE(numout,*) '    Biomass threshold for linear mortality      xminp     =', xminp

         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for microzooplankton'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in microzoo guts  part =', part
         WRITE(numout,*) '    Microzooplankton maximum grazing rate           gmax1 =', gmax1
         WRITE(numout,*) '    Microzooplankton grazing initial slope          aps =', aps
         WRITE(numout,*) '    Microzooplankton specific respiration           zsr1 =', zsr1
         WRITE(numout,*) '    Microzooplankton unassimilated fraction         lambda1 =', lambda1

         WRITE(numout,*) ' Namelist parameters for mesozooplankton'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in mesozoo guts  part2        =', part2
         WRITE(numout,*) '    Mesozooplankton maximum grazing rate           gmax2        =', gmax2
         WRITE(numout,*) '    Mesozooplankton grazing initial slope          apl          =', apl
         WRITE(numout,*) '    Mesozooplankton specific respiration           zsr2         =', zsr2
         WRITE(numout,*) '    Mesozooplankton unassimilated fraction         lambda2      =', lambda2

         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for remineralization'
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

         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for CaCO3 dissolution'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    calcite dissolution rate constant in d^-1           =', kdca
         WRITE(numout,*) '    order of dissolution reaction (not used)            =', nca

      ENDIF

   END SUBROUTINE canoe_nzd_init

   INTEGER FUNCTION canoe_nzd_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_canoe_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to CANOE
      ! ALLOCATE( tab(...) , STAT=trc_sms_canoe_alloc )
      !
      ALLOCATE( denitr(  jpi, jpj, jpk ) , STAT=canoe_nzd_alloc )
      IF( canoe_nzd_alloc /= 0 ) CALL ctl_stop( 'STOP', 'canoe_nzd_alloc : failed to allocate denitr array' )

   END FUNCTION canoe_nzd_alloc

END MODULE canoenzd

