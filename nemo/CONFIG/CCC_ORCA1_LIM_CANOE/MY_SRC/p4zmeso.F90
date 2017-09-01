MODULE p4zmeso
   !!======================================================================
   !!                         ***  MODULE p4zmeso  ***
   !! TOP :   PISCES Compute the sources/sinks for mesozooplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_meso       :   Compute the sources/sinks for mesozooplankton
   !!   p4z_meso_init  :   Initialization of the parameters for mesozooplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zint          !  interpolation and computation of various fields
   USE p4zprod         !  production
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_meso              ! called in p4zbio.F90
   PUBLIC   p4z_meso_init         ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part2      = 0.5_wp          !: part of calcite not dissolved in mesozoo guts (not used)
   REAL(wp), PUBLIC ::  gmax2      = 0.85_wp         !: maximum grazing rate rate
   REAL(wp), PUBLIC ::  apl        = 0.075_wp        !: large zooplankton functional response parameter
   REAL(wp), PUBLIC ::  zsr2       = 0.3_wp          !: specific respiration rate
   REAL(wp), PUBLIC ::  lambda2    = 0.8_wp          !: assimilation efficiency

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmeso.F90 3295 2012-01-30 15:49:07Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_meso( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER  :: ji, jj, jk
      REAL(wp) :: ztn,Tf,lpc,lpn,lpf,chl,grazt,grazp,grazz,R,szc,itfc
      REAL(wp) :: c2n,n2c,c2fe,fe2c,n2fe,fe2n
      REAL(wp) :: cxs,nxs1,fexs1,nxs2,fexs2
      REAL(wp) :: csw1,csw2
      CHARACTER (len=25) :: charout
      REAL(wp) :: zrfact2
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_meso')
      !
      grazing2(:,:,:) = 0.  !: grazing set to zero
      grazing3(:,:,:) = 0.  !: grazing set to zero

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztn = tsn(ji,jj,jk,jp_tem)
               Tf = tgfuncz2(ji,jj,jk)

               lpc = MAX(trn(ji,jj,jk,jpdia),0.)
               lpn = MAX(trn(ji,jj,jk,jpdn),0.)
               lpf = MAX(trn(ji,jj,jk,jpdfe),0.)
               chl = MAX(trn(ji,jj,jk,jpdch),0.)
               szc = MAX(trn(ji,jj,jk,jpzoo),0.)
               itfc=1./(lpc+szc+rtrn)

               c2n=lpc/(lpn+rtrn)
               n2c=lpn/(lpc+rtrn)
               c2fe=lpc/(lpf+rtrn)
               fe2c=lpf/(lpc+rtrn)
               n2fe=lpn/(lpf+rtrn)
               fe2n=lpf/(lpn+rtrn)

! assume grazing hyperbola is determined by total food concentration and the two food types are consumed in proportion to their concentrations (in C units)
               grazt=gmax2*(1.-EXP(-apl*(lpc+szc)))*trn(ji,jj,jk,jpmes)*xstep
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
               R = MAX(zsr2*Tf*trn(ji,jj,jk,jpmes)*xstep-cxs,0.)

               !   Update the arrays TRA which contain the biological sources and sinks
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + R*rr_n2c + nxs1 + nxs2
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - R - cxs
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + R*rr_fe2c + fexs1 + fexs2
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + (R + cxs)*1.E-6
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - grazz
               tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) + lambda2*(grazz+grazp) - R
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - grazp - cxs
               tra(ji,jj,jk,jpdn) = tra(ji,jj,jk,jpdn) - grazp*rr_n2c - (nxs1+nxs2)
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - grazp*rr_fe2c - (fexs1+fexs2)
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - (grazp+cxs)*chl/(lpc+rtrn)
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + (1.-lambda2)*(grazz+grazp)

! I am excluding cxs from grazing diagnostics (represent zooplankton gains rather than phytoplankton losses)
               grazing2(ji,jj,jk) = grazp
               grazing3(ji,jj,jk) = grazz

            END DO
         END DO
      END DO
      !
      IF( ln_diatrc .AND. lk_iomput ) THEN
         zrfact2 = 1.e-3 * rfact2r
!         grazing(:,:,:) = grazing(:,:,:) * zrfact2 * tmask(:,:,:)   ! Total grazing of phyto by zoo
         prodcal(:,:,:) = prodcal(:,:,:) * zrfact2 * tmask(:,:,:)   ! Calcite production
         IF( jnt == nrdttrc ) THEN
            CALL iom_put( "GRAZ2" , grazing2 * zrfact2 * tmask(:,:,:) )  ! Grazing of large phytoplankton
            CALL iom_put( "GRAZ3" , grazing3 * zrfact2 * tmask(:,:,:) )  ! Grazing of microzooplankton
            CALL iom_put( "PCAL" , prodcal  )  ! Calcite production
         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('meso')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_meso')
      !
   END SUBROUTINE p4z_meso

   SUBROUTINE p4z_meso_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the nampismes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismes
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampismes/ part2, gmax2, apl, zsr2, lambda2

      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampismes )


      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' ' 
         WRITE(numout,*) ' Namelist parameters for mesozooplankton, nampismes'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in mesozoo guts  part2        =', part2
         WRITE(numout,*) '    Mesozooplankton maximum grazing rate           gmax2        =', gmax2
         WRITE(numout,*) '    Mesozooplankton grazing initial slope          apl          =', apl
         WRITE(numout,*) '    Mesozooplankton specific respiration           zsr2         =', zsr2
         WRITE(numout,*) '    Mesozooplankton unassimilated fraction         lambda2      =', lambda2
      ENDIF


   END SUBROUTINE p4z_meso_init


#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_meso                    ! Empty routine
   END SUBROUTINE p4z_meso
#endif 

   !!======================================================================
END MODULE  p4zmeso
