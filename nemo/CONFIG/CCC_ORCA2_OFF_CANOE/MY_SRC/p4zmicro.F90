MODULE p4zmicro
   !!======================================================================
   !!                         ***  MODULE p4zmicro  ***
   !! TOP :   PISCES Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
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
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zint          !  interpolation and computation of various fields
   USE p4zprod         !  production
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_micro         ! called in p4zbio.F90
   PUBLIC   p4z_micro_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_micro_alloc    ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part       = 0.5_wp          !: part of calcite not dissolved in microzoo guts (not used)
   REAL(wp), PUBLIC ::  gmax1      = 1.7_wp          !: maximum grazing rate
   REAL(wp), PUBLIC ::  aps        = 0.075_wp        !: small zooplankton functional response parameter
   REAL(wp), PUBLIC ::  zsr1       = 0.3_wp          !: specific respiration rate
   REAL(wp), PUBLIC ::  lambda1    = 0.8_wp          !: assimilation efficiency

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zmicro.F90 3295 2012-01-30 15:49:07Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_micro( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER  :: ji, jj, jk
      REAL(wp) :: ztn,Tf,spc,spn,spf,chl,grazp,R
      REAL(wp) :: c2n,n2c,c2fe,fe2c,n2fe,fe2n
      REAL(wp) :: cxs,nxs1,fexs1,nxs2,fexs2
      REAL(wp) :: csw1,csw2
      REAL(wp) :: zrfact2
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_micro')
      !
      grazing1(:,:,:) = 0.  !: grazing set to zero

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ztn = tsn(ji,jj,jk,jp_tem)
               Tf = tgfuncz(ji,jj,jk) 

               spc = MAX(trn(ji,jj,jk,jpphy),0.)
               spn = MAX(trn(ji,jj,jk,jpnn),0.)
               spf = MAX(trn(ji,jj,jk,jpnfe),0.)
               chl = MAX(trn(ji,jj,jk,jpnch),0.)

               c2n=spc/(spn+rtrn)
               n2c=spn/(spc+rtrn)
               c2fe=spc/(spf+rtrn)
               fe2c=spf/(spc+rtrn)
               n2fe=spn/(spf+rtrn)
               fe2n=spf/(spn+rtrn)

! Micrograzer functional response is determined by phytoplankton C
               grazp=gmax1*(1.-EXP(-aps*spc))*trn(ji,jj,jk,jpzoo)*xstep
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
               R = MAX(zsr1*Tf*trn(ji,jj,jk,jpzoo)*xstep-cxs,0.)

               ! Grazing by microzooplankton
               grazing1(ji,jj,jk) = grazp

               !  Update of the TRA arrays
               !  ------------------------
               !zgrarsig  = zgrarem * sigma1
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + R*rr_n2c + nxs1 + nxs2
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - R - cxs
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + R*rr_fe2c + fexs1 + fexs2
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + (R + cxs)*1.E-6
               !   Update the arrays TRA which contain the biological sources and sinks
               !   --------------------------------------------------------------------
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) + lambda1*grazp - R
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - grazp - cxs
               tra(ji,jj,jk,jpnn) = tra(ji,jj,jk,jpnn) - grazp*rr_n2c - (nxs1+nxs2)
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - grazp*rr_fe2c - (fexs1+fexs2)
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - (grazp+cxs)*chl/(spc+rtrn)
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + (1.-lambda1)*grazp
            END DO
         END DO
      END DO
      !
      IF( ln_diatrc ) THEN
         zrfact2 = 1.e-3 * rfact2r  ! conversion from umol/L/timestep into mol/m3/s
         IF( jnt == nrdttrc ) THEN
          CALL iom_put( "GRAZ1"   , grazing1(:,:,:) * zrfact2 * tmask(:,:,:) )  ! microzooplankton grazing on nanophytoplankton
         ENDIF
      ENDIF

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

      NAMELIST/nampiszoo/ part, gmax1, aps, zsr1, lambda1

      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampiszoo )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for microzooplankton, nampiszoo'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in microzoo guts  part        =', part
         WRITE(numout,*) '    Microzooplankton maximum grazing rate           gmax1       =', gmax1
         WRITE(numout,*) '    Microzooplankton grazing initial slope          aps         =', aps
         WRITE(numout,*) '    Microzooplankton specific respiration           zsr1        =', zsr1
         WRITE(numout,*) '    Microzooplankton unassimilated fraction         lambda1     =', lambda1
      ENDIF

   END SUBROUTINE p4z_micro_init

   INTEGER FUNCTION p4z_micro_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_micro_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( grazing1(jpi,jpj,jpk), STAT=p4z_micro_alloc )
      IF( p4z_micro_alloc /= 0 ) CALL ctl_warn('p4z_micro_alloc : failed to allocate arrays.')

   END FUNCTION p4z_micro_alloc

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
