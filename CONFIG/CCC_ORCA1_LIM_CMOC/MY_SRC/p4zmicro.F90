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
! <CMOC OR 06/17/2014> Code trimming !    USE p4zlim          !  Co-limitations
! <CMOC OR 06/17/2014> Code trimming !     USE p4zsink         !  vertical flux of particulate matter due to sinking
! <CMOC OR 06/17/2014> Code trimming !    USE p4zint          !  interpolation and computation of various fields
   USE p4zprod         !  production
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_micro         ! called in p4zbio.F90
   PUBLIC   p4z_micro_init    ! called in trcsms_pisces.F90
   ! <CMOC OR 06/17/2014> Code trimming !  PUBLIC   p4z_micro_alloc    ! called in trcsms_pisces.F90

   !! * Shared module variables
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  part       = 0.5_wp     !: part of calcite not dissolved in microzoo guts
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  xpref2c    = 0.2_wp     !: microzoo preference for POC 
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  xpref2p    = 1.0_wp     !: microzoo preference for nanophyto
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  xpref2d    = 0.6_wp     !: microzoo preference for diatoms
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  xthreshdia = 1E-8_wp    !: diatoms feeding threshold for microzooplankton 
    REAL(wp), PUBLIC ::  xthreshphy = 2E-7_wp    !: nanophyto threshold for microzooplankton 
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  xthreshpoc = 1E-8_wp    !: poc threshold for microzooplankton 
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  xthresh    = 0._wp      !: feeding threshold for microzooplankton 
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  resrat     = 0.03_wp    !: exsudation rate of microzooplankton
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  mzrat      = 0.0_wp     !: microzooplankton mortality rate 
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  grazrat    = 3.0_wp     !: maximal microzoo grazing rate
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  xkgraz     = 20E-6_wp   !: non assimilated fraction of P by microzoo 
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  unass      = 0.3_wp     !: Efficicency of microzoo growth 
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  sigma1     = 0.6_wp     !: Fraction of microzoo excretion as DOM 
! <CMOC OR 06/17/2014> Code trimming !    REAL(wp), PUBLIC ::  epsher     = 0.3_wp     !: half sturation constant for grazing 1 


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
      INTEGER, INTENT(in) ::   kt ! ocean time step
      INTEGER  :: ji, jj, jk
      ! <CMOC OR 06/17/2014> Code trimming !  REAL(wp) :: zcompadi, zcompaz , zcompaph, zcompapoc
      REAL(wp) :: zcompaph
      ! <CMOC OR 06/17/2014> Code trimming !  REAL(wp) :: zgraze  , zdenom, zdenom2, zncratio
      ! <CMOC OR 06/17/2014> Code trimming !  REAL(wp) :: zfact   , zstep, zfood, zfoodlim
      REAL(wp) :: zstep
      ! <CMOC OR 06/17/2014> Code trimming !  REAL(wp) :: zepshert, zepsherv, zgrarsig, zgraztot, zgraztotf
      ! <CMOC OR 06/17/2014> Code trimming !  REAL(wp) :: zgrarem, zgrafer, zgrapoc, zprcaca, zmortz
      REAL(wp) :: zgrapoc
      ! <CMOC OR 06/17/2014> Code trimming !  REAL(wp) :: zrespz, ztortz, zgrasrat
      ! <CMOC OR 06/17/2014> Code trimming !  REAL(wp) :: zgrazp, zgrazpcmoc, zgrazm, zgrazsd
      REAL(wp) :: zgrazpcmoc
      ! <CMOC OR 06/17/2014> Code trimming !  REAL(wp) :: zgrazmf, zgrazsf, zgrazpf
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_micro')
      !
! <CMOC OR 06/17/2014> Code trimming !       grazing(:,:,:) = 0.  !: grazing set to zero
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! <CMOC OR 06/17/2014> Code trimming !  zcompaz = MAX( ( trn(ji,jj,jk,jpzoo) - 1.e-8 ), 0.e0 )
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               ! <CMOC OR 06/17/2014> Code trimming !  zfact   = zstep * tgfunc2(ji,jj,jk) * zcompaz

               !  Respiration rates of both zooplankton
               !  -------------------------------------
               ! <CMOC OR 06/17/2014> Code trimming !  zrespz = resrat * zfact * trn(ji,jj,jk,jpzoo) / ( 2. * xkmort + trn(ji,jj,jk,jpzoo) )  &
               ! <CMOC OR 06/17/2014> Code trimming !     &   + resrat * zfact * 3. * nitrfac(ji,jj,jk)

               !  Zooplankton mortality. A square function has been selected with
               !  no real reason except that it seems to be more stable and may mimic predation.
               !  ---------------------------------------------------------------
               ! <CMOC OR 06/17/2014> Code trimming !  ! <CMOC OR 06/17/2014> Code trimming !  ztortz = mzrat * 1.e6 * zfact * trn(ji,jj,jk,jpzoo)

               ! <CMOC OR 06/17/2014> Code trimming ! zcompadi  = MIN( MAX( ( trn(ji,jj,jk,jpdia) - xthreshdia ), 0.e0 ), xsizedia )
               zcompaph  = MAX( ( trn(ji,jj,jk,jpphy) - xthreshphy ), 0.e0 )
               ! <CMOC OR 06/17/2014> Code trimming ! zcompapoc = MAX( ( trn(ji,jj,jk,jppoc) - xthreshpoc ), 0.e0 )
               
               !     Microzooplankton grazing
               !     ------------------------
               ! <CMOC OR 06/17/2014> Code trimming ! zfood     = xpref2p * zcompaph + xpref2c * zcompapoc + xpref2d * zcompadi
               ! <CMOC OR 06/17/2014> Code trimming ! zfoodlim  = MAX( 0. , zfood - xthresh )
               ! <CMOC OR 06/17/2014> Code trimming ! zdenom    = zfoodlim / ( xkgraz + zfoodlim )
               ! <CMOC OR 06/17/2014> Code trimming ! zdenom2   = zdenom / ( zfood + rtrn )
               ! <CMOC OR 06/17/2014> Code trimming ! zgraze    = grazrat * zstep * tgfunc2(ji,jj,jk) * trn(ji,jj,jk,jpzoo) 

               ! <CMOC OR 06/17/2014> Code trimming ! zgrazp    = zgraze  * xpref2p * zcompaph  * zdenom2
               zgrazpcmoc= zstep   * rm_cmoc * zcompaph  * trn(ji,jj,jk,jpphy)  / ( kp_cmoc * 1e-6_wp * cnrr_cmoc * kp_cmoc * 1e-6_wp * cnrr_cmoc + trn(ji,jj,jk,jpphy) * trn(ji,jj,jk,jpphy) + rtrn ) * trn(ji,jj,jk,jpzoo) ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones 
! <CMOC OR 02/26/2014> 2 per day changed to 1.05 per day (actual model value - reduce the grazing pressure on phytoplankton )  ! zgraze  * xpref2p * zcompaph  * zdenom2 <CMOC OR 10/30/2013> nanophyto grazing
               ! <CMOC OR 06/17/2014> Code trimming ! zgrazm    = zgraze  * xpref2c * zcompapoc * zdenom2 
               ! <CMOC OR 06/17/2014> Code trimming ! zgrazsd   = zgraze  * xpref2d * zcompadi  * zdenom2 

               ! <CMOC OR 06/17/2014> Code trimming !  zgrazpf   = zgrazp  * trn(ji,jj,jk,jpnfe) / (trn(ji,jj,jk,jpphy) + rtrn)
               ! <CMOC OR 06/17/2014> Code trimming !  zgrazmf   = zgrazm  * trn(ji,jj,jk,jpsfe) / (trn(ji,jj,jk,jppoc) + rtrn)
               ! <CMOC OR 06/17/2014> Code trimming ! zgrazsf   = zgrazsd * trn(ji,jj,jk,jpdfe) / (trn(ji,jj,jk,jpdia) + rtrn)
               !
               ! <CMOC OR 06/17/2014> Code trimming ! zgraztot  = zgrazp  + zgrazm  + zgrazsd 
               ! <CMOC OR 06/17/2014> Code trimming !  zgraztotf = zgrazpf + zgrazsf + zgrazmf 

               ! Grazing by microzooplankton
               ! <CMOC OR 06/17/2014> Code trimming ! grazing(ji,jj,jk) = grazing(ji,jj,jk) + zgraztot

               !    Various remineralization and excretion terms
               !    --------------------------------------------
               ! <CMOC OR 06/17/2014> Code trimming !  zgrasrat  = zgraztotf / ( zgraztot + rtrn )
! <CMOC OR 05/26/2014> trimming PISCES code in p4zprod.F90  zncratio  = ( xpref2p * zcompaph * quotan(ji,jj,jk) &
! <CMOC OR 05/26/2014> trimming PISCES code in p4zprod.F90  &        + xpref2d * zcompadi * quotad(ji,jj,jk) + xpref2c * zcompapoc ) / ( zfood + rtrn )
               ! <CMOC OR 06/17/2014> Code trimming !  zepshert  = 0. ! <CMOC OR 05/26/2014> trimming PISCES code in p4zprod.F90  epsher * MIN( 1., zncratio )
               ! <CMOC OR 06/17/2014> Code trimming !  zepsherv  = zepshert * MIN( 1., zgrasrat / ferat3 )
               ! <CMOC OR 06/17/2014> Code trimming !  zgrafer   = zgraztot * MAX( 0. , ( 1. - unass ) * zgrasrat - ferat3 * zepshert ) 
               ! <CMOC OR 06/17/2014> Code trimming !  zgrarem   = zgraztot * ( 1. - zepsherv - unass )
               zgrapoc   = ( 1._wp - ga_cmoc ) * zgrazpcmoc ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones 
 !zgraztot * unass <CMOC OR 10/02/2013> detritus input from microzooplankton

               !  Update of the TRA arrays
               !  ------------------------
               ! <CMOC OR 06/17/2014> Code trimming !  zgrarsig  = zgrarem * sigma1
               ! <CMOC OR 06/17/2014> Code trimming !  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarsig
               ! <CMOC OR 06/17/2014> Code trimming !  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgrarsig
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + mzn_cmoc * zstep * trn(ji,jj,jk,jpzoo) ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones 
 !<CMOC OR 10/02/2013> contribution of microzooplankton exudation (has to be added) to dissolved nitrogen 
               ! <CMOC OR 06/17/2014> Code trimming !  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgrarem - zgrarsig
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - mzn_cmoc * zstep * trn(ji,jj,jk,jpzoo) ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones 
 !* 1.3_wp !- o2ut * zgrarsig <CMOC OR 11/14/2013>
               ! <CMOC OR 06/17/2014> Code trimming !  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgrafer
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zgrapoc  ! <CMOC OR 12/11/2013> this term is set back to its original CMOC form ! <    CMOC OR 12/09/2013> grazing term set to 0 
               ! <CMOC OR 06/17/2014> Code trimming !  tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zgraztotf * unass
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + mzn_cmoc * zstep * trn(ji,jj,jk,jpzoo)  ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones 
 !+ zgrarsig <CMOC OR 10/29/2013> grazing effect on DIC
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - mzn_cmoc * zstep * trn(ji,jj,jk,jpzoo) * ncrr_cmoc ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones 
 !+ rno3 * zgrarsig <CMOC OR 10/29/2013> grazing effect on TA
! <CMOC OR 06/17/2014> Code trimming ! #if defined key_kriest
! <CMOC OR 06/17/2014> Code trimming !                tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + zgrapoc * xkr_ddiat
! <CMOC OR 06/17/2014> Code trimming ! #endif
               !   Update the arrays TRA which contain the biological sources and sinks
               !   --------------------------------------------------------------------
               ! <CMOC OR 06/17/2014> Code trimming !  zmortz = ztortz + zrespz
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) +  ga_cmoc * zgrazpcmoc  - ( mzn_cmoc + mzd_cmoc ) * zstep * trn(ji,jj,jk,jpzoo) - mz2_cmoc * ncrr_cmoc * 1.0e3_wp * zstep * trn(ji,jj,jk,jpzoo) * trn(ji,jj,jk,jpzoo)  ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones 
 ! <CMOC OR 12/05/2013> fixed volume in 0.1_wp parameter (from m^-3 to L^-1)  !- zmortz + zepsherv * zgraztot <CMOC OR 10/02/2013> zooplankton mortality and growth
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgrazpcmoc !- zgrazp <CMOC OR 10/02/2013> nanophyto grazing
               ! <CMOC OR 06/17/2014> Code trimming ! tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazsd
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgrazpcmoc * trn(ji,jj,jk,jpnch)/(trn(ji,jj,jk,jpphy)+rtrn) !- zgrazp * trn(ji,jj,jk,jpnch)/(trn(ji,jj,jk,jpphy)+rtrn) <CMOC OR 10/02/2013> nanophyto grazing, chlorophyll removal
               ! <CMOC OR 06/17/2014> Code trimming ! tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zgrazsd * trn(ji,jj,jk,jpdch)/(trn(ji,jj,jk,jpdia)+rtrn)
               ! <CMOC OR 06/17/2014> Code trimming ! tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zgrazsd * trn(ji,jj,jk,jpdsi)/(trn(ji,jj,jk,jpdia)+rtrn)
               ! <CMOC OR 06/17/2014> Code trimming ! tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zgrazsd * trn(ji,jj,jk,jpdsi)/(trn(ji,jj,jk,jpdia)+rtrn)
               ! <CMOC OR 06/17/2014> Code trimming ! tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgrazpf
               ! <CMOC OR 06/17/2014> Code trimming ! tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazsf
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + mzd_cmoc * zstep * trn(ji,jj,jk,jpzoo) + mz2_cmoc * ncrr_cmoc * 1.0e3_wp * zstep * trn(ji,jj,jk,jpzoo) * trn(ji,jj,jk,jpzoo) ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones 
 ! <CMOC OR 12/20/2013> mzn=0.2 was used instead of mzd=0.05 for the linear coefficient (first term on the right of tra ) ! <CMOC OR 12/05/2013> fixed volume in 0.1_wp parameter (from m^-3 to L^-1)  !+ zmortz - zgrazm <CMOC OR 10/02/2013> contribution of microzooplankton mortality to detritus pool
               ! <CMOC OR 06/17/2014> Code trimming ! tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ferat3 * zmortz - zgrazmf
               ! <CMOC OR 06/17/2014> Code trimming ! zprcaca = xfracal(ji,jj,jk) * zgrazpcmoc ! <CMOC OR 01/07/2014> ! PISCES calcification ! * zgrazp
               !
               ! calcite production
               ! <CMOC OR 06/17/2014> Code trimming ! prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
               ! <CMOC OR 06/17/2014> Code trimming ! zprcaca = part * zprcaca
               ! <CMOC OR 06/17/2014> Code trimming !  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) ! <CMOC OR 02/19/2014> - zprcaca      ! <CMOC OR 01/07/2014> PISCES calcification
               ! <CMOC OR 06/17/2014> Code trimming !  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) ! <CMOC OR 02/19/2014> - 2. * zprcaca ! <CMOC OR 01/07/2014> PISCES calcification
               ! <CMOC OR 06/17/2014> Code trimming ! tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
! <CMOC OR 06/17/2014> Code trimming ! #if defined key_kriest
! <CMOC OR 06/17/2014> Code trimming !                tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) + ( zmortz - zgrazm ) * xkr_ddiat
! <CMOC OR 06/17/2014> Code trimming ! #endif
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

      ! <CMOC OR 06/17/2014> Code trimming !  NAMELIST/nampiszoo/ part, grazrat, resrat, mzrat, xpref2c, xpref2p, &
      ! <CMOC OR 06/17/2014> Code trimming !     &                xpref2d,  xthreshdia,  xthreshphy,  xthreshpoc, &
      ! <CMOC OR 06/17/2014> Code trimming !     &                xthresh, xkgraz, epsher, sigma1, unass
      NAMELIST/nampiszoo/ xthreshphy

      ! <CMOC OR 03/08/2014> CMOC namelist
      NAMELIST/namcmoczoo/ rm_cmoc, kp_cmoc, ga_cmoc, mzn_cmoc, mzd_cmoc, mz2_cmoc
      ! <CMOC OR 03/08/2014> CMOC namelist end 
      !!----------------------------------------------------------------------
      
      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampiszoo )

      REWIND( numcmoc )                     ! <CMOC OR 03/10/2014> ! read numcmoc, cmoczoo
      READ  ( numcmoc, namcmoczoo )
      
      
      IF(lwp) THEN                         ! control print
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

! <CMOC OR 06/17/2014> Code trimming !    INTEGER FUNCTION p4z_micro_alloc()
!      !!----------------------------------------------------------------------
!      !!                     ***  ROUTINE p4z_micro_alloc  ***
!      !!----------------------------------------------------------------------
!      ALLOCATE( grazing(jpi,jpj,jpk), STAT=p4z_micro_alloc )
!      IF( p4z_micro_alloc /= 0 ) CALL ctl_warn('p4z_micro_alloc : failed to allocate arrays.')
!
! <CMOC OR 06/17/2014> Code trimming !    END FUNCTION p4z_micro_alloc

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
