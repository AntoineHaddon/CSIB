MODULE canoeprod
   !!======================================================================
   !!                         ***  MODULE canoeprod  ***
   !! TOP :  Growth Rate of the two phytoplanktons groups 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!          CanOE2  !  2022-23  (J. Christian, O. Riche) CanOE in NEMO4 phytoplankton growth and primary production
   !!                  !                                    use trc_opt and par_3bands
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_canoe'?                                     CanOE bio-model
   !!----------------------------------------------------------------------
   !!   canoe_prod       :   Compute the growth Rate of the two phytoplanktons groups
   !!   canoe_prod_init  :   Initialization of the parameters for growth
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   
   ! O Riche Sept 13th 2022
   ! sms_top needed?
   USE sms_top_canbgc     ! TOP Source Minus Sink variables
   USE sms_canoe          ! CanOE specific parameters declaration

   USE trc_closea_canbgc  ! tmask_bgc_closea
   USE trcopt_canbgc      ! PAR attenuation
   !
   USE canoetemp          ! CanOE temperature dependencies module
   !
   ! access par_1band array and requires trcsms_cmoc to call trc_opt_1band
   ! to update par_1band every time step
   
   USE prtctl          !  print control for debugging
   USE lib_mpp         !  ctl_stop on failed mem allocate check
   USE lib_fortran     !  access glob_sum function
   USE iom             !  I/O manager

   ! timing modules
   USE in_out_manager  ! nn_timing integer
   USE timing          ! *_timing subroutines

   IMPLICIT NONE
   PRIVATE

   PUBLIC   canoe_prod         ! called in trcsms_canoe.F90
   PUBLIC   canoe_prod_init    ! called in trcini_canoe.F90
   PUBLIC   trc_n2fx_canoe
   PUBLIC   trc_n2fx_init_canoe

   ! CanOE PP parameters
   ! these are hardwired parameters
   REAL(wp), PUBLIC ::  mw_c       = 12._wp            !:
   REAL(wp), PUBLIC ::  mw_n       = 14._wp            !:
   REAL(wp), PUBLIC ::  mw_fe      = 55.845_wp         !:
   !! these are input parameters in namelist_pisces
   REAL(wp), SAVE, PUBLIC ::  QNmax1     = 0.172_wp          !: Small phytoplankton max N quota
   REAL(wp), SAVE, PUBLIC ::  QNmin1     = 0.04_wp           !: Small phytoplankton min N quota
   REAL(wp), SAVE, PUBLIC ::  QNmax2     = 0.172_wp          !: Large phytoplankton max N quota
   REAL(wp), SAVE, PUBLIC ::  QNmin2     = 0.04_wp           !: Large phytoplankton min N quota
   REAL(wp), SAVE, PUBLIC ::  VCNref     = 0.6_wp            !: Reference rate of N uptake
   REAL(wp), SAVE, PUBLIC ::  QFemax1    = 93.075_wp         !: Small phytoplankton max Fe quota
   REAL(wp), SAVE, PUBLIC ::  QFemin1    = 4.65_wp           !: Small phytoplankton min Fe quota
   REAL(wp), SAVE, PUBLIC ::  QFemax2    = 69.8063_wp        !: Large phytoplankton max Fe quota
   REAL(wp), SAVE, PUBLIC ::  QFemin2    = 4.65_wp           !: Large phytoplankton min Fe quota
   REAL(wp), SAVE, PUBLIC ::  VCFref     = 79._wp            !: Reference rate of Fe uptake
   REAL(wp), SAVE, PUBLIC ::  PCref      = 3._wp             !: Reference rate of photosynthesis
   REAL(wp), SAVE, PUBLIC ::  alphachl   = 1.08_wp           !: Initial slope of P-E curve
   REAL(wp), SAVE, PUBLIC ::  kn1        = 0.1_wp            !: Small P half-saturation for NO3 uptake
   REAL(wp), SAVE, PUBLIC ::  ka1        = 0.05_wp           !: Small P half-saturation for NH4 uptake
   REAL(wp), SAVE, PUBLIC ::  kf1        = 100._wp           !: Small P half-saturation for Fe uptake
   REAL(wp), SAVE, PUBLIC ::  kn2        = 0.5_wp            !: Large P half-saturation for NO3 uptake
   REAL(wp), SAVE, PUBLIC ::  ka2        = 0.05_wp           !: Large P half-saturation for NH4 uptake
   REAL(wp), SAVE, PUBLIC ::  kf2        = 200._wp           !: Large P half-saturation for Fe uptake
   REAL(wp), SAVE, PUBLIC ::  thetamax   = 0.18_wp           !: Maximum chlorophyll/nitrogen ratio
   REAL(wp), SAVE, PUBLIC ::  eta        = 2._wp             !: Metabolic cost of biosynthesis
   REAL(wp), SAVE, PUBLIC ::  kexh       = 1.7_wp            !: exhudation of excess intracellular C
   ! CanOE N2-fixation parameters (namelist namcanoenfx)
   REAL(wp), SAVE, PUBLIC ::  nitrfix    = 2.25E-2_wp        !: Reference rate of dinitrogen fixation
   REAL(wp), SAVE, PUBLIC ::  kni        = 0.1_wp            !: DNF N inhibition parameter
   REAL(wp), SAVE, PUBLIC ::  diazolight = 50._wp            !: DNF irradiance dependence parameter
   REAL(wp), SAVE, PUBLIC ::  concfediaz = 100._wp           !: DNF iron concentration dependence parameter
   INTEGER, SAVE, PUBLIC ::   jk_max_dnf = 25                !: layer index for max depth of nitrogen fixation

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prmax    !: optimal production = f(temperature)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotan   !: proxy of N quota in Nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotad   !: proxy of N quota in diatomee
   
   REAL(wp) :: tpp                    !: Total primary production

   !!* Substitution
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: canoeprod.F90 3773 2013-02-07 11:06:58Z cbricaud $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE canoe_prod( kt , jnt , Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE canoe_prod  ***
      !!
      !! ** Purpose :   Compute the phytoplankton production depending on
      !!                light, temperature and nutrient availability
      !!
      !! ** Method  : - forward time integration (Euler or Leapfrog)
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt, jnt
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zfact, znanotot, zdiattot, zconctemp, zconctemp2
      REAL(wp) :: zratio, zmax, ztn, zadap
      REAL(wp) :: zlim, zprod, zproreg, zproreg2
      REAL(wp) :: zmxltst, zmxlday, zmaxday
      REAL(wp) :: zrum, zcodel, zargu, zval
      REAL(wp) :: zrfact2
! some local variables required by the revised model
      REAL(wp) :: phyc,phyn,phyfe,chl,Ni,Na,Fe,Tf,QN,qndep,qfedep,VCNmax,Alim,Nlim,VCN,QFe,VCFmax,VCF
      REAL(wp) :: PCmax,thetac,PCphot,rhochl,ei,xsphsyn, mwr_n2c, imw_n
      CHARACTER (len=25) :: charout
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zprdia, zprbio, zprdch, zprnch, zysopt   
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zprorca, zprorcad, zprofed, zprofen, zpronew, zpronewd
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zprocn, zprocd, zpronn, zprond
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('canoe_prod')
      !
      !
      IF( lwp ) THEN
        WRITE(numout,*)
        WRITE(numout,*), 'canoe_prod: compute phytoplankton and chlorophyl production'
        WRITE(numout,*), '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*)
        CALL FLUSH(numout)
      END IF
      !  Allocate temporary workspace
      ALLOCATE( zprdia(  jpi, jpj, jpk ),  zprbio(   jpi, jpj, jpk ) )
      ALLOCATE( zprdch(  jpi, jpj, jpk ),  zprnch(   jpi, jpj, jpk ) )
      ALLOCATE( zysopt(  jpi, jpj, jpk )                             ) 
      ALLOCATE( zprorca( jpi, jpj, jpk ),  zprorcad( jpi, jpj, jpk ) )
      ALLOCATE( zprofed( jpi, jpj, jpk ),  zprofen(  jpi, jpj, jpk ) )
      ALLOCATE( zpronew( jpi, jpj, jpk ),  zpronewd( jpi, jpj, jpk ) )
      ALLOCATE( zprocn(  jpi, jpj, jpk ),  zprocd(   jpi, jpj, jpk ) )
      ALLOCATE( zpronn(  jpi, jpj, jpk ),  zprond(   jpi, jpj, jpk ) ) 
      !
      zprorca (:,:,:) = 0._wp
      zprorcad(:,:,:) = 0._wp
      zprofed (:,:,:) = 0._wp
      zprofen (:,:,:) = 0._wp
      zprochln(:,:,:) = 0._wp
      zprochld(:,:,:) = 0._wp
      zpronew (:,:,:) = 0._wp
      zpronewd(:,:,:) = 0._wp
      xlimdn  (:,:,:) = 0._wp
      xlimdfe0(:,:,:) = 0._wp
      xlimnn  (:,:,:) = 0._wp
      xlimdn  (:,:,:) = 0._wp
      zprdia  (:,:,:) = 0._wp
      zprbio  (:,:,:) = 0._wp
      zprdch  (:,:,:) = 0._wp
      zprnch  (:,:,:) = 0._wp
      zysopt  (:,:,:) = 0._wp
      zprocn  (:,:,:) = 0._wp
      zprocd  (:,:,:) = 0._wp
      zpronn  (:,:,:) = 0._wp
      zprond  (:,:,:) = 0._wp
! precalculate some constant terms to minimize divisions
      mwr_n2c = mw_n/mw_c
      imw_n   =   1./mw_n
      WRITE(numout,*) "start loops"
      ! Computation of the various production terms 
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( par_3bands(ji,jj,jk) > 1.E-3 ) THEN
                      ztn  = ts(ji,jj,jk,jp_tem,Kmm)
                      phyc  = MAX(tr(ji,jj,jk,jrphy, Kbb),0.)*mw_c
                      phyn  = MAX(tr(ji,jj,jk,jrnn, Kbb) ,0.)*mw_n
                      phyfe = MAX(tr(ji,jj,jk,jrnfe, Kbb),0.)*mw_fe
                      chl   = MAX(tr(ji,jj,jk,jrnch, Kbb),0.)
                      Ni    = MAX(tr(ji,jj,jk,jqno3, Kbb),0.)
                      Na    = MAX(tr(ji,jj,jk,jrnh4, Kbb),0.)
                      Fe    = MAX(tr(ji,jj,jk,jrfer, Kbb),0.)                        ! Fe variables are in nmol m^-3, others in mmol m^-3
                      ei    = par_3bands(ji,jj,jk)*4.15                        ! convert to umol m^-2 s^-1

! this is modified from ~/mexfiles/vrm/fwd/bsource_vrm.f via chemo_2P2Z_gmk.F
! small phytoplankton

                      Tf = tgfuncp0(ji,jj,jk)
                      QN = MIN(QNmax1,phyn/(phyc+rtrn))
                      QN = MAX(QNmin1,QN)
                      qndep = MAX((QNmax1-QN)/(QNmax1-QNmin1),0.)                  ! in principle this should be nonegative but if roundoff makes it even slightly negative the exponent could go NaN
                      VCNmax = VCNref*Tf*qndep**0.05
                      Alim  =  Na/(ka1+Na)
                      Nlim  =  Ni/(kn1+Ni)
                      VCN = VCNmax*((1.-Alim)*Nlim+Alim)

                      QFe = MIN(QFemax1,phyfe/(phyc+rtrn))
                      QFe = MAX(QFemin1,QFe)
                      qfedep = MAX((QFemax1-QFe)/(QFemax1-QFemin1),0.) 
                      VCFmax = VCFref*Tf*qfedep**0.05
                      VCF = VCFmax*Fe/(kf1+Fe)

                      PCmax = PCref*Tf*MIN((QFe-QFemin1+rtrn*1.e6)/(QFemax1-QFemin1),(QN-QNmin1+rtrn)/(QNmax1-QNmin1))

                      PCmax = MAX(PCmax,1.0e-10)
                      thetac = MAX(chl/(phyc+rtrn),0.001)
                      PCphot = PCmax*(1.-EXP(-alphachl*ei*thetac/PCmax))
                      rhochl = thetamax*(PCphot/(alphachl*thetac*MAX(ei,0.001)))

! calculate excess intracellular C for exhudation
                      xsphsyn=(phyc/(phyn+rtrn)*mwr_n2c-rr_c2n)*phyn*imw_n
                      xsphsyn=MAX(xsphsyn,0.)

                      zprocn(ji,jj,jk) = (PCphot-eta*VCN)*tr(ji,jj,jk,jrphy, Kbb)*xstepb-kexh*xsphsyn*xstepb ! C production rate (in molar units)
                      zpronn(ji,jj,jk) = VCN/QN*tr(ji,jj,jk,jrnn, Kbb)*xstepb                                ! N uptake rate
                      zprofen(ji,jj,jk) = VCF/QFe*tr(ji,jj,jk,jrnfe, Kbb)*xstepb                             ! Fe uptake rate
                      zprochln(ji,jj,jk) = rhochl*VCN/thetac*tr(ji,jj,jk,jrnch, Kbb)*xstepb                  ! Chl production rate
                      zpronew(ji,jj,jk) = zpronn(ji,jj,jk)*(1.-Alim)*Nlim/(Alim+(1.-Alim)*Nlim+rtrn)     ! NO3 uptake
                      xlimnn(ji,jj,jk)   = 1.-qndep 
                      xlimnfe0(ji,jj,jk) = 1.-qfedep

! large phytoplankton

                      phyc = MAX(tr(ji,jj,jk,jrdia, Kbb),0.)*mw_c
                      phyn = MAX(tr(ji,jj,jk,jrdn, Kbb),0.)*mw_n
                      phyfe= MAX(tr(ji,jj,jk,jrdfe, Kbb),0.)*mw_fe
                      chl =  MAX(tr(ji,jj,jk,jrdch, Kbb),0.)

                      QN = MIN(QNmax2,phyn/(phyc+rtrn))
                      QN = MAX(QNmin2,QN)
                      qndep = MAX((QNmax2-QN)/(QNmax2-QNmin2),0.) 
                      VCNmax = VCNref*Tf*qndep**0.05
                      Alim  =  Na/(ka2+Na)
                      Nlim  =  Ni/(kn2+Ni)
                      VCN = VCNmax*((1.-Alim)*Nlim+Alim)

                      QFe = MIN(QFemax2,phyfe/(phyc+rtrn))
                      QFe = MAX(QFemin2,QFe)
                      qfedep = MAX((QFemax2-QFe)/(QFemax2-QFemin2),0.) 
                      VCFmax = VCFref*Tf*qfedep**0.05
                      VCF = VCFmax*Fe/(kf2+Fe)

                      PCmax = PCref*Tf*MIN((QFe-QFemin2+rtrn)/(QFemax2-QFemin2),(QN-QNmin2+rtrn)/(QNmax2-QNmin2))

                      PCmax = MAX(PCmax,1.0e-10)
                      thetac = MAX(chl/(phyc+rtrn),0.001)
                      PCphot = PCmax*(1.-EXP(-alphachl*ei*thetac/PCmax))
                      rhochl = thetamax*(PCphot/(alphachl*thetac*MAX(ei,0.001)))

                      xsphsyn=(phyc/(phyn+rtrn)*mwr_n2c-rr_c2n)*phyn*imw_n
                      xsphsyn=MAX(xsphsyn,0.)

                      zprocd(ji,jj,jk) = (PCphot-eta*VCN)*tr(ji,jj,jk,jrdia, Kbb)*xstepb-kexh*xsphsyn*xstepb       ! C production rate (in molar units)
                      zprond(ji,jj,jk) = VCN/QN*tr(ji,jj,jk,jrdn, Kbb)*xstepb                                     ! N uptake rate
                      zprofed(ji,jj,jk) = VCF/QFe*tr(ji,jj,jk,jrdfe, Kbb)*xstepb                                  ! Fe uptake rate
                      zprochld(ji,jj,jk) = rhochl*VCN/thetac*tr(ji,jj,jk,jrdch, Kbb)*xstepb                       ! Chl production rate
                      zpronewd(ji,jj,jk) = zprond(ji,jj,jk)*(1.-Alim)*Nlim/(Alim+(1.-Alim)*Nlim+rtrn)        ! NO3 uptake
                      xlimdn(ji,jj,jk)   = 1.-qndep 
                      xlimdfe0(ji,jj,jk) = 1.-qfedep 

               ENDIF
            END DO
         END DO
      END DO
      WRITE(numout,*) "xlimdfe0:",xlimdfe0
      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1 ,jpi
              zproreg  = zpronn(ji,jj,jk) - zpronew(ji,jj,jk)
              zproreg2 = zprond(ji,jj,jk) - zpronewd(ji,jj,jk)
              tr(ji,jj,jk,jqno3, Krhs) = tr(ji,jj,jk,jqno3, Krhs) - zpronew(ji,jj,jk) 
              tr(ji,jj,jk,jqno3, Krhs) = tr(ji,jj,jk,jqno3, Krhs) - zpronewd(ji,jj,jk)
              tr(ji,jj,jk,jrnh4, Krhs) = tr(ji,jj,jk,jrnh4, Krhs) - zproreg                       
              tr(ji,jj,jk,jrnh4, Krhs) = tr(ji,jj,jk,jrnh4, Krhs) - zproreg2
              tr(ji,jj,jk,jrphy, Krhs) = tr(ji,jj,jk,jrphy, Krhs) + zprocn(ji,jj,jk)
              tr(ji,jj,jk,jrnn, Krhs)  = tr(ji,jj,jk,jrnn, Krhs)  + zpronn(ji,jj,jk)
              tr(ji,jj,jk,jrnch, Krhs) = tr(ji,jj,jk,jrnch, Krhs) + zprochln(ji,jj,jk) 
              tr(ji,jj,jk,jrnfe, Krhs) = tr(ji,jj,jk,jrnfe, Krhs) + zprofen(ji,jj,jk)
              tr(ji,jj,jk,jrdia, Krhs) = tr(ji,jj,jk,jrdia, Krhs) + zprocd(ji,jj,jk) 
              tr(ji,jj,jk,jrdn, Krhs)  = tr(ji,jj,jk,jrdn, Krhs)  + zprond(ji,jj,jk) 
              tr(ji,jj,jk,jrdch, Krhs) = tr(ji,jj,jk,jrdch, Krhs) + zprochld(ji,jj,jk) 
              tr(ji,jj,jk,jrdfe, Krhs) = tr(ji,jj,jk,jrdfe, Krhs) + zprofed(ji,jj,jk) 
! O2 production equals DIC reduction + an additional nitrate term based on Laws 1991; this term is set to conserve O2 globally at steady state, i.e. 0.301887 = 2/rr_c2n where 2 mol O2 / mol N is the O2 sink to nitrification
              tr(ji,jj,jk,jqoxy, Krhs) = tr(ji,jj,jk,jqoxy, Krhs) + zprocn(ji,jj,jk)                       
              tr(ji,jj,jk,jqoxy, Krhs) = tr(ji,jj,jk,jqoxy, Krhs) + zprocd(ji,jj,jk)  
              tr(ji,jj,jk,jqoxy, Krhs) = tr(ji,jj,jk,jqoxy, Krhs) + 0.301887 * zpronew(ji,jj,jk) * rr_c2n 
              tr(ji,jj,jk,jqoxy, Krhs) = tr(ji,jj,jk,jqoxy, Krhs) + 0.301887 * zpronewd(ji,jj,jk) * rr_c2n   
              tr(ji,jj,jk,jrfer, Krhs) = tr(ji,jj,jk,jrfer, Krhs) - zprofen(ji,jj,jk) 
              tr(ji,jj,jk,jrfer, Krhs) = tr(ji,jj,jk,jrfer, Krhs) - zprofed(ji,jj,jk)
              tr(ji,jj,jk,jqdic, Krhs) = tr(ji,jj,jk,jqdic, Krhs) - zprocn(ji,jj,jk)*1.E-6
              tr(ji,jj,jk,jqdic, Krhs) = tr(ji,jj,jk,jqdic, Krhs) - zprocd(ji,jj,jk)*1.E-6
              tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) + zpronew(ji,jj,jk)*1.E-6 
              tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) + zpronewd(ji,jj,jk)*1.E-6 
              tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) - zproreg*1.E-6 
              tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) - zproreg2*1.E-6
          END DO
        END DO
     END DO
     WRITE(numout,*) "done loops"
     ! Total primary production per year
     tpp = tpp + glob_sum( 'canoe_prod' , ( zprorca(:,:,:) + zprorcad(:,:,:) ) * cvol(:,:,:) )

     IF( kt == nitend .AND. jnt == qnrdttrc ) THEN
       WRITE(numout,*) 'Total PP (Gtc) :'
       WRITE(numout,*) '-------------------- : ',tpp * 12. / 1.E12
       WRITE(numout,*) 
     ENDIF
     WRITE(numout,*) "before iom_puts"
      !
     zrfact2 = 1.e-3 * qfact2r  ! conversion from umol/L/timestep into mol/m3/s
     IF( jnt == qnrdttrc ) THEN
       CALL iom_put( "PPPHY"   , zprocn (:,:,:)  * zrfact2 * tmask_bgc_closea(:,:,:) )  ! primary production by nanophyto
       CALL iom_put( "PPPHY2"  , zprocd (:,:,:)  * zrfact2 * tmask_bgc_closea(:,:,:) )  ! primary production by diatom
       CALL iom_put( "PPNEWN"  , zpronew (:,:,:) * zrfact2 * tmask_bgc_closea(:,:,:) )  ! new primary production by nanophyto
       CALL iom_put( "PPNEWD"  , zpronewd(:,:,:) * zrfact2 * tmask_bgc_closea(:,:,:) )  ! new primary production by diatom
       CALL iom_put( "PFeD"    , zprofed (:,:,:) * zrfact2 * tmask_bgc_closea(:,:,:) )  ! biogenic iron production by diatom
       CALL iom_put( "PFeN"    , zprofen (:,:,:) * zrfact2 * tmask_bgc_closea(:,:,:) )  ! biogenic iron production by nanophyto
       CALL iom_put( "LNN"     , xlimnn  (:,:,:)           * tmask_bgc_closea(:,:,:) )  ! Nitrogen limitation term
       CALL iom_put( "LDN"     , xlimdn  (:,:,:)           * tmask_bgc_closea(:,:,:) )  ! Nitrogen limitation term
       !CALL iom_put( "LNFe"    , xlimnfe0 (:,:,:)          * tmask_bgc_closea(:,:,:) )  ! Iron limitation term
       !CALL iom_put( "LDFe"    , xlimdfe0 (:,:,:)          * tmask_bgc_closea(:,:,:) )  ! Iron limitation term
       CALL iom_put( "PAR"     , par_3bands (:,:,:)        * tmask_bgc_closea(:,:,:) )  ! Irradiance
     ENDIF
     !
     WRITE(numout,*) "Done iom_puts"
     IF( sn_cfctl%l_prttrc )   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('prod')")
        CALL prt_ctl_info(charout, cdcomp = 'top')
        CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask, clinfo=ctrcnm)
     ENDIF
     !
     DEALLOCATE( zprdia,  zprbio,   zprdch,  zprnch,  zysopt            ) 
     DEALLOCATE( zprorca, zprorcad, zprofed, zprofen, zpronew, zpronewd )
     DEALLOCATE( zprocn,  zprocd,   zpronn,  zprond                     ) 
     !
     IF( ln_timing )  CALL timing_stop('canoe_prod')
     !
   END SUBROUTINE canoe_prod


   SUBROUTINE canoe_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE canoe_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      !
      INTEGER ::   ios, ierr ! Local integers
      NAMELIST/namcanoeprod/ QNmax1, QNmin1, QNmax2, QNmin2, VCNref, QFemax1,     &
         &                   QFemin1, QFemax2, QFemin2, VCFref, PCref, alphachl,  &
         &                   kn1, ka1, kf1, kn2, ka2, kf2, thetamax, eta, kexh
      !
      !!----------------------------------------------------------------------

      REWIND( numnatp_refb )              ! Namelist namcanoeprod in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcanoeprod, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcanoeprod in reference namelist_canoe' )
      REWIND( numnatp_cfgb )              ! Namelist namcanoeprod in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcanoeprod, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcanoeprod in configuration namelist_canoe' )

      IF(lwm) WRITE( numonpb, namcanoeprod )   
    
      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' canoe_prod_init: namcanoeprod'         
         WRITE(numout,*) ' Namelist parameters for phytoplankton growth, nampisprod'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Small phytoplankton max N quota          QNmax1       =', QNmax1
         WRITE(numout,*) '    Small phytoplankton min N quota          QNmin1       =', QNmin1
         WRITE(numout,*) '    Large phytoplankton max N quota          QNmax2       =', QNmax2
         WRITE(numout,*) '    Large phytoplankton min N quota          QNmin2       =', QNmin2
         WRITE(numout,*) '    Reference rate of N uptake               VCNref       =', VCNref
         WRITE(numout,*) '    Small phytoplankton max Fe quota         QFemax1      =', QFemax1
         WRITE(numout,*) '    Small phytoplankton min Fe quota         QFemin1      =', QFemin1
         WRITE(numout,*) '    Large phytoplankton max Fe quota         QFemax2      =', QFemax2
         WRITE(numout,*) '    Large phytoplankton min Fe quota         QFemin2      =', QFemin2
         WRITE(numout,*) '    Reference rate of Fe uptake              VCFref       =', VCFref
         WRITE(numout,*) '    Reference rate of photosynthesis         PCref        =', PCref
         WRITE(numout,*) '    Initial slope of P-E curve               alphachl     =', alphachl
         WRITE(numout,*) '    Small P half-saturation for NO3 uptake   kn1          =', kn1
         WRITE(numout,*) '    Small P half-saturation for NH4 uptake   ka1          =', ka1
         WRITE(numout,*) '    Small P half-saturation for Fe uptake    kf1          =', kf1
         WRITE(numout,*) '    Large P half-saturation for NO3 uptake   kn2          =', kn2
         WRITE(numout,*) '    Large P half-saturation for NH4 uptake   ka2          =', ka2
         WRITE(numout,*) '    Large P half-saturation for Fe uptake    kf2          =', kf2
         WRITE(numout,*) '    Maximum chlorophyll/carbon ratio         thetamax     =', thetamax
         WRITE(numout,*) '    Metabolic cost of biosynthesis           eta          =', eta
         WRITE(numout,*) '    Exudation rate of excess intracellular C kexh         =', kexh
      ENDIF
      !
      tpp = 0._wp
      !
      ! Allocate mem to limitation functions
      ALLOCATE( xlimnfe0(jpi,jpj,jpk), xlimdfe0(jpi,jpj,jpk),       &
         &      xlimnn(jpi,jpj,jpk), xlimdn(jpi,jpj,jpk), zn2fix(jpi,jpj,jpk),       &
         &      STAT=ierr )
      !
      IF( ierr /= 0 ) CALL ctl_stop( 'STOP', 'canoe_prod_init : failed to allocate xlim* arrays' )
      !
      ! Allocate mem to chl-a production arrays
      ALLOCATE( zprochln(jpi,jpj,jpk) , zprochld(jpi,jpj,jpk) , STAT=ierr )
      !
      IF( ierr /= 0 ) CALL ctl_stop( 'STOP', 'canoe_prod_init : failed to allocate zprochl* arrays' )
      !
   END SUBROUTINE canoe_prod_init

   SUBROUTINE trc_n2fx_canoe( kt, jnt, Kbb, Kmm, Krhs )
      ! compute N2 fixation 
      !REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) :: zpar  ! any PAR array
      !
      !LOGICAL, OPTIONAL, INTENT(in) :: write_rhs_flag   ! 
      !LOGICAL                       :: write_rhs_flag0  ! 
      !
      INTEGER                       :: ji, jj, jk      ! nested loop indices
      INTEGER, INTENT(in)           :: kt, jnt ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices

      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zn2fixtot
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: znitrpot, zwork
      REAL(wp)   :: zrtn, zlim, zfact 
      !
      ALLOCATE( znitrpot(jpi, jpj, jpk), zwork(jpi, jpj, jpk) )
      ALLOCATE( zn2fixtot(jpi, jpj  ) )
      ! Nitrogen fixation
      ! ----------------------------------------------------------

      ! <CMOC code OR 10/15/2015> Initialization of CMOC arrays
      zn2fix   (:,:,:) = 0._wp
      zn2fixtot(:,:)   = 0._wp
      !
      ! hardwiring index for now because jk_eud_cmoc is in a namelist only read in CMOC runs
      DO jk = 1, jk_max_dnf
         DO jj = 1, jpj
            DO ji = 1, jpi
                   ! this is copied from CanESM5 p4z_sed
                   zlim = kni / (kni + tr(ji,jj,jk,jqno3, Kbb) + tr(ji,jj,jk,jrnh4, Kbb))        ! CMOC function for NO3-inhibition
                   znitrpot(ji,jj,jk) =  MAX( 0.e0, ( tgfuncp0(ji,jj,jk) - 0.773 ) ) * 1.962   &
                   &                 *  zlim * tr(ji,jj,jk,jrfer, Kbb) / ( concfediaz + tr(ji,jj,jk,jrfer, Kbb) ) &
                   &                 * ( 1.- EXP( -par_3bands(ji,jj,jk) / diazolight ) ) &
                   &                 * tmask_bgc_closea(ji,jj,jk)        ! open ocean / land mask
                   !&                 * oomask(ji,jj) * tmask_bgc_closea(ji,jj,jk)        ! open ocean / land mask
                   !
                   ! diagnostic output
                   zn2fix(ji,jj,jk) = znitrpot(ji,jj,jk) * nitrfix
                   ! column integral in mmol m^-2 d^-1
                   zn2fixtot(ji,jj) = zn2fixtot(ji,jj) + znitrpot(ji,jj,jk) * nitrfix * e3t(ji,jj,jk, Kmm)     
            END DO
         END DO
      END DO
      !     --------------------------------------------------------------------
      !     Update the arrays TRA which contain the biological sources and sinks
      !     --------------------------------------------------------------------
      !
      DO jk = 1, jk_max_dnf
         DO jj = 1, jpj
            DO ji = 1, jpi
                   zfact = znitrpot(ji,jj,jk) * nitrfix * xstepb
                   tr(ji,jj,jk,jrnh4, Krhs) = tr(ji,jj,jk,jrnh4, Krhs) + zfact
                   tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) + 1.e-6 * zfact
            END DO
         END DO
      END DO

      ! area-weighted total in TgN/y (note that this is the total for current tile only)
      WRITE(numout,*) 'DNF sum:', SUM(zn2fixtot(:,:)*e1t(:,:)*e2t(:,:))*365.*14.*1.e-15

      IF( lk_iomput ) THEN
         IF( jnt == qnrdttrc ) THEN
            ! nitrogen fixation in molN m^-3 s^-1 
            zwork(:,:,:)  =  zn2fix(:,:,:) * 0.001/rday * tmask_bgc_closea(:,:,:)
            CALL iom_put( "Nfix"   , zwork )
       ENDIF
      ENDIF
      !
      DEALLOCATE(znitrpot, zn2fixtot, zwork)
      !
   END SUBROUTINE trc_n2fx_canoe

   SUBROUTINE trc_n2fx_init_canoe
      !
      INTEGER ::   ios  
      ! 
      NAMELIST/namcanoenfx/ kni, concfediaz, nitrfix, diazolight, jk_max_dnf
      REWIND( numnatp_refb )              ! Namelist namcanoenfx in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcanoenfx, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcanoenfx in reference namelist_cmoc' )
      REWIND( numnatp_cfgb )              ! Namelist namcanoenfx in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcanoenfx, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcanoenfx in configuration namelist_cmoc' )
      !
      IF(lwm) WRITE( numonpb, namcanoenfx )      
      !
      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for dinitrogen fixation , namcanoenfx'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Reference rate of dinitrogen fixation          nitrfix =',    nitrfix
         WRITE(numout,*) '    DNF N inhibition parameter                     kni =', kni
         WRITE(numout,*) '    DNF irradiance dependence parameter            diazolight =', diazolight
         WRITE(numout,*) '    DNF iron concentration dependence parameter    concfediaz =', concfediaz
         WRITE(numout,*) '    layer index for max depth of nitrogen fixation jk_max_dnf =', jk_max_dnf
         WRITE(numout,*) ' '
      END IF   
      !
   END SUBROUTINE trc_n2fx_init_canoe

END MODULE  canoeprod

