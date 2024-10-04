MODULE cmocprod
   !!======================================================================
   !!                         ***  MODULE cmocprod  ***
   !! TOP :  Growth Rate of the two phytoplanktons groups 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!          CMOC 1  !  2013-2015(O. Riche) phytoplankton growth and primary production
   !!          CMOC 2  !  2022-23  (O. Riche) NEMO4 phytoplankton growth and primary production
   !!                  !            use trc_opt_1band and par_1band  
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_cmoc'?                                       CMOC bio-model
   !!----------------------------------------------------------------------
   !!   cmoc_prod       :   Compute the growth Rate of the two phytoplanktons groups
   !!   cmoc_prod_init  :   Initialization of the parameters for growth
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   
   ! O Riche Sept 13th 2022
   ! sms_top needed?
   USE sms_top_canbgc     !  TOP Source Minus Sink variables
   USE sms_cmoc           !  CMOC specific parameters declaration

   USE trc_closea_canbgc  !  tmask_bgc_closea
   USE trcopt_canbgc
   ! access par_1band array and requires trcsms_cmoc to call trc_opt_1band
   ! to update par_1band every time step
   
   USE prtctl      !  print control for debugging
   USE iom             !  I/O manager

   ! timing modules
   USE in_out_manager  ! nn_timing integer
   USE timing          ! *_timing subroutines

   IMPLICIT NONE
   PRIVATE

   PUBLIC   cmoc_prod         ! called in trcsms_cmoc.F90
   PUBLIC   cmoc_prod_init    ! called in trcini_cmoc.F90

   REAL(wp), SAVE :: r1_rday                !: 1 / rday
   

   !!* Substitution
!#  include "top_substitute.h90"
!#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: cmocprod.F90 3773 2013-02-07 11:06:58Z cbricaud $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE cmoc_prod( kt , jnt , Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE cmoc_prod  ***
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
      INTEGER  ::   ji, jj, jk
      INTEGER  ::   jn
      REAL(wp) ::   zfact
      REAL(wp) ::   ztn, zadap
      REAL(wp) ::   zprod
      REAL(wp) ::   zpislopen, ztheta 
      REAL(wp) ::   zrfact2
      CHARACTER (len=25) :: charout
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zpislopead, zprbio, zprnch
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zprorca, zprochln
      ! <CMOC code OR 10/20/2015> nitrogen and light limitation functions
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zlimn,zliml
      ! <CMOC code OR 10/30/2015> etot is replaced by zetot = qsr * 0.43 and CMOC light attenuation
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zetot
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('cmoc_prod')
      !
      !
      IF( lwp ) THEN
        WRITE(numout,*)
        WRITE(numout,*), 'cmoc_prod: compute phytoplankton and chlorophyl production'
        WRITE(numout,*), '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*)
        CALL FLUSH(numout)
      END IF
      !  Allocate temporary workspace
      ALLOCATE( zpislopead(jpi, jpj, jpk), zprbio(jpi, jpj, jpk), zprnch(jpi, jpj, jpk)  )
      ALLOCATE( zprorca(jpi, jpj, jpk), zprochln(jpi, jpj, jpk) )
      ALLOCATE( zlimn(jpi, jpj, jpk), zliml(jpi, jpj, jpk)      )
      !
      ! <CMOC code OR 10/30/2015> etot is replaced by zetot = qsr * 0.43 and CMOC light attenuation
      ALLOCATE( zetot(jpi, jpj, jpk) )
      !
      zetot   (:,:,:) = 0._wp
      zprorca (:,:,:) = 0._wp
      zprochln(:,:,:) = 0._wp
      zpislopead(:,:,:) = 0._wp
      zprbio  (:,:,:) = 0._wp
      zprnch  (:,:,:) = 0._wp
      zlimn   (:,:,:) = 1._wp
      zliml   (:,:,:) = 1._wp
      !
      DO jk = 1, jpkm1
      !
        DO jj = 1, jpj
        !
           DO ji = 1, jpi
            ! <CMOC code OR 10/30/2015> etot is replaced by zetot = qsr * 0.43 and CMOC light attenuation
            ! zetot(ji,jj,jk) = qsr(ji,jj) * 0.43_wp & 
            ! !
            ! &               * exp ( - ( (0.04 + 0.03 * tr(ji,jj,1,jqnch,Kbb) * 1e6_wp) * gdept(ji,jj,jk,Kmm) ) )
            !
            ! O Riche Sept 13th 2022
            ! use trc_opt_1band
            zetot(ji,jj,jk) = par_1band(ji,jj,jk)
              !
              !
              ! <CMOC code OR 10/20/2015>
              ! photosynthetic phytoplankton growth rate
              ! -------------------------
              ! 
              ! original CMOC condition for PAR
              IF( zetot(ji,jj,jk) > 1.E-3 ) THEN
                ztn    = ts(ji,jj,jk,jp_tem,Kmm) + 273.15_wp
                ! ep_cmoc is in kJ mol^-1 and 8.31 is the ideal gas constant in J mol^-1 K^-1
                zadap  = ep_cmoc * 1.e3_wp / 8.31_wp * ( 1._wp / ( ztn + rtrn ) - 1._wp / ( tvm_cmoc + 273.15_wp) )
                zfact  = EXP ( -zadap )
                ! zfact is the Arrhenius function, vm_cmoc the growth rate at 30oC in d^-1
                ! zpislopead is the photosynthetic growth in s^-1
                zpislopead (ji,jj,jk) = vm_cmoc * r1_rday * zfact
                !
                ! phytoplankton photoacclimation used in light limitation
                ! tr(...,jqnch,Kmm) / tr(...,jqphy,Kmm) / 12. is theta in gChl per gC
                ! ztheta is set to a maximum of thm_cmoc so as to prevent appearance of light-saturation in case zetot is small but trn(ji,jj,jk,jqphy) is 0
                ! tr(ji,jj,jk,jqphy,Kmm) is 0
                ztheta = MIN(thm_cmoc,tr(ji,jj,jk,jqnch,Kbb)/(tr(ji,jj,jk,jqphy,Kbb)*12._wp+rtrn))
                zpislopen =  MAX( achl_cmoc * ztheta / ( zpislopead(ji,jj,jk) * rday  + rtrn ), 0.)
                ! zpislopead * rday is growth rate in d^-1 at temperature ToC as achl_cmoc is in d^-1
                !
                ! limitation functions
                ! --------------------
                ! light
                zliml (ji,jj,jk) = 1.- EXP( -zpislopen  * zetot(ji,jj,jk) )
                ! DIN
                zlimn (ji,jj,jk) = tr(ji,jj,jk,jqno3,Kbb) / ( kn_cmoc * 1e-6_wp * cnrr_cmoc + tr(ji,jj,jk,jqno3,Kbb)+ rtrn )
                zlimn(ji,jj,jk) = MAX(zlimn(ji,jj,jk),0.)
                zlimn(ji,jj,jk) = MIN(zlimn(ji,jj,jk),1.)
                ! iron is a constant and prescribed mask (xlimnfecmoc) see Zahariev et al 2008
                ! update growth rate
                zprbio(ji,jj,jk) = zpislopead(ji,jj,jk) * min ( zliml(ji,jj,jk) , zlimn(ji,jj,jk) , xlimnfecmoc(ji,jj) ) 
                !  Computation of balanced chlorophyll based on balanced chlorophyll to carbon ratio; unit is gChl per molC
                !  see Zahariev et al 2008 and Zahariev Environment Canada report (Canadian Model of Ocean Carbon v1.0)
                !  balanced is defined as chlorophyll to carbon ratio in steady-state (Geider et al. 1996-1997)
                !  p.40 Eq. 4.65 (note that in the report phytoplankton currency is N not C).
                !  12._wp (gC molC^-1) to convert tr(...,jqphy,Kmm) from moles to grams in the tr(...,jqchn, Krhs) equations.
                !  zprnch must be in gchl L^-1 per molC L^-1.
                zprnch(ji,jj,jk) = 12._wp * thm_cmoc  * 2._wp    *  zpislopead(ji,jj,jk)  /  &
                &                 ( 2._wp * zpislopead(ji,jj,jk) +                           &
                &                  achl_cmoc * thm_cmoc  * zetot(ji,jj,jk) * r1_rday + rtrn )
                !
              ENDIF
              !
          END DO
        END DO
      END DO
      !
      DO jk = 1, jpkm1
        !
        DO jj = 1, jpj
          !
          DO ji = 1, jpi
            !
            IF( zetot(ji,jj,jk) > 1.E-3 ) THEN
              ! Prognostic phytoplankton and chlorophyll tendencies
              ! ---------------------------------------------------
              !
              ! phytoplankton production term over a time step
              ! zprbio is photosynthetic growth rate in s^-1 (only)
              zprorca(ji,jj,jk) =  zprbio(ji,jj,jk)  * tr(ji,jj,jk,jqphy,Kbb) * qfact2
              ! chlorophyll production term   over a time step
              zprod =              zprbio(ji,jj,jk)  * tr(ji,jj,jk,jqnch,Kbb) * qfact2
              ! nudge chlorophyll back to balanced growth, Zahariev et al 2008
              zprochln(ji,jj,jk) = zprod + (MAX(zprnch(ji,jj,jk)*tr(ji,jj,jk,jqphy,Kbb),0.) - &
              &                             MAX(tr(ji,jj,jk,jqnch,Kbb),0.)                       &
              &                            ) * itau_cmoc * r1_rday * qfact2                
              !
            ENDIF
            !
          END DO
        END DO
      END DO
      !
      !   Update the arrays TRA which contain the biological sources and sinks
      !   --------------------------------------------------------------------
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1 ,jpi
            !
            IF( zetot(ji,jj,jk) > 1.E-3 ) THEN
              tr(ji,jj,jk,jqno3, Krhs) = tr(ji,jj,jk,jqno3, Krhs) - zprorca(ji,jj,jk)
              tr(ji,jj,jk,jqphy, Krhs) = tr(ji,jj,jk,jqphy, Krhs) + zprorca(ji,jj,jk)
              tr(ji,jj,jk,jqnch, Krhs) = tr(ji,jj,jk,jqnch, Krhs) + zprochln(ji,jj,jk)
              tr(ji,jj,jk,jqoxy, Krhs) = tr(ji,jj,jk,jqoxy, Krhs) + zprorca(ji,jj,jk)
              tr(ji,jj,jk,jqdic, Krhs) = tr(ji,jj,jk,jqdic, Krhs) - zprorca(ji,jj,jk)
              tr(ji,jj,jk,jqtal, Krhs) = tr(ji,jj,jk,jqtal, Krhs) + ncrr_cmoc * zprorca(ji,jj,jk)
              ! O Riche Sept 14th can be uncommented or moved to TOP
              ! tr(ji,jj,jk,jqdnt, Krhs) = tr(ji,jj,jk,jqdnt, Krhs) - zprorca(ji,jj,jk)
              !
            ENDIF
          END DO
        END DO
     END DO
     !
     ! O Riche Sept 14th 2022
     ! Can be uncommented when diagnostics below
     ! have been added to xml definition files
     zrfact2 = 1.e3 * qfact2r  ! conversion from mol L^-1 timestep^-1 into mol m^-3 s^-1
     
     IF( lk_iomput ) THEN
       IF( jnt == qnrdttrc ) THEN
          CALL iom_put( "PPPHY"   , zprorca (:,:,:)   * tmask_bgc_closea(:,:,:) )
          CALL iom_put( "chlaP"   , zprochln(:,:,:)   * tmask_bgc_closea(:,:,:) )
          CALL iom_put( "photor"  , zprbio(:,:,:)     * tmask_bgc_closea(:,:,:) )
          CALL iom_put( "Mumax"   , zpislopead(:,:,:) * tmask_bgc_closea(:,:,:) )
          CALL iom_put( "thchl2C" , zprnch  (:,:,:)   * tmask_bgc_closea(:,:,:) )
          CALL iom_put( "LNnut"   , zlimn   (:,:,:)   * tmask_bgc_closea(:,:,:) )
          CALL iom_put( "LNFe"    , xlimnfecmoc (:,:) * tmask_bgc_closea(:,:,1) )
          CALL iom_put( "LNlight" , zliml   (:,:,:)   * tmask_bgc_closea(:,:,:) )
          CALL iom_put( "PARCMOC" , zetot   (:,:,:)   * tmask_bgc_closea(:,:,:) )
       ENDIF
      ENDIF
      
      IF( sn_cfctl%l_prttrc )   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_info(charout, cdcomp = 'top')
         CALL prt_ctl(tab4d_1=tr(:,:,:,:, Krhs), mask1=tmask_bgc_closea, clinfo=ctrcnm)
      ENDIF
      !
      DEALLOCATE( zpislopead, zprbio, zprnch )
      DEALLOCATE( zprorca, zprochln          )
      DEALLOCATE( zlimn, zliml               )
      ! <CMOC code OR 10/30/2015> etot is replaced by zetot = qsr * 0.43 and CMOC light attenuation
      DEALLOCATE( zetot                      )
      !
      IF( ln_timing )  CALL timing_stop('cmoc_prod')
      !
   END SUBROUTINE cmoc_prod


   SUBROUTINE cmoc_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE cmoc_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      !
      USE trcsrc_canbgc
      ! <CMOC code OR 09/13/2022> add extra integer to handle NEMO4-style namelists
      INTEGER ::   ios       ! Local integer
      ! <CMOC code OR 10/20/2015> CMOC namelist
      NAMELIST/namcmocphy/ achl_cmoc, thm_cmoc, tau_cmoc, itau_cmoc, ep_cmoc,     &
         &                 tvm_cmoc, vm_cmoc, kn_cmoc
      ! <CMOC code OR 10/20/2015> CMOC namelist end 
      !!----------------------------------------------------------------------
      NAMELIST/namcmocrr/ cnrr_cmoc, ncrr_cmoc

      REWIND( numnatp_refb )              ! Namelist namcmocphy in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcmocphy, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcmocphy in reference namelist_cmoc' )
      REWIND( numnatp_cfgb )              ! Namelist namcmocphy in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcmocphy, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcmocphy in configuration namelist_cmoc' )

      IF(lwm) WRITE( numonpb, namcmocphy )
      

      REWIND( numnatp_refb )              ! Namelist namcmocrr in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcmocrr, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcmocrr in reference namelist_cmoc' )
      REWIND( numnatp_cfgb )              ! Namelist namcmocrr in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcmocrr, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcmocrr in configuration namelist_cmoc' )
      IF(lwm) WRITE( numonpb, namcmocrr )
      
      IF(lwp) THEN                         ! control print

        WRITE(numout,*) ' Namelist parameters for phytoplankton growth, namcmocphy'
        WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*) '    Initial slope of the P-I curve            achl_cmoc    =', achl_cmoc
        WRITE(numout,*) '    Maximum chlorophyll relaxation time       thm_cmoc     =', thm_cmoc
        WRITE(numout,*) '    Chlorophyll relaxation time               tau_cmoc     =', tau_cmoc
        WRITE(numout,*) '    Inverse of chlorophyll relaxation time    itau_cmoc    =', itau_cmoc
        WRITE(numout,*) '    Activation energy for growth              ep_cmoc      =', ep_cmoc
        WRITE(numout,*) '    Reference ocean temperature               tvm_cmoc     =', tvm_cmoc
        WRITE(numout,*) '    Reference maximum photosynth. rate        vm_cmoc      =', vm_cmoc
        WRITE(numout,*) '    Half-saturation constant for N            kn_cmoc      =', kn_cmoc
        WRITE(numout,*) ' '
        WRITE(numout,*) ' Namelist parameters for phytoplankton growth, namcmocrr'
        WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*) '    Classic Redfield C:N ratio                cnrr_cmoc    =', cnrr_cmoc
        WRITE(numout,*) '    Classic Redfield N:C ratio                ncrr_cmoc    =', ncrr_cmoc

      ENDIF
      !
      r1_rday   = 1._wp / rday 
      !
      ! initialize iron mask with IC file saved in the 2d src arrays stack
      CALL trc_src2d( nittrc000, js2d_femask )       ! 1st time step / nittrc000 same nit000 in ocean physics
      xlimnfecmoc(:,:) = src2d_dta(:,:,js2d_femask)
      !
   END SUBROUTINE cmoc_prod_init

END MODULE  cmocprod
