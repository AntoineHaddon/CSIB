MODULE p4zprod
   !!======================================================================
   !!                         ***  MODULE p4zprod  ***
   !! TOP :  Growth Rate of the two phytoplanktons groups 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
  !!          CMOC 1  !  2013-2015(O. Riche) phytoplankton growth and primary production
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_prod       :   Compute the growth Rate of the two phytoplanktons groups
   !!   p4z_prod_init  :   Initialization of the parameters for growth
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_prod         ! called in p4zbio.F90
   PUBLIC   p4z_prod_init    ! called in trcsms_pisces.F90
   REAL(wp) :: r1_rday                !: 1 / rday


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zprod.F90 3773 2013-02-07 11:06:58Z cbricaud $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_prod( kt , jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod  ***
      !!
      !! ** Purpose :   Compute the phytoplankton production depending on
      !!              light, temperature and nutrient availability
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt, jnt
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zfact
      REAL(wp) ::   ztn, zadap
      REAL(wp) ::   zprod
      REAL(wp) ::   zpislopen, ztheta 
      REAL(wp) ::   zrfact2
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zpislopead, zprbio, zprnch
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zprorca, zprochln
      ! <CMOC code OR 10/20/2015> nitrogen and light limitation functions
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zlimn,zliml
      ! <CMOC code OR 10/30/2015> etot is replaced by zetot = qsr * 0.43 and CMOC light attenuation
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zetot
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_prod')
      !
      !  Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj, jpk, zpislopead, zprbio, zprnch  )
      CALL wrk_alloc( jpi, jpj, jpk, zprorca, zprochln )
      CALL wrk_alloc( jpi, jpj, jpk, zlimn, zliml      )
      !
      ! <CMOC code OR 10/30/2015> etot is replaced by zetot = qsr * 0.43 and CMOC light attenuation
      CALL wrk_alloc( jpi, jpj, jpk, zetot             )
      !
      zetot   (:,:,:) = 0._wp
      zprorca (:,:,:) = 0._wp
      zprochln(:,:,:) = 0._wp
      zprbio  (:,:,:) = 0._wp
      zprnch  (:,:,:) = 0._wp
      zlimn   (:,:,:) = 0._wp
      zliml   (:,:,:) = 0._wp

!CDIR NOVERRCHK
         DO jk = 1, jpkm1
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi

		  ! <CMOC code OR 10/30/2015> etot is replaced by zetot = qsr * 0.43 and CMOC light attenuation
		  zetot(ji,jj,jk) = qsr(ji,jj) * 0.43_wp & 
		  !
		  &               * exp ( - ( (0.04 + 0.03 * trn(ji,jj,1,jpnch) * 1e6_wp) * fsdept(ji,jj,jk) ) )
		  !
		  !
                  ! <CMOC code OR 10/20/2015>
                  ! photosynthetic phytoplankton growth rate
                  ! -------------------------
                  
                  ! original PISCES condition for PAR
                  IF( zetot(ji,jj,jk) > 1.E-3 ) THEN
                      ztn    = tsn(ji,jj,jk,jp_tem) + 273.15_wp
                      ! ep_cmoc is in kJ mol^-1 and 8.31 is the ideal gas constant in J mol^-1 K^-1
                      zadap  = ep_cmoc * 1e+3_wp / 8.31_wp * ( 1._wp / ( ztn + rtrn ) - 1._wp / ( tvm_cmoc + 273.15_wp) )
                      zfact  = EXP ( -zadap )
                      ! zfact is the Arrhenius function, vm_cmoc the growth rate at 30oC in d^-1
                      ! zpislopead is the photosynthetic growth in s^-1
                      zpislopead (ji,jj,jk) = vm_cmoc * r1_rday * zfact
                      
                      ! phytoplankton photoacclimation used in light limitation
                      ! trn(...,jpnchl) / trn(...,jpphy) / 12. is theta in gChl per gC
                      ! ztheta is set to a maximum of thm_cmoc so as prevent appearance of light-saturation in case when zetot is small but trn(ji,jj,jk,jpphy) is 0
                      ztheta = MIN(thm_cmoc,trn(ji,jj,jk,jpnch)/(trn(ji,jj,jk,jpphy)*12._wp+rtrn))
                      zpislopen =  achl_cmoc * ztheta / ( zpislopead(ji,jj,jk) * rday  + rtrn )
                      ! zpislopead * rday is growth rate in d^-1 at temperature ToC as achl_cmoc is in d^-1

                      ! limitation functions
                      ! --------------------
                      ! light
                      zliml (ji,jj,jk) = 1.- EXP( -zpislopen  * zetot(ji,jj,jk) )
                      ! DIN
                      zlimn (ji,jj,jk) = trn(ji,jj,jk,jpno3) / ( kn_cmoc * 1e-6_wp * cnrr_cmoc + trn(ji,jj,jk,jpno3)+ rtrn )
                      ! iron is a constant and prescribed mask (xlimnfecmoc) see Zahariev et al 2008
                      ! update growth rate
                      zprbio(ji,jj,jk) = zpislopead(ji,jj,jk) * min ( zliml(ji,jj,jk) , zlimn(ji,jj,jk) , xlimnfecmoc(ji,jj) ) 
                      !  Computation of balanced chlorophyll based on balanced chlorophyll to carbon ratio; unit is gChl per molC
                      !  see Zahariev et al 2008 and Zahariev Environment Canada report (Canadian Model of Ocean Carbon v1.0)
                      !  balanced is defined as chlorophyll to carbon ratio in steady-state (Geider et al. 1996-1997)
                      !  p.40 Eq. 4.65 (note that in the report phytoplankton currency is N not C).
                      !  12._wp (gC molC^-1) to convert trn(...,jpphy) from moles to grams in the tra(...,jpchn) equations.
                      !  zprnch must be in gchl L^-1 per molC L^-1.
                      zprnch(ji,jj,jk) = 12._wp * thm_cmoc  * 2._wp    *  zpislopead(ji,jj,jk)  /  &
                      &                 ( 2._wp * zpislopead(ji,jj,jk) +                           &
                      &                  achl_cmoc * thm_cmoc  * zetot(ji,jj,jk) * r1_rday + rtrn )

                  ENDIF
               END DO
            END DO
         END DO

!CDIR NOVERRCHK
      DO jk = 1, jpkm1
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               IF( zetot(ji,jj,jk) > 1.E-3 ) THEN
               
                  ! Prognostic phytoplankton and chlorophyll tendencies
                  ! ---------------------------------------------------
                  !
                  ! phytoplankton production term over a time step
		  ! zprbio is photosynthetic growth rate in s^-1 (only)
                  zprorca(ji,jj,jk) =  zprbio(ji,jj,jk)  * trn(ji,jj,jk,jpphy) * rfact2
                  ! chlorophyll production term   over a time step
		  zprod =              zprbio(ji,jj,jk)  * trn(ji,jj,jk,jpnch) * rfact2
		  ! nudge chlorophyll back to balanced growth, Zahariev et al 2008
		  zprochln(ji,jj,jk) = zprod + (zprnch (ji,jj,jk) * trn(ji,jj,jk,jpphy) - &
                  &                             trn(ji,jj,jk,jpnch)                       &
                  &                            ) * itau_cmoc * r1_rday * rfact2
                     
                  ENDIF
               END DO
            END DO
         END DO

      !   Update the arrays TRA which contain the biological sources and sinks
      !   --------------------------------------------------------------------
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1 ,jpi
              
              tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zprorca(ji,jj,jk)
              tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprorca(ji,jj,jk)
              tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) + zprochln(ji,jj,jk)
              tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + zprorca(ji,jj,jk)
              tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprorca(ji,jj,jk)
              tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + ncrr_cmoc * zprorca(ji,jj,jk)
              tra(ji,jj,jk,jpdnt) = tra(ji,jj,jk,jpdnt) - zprorca(ji,jj,jk)
              
          END DO
        END DO
     END DO
         !
         zrfact2 = 1.e3 * rfact2r  ! conversion from mol L^-1 timestep^-1 into mol m^-3 s^-1
         IF( lk_iomput ) THEN
           IF( jnt == nrdttrc ) THEN
              CALL iom_put( "PPPHY"   , zprorca (:,:,:) * zrfact2 * tmask(:,:,:) )
              CALL iom_put( "Mumax"   , zpislopead  (:,:,:) * rday * tmask(:,:,:) )
              CALL iom_put( "LNnut"   , zlimn   (:,:,:) * tmask(:,:,:) )
              CALL iom_put( "LNFe"    , xlimnfecmoc (:,:) * tmask(:,:,1) )
              CALL iom_put( "LNlight" , zliml   (:,:,:) * tmask(:,:,:) )
              
           ENDIF
           
          ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zpislopead, zprbio, zprnch )
      CALL wrk_dealloc( jpi, jpj, jpk, zprorca, zprochln          )
      CALL wrk_dealloc( jpi, jpj, jpk, zlimn, zliml               )
      ! <CMOC code OR 10/30/2015> etot is replaced by zetot = qsr * 0.43 and CMOC light attenuation
      CALL wrk_dealloc( jpi, jpj, jpk, zetot                      )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_prod')
      !
   END SUBROUTINE p4z_prod


   SUBROUTINE p4z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      !
      ! <CMOC code OR 10/20/2015> CMOC namelist
      NAMELIST/namcmocphy/ achl_cmoc, thm_cmoc, tau_cmoc, itau_cmoc, ep_cmoc,     &
         &                 tvm_cmoc, vm_cmoc, kn_cmoc
      ! <CMOC code OR 10/20/2015> CMOC namelist end 
         !!----------------------------------------------------------------------

      REWIND( numcmoc )
      READ  ( numcmoc, namcmocphy )
      
      IF(lwp) THEN                         ! control print

         WRITE(numout,*) ' Namelist parameters for phytoplankton growth, namcmocphy'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Initial slope of the P-I curve            achl_cmoc    =', achl_cmoc
         WRITE(numout,*) '    Maximum chlorophyll relaxation time       thm_cmoc     =', thm_cmoc
         WRITE(numout,*) '    Chlorophyll relaxation time               tau_cmoc     =', tau_cmoc
         WRITE(numout,*) '    Inverse of chlorophyll relaxation time    itau_cmoc    =', itau_cmoc
         WRITE(numout,*) '    Activation energy for growth              ep_cmoc      =', ep_cmoc
         WRITE(numout,*) '    Reference ocean temperature               tvm_cmoc     =', tvm_cmoc
         WRITE(numout,*) '    Reference maximum photosynth. rate        vm_cmoc      =', vm_cmoc
         WRITE(numout,*) '    Half-saturation constant for N            kn_cmoc      =', kn_cmoc
         WRITE(numout,*) ' '

      ENDIF
      !
      r1_rday   = 1._wp / rday 
      !
   END SUBROUTINE p4z_prod_init


#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_prod                    ! Empty routine
   END SUBROUTINE p4z_prod
#endif 

   !!======================================================================
END MODULE  p4zprod
