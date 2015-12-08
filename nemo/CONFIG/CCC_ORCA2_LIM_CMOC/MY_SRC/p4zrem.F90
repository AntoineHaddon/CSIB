MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   PISCES Compute remineralization/scavenging of organic compounds
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!          CMOC 1  !  2013-2015(O. Riche) remineralization, PIC export and nitrogen fixation, based on Zahariev et al 2008
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_top'       and                                      TOP models
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/scavenging of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zopt          !  optical model
   USE prtctl_trc      !  print control for debugging
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
   USE iom             !
   USE fldread         !
   USE p4zsink         !
   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p4z_rem_init    ! called in trcsms_pisces.F90

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zrem.F90 3558 2012-11-14 19:15:05Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_rem( kt, jnt )! <CMOC code OR 10/15/2015> jnt is to be known to prevent iom_put('Nfix') to cause an error (see below) 
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, jnt ! <CMOC code OR 10/15/2015> add the time split index jnt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      CHARACTER (len=25) :: charout
      ! <CMOC code OR 10/15/2015> arrays for total water column remineralisation, total euphotic zone nitrogen fixation, temporary array for DNF diagnostics, pon flux (euphotic zone bottom) for PIC burial diagnostics, PIC flux at the bottom, bottom POC
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zredettot, zn2fixtot, zwork, zfpon, zbpon, zbpoc
      ! <CMOC code OR 10/15/2015> arrays for depth-dependent rates, zJNd is used to compute the balance between denitrification and nitrogen fixation
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zredet,    zn2fix,   zJNd
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_rem')
      !
      CALL wrk_alloc( jpi, jpj,      zredettot, zn2fixtot, zwork , zfpon, zbpon, zbpoc   )       
      CALL wrk_alloc( jpi, jpj, jpk, zredet,    zn2fix,    zJNd          )                       

      ! <CMOC code OR 10/15/2015> Initialization of CMOC arrays
       zredet   (:,:,:) = 0._wp
       zredettot(:,:)   = 0._wp
       zn2fix   (:,:,:) = 0._wp
       zn2fixtot(:,:)   = 0._wp

       zJNd     (:,:,:) = 0._wp
       zwork    (:,:)   = 0._wp
       zfpon    (:,:)   = 0._wp
       zbpon    (:,:)   = 0._wp
       zbpoc    (:,:)   = 0._wp

      ! <CMOC code OR 10/15/2015> remineralisation rate of detritus
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               !
               ! one way to think of the global variable "xstep" is the time step in days
               zredet (ji,jj,jk) = reref_cmoc * xstep &
               !
               &                 * exp ( -ed_cmoc * 1e3_wp / 8.31_wp * ( 1._wp / ( tsn(ji,jj,jk,jp_tem) + 273.15_wp + rtrn ) - 1._wp / ( tvm_cmoc + 273.15_wp ) ) )
               !
            END DO
         END DO
      END DO

      ! <CMOC code OR 10/15/2015> N2-fixation and denitrification implementation
      ! <CMOC code OR 10/15/2015> remineralization integration over the water column, from 111 m (k level = 12), excluded, down to the bottom.
      DO jk = 12, jpkm1
   
                zredettot(:,:) = zredettot(:,:) + zredet(:,:,jk) * trn(:,:,jk,jppoc)  * fse3t(:,:,jk)
            
      END DO
      ! <CMOC code OR 10/15/2015> N2 fixation integration over the surface layer, between k level 1 and k level 11
      DO jk = 1, 11
   
                zn2fix (:,:,jk) = pnf_cmoc * cnrr_cmoc * 1e-12_wp / 3600._wp * rfact2 & ! reference rate
                !
                &                 * kn_cmoc * 1e-6_wp / ( kn_cmoc * 1e-6_wp + trn(:,:,jk,jpno3) + rtrn) & ! N inhibition
                !
                &                 * qsr(:,:)*0.43_wp * exp ( - ( (0.04 + 0.03 * trn(:,:,1,jpnch) * 1e6_wp) * fsdept(:,:,jk) ) ) / inf_cmoc & ! ligh sensitivity
                !
                &                 * ( max(tsn(:,:,jk,jp_tem), tnfmi_cmoc ) - tnfmi_cmoc ) / ( tnfMa_cmoc - tnfmi_cmoc ) & ! temperature dependence
                !
                &                 * ( phinf_cmoc * exp( 1._wp ) * anf_cmoc * fsdept(:,:,jk) * exp ( -anf_cmoc * fsdept(:,:,jk) ) + phi0_cmoc ) ! diazotroph abundance dependence
                !
                zn2fixtot(:,:) = zn2fixtot(:,:) + zn2fix(:,:,jk) * fse3t(:,:,jk) ! total nitrogen fixation on the current 1/4 time step
 
      END DO 

      ! <CMOC code OR 10/15/2015> In the 2 loops below, compute the CMOC term of N2-fixation/denitrification and scale denitrification according to N2-fixation to have a balance between them at each grid point
      ! ---------------------------------------------------------------------
      DO jk = 1, 11 
      !
            zJNd(:,:,jk) =  zn2fix(:,:,jk)
      !
      END DO 

      DO jk = 12, jpkm1
      !
            zJNd(:,:,jk) = -zredet(:,:,jk) * trn(:,:,jk,jppoc) * zn2fixtot(:,:) / (zredettot(:,:) + rtrn)
      !
      END DO
      
      ! print mean trends (used for debugging)
       IF(ln_ctl)   THEN
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! <CMOC code OR 10/15/2015> POC remineralization
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zredet (ji,jj,jk) *  trn(ji,jj,jk,jppoc) 

            END DO
         END DO
      END DO
      
      ! print mean trends (used for debugging)
       IF(ln_ctl)   THEN
         WRITE(charout, FMT="('rem3')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      !     Calcite flux <CMOC code OR 10/15/2015>
      !     --------------------------------------------------------------------
      ! <CMOC code OR 10/15/2015> rain ratio at level jk, temperature is temperature in the 1st layer 
         xrcico(:,:) = rmcico_cmoc * exp ( aci_cmoc * ( tsn(:,:,1,jp_tem) - trcico_cmoc ) ) / ( 1._wp + exp ( aci_cmoc * ( tsn(:,:,1,jp_tem) - trcico_cmoc ) ) + rtrn )

      DO jk = 1, jpkm1
      !     --------------------------------------------------------------------
      !     Update the arrays TRA which contain the biological sources and sinks
      !     --------------------------------------------------------------------
         tra(:,:,jk,jpno3) = tra(:,:,jk,jpno3) + zredet (:,:,jk) *    trn(:,:,jk,jppoc) & 
         &                   + zJNd(:,:,jk)
         tra(:,:,jk,jpoxy) = tra(:,:,jk,jpoxy) - zredet (:,:,jk) *    trn(:,:,jk,jppoc)
         tra(:,:,jk,jpdic) = tra(:,:,jk,jpdic) + zredet (:,:,jk) *    trn(:,:,jk,jppoc) 
         tra(:,:,jk,jptal) = tra(:,:,jk,jptal) - zredet (:,:,jk) *    trn(:,:,jk,jppoc) * ncrr_cmoc
      END DO

       DO ji = 1, jpi
        DO jj = 1, jpj

	! <CMOC code OR 10/15/2015> Alkalinity and DIC balance due to PIC
	! <CMOC code OR 10/15/2015> open ocean criterion based on ocean bottom depth
	! k = 15 equivalent to z about 150 m
        IF ( 15 <= mbkt(ji,jj) ) THEN

         ! <CMOC code OR 10/15/2015> POC flux at the bottom of the euphotic zone taken as k = 11 as in the rest of the code  
         zfpon(ji,jj) = xrcico(ji,jj) * wsbio3(ji,jj,11) * xstep * trn(ji,jj,11,jppoc) ! PIC export at the bottom of the euphotic zone based on Zahariev et al 2008 p.59
         ! Diagnostic terms follow
         zbpon(ji,jj) = exp(-1* (fsdepw(ji,jj,mbkt(ji,jj)+1)-fsdepw(ji,jj,11)) /dci_cmoc) ! fraction of PIC left at the bottom of the water column (used by BUCALC diagnostics)
         zbpoc(ji,jj) = trn(ji,jj,mbkt(ji,jj),jppoc)   ! POC at the bottom at the water column (used by BUPOC diagnostics)
         ! <CMOC code OR 10/15/2015> fsdepw is the w grid, according to Nemo 23 the depth of the top and bottom of the layer where a tracer at level k is located, k is the top and k+1 is the bottom.
         
         ! <CMOC code OR 10/15/2015> Calcification in the euphotic zone
         ! PIC export below the euphotic zone is balanced by calcification and the associated loss of DIC and alkalinity at the surface, evenly distributed over the euphotic zone
         DO jk =1, 11
	    tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) -     zfpon(ji,jj) * ( 1 - zbpon(ji,jj) ) * ideup_cmoc
	    tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2 * zfpon(ji,jj) * ( 1 - zbpon(ji,jj) ) * ideup_cmoc
         END DO
         ! <CMOC code OR 10/15/2015> ( 1 - zbpon ) is the fraction of PIC export that went into the water column below the euphotic zone
         ! the extra term - exp(), which is in fact + zfpon * exp(), is a balance term to enforce long-term equilibrium between surface and bottom conditions of DIC and Akalinity  

         ! <CMOC code OR 10/15/2015> below the euphotic zone
         DO jk = 12, mbkt(ji,jj)
	    ! <CMOC code OR 10/15/2015> Calcite Dissolution, source of Alkalinity/DIC
	    ! source/sink is calculated as the finite difference between 2 consecutive depth levels
	    tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) +     zfpon(ji,jj) * ( exp(-1._wp*(fsdepw(ji,jj,jk)-fsdepw(ji,jj,11))/dci_cmoc) - exp(-1._wp*(fsdepw(ji,jj,jk+1)-fsdepw(ji,jj,11))/dci_cmoc) ) / fse3w(ji,jj,jk)
	    tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 2 * zfpon(ji,jj) * ( exp(-1._wp*(fsdepw(ji,jj,jk)-fsdepw(ji,jj,11))/dci_cmoc) - exp(-1._wp*(fsdepw(ji,jj,jk+1)-fsdepw(ji,jj,11))/dci_cmoc) ) / fse3w(ji,jj,jk)
	    ! <CMOC code OR 10/15/2015> Below the euphotic zone calcite dissolution dominates; the POC flux varies as zfpon*exp(-(z-110)/2700), zfpon being the flux at the bottom of the euphotic zone 
         END DO

        ENDIF

        END DO
       END DO

      ! print mean trends (used for debugging)
      IF(ln_ctl)   THEN
         WRITE(charout, FMT="('rem6')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      IF( ln_diatrc ) THEN  
        IF( lk_iomput ) THEN
        
	! <CMOC code OR 10/15/2015>
         IF( jnt == nrdttrc ) THEN
              zwork(:,:)  =  zn2fixtot(:,:) * ncrr_cmoc * 1.e+3_wp * rfact2r * tmask(:,:,1) ! <CMOC code OR 10/15/2015> 1.e+3_wp is to convert from L^-1 to m^-3 (left in the sum line #119); the diagnostics has to be rescaled to per second by dividing by rfact2.
              CALL iom_put( "Nfix"   , zwork )                                         ! nitrogen fixation in molN m^-2 s^-1 
              CALL iom_put( "BUPOC"  , wsbio3(:,:,11) /rday * zbpoc(:,:) * 1e+3_wp  )  ! POC burial flux
              CALL iom_put( "BUCALC" , zfpon(:,:) * 1e+3_wp / rfact2 * zbpon(:,:)  )   ! PIC burial flux

         ENDIF
        ENDIF
      ENDIF

      CALL wrk_dealloc( jpi, jpj,      zredettot, zn2fixtot, zwork, zfpon, zbpon, zbpoc  )
      CALL wrk_dealloc( jpi, jpj, jpk, zredet,    zn2fix,   zJNd         )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_rem')
      !
   END SUBROUTINE p4z_rem


   SUBROUTINE p4z_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------

      ! <CMOC code OR 10/15/2015> CMOC namelist
      NAMELIST/namcmocpoc/ ed_cmoc, reref_cmoc
      NAMELIST/namcmoccal/ rmcico_cmoc, trcico_cmoc, aci_cmoc, dci_cmoc
      NAMELIST/namcmocnfx/ phinf_cmoc, phi0_cmoc, anf_cmoc, pnf_cmoc, inf_cmoc, tnfMa_cmoc, tnfmi_cmoc
      NAMELIST/namcmocdeu/  deup_cmoc, ideup_cmoc

      REWIND( numcmoc )            
      READ  ( numcmoc, namcmocpoc )
      REWIND( numcmoc )            
      READ  ( numcmoc, namcmoccal )
      REWIND( numcmoc )            
      READ  ( numcmoc, namcmocnfx )
      REWIND( numcmoc )            
      READ  ( numcmoc, namcmocdeu )
      ! <CMOC code OR 10/15/2015> CMOC namelist end 
      !!----------------------------------------------------------------------
      
      ! control print
      IF(lwp) THEN

         WRITE(numout,*) ' Namelist parameters for remineralization, namcmocpoc'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Remineralisation rate of POC              reref_cmoc=', reref_cmoc
         WRITE(numout,*) '    Activation energy for remineralization    ed_cmoc   =', ed_cmoc   
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for calcite export  , namcmoccal'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Maximum rain ratio                      rmcico_cmoc =', rmcico_cmoc
         WRITE(numout,*) '    Rain ratio half-point temperature       trcico_cmoc =', trcico_cmoc
         WRITE(numout,*) '    Rain ratio scaling factor                  aci_cmoc =',    aci_cmoc
         WRITE(numout,*) '    CaCO3 redissolution depth scale            dci_cmoc =',    dci_cmoc
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for dinitrogen fix. , namcmocnfx'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Maximum ref. diazotroph concentration    phinf_cmoc =',  phinf_cmoc
         WRITE(numout,*) '    Surface ref. diazotroph concentration     phi0_cmoc =',   phi0_cmoc
         WRITE(numout,*) '    Inverse depth of diazotroph conc.max.      anf_cmoc =',    anf_cmoc
         WRITE(numout,*) '    Maximum ref. rate of dinitrogen fix.       pnf_cmoc =',    pnf_cmoc
         WRITE(numout,*) '    Maximum ref. dinitrogen fix surf. irr.     inf_cmoc =',    inf_cmoc
         WRITE(numout,*) '    Maximum ref. dinitrogen fix SST          tnfMa_cmoc =',  tnfMa_cmoc
         WRITE(numout,*) '    Minimum ref. dinitrogen fix SST          tnfmi_cmoc =',  tnfmi_cmoc
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters, Euphotic zone depth , namcmocdeu'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Depth of euphotic zone                    deup_cmoc =',   deup_cmoc
         WRITE(numout,*) '    Inverse of the depth of euphotic zone    ideup_cmoc =',  ideup_cmoc
      ENDIF
      !
   END SUBROUTINE p4z_rem_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_rem                    ! Empty routine
   END SUBROUTINE p4z_rem
#endif 

   !!======================================================================
END MODULE p4zrem
