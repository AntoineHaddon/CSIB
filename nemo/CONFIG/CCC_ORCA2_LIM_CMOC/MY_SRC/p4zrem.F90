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

   SUBROUTINE p4z_rem( kt, jnt ) ! <CMOC code OR 10/15/2015> jnt is to be i
                                 ! known to prevent iom_put('Nfix') to cause 
                                 ! an error (see below) 
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, jnt ! <CMOC code OR 10/15/2015> add the time 
                                       ! split index jnt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      CHARACTER (len=25) :: charout
      ! <CMOC code OR 10/15/2015> arrays for total water column remineralisation, 
      ! total euphotic zone nitrogen fixation, temporary array for DNF diagnostics, 
      ! pon flux (euphotic zone bottom) for PIC burial diagnostics, PIC flux at the 
      ! bottom, bottom POC
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zredettot, zn2fixtot, zwork, zfpon, zbpon
      REAL(wp), POINTER, DIMENSION(:,:  ) :: zbpoc, zdenittot 
      ! <CMOC code OR 10/15/2015> arrays for depth-dependent rates, zJNd is used to 
      !compute the balance between denitrification and nitrogen fixation
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zredet,    zn2fix,   zJNd
      REAL(wp) :: zdeup, zideup  ! <CMOC code OR 01/23/2016>  scale of the euphotic 
                                 ! defined according to jk_eud_cmoc
      REAL(wp), POINTER, DIMENSION(:) :: zdepw        ! computation of depths between t-grid cells.
      REAL(wp), POINTER, DIMENSION(:) :: zcalflxexp   ! exponential decay of calcite flux with depth
      REAL(wp) :: zcaldiv, zcalbotflx                 ! divergence of the calcite flux, bottom flx

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_rem')
      !
      CALL wrk_alloc( jpi, jpj,      zredettot, zn2fixtot, zwork , zfpon, zbpon, zbpoc, zdenittot) 
      CALL wrk_alloc( jpi, jpj, jpk, zredet,    zn2fix,    zJNd          )                       
      CALL wrk_alloc( jpk, zdepw, zcalflxexp )

      ! <CMOC code OR 10/15/2015> Initialization of CMOC arrays
       zredet   (:,:,:) = 0._wp
       zredettot(:,:)   = 0._wp
       zn2fix   (:,:,:) = 0._wp
       zn2fixtot(:,:)   = 0._wp
       ! <CMOC code OR 12/11/2015> Total denitrification diagnostics
       zdenittot(:,:)   = 0._wp

       zJNd     (:,:,:) = 0._wp
       zwork    (:,:)   = 0._wp
       zfpon    (:,:)   = 0._wp
       zbpon    (:,:)   = 0._wp
       zbpoc    (:,:)   = 0._wp

      ! <CMOC code OR 10/15/2015> remineralisation rate of detritus
      ! <CMOC code OR 12/17/2015> imposing the seafloor bathymetry 
      !                           as the bottom of the water column 
      !                           instead of the maximum level on 
      !                           the vertical grid
 
      DO ji = 1, jpi
          DO jj = 1, jpj

            DO jk = jk_eud_cmoc, mbkt(ji,jj) ! <CMOC code OR 01/23/2016> introduce 
                                             ! the jk_eud_cmoc index 
               ! one way to think of the global variable "xstep" is the time step in days
               zredet (ji,jj,jk) = reref_cmoc * xstep &
               &                 * exp ( -ed_cmoc * 1e3_wp / 8.31_wp *        &
               &                ( 1._wp / ( tsn(ji,jj,jk,jp_tem) + 273.15_wp  &
               &                 + rtrn ) - 1._wp / ( tvm_cmoc + 273.15_wp )  &
               &                )      )
            END DO

      ! <CMOC code OR 10/15/2015> N2-fixation and denitrification implementation
      ! <CMOC code OR 10/15/2015> remineralization integration over the water column, 
      !                           from 111 m (k level = 12), included, down to the bottom.
      ! <CMOC code OR 12/17/2015> imposing the seafloor bathymetry as the bottom i
      !                           of the water column instead of the maximum level 
      !                           on the vertical grid 
            DO jk = jk_eud_cmoc+1, mbkt(ji,jj)  
                zredettot(ji,jj) = zredettot(ji,jj) + zredet(ji,jj,jk)  &
                &                    *   trn(ji,jj,jk,jppoc)            &
                &                    * fse3t(ji,jj,jk)                  &
                &                    * tmask(ji,jj,jk)
            END DO

      ! <CMOC code OR 10/15/2015> N2 fixation integration over the surface layer, between k level 1 and k level 11

            ! <CMOC code OR 01/23/2016> Move DNF inside the conditional 
            ! branching for open ocean criterion (nk_bal_cmoc = 15)
            IF ( nk_bal_cmoc <= mbkt(ji,jj) ) THEN
            
            DO jk = 1, jk_eud_cmoc
   
                zn2fix (ji,jj,jk) = pnf_cmoc * cnrr_cmoc * 1e-12_wp / 3600._wp * rfact2 & ! reference rate
                !
                &                 * kn_cmoc * 1e-6_wp / ( kn_cmoc * 1e-6_wp                  &
                                  + trn(ji,jj,jk,jpno3) + rtrn) & ! N inhibition
                !
                &                 * qsr(ji,jj)*0.43_wp * exp ( - ( (0.04 + 0.03              &
                &                 * trn(ji,jj,1,jpnch) * 1e6_wp) * fsdept(ji,jj,jk) ) )      &
                &                 / inf_cmoc                                                 & ! ligh sensitivity
                !
                &                 * ( max(tsn(ji,jj,jk,jp_tem), tnfmi_cmoc ) - tnfmi_cmoc )   &
                &                 / ( tnfMa_cmoc - tnfmi_cmoc ) & ! temperature dependence
                !
                &                 * ( phinf_cmoc * exp( 1._wp ) * anf_cmoc * fsdept(ji,jj,jk) &
                &                 * exp ( -anf_cmoc * fsdept(ji,jj,jk) ) + phi0_cmoc ) ! diazotroph abundance dependence
                !
                ! total nitrogen fixation on the current 1/4 time step
                zn2fixtot(ji,jj) = zn2fixtot(ji,jj) + zn2fix(ji,jj,jk) * fse3t(ji,jj,jk)      &
                &                                    *  tmask(ji,jj,jk) 
 
   
      ! <CMOC code OR 10/15/2015> In the 2 loops below, compute the CMOC term of  
      ! N2-fixation/denitrification and scale denitrification according to 
      ! N2-fixation to have a balance between them at each grid point
      ! ---------------------------------------------------------------------
      
      ! <CMOC code OR 12/17/2015> ! imposing an open ocean criterion on DNF - denit term
      !
                zJNd(ji,jj,jk) =  zn2fix(ji,jj,jk)
      !
             END DO 

             DO jk = jk_eud_cmoc+1 , mbkt(ji,jj) 
                zJNd(ji,jj,jk)  = -zredet(ji,jj,jk) * trn(ji,jj,jk,jppoc) *   &
                &               zn2fixtot(ji,jj) / (zredettot(ji,jj) + rtrn)
! <CMOC code OR 12/11/2015> Total denitrification (negative at this stage); NOTE: since denitrification is 
! over the whole water column I apply the land mask at each depth instead of once using the surface 
! (as in the case of DNF diagnostics, see farther below)j) + rtrn)
                zdenittot(ji,jj) = zdenittot(ji,jj) +                         &
                &                       zJNd(ji,jj,jk) * fse3t(ji,jj,jk) *    &
                &                                        tmask(ji,jj,jk)      
             END DO
            ENDIF
          END DO 
      END DO      
 
      !     --------------------------------------------------------------------
      !     Update the arrays TRA which contain the biological sources and sinks
      !     --------------------------------------------------------------------
      DO jk = 1, jpkm1
         tra(:,:,jk,jppoc) = tra(:,:,jk,jppoc) - zredet(:,:,jk) * trn(:,:,jk,jppoc) 
         tra(:,:,jk,jpno3) = tra(:,:,jk,jpno3) + zredet(:,:,jk) * trn(:,:,jk,jppoc) & 
         &                                     +   zJNd(:,:,jk)
         tra(:,:,jk,jpoxy) = tra(:,:,jk,jpoxy) - zredet(:,:,jk) * trn(:,:,jk,jppoc)
         tra(:,:,jk,jpdic) = tra(:,:,jk,jpdic) + zredet(:,:,jk) * trn(:,:,jk,jppoc) 
         tra(:,:,jk,jptal) = tra(:,:,jk,jptal) - zredet(:,:,jk) * trn(:,:,jk,jppoc) * ncrr_cmoc
      END DO

      !     Calcite flux
      !     --------------------------------------------------------------------
      ! <CMOC code OR 10/15/2015> rain ratio at level jk, temperature is temperature in the 1st layer 
         xrcico(:,:) = rmcico_cmoc * exp(aci_cmoc * ( tsn(:,:,1,jp_tem)       &
         &                                - trcico_cmoc ) ) /                 &
         &                              (1._wp + exp(aci_cmoc *               &
         &                                (tsn(:,:,1,jp_tem) - trcico_cmoc))  &
         &                               + rtrn                               &
         &                              )

       DO ji = 1, jpi
          DO jj = 1, jpj
             ! Open ocean criterion based on ocean bottom depth; generally nk_bal_cmoc=15
             ! Use this to avoid calcite flux over areas shallower than the euphotic zone.
             IF ( nk_bal_cmoc <= mbkt(ji,jj) ) THEN
                ! Compute the depth of the euphotic zone. Note the sum of e3t is nearly, but not exactly
                ! equal to depw depths.
                zdeup = 0._wp
                DO jk =1, jk_eud_cmoc
                   zdeup = zdeup + fse3t(ji,jj,jk)
                ENDDO
                ! Inverse depth, for computing averages.
                zideup= 1._wp / zdeup

                ! Calculate the center points between t-grid. Should be the w-grid points, but there is
                ! a small offset.
                zdepw(:) = 0._wp
                DO jk = 2, jpk
                   zdepw(jk) = zdepw(jk-1) + fse3t(ji,jj,jk-1)
                ENDDO
                !WRITE(numout,*) 'depth comp', zdeup, zdepw(jk_eud_cmoc+1), fsdepw(ji,jj,jk_eud_cmoc+1)

                ! PIC export at the bottom of the euphotic zone based on Zahariev et al 2008 p.59
                zfpon(ji,jj) = xrcico(ji,jj) * wsbio3(ji,jj,jk_eud_cmoc) * xstep * trn(ji,jj,jk_eud_cmoc,jppoc) 

                ! Exponential decay of calcite flux with depth. 
                zcalflxexp(:) = 0._wp
                DO jk = jk_eud_cmoc+1, mbkt(ji,jj)+1
                    zcalflxexp(jk) = zfpon(ji,jj) * exp(-1._wp*(zdepw(jk)-zdeup) / dci_cmoc)
                ENDDO

                ! Bottom PIC flux into sediments, which is removed from the deepest layer and
                ! added back to the surface (below).
                zcalbotflx = zcalflxexp( mbkt(ji,jj)+1 )

                ! Set the bottom boundary condition on the calcite flux to zero (no flux to sediment),
                ! and deal with the sediment flux separately from the sinking in the sections below.
                zcalflxexp( mbkt(ji,jj)+1 ) = 0._wp
               
                ! Over the levels of the euphotic zone, remove the euphotic-zone averaged
                ! PIC flux (mol/m3) from each level. 
                DO jk =1, jk_eud_cmoc
                    tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) -                       &
                   &                         zfpon(ji,jj) * zideup

                    tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) -                       &
                   &                     2 * zfpon(ji,jj) * zideup
                END DO

                ! Below the euphotic zone; compute the divergence of the PIC flux
                ! and distribute it over the t-cell. No sinking flux through the bottom here.
                DO jk = jk_eud_cmoc+1, mbkt(ji,jj)
                   zcaldiv =  ( zcalflxexp(jk) - zcalflxexp(jk+1) ) / fse3t(ji,jj,jk)

                   tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) +         zcaldiv                      
                   tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 2._wp * zcaldiv                      
                END DO

                ! Do the bottom sedimentation of calcite. The sedimenting flux is added back
                ! to the surface layer (psuedo "river flux") for conservation.
                !tra(ji,jj,mbkt(ji,jj),jptal) = tra(ji,jj,mbkt(ji,jj),jptal) - zcalbotflx / fse3t(ji,jj,mbkt(ji,jj))
                !tra(ji,jj,1,jptal) = tra(ji,jj,1,jptal)  + zcalbotflx / fse3t(ji,jj, 1)
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
              ! <CMOC code OR 10/15/2015> 1.e+3_wp is to convert from L^-1 to m^-3
              !  (left in the sum line #119); the diagnostics has to be rescaled 
              ! to per second by dividing by rfact2.
              zwork(:,:)  =  zn2fixtot(:,:) * ncrr_cmoc * 1.e+3_wp * rfact2r * tmask(:,:,1)
              ! nitrogen fixation in molN m^-2 s^-1 
              CALL iom_put( "Nfix"   , zwork )         
              ! <CMOC code OR 12/11/2015> 1.e+3_wp is to convert from L^-1 to 
              ! m^-3 (left in the sum line #119); the diagnostics has to be 
              ! rescaled to per second by dividing by rfact2; NOTE: land mask 
              ! already taken into account
              zwork(:,:)  = -zdenittot(:,:) * ncrr_cmoc * 1.e+3_wp * rfact2r                
              CALL iom_put( "Denit"  , zwork )                                         ! denitrification in molN m^-2 s^-1 
              ! <CMOC code OR 12/11/2015> denitrification ! CALL iom_put( "BUPOC"  , wsbio3(:,:,11) /rday * zbpoc(:,:) * 1e+3_wp  )  ! POC burial flux
              ! <CMOC code OR 12/11/2015> denitrification ! CALL iom_put( "BUCALC" , zfpon(:,:) * 1e+3_wp * rfact2r * zbpon(:,:)  )  ! <CMOC code OR 12/11/2015> *rfact2r replaces /rfact2 ! PIC burial flux

         ENDIF
        ENDIF
      ENDIF

      CALL wrk_dealloc( jpi, jpj,      zredettot, zn2fixtot, zwork, zfpon, zbpon, zbpoc, zdenittot  ) ! <CMOC code OR 12/11/2015> Total denitrification
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
      NAMELIST/namcmocdeu/ jk_eud_cmoc, nk_bal_cmoc

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
         WRITE(numout,*) '    Index of the bottom of the euphotic zone jk_eud_cmoc =', jk_eud_cmoc
         WRITE(numout,*) '    Open ocean, min. number of vert. layers  nk_bal_cmoc =', nk_bal_cmoc

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
