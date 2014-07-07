MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   PISCES Compute remineralization/scavenging of organic compounds
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_top'       and                                      TOP models
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/scavenging of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!   p4z_rem_alloc :  Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zopt          !  optical model
   USE p4zche          !  chemical model
   USE p4zprod         !  Growth rate of the 2 phyto groups
! <CMOC OR 05/21/2014> Removal of p4zmeso module    USE p4zmeso         !  Sources and sinks of mesozooplankton
   USE p4zint          !  interpolation and computation of various fields
   USE prtctl_trc      !  print control for debugging
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
   USE iom             ! <CMOC OR 01/16/2014> iom input used by diagnostic file, in particular to import total N2-fixation
   USE fldread         !  <CMOC OR 01/16/2014> time interpolation 
   USE p4zsink         ! <CMOC OR 02/19/2014> sinking flux/FPON
   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p4z_rem_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_rem_alloc

   !! * Shared module variables
   REAL(wp), PUBLIC ::  xremik    = 0.3_wp     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::  xremip    = 0.025_wp   !: remineralisation rate of DOC
   REAL(wp), PUBLIC ::  nitrif    = 0.05_wp    !: NH4 nitrification rate 
   REAL(wp), PUBLIC ::  xsirem    = 0.003_wp   !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::  xsiremlab = 0.025_wp   !: fast remineralisation rate of POC 
   REAL(wp), PUBLIC ::  xsilab    = 0.31_wp    !: fraction of labile biogenic silica 
   REAL(wp), PUBLIC ::  xlam1     = 0.005_wp   !: scavenging rate of Iron 
   REAL(wp), PUBLIC ::  oxymin    = 1.e-6_wp   !: halk saturation constant for anoxia 
   REAL(wp), PUBLIC ::  ligand    = 0.6E-9_wp  !: ligand concentration in the ocean 


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr     !: denitrification array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitnh4   !: -    -    -    -   -


   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zrem.F90 3558 2012-11-14 19:15:05Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_rem( kt, jnt )    ! <CMOC OR 01/20/2014> jnt is to be known to prevent iom_put('Nfix') to cause an error (see below) 
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, jnt ! <CMOC OR 01/20/2014> add the time splitting index jnt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zremip, zremik , zlam1b, zdepbac2
      REAL(wp) ::   zkeq  , zfeequi, zsiremin, zfesatur
      REAL(wp) ::   zsatur, zsatur2, znusil, zdep, zfactdep
      REAL(wp) ::   zbactfer, zorem, zofer ! <CMOC OR 05/05/2014> Removal of GOC tracer ! zorem2, zofer
      REAL(wp) ::   zosil, zdenom1, zscave, zaggdfe, zcoag
#if ! defined key_kriest
      REAL(wp) ::   zofer2, zdenom ! <CMOC OR 05/05/2014> Removal of GOC tracer ! , zdenom2
#endif
      REAL(wp) ::   zlamfac, zonitr, zstep
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:  ) :: ztempbac, zredettot, zn2fixtot, zwork, zfpon, zbpon, zbpoc   ! <CMOC OR 03/14/2014> Bottom POC ! <CMOC OR 03/13/2014> add Burial of PIC ! <CMOC OR 02/25/2014> zfpon euphotic zone botton POC flux ! <CMOC OR 01/15/2014> zwork added for N2 fixation diagnosis ! <CMOC OR 12/11/2013> total water column nitrogen fixation and detritus remineralization 
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zdepbac, zolimi, zolimi2, zredet,    zn2fix,   zJNd          ! <CMOC OR 12/11/2013> zn2fix rate of nitrogen fixation ! <CMOC OR 10/02/2013> zredet rate of remineralization of detritus
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_rem')
      !
      ! Allocate temporary workspace
      CALL wrk_alloc( jpi, jpj,      ztempbac, zredettot, zn2fixtot, zwork , zfpon, zbpon, zbpoc   )        ! <CMOC OR 03/14/2014> Bottom POC ! <CMOC OR 02/25/2014> ! <CMOC OR 01/15/2014> ! <CMOC OR 12/11/2013> ! <CMOC OR 10/02/2013>
      CALL wrk_alloc( jpi, jpj, jpk, zdepbac, zolimi, zolimi2, zredet,    zn2fix,    zJNd          )        ! <CMOC OR 12/11/2013> zJNd: N2-fixation flux ! <CMOC OR 10/02/2013> zredet rate of remineralization of detritus

       ! Initialisation of temprary arrys
       zdepbac  (:,:,:) = 0._wp
       zolimi   (:,:,:) = 0._wp
       zolimi2  (:,:,:) = 0._wp
       ztempbac (:,:)   = 0._wp
       zredet   (:,:,:) = 0._wp !<CMOC OR 11/30/2013> zredet is initialized to zero 
       zredettot(:,:)   = 0._wp !<CMOC OR 12/11/2013> zredettot initialized to zero
       zn2fix   (:,:,:) = 0._wp !<CMOC OR 12/11/2013> zn2fix    initialized to zero
       zn2fixtot(:,:)   = 0._wp !<CMOC OR 12/11/2013> znefixtot initialized to zero
       IF(jnt == 1 ) THEN ! On time-split 1 re-initialize xn2fixdia N2 fixation diagnostic variable ! <CMOC OR 04/28/2014> fix N2-fixation diagnostics
          xn2fixdia(:,:)   = 0._wp !<CMOC OR 12/11/2013> znefixtot initialized to zero
       ENDIF
       zJNd     (:,:,:) = 0._wp !<CMOC OR 12/11/2013> zJNd      initialized to zero
       zwork    (:,:)   = 0._wp !<CMOC OR 01/15/2014>
       zfpon    (:,:)   = 0._wp !<CMOC OR 02/25/2014>
       zbpon    (:,:)   = 0._wp !<CMOC OR 03/13/2014> Burial of PIC
       zbpoc    (:,:)   = 0._wp !<CMOC OR 03/14/2014> Bottom POC    
      !  Computation of the mean phytoplankton concentration as
      !  a crude estimate of the bacterial biomass
      !   --------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep = MAX( hmld(ji,jj), heup(ji,jj) )
               IF( fsdept(ji,jj,jk) < zdep ) THEN
! <CMOC OR 05/21/2014> Removal of p4zmeso module                   zdepbac(ji,jj,jk) = MIN( 0.7 * ( trn(ji,jj,jk,jpzoo) + 2.* trn(ji,jj,jk,jpmes) ), 4.e-6 )
                  zdepbac(ji,jj,jk) = MIN( 0.7 * trn(ji,jj,jk,jpzoo) , 4.e-6 )
                  ztempbac(ji,jj)   = zdepbac(ji,jj,jk)
               ELSE
                  zdepbac(ji,jj,jk) = MIN( 1., zdep / fsdept(ji,jj,jk) ) * ztempbac(ji,jj)
               ENDIF
            END DO
         END DO
      END DO

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! denitrification factor computed from O2 levels
               nitrfac(ji,jj,jk) = MAX(  0.e0, 0.4 * ( 6.e-6  - trn(ji,jj,jk,jpoxy) )    &
                  &                                / ( oxymin + trn(ji,jj,jk,jpoxy) )  )
               nitrfac(ji,jj,jk) = MIN( 1., nitrfac(ji,jj,jk) )
            END DO
         END DO
      END DO

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               ! DOC ammonification. Depends on depth, phytoplankton biomass
               !     and a limitation term which is supposed to be a parameterization
               !     of the bacterial activity. 
               zremik = xremik * zstep / 1.e-6 * xlimbac(ji,jj,jk) * zdepbac(ji,jj,jk) 
               zremik = MAX( zremik, 2.e-4 * xstep )
               !     Ammonification in oxic waters with oxygen consumption
               !     -----------------------------------------------------
               zolimi (ji,jj,jk) = zremik * ( 1.- nitrfac(ji,jj,jk) ) * trn(ji,jj,jk,jpdoc) 
               zolimi2(ji,jj,jk) = MIN( ( trn(ji,jj,jk,jpoxy) - rtrn ) / o2ut, zolimi(ji,jj,jk) ) 
               !     Ammonification in suboxic waters with denitrification
               !     -------------------------------------------------------
               denitr(ji,jj,jk)  = MIN(  ( trn(ji,jj,jk,jpno3) - rtrn ) / rdenit,   &
                  &                     zremik * nitrfac(ji,jj,jk) * trn(ji,jj,jk,jpdoc)  )
               !
               zolimi (ji,jj,jk) = MAX( 0.e0, zolimi (ji,jj,jk) )
               zolimi2(ji,jj,jk) = MAX( 0.e0, zolimi2(ji,jj,jk) )
               denitr (ji,jj,jk) = MAX( 0.e0, denitr (ji,jj,jk) )
               !
               zredet (ji,jj,jk) = reref_cmoc * zstep * exp ( -ed_cmoc * 1e3_wp / 8.31_wp * ( 1._wp / ( tsn(ji,jj,jk,jp_tem) + 273.15_wp + rtrn ) - 1._wp / ( tvm_cmoc + 273.15_wp ) ) ) ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones ! <CMOC OR 10/03/2013> CMOC remineralisation rate of detritus remineralization
               !
            END DO
         END DO
      END DO

      ! <CMOC OR 12/11/2013> CMOC N2-fixation and denitrification implementation

      DO jk = 12, jpkm1 ! <CMOC OR 12/11/2013> remineralization integration over the water column, from 100 m (k level = 11), excluded, down to the bottom.
   
                zredettot(:,:) = zredettot(:,:) + zredet(:,:,jk) * trn(:,:,jk,jppoc) * fse3t(:,:,jk) ! <CMOC OR 12/11/2013> compute the total remineralization over the water column below 100 m, i.e. k>=11
            
      END DO

      DO jk = 1, 11 ! <CMOC OR 12/11/2013> N2 fixation integration over the surface layer, between k level 1 and k level 11
   
                zn2fix (:,:,jk) = pnf_cmoc * cnrr_cmoc * 1e-12_wp * 1._wp / 3600._wp * rfact2 & ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones ! <CMOC OR 12/11/2013> pNfref x biology timestep (1440 s)
                !
                &                 * kn_cmoc * 1e-6_wp / ( kn_cmoc * 1e-6_wp + trn(:,:,jk,jpno3) + rtrn) & ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones ! <CMOC OR 12/11/2013> KN/(KN+N)
                !
                &                 * etot(:,:,jk) / inf_cmoc & ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones ! <CMOC OR 12/11/2013>  Ipar(0)/INfmax * exp ( -( kw + kchl * chl ) * z) ! <CMOC OR 12/18/2013> Ipar(0)/INfmax * exp ( -( kw + kchl * chl ) * z) replaced by Ipar_Pisces(z)/INfmax
                !
                ! <CMOC OR 12/18/2013>  &                 * ( 0.04_wp + 0.03_wp * trn(ji,jj,jk,jpnch) * 1.e6_wp ) ) & ! <CMOC OR 12/11/2013> CMOC light extinction coefficient
                !
                &                 * ( max(tsn(:,:,jk,jp_tem), tnfmi_cmoc ) - tnfmi_cmoc ) / ( tnfMa_cmoc - tnfmi_cmoc ) & ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones ! <CMOC OR 02/18/2014> Set the temperature dependence to be 0 if temperature is equal to or smaller than 20oC  ! <CMOC OR 12/11/2013> (T-TNfmin)/(TNfmax-TNfmin)
                !
                &                 * ( phinf_cmoc * exp( 1._wp ) * anf_cmoc * fsdept(:,:,jk) * exp ( -anf_cmoc * fsdept(:,:,jk) ) + phi0_cmoc ) ! <CMOC OR 03/10/2014> "softwired" parameters replace hardwired ones ! <CMOC OR 12/11/2013> trichomes concentration ! <CMOC OR 12/11/2013> CMOC nitrogen fixation rate 
                !
                zn2fixtot(:,:) = zn2fixtot(:,:) + zn2fix(:,:,jk) * fse3t(:,:,jk) ! <CMOC OR 12/11/2013> compute the total nitrogen fixation rate over the water column down to 100 m, i.e. k<=10
 
      END DO 

      ! <CMOC OR 12/11/2013> In the 2 loops below, compute the CMOC term of N2-fixation/Denitrification 
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


      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               !    NH4 nitrification to NO3. Ceased for oxygen concentrations
               !    below 2 umol/L. Inhibited at strong light 
               !    ----------------------------------------------------------
               zonitr  =nitrif * zstep * trn(ji,jj,jk,jpnh4) / ( 1.+ emoy(ji,jj,jk) ) * ( 1.- nitrfac(ji,jj,jk) ) 
               denitnh4(ji,jj,jk) = nitrif * zstep * trn(ji,jj,jk,jpnh4) * nitrfac(ji,jj,jk) 
               !   Update of the tracers trends
               !   ----------------------------
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) !- zonitr - denitnh4(ji,jj,jk) <CMOC OR 11/08/2013> Testing if the source of NO3 is linked to NH4
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) ! + zonitr - rdenita * denitnh4(ji,jj,jk) <CMOC OR 10/02/2013> There is no nitrification nor anammox reaction in CMOC
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) !- o2nit * zonitr <CMOC OR 10/04/2013> No nitrification in CMOC for the moment
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) !- 2 * rno3 * zonitr + rno3 * ( rdenita - 1. ) * denitnh4(ji,jj,jk) <CMOC OR 10/02/2013> TA is part of CMOC chemistry component 
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               !    Bacterial uptake of iron. No iron is available in DOC. So
               !    Bacteries are obliged to take up iron from the water. Some
               !    studies (especially at Papa) have shown this uptake to be significant
               !    ----------------------------------------------------------
               zdepbac2 = zdepbac(ji,jj,jk) * zdepbac(ji,jj,jk)
! <CMOC OR 05/26/2014> trimming PISCES code in p4zprod.F90                 zbactfer = 20.e-6 * rfact2 * prmax(ji,jj,jk)                                 &
! <CMOC OR 05/26/2014> trimming PISCES code in p4zprod.F90     &              * trn(ji,jj,jk,jpfer) / ( 5E-10 + trn(ji,jj,jk,jpfer) )    &
! <CMOC OR 05/21/2014>  Removal of p4zmeso !              &              * zdepbac2 / ( xkgraz2 + zdepbac(ji,jj,jk) )               &
! <CMOC OR 05/26/2014> trimming PISCES code in p4zprod.F90     &              * zdepbac2 / (  zdepbac(ji,jj,jk) )               &
! <CMOC OR 05/26/2014> trimming PISCES code in p4zprod.F90     &              * ( 0.5 + SIGN( 0.5, trn(ji,jj,jk,jpfer) -2.e-11 )  )

               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) ! <CMOC OR 05/26/2014> trimming PISCES code in p4zprod.F90 - zbactfer
#if defined key_kriest
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) ! <CMOC OR 05/26/2014> trimming PISCES code in p4zprod.F90+ zbactfer
#else
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) ! <CMOC OR 05/26/2014> trimming PISCES code in p4zprod.F90 + zbactfer
#endif
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem2')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               !    POC disaggregation by turbulence and bacterial activity. 
               !    -------------------------------------------------------------
               zremip = xremip * zstep * tgfunc(ji,jj,jk) * ( 1.- 0.7 * nitrfac(ji,jj,jk) ) 

               !    POC disaggregation rate is reduced in anoxic zone as shown by
               !    sediment traps data. In oxic area, the exponent of the martin s
               !    law is around -0.87. In anoxic zone, it is around -0.35. This
               !    means a disaggregation constant about 0.5 the value in oxic zones
               !    -----------------------------------------------------------------
               zorem  = zremip * trn(ji,jj,jk,jppoc)
               zofer  = zremip * trn(ji,jj,jk,jpsfe)
#if ! defined key_kriest
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! zorem2 = zremip * trn(ji,jj,jk,jpgoc)
               zofer2 = zremip * trn(ji,jj,jk,jpbfe)
! <CMOC OR 05/05/2014> Removal of GOC tracer ! #else
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! zorem2 = zremip * trn(ji,jj,jk,jpnum)
#endif

               !  Update the appropriate tracers trends
               !  -------------------------------------

               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zorem
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zofer
#if defined key_kriest
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zorem
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! tra(ji,jj,jk,jpnum) = tra(ji,jj,jk,jpnum) - zorem2
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zofer
#else
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zredet (ji,jj,jk) *  trn(ji,jj,jk,jppoc) ! <CMOC OR 12/10/2013> ! missing POC removal found by Neil + zorem2 - zorem <CMOC OR 10/02/2013> POC is D, detritus, in CMOC
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) - zorem2
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zofer2 - zofer
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) - zofer2
#endif

            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem3')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               !     Remineralization rate of BSi depedant on T and saturation
               !     ---------------------------------------------------------
               zsatur   = ( sio3eq(ji,jj,jk) - trn(ji,jj,jk,jpsil) ) / ( sio3eq(ji,jj,jk) + rtrn )
               zsatur   = MAX( rtrn, zsatur )
               zsatur2  = zsatur * ( 1. + tsn(ji,jj,jk,jp_tem) / 400.)**4
               znusil   = 0.225  * ( 1. + tsn(ji,jj,jk,jp_tem) / 15.) * zsatur + 0.775 * zsatur2**9.25
               zdep     = MAX( hmld(ji,jj), heup(ji,jj) ) 
               zdep     = MAX( 0., fsdept(ji,jj,jk) - zdep )
               zfactdep = xsilab * EXP(-( xsiremlab - xsirem ) * zdep / wsbio2 )
               zsiremin = ( xsiremlab * zfactdep + xsirem * ( 1. - zfactdep ) ) * zstep * znusil
               zosil    = zsiremin * trn(ji,jj,jk,jpgsi)
               !
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) - zosil
               tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) + zosil
               !
            END DO
         END DO
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem4')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      zfesatur = ligand
!CDIR NOVERRCHK
      DO jk = 1, jpkm1
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zstep   = xstep
# if defined key_degrad
               zstep = zstep * facvol(ji,jj,jk)
# endif
               !  Compute de different ratios for scavenging of iron
               !  --------------------------------------------------

#if  defined key_kriest
               zdenom1 = trn(ji,jj,jk,jppoc) / &
           &           ( trn(ji,jj,jk,jppoc) + trn(ji,jj,jk,jpgsi) + trn(ji,jj,jk,jpcal) + rtrn )
#else
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! zdenom = 1. / ( trn(ji,jj,jk,jppoc) + trn(ji,jj,jk,jpgoc) + trn(ji,jj,jk,jpgsi) + trn(ji,jj,jk,jpcal) + rtrn )
               zdenom = 1. / ( trn(ji,jj,jk,jppoc) + trn(ji,jj,jk,jpgsi) + trn(ji,jj,jk,jpcal) + rtrn ) ! <CMOC OR 05/05/2014> Removal of GOC tracer ! 
               zdenom1 = trn(ji,jj,jk,jppoc) * zdenom
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! zdenom2 = trn(ji,jj,jk,jpgoc) * zdenom
#endif
               !  scavenging rate of iron. this scavenging rate depends on the load in particles
               !  on which they are adsorbed. The  parameterization has been taken from studies on Th
               !     ------------------------------------------------------------
               zkeq = fekeq(ji,jj,jk)
               zfeequi = ( -( 1. + zfesatur * zkeq - zkeq * trn(ji,jj,jk,jpfer) )               &
                  &        + SQRT( ( 1. + zfesatur * zkeq - zkeq * trn(ji,jj,jk,jpfer) )**2       &
                  &               + 4. * trn(ji,jj,jk,jpfer) * zkeq) ) / ( 2. * zkeq )

#if defined key_kriest
               zlam1b = 3.e-5 + xlam1 * (  trn(ji,jj,jk,jppoc)                   &
                  &                      + trn(ji,jj,jk,jpcal) + trn(ji,jj,jk,jpgsi)  ) * 1.e6
#else
               zlam1b = 3.e-5 + xlam1 * (  trn(ji,jj,jk,jppoc) & ! <CMOC OR 05/05/2014> Removal of GOC tracer ! + trn(ji,jj,jk,jpgoc)   &
                  &                      + trn(ji,jj,jk,jpcal) + trn(ji,jj,jk,jpgsi)  ) * 1.e6
#endif
               zscave = zfeequi * zlam1b * zstep

               !  Increased scavenging for very high iron concentrations
               !  found near the coasts due to increased lithogenic particles
               !  and let say it is unknown processes (precipitation, ...)
               !  -----------------------------------------------------------
               zlam1b  = xlam1 * MAX( 0.e0, ( trn(ji,jj,jk,jpfer) * 1.e9 - 1. ) )
               zcoag   = zfeequi * zlam1b * zstep
               zlamfac = MAX( 0.e0, ( gphit(ji,jj) + 55.) / 30. )
               zlamfac = MIN( 1.  , zlamfac )
               zdep    =  MIN(1., 1000. / fsdept(ji,jj,jk) )
#if ! defined key_kriest
               zlam1b = (  80.* ( trn(ji,jj,jk,jpdoc) + 35.e-6 )                           &
                  &     + 698.*   trn(ji,jj,jk,jppoc) ) & ! <CMOC OR 05/05/2014> Removal of GOC tracer ! + 1.05e4 * trn(ji,jj,jk,jpgoc)  )    &
                  &   * xdiss(ji,jj,jk) + 1E-4 * ( 1. - zlamfac ) * zdep
#else
               zlam1b = (  80.* (trn(ji,jj,jk,jpdoc) + 35E-6)              &
                  &     + 698.*  trn(ji,jj,jk,jppoc)  )                    &
                  &   * xdiss(ji,jj,jk) + 1E-4 * ( 1. - zlamfac ) * zdep
#endif
               zaggdfe = zlam1b * zstep * 0.5 * ( trn(ji,jj,jk,jpfer) - zfeequi )
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zscave - zaggdfe - zcoag
#if defined key_kriest
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zscave * zdenom1
#else
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zscave * zdenom1
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zscave * zdenom2
#endif
            END DO
         END DO
      END DO
      !

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem5')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      !     Calcite flux <CMOC OR 02/19/2014>
      !     --------------------------------------------------------------------
         xrcico(:,:) = rmcico_cmoc * exp ( aci_cmoc * ( tsn(:,:,1,jp_tem) - trcico_cmoc ) ) / ( 1._wp + exp ( aci_cmoc * ( tsn(:,:,1,jp_tem) - trcico_cmoc ) ) + rtrn )
         ! <CMOC OR 02/19/2014> rain ratio at level jk  

      DO jk = 1, jpkm1
      !     --------------------------------------------------------------------
      !     Update the arrays TRA which contain the biological sources and sinks
      !     --------------------------------------------------------------------
         tra(:,:,jk,jppo4) = tra(:,:,jk,jppo4) + zolimi (:,:,jk) + denitr(:,:,jk)
         tra(:,:,jk,jpnh4) = tra(:,:,jk,jpnh4) ! + zolimi (:,:,jk) + denitr(:,:,jk) <CMOC OR 11/08/2013> Testing if NH4 terms are the sources of the high concentrations of NO3/DIN in CMOC
         tra(:,:,jk,jpno3) = tra(:,:,jk,jpno3) + zredet (:,:,jk) *    trn(:,:,jk,jppoc) & !- denitr (:,:,jk) * rdenit <CMOC OR 10/02/2013> + re x D
         &                + zJNd(:,:,jk) ! <CMOC OR 12/11/2013> CMOC N2-fixation/denitrification rate
         tra(:,:,jk,jpdoc) = tra(:,:,jk,jpdoc) - zolimi (:,:,jk) - denitr(:,:,jk)
         tra(:,:,jk,jpoxy) = tra(:,:,jk,jpoxy) - zredet (:,:,jk) *    trn(:,:,jk,jppoc) !* 1.3_wp !- zolimi2(:,:,jk) * o2ut <CMOC OR 11/14/2013> remineralization used up oxygen in CMOC
         tra(:,:,jk,jpdic) = tra(:,:,jk,jpdic) + zredet (:,:,jk) *    trn(:,:,jk,jppoc)          !+ zolimi (:,:,jk) + denitr(:,:,jk) <CMOC OR 10/29/2013>
         tra(:,:,jk,jptal) = tra(:,:,jk,jptal) - zredet (:,:,jk) *    trn(:,:,jk,jppoc) * ncrr_cmoc !+ rno3 * ( zolimi(:,:,jk) + ( rdenit + 1.) * denitr(:,:,jk) ) <CMOC OR 10/29/2013>
      END DO

       DO ji = 1, jpi
        DO jj = 1, jpj

         ! <CMOC OR 03/12/2014> alk/dic balance only in the open ocean 
        IF ( 15 <= mbkt(ji,jj) ) THEN   ! <CMOC OR 03/13/2014> open ocean criterion 

         ! POC flux at the bottom of the euphotic zone <CMOC OR 02/25/2014> 
         zfpon(ji,jj) = xrcico(ji,jj) * wsbio3(ji,jj,11) * zstep * trn(ji,jj,11,jppoc) ! <CMOC OR 02/25/2014> zstep = x    step = rfact2 / rday
         !            rain ratio     * sinking rate   * dt/day in s * POC @ 100 m / 100 m   * evenly spread over the 11 layers of the euphotic zone 
         ! <CMOC OR 03/13/2014> Exponential dependency of the PIC burial 
         zbpon(ji,jj) = exp(-1* (fsdepw(ji,jj,mbkt(ji,jj)+1)-fsdepw(ji,jj,11)) /dci_cmoc)
         zbpoc(ji,jj) = trn(ji,jj,mbkt(ji,jj),jppoc)   ! <CMOC OR 03/14/2014> Bottom POC      

         DO jk =1, 11              ! in the euphotic zone

         ! <CMOC OR 02/19/2014> Calcite Precipitation, Sink of Alkalinity/DIC 
         tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) -     zfpon(ji,jj) * ( 1 - zbpon(ji,jj) ) * ideup_cmoc ! <CMOC OR 04/25/2014> the source associated with the zfpon flux should be zfpon/ideup_cmoc not zfpon/ideup_cmoc/11 ! <CMOC OR 03/12/2014> bathymetry dependence fixed ! <CMOC OR 02/24/2014> ! <CMOC OR 02/    19/2014>
         tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2 * zfpon(ji,jj) * ( 1 - zbpon(ji,jj) ) * ideup_cmoc ! <CMOC OR 04/25/2014> the source associated with the zfpon flux should be zfpon/ideup_cmoc not zfpon/ideup_cmoc/11 ! <CMOC OR 03/12/2014> bathymetry dependence fixed ! <CMOC OR 02/24/2014> ! <CMOC OR 02/    19/2014>
         END DO

         ! <CMOC OR 02/25/2014> the extra term - exp(), which is in fact + zfpon * exp(), is a balance term to enforce long-term equilibrium between surface and bottom conditions of DIC and Akalinity  

         DO jk = 12, mbkt(ji,jj)          ! <CMOC OR 03/12/2014> ! below the euphotic zone
         ! <CMOC OR 02/19/2014> Calcite Dissolution, Source of Alkalinity/DIC 
         ! <CMOC OR 03/13/2014> Burial of PIC 
         tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) +     zfpon(ji,jj) * ( exp(-1._wp*(fsdepw(ji,jj,jk)-fsdepw(ji,jj,11))/dci_cmoc) - exp(-1._wp*(fsdepw(ji,jj,jk+1)-fsdepw(ji,jj,11))/dci_cmoc) ) / fse3w(ji,jj    ,jk)
         tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 2 * zfpon(ji,jj) * ( exp(-1._wp*(fsdepw(ji,jj,jk)-fsdepw(ji,jj,11))/dci_cmoc) - exp(-1._wp*(fsdepw(ji,jj,jk+1)-fsdepw(ji,jj,11))/dci_cmoc) ) / fse3w(ji,jj    ,jk)
         ! <CMOC OR 02/19/2014>                          Particulate Inorganic Carbon Export at 100 m         / 100 m   *      dz       * below-euphotic-layer thickness 
         ! <CMOC OR 02/25/2014> Below the euphotic zone calcite dissolution dominates; the POC flux varies as zfpon*exp(-(z-105)/2700), zfpon being the flux at the bottom of the euphotic zone 
         END DO

        ENDIF

        END DO
       END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem6')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      ! <CMOC OR 04/28/2014> fix N2-fixation diagnostics
      xn2fixdia(:,:) = xn2fixdia(:,:) + zn2fixtot(:,:) ! add the intermediate time-splits to get the N2-fixation rate over a time-step
      
      IF( ln_diatrc ) THEN                      ! <CMOC OR 01/15/2014> CMOC diagnosis of N2 Fixation
        IF( lk_iomput ) THEN 
         IF( jnt == nrdttrc ) THEN              ! <CMOC OR 01/20/2014>
              zwork(:,:)  =  xn2fixdia(:,:) * ncrr_cmoc * 1.e+3_wp * rfact2r / 4 * tmask(:,:,1) ! <CMOC OR 05/13/2014> fix the time step issue with adding all the n2 fixation in the diagnostic variable ! <CMOC OR 04/28/2014> fix N2-fixation diagnostics ! convert into mol N m^-2 s^-1
              CALL iom_put( "Nfix"   , zwork )  ! nitrogen fixation in molN m^-2 s^-1 
              CALL iom_put( "BUPOC"  , wsbio3(:,:,11) /rday * zbpoc(:,:) * 1e+3_wp  )  ! <CMOC OR 03/12/2014> POC burial flux
              CALL iom_put( "BUCALC" , zfpon(:,:) * 1e+3_wp / rfact2 * zbpon(:,:)  )  ! <CMOC OR 03/12/2014> PIC burial flux

         ENDIF                                  ! <CMOC OR 01/20/2014>
        ENDIF
      ENDIF                                     ! <CMOC OR 01/15/2014> End Diagnosis Block


      !
      CALL wrk_dealloc( jpi, jpj,      ztempbac, zredettot, zn2fixtot, zwork, zfpon, zbpon, zbpoc  ) ! <CMOC OR 03/14/2014> Bottom POC ! <CMOC OR 03/13/2014> add Burial of PIC ! <CMOC OR 02/25/2014> ! <CMOC OR 01/15/2014> !<CMOC OR 12/11/2013> !<CMOC OR 10/02/2013> 
      CALL wrk_dealloc( jpi, jpj, jpk, zdepbac, zolimi, zolimi2, zredet,    zn2fix,   zJNd         ) ! <CMOC OR 12/11/2013> !<CMOC OR 10/02/2013> remineralization rate of detritus
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
      NAMELIST/nampisrem/ xremik, xremip, nitrif, xsirem, xsiremlab, xsilab,   &
      &                   xlam1, oxymin, ligand 

      ! <CMOC OR 03/08/2014> CMOC namelist
      NAMELIST/namcmocpoc/ ed_cmoc, reref_cmoc
      NAMELIST/namcmoccal/ rmcico_cmoc, trcico_cmoc, aci_cmoc, dci_cmoc
      NAMELIST/namcmocnfx/ phinf_cmoc, phi0_cmoc, anf_cmoc, pnf_cmoc, inf_cmoc, tnfMa_cmoc, tnfmi_cmoc
      NAMELIST/namcmocdeu/  deup_cmoc, ideup_cmoc
      ! <CMOC OR 03/08/2014> CMOC namelist end 
      !!----------------------------------------------------------------------

      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampisrem )

      REWIND( numcmoc )                     ! <CMOC OR 03/10/2014>  ! read numcmoc, cmocpoc
      READ  ( numcmoc, namcmocpoc )

      REWIND( numcmoc )                     ! <CMOC OR 03/10/2014>  ! read numcmoc, cmoccal
      READ  ( numcmoc, namcmoccal )

      REWIND( numcmoc )                     ! <CMOC OR 03/10/2014>  ! read numcmoc, cmocnfx
      READ  ( numcmoc, namcmocnfx )

      REWIND( numcmoc )                     ! <CMOC OR 03/10/2014>  ! read numcmoc, cmocdeu
      READ  ( numcmoc, namcmocdeu )
      
      
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for remineralization, nampisrem'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    remineralisation rate of POC              xremip    =', xremip
         WRITE(numout,*) '    remineralization rate of DOC              xremik    =', xremik
         WRITE(numout,*) '    remineralization rate of Si               xsirem    =', xsirem
         WRITE(numout,*) '    fast remineralization rate of Si          xsiremlab =', xsiremlab
         WRITE(numout,*) '    fraction of labile biogenic silica        xsilab    =', xsilab
         WRITE(numout,*) '    scavenging rate of Iron                   xlam1     =', xlam1
         WRITE(numout,*) '    NH4 nitrification rate                    nitrif    =', nitrif
         WRITE(numout,*) '    halk saturation constant for anoxia       oxymin    =', oxymin
         WRITE(numout,*) '    ligand concentration in the ocean         ligand    =', ligand
         WRITE(numout,*) ' '
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
      nitrfac (:,:,:) = 0._wp
      denitr  (:,:,:) = 0._wp
      denitnh4(:,:,:) = 0._wp
      !
   END SUBROUTINE p4z_rem_init


   INTEGER FUNCTION p4z_rem_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( denitr(jpi,jpj,jpk), denitnh4(jpi,jpj,jpk), STAT=p4z_rem_alloc )
      !
      IF( p4z_rem_alloc /= 0 )   CALL ctl_warn('p4z_rem_alloc: failed to allocate arrays')
      !
   END FUNCTION p4z_rem_alloc

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
