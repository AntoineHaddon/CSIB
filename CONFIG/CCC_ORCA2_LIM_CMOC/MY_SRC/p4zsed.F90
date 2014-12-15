MODULE p4zsed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :   PISCES Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) USE of fldread
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_sed        :  Compute loss of organic matter in the sediments
   !!   p4z_sbc        :  Read and interpolate time-varying nutrients fluxes
   !!   p4z_sed_init   :  Initialization of p4z_sed
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  USE p4zopt          !  optical model
 ! <CMOC OR 07/15/2014> Revert for consistency !    USE p4zlim          !  Co-limitations of differents nutrients
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  USE p4zrem          !  Remineralisation of organic matter
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  USE p4zint          !  interpolation and computation of various fields
   USE iom             !  I/O manager
   USE fldread         !  time interpolation
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sed   
   PUBLIC   p4z_sed_init   
   PUBLIC   p4z_sed_alloc

   !! * Shared module variables
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  LOGICAL  :: ln_dust     = .FALSE.    !: boolean for dust input from the atmosphere
   LOGICAL  :: ln_river    = .FALSE.    !: boolean for river input of nutrients
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  LOGICAL  :: ln_ndepo    = .FALSE.    !: boolean for atmospheric deposition of N
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  LOGICAL  :: ln_ironsed  = .FALSE.    !: boolean for Fe input from sediments

   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp) :: sedfeinput  = 1.E-9_wp   !: Coastal release of Iron
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp) :: dustsolub   = 0.014_wp   !: Solubility of the dust
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp) :: wdust       = 2.0_wp     !: Sinking speed of the dust 
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp) :: nitrfix     = 1E-7_wp    !: Nitrogen fixation rate   
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp) :: diazolight  = 50._wp     !: Nitrogen fixation sensitivty to light 
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp) :: concfediaz  = 1.E-10_wp  !: Fe half-saturation Cste for diazotrophs 


   !! * Module variables
   REAL(wp) :: ryyss                  !: number of seconds per year 
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp) :: r1_ryyss                 !: inverse of ryyss
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp) :: rmtss                  !: number of seconds per month
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp) :: r1_rday                  !: inverse of rday
   LOGICAL  :: ll_sbc

   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_dust      ! structure of input dust
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_riverdic  ! structure of input riverdic
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_riverdoc  ! structure of input riverdoc
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_ndepo     ! structure of input nitrogen deposition
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_ironsed   ! structure of input iron from sediment
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_fmsk      ! structure of input iron limitation mask ! <CMOC OR 01/21/2014>

   INTEGER , PARAMETER :: nbtimes = 365  !: maximum number of times record in a file
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  INTEGER  :: ntimes_dust, ntimes_riv, ntimes_ndep ! <CMOC OR 01/22/2014>, ntimes_fmsk  ! <CMOC OR 01/22/2014> add ntimes_fmsk     ! number of time steps in a file
   INTEGER  :: ntimes_riv ! <CMOC OR 01/22/2014>, ntimes_fmsk  ! <CMOC OR 01/22/2014> add ntimes_fmsk     ! number of time steps in a file

   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:) :: dust      !: dust fields
   REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:) :: rivinp, cotdep    !: river input fields
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp), ALLOCATABLE, SAVE,   DIMENSION(:,:) :: nitdep    !: atmospheric N deposition 
   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: ironsed   !: Coastal supply of iron

   ! <CMOC OR 06/02/2014 Streamlining p4zsed >  REAL(wp) :: sumdepsi, rivalkinput, rivpo4input, nitdepinput

   !!* Substitution
#  include "top_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Header:$ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS


   SUBROUTINE p4z_sed( kt, jnt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed  ***
      !!
      !! ** Purpose :   Compute loss of organic matter in the sediments. This
      !!              is by no way a sediment model. The loss is simply 
      !!              computed to balance the inout from rivers and dust
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, jnt ! ocean time step
      INTEGER  ::   ji, jj, jk, ikt
! <CMOC OR 06/02/2014 Streamlining p4zsed >  #if ! defined key_sed
! <CMOC OR 06/02/2014 Streamlining p4zsed >        REAL(wp) ::   zsumsedsi, zsumsedpo4, zsumsedcal
! <CMOC OR 06/02/2014 Streamlining p4zsed >        REAL(wp) ::   zrivalk, zrivsil, zrivpo4
! <CMOC OR 06/02/2014 Streamlining p4zsed >  #endif
! <CMOC OR 06/02/2014 Streamlining p4zsed >        REAL(wp) ::   zdenitot, znitrpottot, zlim, zfact, zfactcal
! <CMOC OR 06/02/2014 Streamlining p4zsed >        REAL(wp) ::   zsiloss, zcaloss, zwsbio3, zwsbio4, zwscal, zdep
      REAL(wp) ::   zwsbio3, zwsbio4, zwscal, zdep
      CHARACTER (len=25) :: charout
! <CMOC OR 06/02/2014 Streamlining p4zsed >        REAL(wp), POINTER, DIMENSION(:,:  ) :: zsidep, zwork1, zwork2, zwork3
! <CMOC OR 06/02/2014 Streamlining p4zsed >        REAL(wp), POINTER, DIMENSION(:,:,:) :: znitrpot, zirondep
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sed')
      !
      ! Allocate temporary workspace
! <CMOC OR 06/02/2014 Streamlining p4zsed >        CALL wrk_alloc( jpi, jpj,      zsidep, zwork1, zwork2, zwork3 )
! <CMOC OR 06/02/2014 Streamlining p4zsed >        CALL wrk_alloc( jpi, jpj, jpk, znitrpot, zirondep             )

      ! Iron limitation mask <CMOC OR 01/21/2014>
      ! -----------------------------------------
      !IF( kt == nit000 .AND. jnt == 1 ) THEN
      !   CALL fld_read( kt, 1, sf_fmsk )
      !   xlimnfecmoc(:,:) = sf_fmsk(1)%fnow(:,:) ! No mask
      ! <CMOC OR 01/21/2014>
      !END IF    

      IF( jnt == 1 .AND. ll_sbc ) CALL p4z_sbc( kt )

      ! <CMOC OR 06/02/2014 Streamlining p4zsed > zirondep(:,:,:) = 0.e0          ! Initialisation of variables USEd to compute deposition
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > zsidep  (:,:)   = 0.e0

      ! Iron and Si deposition at the surface
      ! -------------------------------------
! <CMOC OR 06/02/2014 Streamlining p4zsed >       DO jj = 1, jpj
!         DO ji = 1, jpi
!            zdep  = rfact2 / fse3t(ji,jj,1)
!            zirondep(ji,jj,1) = ( dustsolub * dust(ji,jj) / ( 55.85 * rmtss ) + 3.e-10 * r1_ryyss ) * zdep
!            zsidep  (ji,jj)   = 8.8 * 0.075 * dust(ji,jj) * zdep / ( 28.1 * rmtss )
!         END DO
!      END DO
!
!      ! Iron solubilization of particles in the water column
!      ! ----------------------------------------------------
!      DO jk = 2, jpkm1
!         zirondep(:,:,jk) = dust(:,:) / ( wdust * 55.85 * rmtss ) * rfact2 * 1.e-4 * EXP( -fsdept(:,:,jk) / 1000. )
! <CMOC OR 06/02/2014 Streamlining p4zsed >       END DO

      ! Add the external input of nutrients, carbon and alkalinity
      ! ----------------------------------------------------------
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > trn(:,:,1,jppo4) = trn(:,:,1,jppo4) + rivinp(:,:) * rfact2 
      !trn(:,:,1,jpno3) = trn(:,:,1,jpno3) + (rivinp(:,:) + nitdep(:,:)) * rfact2
      trn(:,:,1,jpno3) = trn(:,:,1,jpno3) + rivinp(:,:) * rfact2
      !trn(:,:,1,jpfer) = trn(:,:,1,jpfer) + rivinp(:,:) * 3.e-5 * rfact2
      !trn(:,:,1,jpsil) = trn(:,:,1,jpsil) + zsidep (:,:) + cotdep(:,:)   * rfact2 / 6.
      trn(:,:,1,jpdic) = trn(:,:,1,jpdic) + rivinp(:,:) * 2.631 * rfact2
      !trn(:,:,1,jptal) = trn(:,:,1,jptal) + (cotdep(:,:) - rno3*(rivinp(:,:) +  nitdep(:,:) ) ) * rfact2
      trn(:,:,1,jptal) = trn(:,:,1,jptal) + (cotdep(:,:) - rno3*rivinp(:,:) ) * rfact2
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > 

      ! Add the external input of iron which is 3D distributed
      ! (dust, river and sediment mobilization)
      ! ------------------------------------------------------
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > DO jk = 1, jpkm1
      ! <CMOC OR 06/02/2014 Streamlining p4zsed >    trn(:,:,jk,jpfer) = trn(:,:,jk,jpfer) + zirondep(:,:,jk) + ironsed(:,:,jk) * rfact2
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > END DO

! <CMOC OR 06/02/2014 Streamlining p4zsed > #if ! defined key_sed
!      ! Loss of biogenic silicon, Caco3 organic carbon in the sediments. 
!      ! First, the total loss is computed.
!      ! The factor for calcite comes from the alkalinity effect
!      ! -------------------------------------------------------------
!      DO jj = 1, jpj
!         DO ji = 1, jpi
!            ikt = mbkt(ji,jj) 
!# if defined key_kriest
!            zwork1(ji,jj) = trn(ji,jj,ikt,jpgsi) * wscal (ji,jj,ikt)
!            zwork2(ji,jj) = trn(ji,jj,ikt,jppoc) * wsbio3(ji,jj,ikt)
!# else
!            zwork1(ji,jj) = trn(ji,jj,ikt,jpgsi) * wsbio4(ji,jj,ikt)
!            ! <CMOC OR 05/05/2014> Removal of GOC tracer ! zwork2(ji,jj) = trn(ji,jj,ikt,jpgoc) * wsbio4(ji,jj,ikt) + trn(ji,jj,ikt,jppoc) * wsbio3(ji,jj,ikt) 
!            zwork2(ji,jj) = trn(ji,jj,ikt,jppoc) * wsbio3(ji,jj,ikt) ! <CMOC OR 05/05/2014> Removal of GOC tracer ! 
!# endif
!            ! For calcite, burial efficiency is made a function of saturation
!            zfactcal      = MIN( excess(ji,jj,ikt), 0.2 )
!            zfactcal      = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )
!            zwork3(ji,jj) = trn(ji,jj,ikt,jpcal) * wscal (ji,jj,ikt) * 2.e0 * zfactcal
!         END DO
!      END DO
!      zsumsedsi  = glob_sum( zwork1(:,:) * e1e2t(:,:) ) * r1_rday
!      zsumsedpo4 = glob_sum( zwork2(:,:) * e1e2t(:,:) ) * r1_rday
!      zsumsedcal = glob_sum( zwork3(:,:) * e1e2t(:,:) ) * r1_rday
! <CMOC OR 06/02/2014 Streamlining p4zsed > #endif

      ! THEN this loss is scaled at each bottom grid cell for
      ! equilibrating the total budget of silica in the ocean.
      ! Thus, the amount of silica lost in the sediments equal
      ! the supply at the surface (dust+rivers)
      ! ------------------------------------------------------
! <CMOC OR 06/02/2014 Streamlining p4zsed > #if ! defined key_sed
!      zrivsil =  1._wp - ( sumdepsi + rivalkinput * r1_ryyss / 6. ) / zsumsedsi 
!      zrivpo4 =  1._wp - ( rivpo4input * r1_ryyss ) / zsumsedpo4 
!#endif
!
!      DO jj = 1, jpj
!         DO ji = 1, jpi
!            ikt  = mbkt(ji,jj)
!            zdep = xstep / fse3t(ji,jj,ikt)
!            zwsbio4 = wsbio4(ji,jj,ikt) * zdep
!            zwscal  = wscal (ji,jj,ikt) * zdep
!! <CMOC OR 06/02/2014 Streamlining p4zsed > # if defined key_kriest
!! <CMOC OR 06/02/2014 Streamlining p4zsed >             zsiloss = trn(ji,jj,ikt,jpgsi) * zwsbio4
!! <CMOC OR 06/02/2014 Streamlining p4zsed > # else
!            zsiloss = trn(ji,jj,ikt,jpgsi) * zwscal
!! <CMOC OR 06/02/2014 Streamlining p4zsed > # endif
!            zcaloss = trn(ji,jj,ikt,jpcal) * zwscal
!            !
!            trn(ji,jj,ikt,jpgsi) = trn(ji,jj,ikt,jpgsi) - zsiloss
!            trn(ji,jj,ikt,jpcal) = trn(ji,jj,ikt,jpcal) - zcaloss
!! <CMOC OR 06/02/2014 Streamlining p4zsed > #if ! defined key_sed
!            trn(ji,jj,ikt,jpsil) = trn(ji,jj,ikt,jpsil) + zsiloss * zrivsil 
!            zfactcal = MIN( excess(ji,jj,ikt), 0.2 )
!            zfactcal = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )
!            zrivalk  =  1._wp - ( rivalkinput * r1_ryyss ) * zfactcal / zsumsedcal 
!            trn(ji,jj,ikt,jptal) =  trn(ji,jj,ikt,jptal) ! + zcaloss * zrivalk * 2.0 !<CMOC OR 11/29/2013> 
!            trn(ji,jj,ikt,jpdic) =  trn(ji,jj,ikt,jpdic) !+ zcaloss * zrivalk       !<CMOC OR 11/29/2013> 
!! <CMOC OR 06/02/2014 Streamlining p4zsed > #endif
!         END DO
! <CMOC OR 06/02/2014 Streamlining p4zsed >       END DO

      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstep / fse3t(ji,jj,ikt)
            zwsbio4 = wsbio4(ji,jj,ikt) * zdep
            zwsbio3 = wsbio3(ji,jj,ikt) * zdep
! <CMOC OR 06/02/2014 Streamlining p4zsed > # if ! defined key_kriest
            ! <CMOC OR 05/05/2014> Removal of GOC tracer ! trn(ji,jj,ikt,jpgoc) = trn(ji,jj,ikt,jpgoc) - trn(ji,jj,ikt,jpgoc) * zwsbio4
            trn(ji,jj,ikt,jppoc) = trn(ji,jj,ikt,jppoc) - trn(ji,jj,ikt,jppoc) * zwsbio3
            ! <CMOC OR 06/02/2014 Streamlining p4zsed > trn(ji,jj,ikt,jpbfe) = trn(ji,jj,ikt,jpbfe) - trn(ji,jj,ikt,jpbfe) * zwsbio4
            ! <CMOC OR 06/02/2014 Streamlining p4zsed > trn(ji,jj,ikt,jpsfe) = trn(ji,jj,ikt,jpsfe) - trn(ji,jj,ikt,jpsfe) * zwsbio3
! <CMOC OR 06/02/2014 Streamlining p4zsed > #if ! defined key_sed
            ! <CMOC OR 06/02/2014 Streamlining p4zsed > trn(ji,jj,ikt,jpdoc) = trn(ji,jj,ikt,jpdoc) &
               ! <CMOC OR 05/05/2014> Removal of GOC tracer ! &               + ( trn(ji,jj,ikt,jpgoc) * zwsbio4 + trn(ji,jj,ikt,jppoc) * zwsbio3 ) * zrivpo4
               ! <CMOC OR 06/02/2014 Streamlining p4zsed > &               + ( trn(ji,jj,ikt,jppoc) * zwsbio3 ) * zrivpo4 ! <CMOC OR 05/05/2014> Removal of GOC tracer ! 
            trn(ji,jj,ikt,jpdic) = trn(ji,jj,ikt,jpdic) &
               &               + trn(ji,jj,ikt,jppoc) * zwsbio3 ! <CMOC OR 01/15/2014> instantaneously remineralize bottom sinking POC into DIC, DIN and take up O2  
            trn(ji,jj,ikt,jpno3) = trn(ji,jj,ikt,jpno3) &
               &               + trn(ji,jj,ikt,jppoc) * zwsbio3 ! <CMOC OR 01/15/2014> instantaneously remineralize bottom sinking POC into DIC, DIN and take up O2
            trn(ji,jj,ikt,jpoxy) = trn(ji,jj,ikt,jpoxy) &
               &               - trn(ji,jj,ikt,jppoc) * zwsbio3 ! <CMOC OR 01/15/2014> instantaneously remineralize bottom sinking POC into DIC, DIN and take up O2 
! <CMOC OR 06/02/2014 Streamlining p4zsed > #endif

! <CMOC OR 06/02/2014 Streamlining p4zsed > # else
!            trn(ji,jj,ikt,jpnum) = trn(ji,jj,ikt,jpnum) - trn(ji,jj,ikt,jpnum) * zwsbio4
!            trn(ji,jj,ikt,jppoc) = trn(ji,jj,ikt,jppoc) - trn(ji,jj,ikt,jppoc) * zwsbio3
!            trn(ji,jj,ikt,jpsfe) = trn(ji,jj,ikt,jpsfe) - trn(ji,jj,ikt,jpsfe) * zwsbio3
! <CMOC OR 06/02/2014 Streamlining p4zsed > #if ! defined key_sed
!            trn(ji,jj,ikt,jpdoc) = trn(ji,jj,ikt,jpdoc) &
!               &               + ( trn(ji,jj,ikt,jpnum) * zwsbio4 + trn(ji,jj,ikt,jppoc) * zwsbio3 ) * zrivpo4
!#endif
!
! <CMOC OR 06/02/2014 Streamlining p4zsed > # endif
         END DO
      END DO


      ! Nitrogen fixation (simple parameterization). The total gain
      ! from nitrogen fixation is scaled to balance the loss by 
      ! denitrification
      ! -------------------------------------------------------------

      ! <CMOC OR 06/02/2014 Streamlining p4zsed > zdenitot = glob_sum(  ( denitr(:,:,:) * rdenit + denitnh4(:,:,:) * rdenita ) * cvol(:,:,:) ) 
!
!      ! Potential nitrogen fixation dependant on temperature and iron
!      ! -------------------------------------------------------------
!
!!CDIR NOVERRCHK
!      DO jk = 1, jpk
!!CDIR NOVERRCHK
!         DO jj = 1, jpj
!!CDIR NOVERRCHK
!            DO ji = 1, jpi
!               zlim = ( 1.- xnanono3(ji,jj,jk) - xnanonh4(ji,jj,jk) )
!               IF( zlim <= 0.2 )   zlim = 0.01
!#if defined key_degrad
!               zfact = zlim * rfact2 * facvol(ji,jj,jk)
!#else
!               zfact = zlim * rfact2 
!#endif
!               znitrpot(ji,jj,jk) =  MAX( 0.e0, ( 0.6 * tgfunc(ji,jj,jk) - 2.15 ) * r1_rday )   &
!                 &                 * 1._wp                                                      & ! <CMOC OR 01/21/2014> ! <CMOC OR 11/27/2013> no iron limitation for the time being &                 *  zfact * trn(ji,jj,jk,jpfer) / ( concfediaz + trn(ji,jj,jk,jpfer) ) &
!                 &                 * 0.1E-6_wp * 6.6_wp / ( 0.1E-6_wp * 6.6_wp + trn(ji,jj,jk,jpno3) + rtrn )  & ! <CMOC OR 11/27/2013> introduced a nitrogen limitation coefficient similar to CMOC coefficient 
!                 &                 * ( 1.- EXP( -etot(ji,jj,jk) / diazolight ) ) * rfact2 ! <CMOC OR 11/27/2013> original line is multiplied by <rfact2> to disappearance of zfact with the line removed two rows above 
!            END DO
!         END DO 
!      END DO
!
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > znitrpottot = glob_sum( znitrpot(:,:,:) * cvol(:,:,:) )

! <CMOC OR 06/02/2014 Streamlining p4zsed >       ! Nitrogen change due to nitrogen fixation
!      ! ----------------------------------------
!      DO jk = 1, jpk
!         DO jj = 1, jpj
!            DO ji = 1, jpi
!               zfact = znitrpot(ji,jj,jk) * nitrfix
!               trn(ji,jj,jk,jpnh4) = trn(ji,jj,jk,jpnh4) !+ zfact <CMOC OR 11/08/2013> remove N2 fixation to test if it is the external source of NO3 
!               trn(ji,jj,jk,jptal) = trn(ji,jj,jk,jptal) !+ rno3 * zfact <CMOC OR 10/04/2013> N2 fixation disabled for now
!               trn(ji,jj,jk,jpoxy) = trn(ji,jj,jk,jpoxy) !+ zfact   * o2nit <CMOC OR 10/04/2013> N2 fixation disabled for now
!               trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4) + 30. / 46. * zfact
!           !    trn(ji,jj,jk,jppo4) = trn(ji,jj,jk,jppo4) + zfact
!           END DO
!         END DO 
!      END DO
!      !
!      IF( ln_diatrc ) THEN
!         zfact = 1.e+3 * rfact2r
!         IF( lk_iomput ) THEN
!            zwork1(:,:)  =  ( zirondep(:,:,1) + ironsed(:,:,1) * rfact2 ) * zfact * fse3t(:,:,1) * tmask(:,:,1) 
!          ! <CMOC OR 01/15/2014>  zwork2(:,:)  =    znitrpot(:,:,1) * nitrfix                   * zfact * fse3t(:,:,1) * tmask(:,:,1)
!            IF( jnt == nrdttrc ) THEN
!               CALL iom_put( "Irondep", zwork1  )  ! surface downward net flux of iron
!          ! <CMOC OR 01/15/2014>     CALL iom_put( "Nfix"   , zwork2 )  ! nitrogen fixation at surface
!            ENDIF
!         ELSE
!            trc2d(:,:,jp_pcs0_2d + 11) = zirondep(:,:,1)           * zfact * fse3t(:,:,1) * tmask(:,:,1)
!          ! <CMO OR 01/15/2014>  trc2d(:,:,jp_pcs0_2d + 12) = znitrpot(:,:,1) * nitrfix * zfact * fse3t(:,:,1) * tmask(:,:,1)
!         ENDIF
!      ENDIF
!      !
!      IF(ln_ctl) THEN  ! print mean trends (USEd for debugging)
!         WRITE(charout, fmt="('sed ')")
!         CALL prt_ctl_trc_info(charout)
!         CALL prt_ctl_trc(tab4d=trn, mask=tmask, clinfo=ctrcnm)
!      ENDIF
!      !
!      CALL wrk_dealloc( jpi, jpj,      zsidep, zwork1, zwork2, zwork3 )
! <CMOC OR 06/02/2014 Streamlining p4zsed >       CALL wrk_dealloc( jpi, jpj, jpk, znitrpot, zirondep             )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sed')
      !
   END SUBROUTINE p4z_sed

!! <CMOC OR 06/02/2014 Streamlining p4zsed >
   SUBROUTINE p4z_sbc( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_sbc  ***
      !!
      !! ** purpose :   read and interpolate the external sources of 
      !!                nutrients
      !!
      !! ** method  :   read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * local declarations
      INTEGER  :: ji,jj 
      REAL(wp) :: zcoef
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sbc')
      !
      ! Compute dust at nit000 or only if there is more than 1 time record in dust file
!      IF( ln_dust ) THEN
!         IF( kt == nit000 .OR. ( kt /= nit000 .AND. ntimes_dust > 1 ) ) THEN
!            CALL fld_read( kt, 1, sf_dust )
!            dust(:,:) = sf_dust(1)%fnow(:,:,1)
!         ENDIF
!      ENDIF
!
      ! N/P and Si releases due to coastal rivers
      ! Compute river at nit000 or only if there is more than 1 time record in river file
      ! -----------------------------------------
      IF( ln_river ) THEN
         IF( kt == nit000 .OR. ( kt /= nit000 .AND. ntimes_riv > 1 ) ) THEN
            CALL fld_read( kt, 1, sf_riverdic )
            CALL fld_read( kt, 1, sf_riverdoc )
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zcoef = ryyss * cvol(ji,jj,1) 
                  cotdep(ji,jj) =   sf_riverdic(1)%fnow(ji,jj,1)                                  * 1E9 / ( 12. * zcoef + rtrn )
                  rivinp(ji,jj) = ( sf_riverdic(1)%fnow(ji,jj,1) + sf_riverdoc(1)%fnow(ji,jj,1) ) * 1E9 / ( 31.6* zcoef + rtrn )
               END DO
            END DO
         ENDIF
      ENDIF

!      ! Compute N deposition at nit000 or only if there is more than 1 time record in N deposition file
!      IF( ln_ndepo ) THEN
!         IF( kt == nit000 .OR. ( kt /= nit000 .AND. ntimes_ndep > 1 ) ) THEN
!            CALL fld_read( kt, 1, sf_ndepo )
!            DO jj = 1, jpj
!               DO ji = 1, jpi
!                  nitdep(ji,jj) = 7.6 * sf_ndepo(1)%fnow(ji,jj,1) / ( 14E6 * ryyss * fse3t(ji,jj,1) + rtrn )
!               END DO
!            END DO
!         ENDIF
!      ENDIF
!      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sbc')
!      !
   END SUBROUTINE p4z_sbc
! <CMOC OR 06/02/2014 Streamlining p4zsed > 

   SUBROUTINE p4z_sed_init

      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_sed_init  ***
      !!
      !! ** purpose :   initialization of the external sources of nutrients
      !!
      !! ** method  :   read the files and compute the budget
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER  :: ji, jj, jk, jm
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > INTEGER  :: numdust, numriv, numiron, numdepo, numfmsk ! <CMOC OR 01/21/2014> add numfmsk (iron limitation mask)
      INTEGER  :: numfmsk, numriv ! <CMOC OR 01/21/2014> add numfmsk (iron limitation mask)
      INTEGER  :: ierr, ierr1, ierr2 ! <CMOC OR 06/02/2014 Streamlining p4zsed >  , ierr1, ierr2, ierr3
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > REAL(wp) :: zexpide, zdenitide, zmaskt
      REAL(wp), DIMENSION(nbtimes) :: zsteps                 ! times records
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zdust, zndepo, zriverdic, zriverdoc, zcmask ! <CMOC OR 01/22/2014>, zfmask ! <OR CMOC 01/22/2014> add zfmask for iron limitation mask
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zriverdic, zriverdoc ! <CMOC OR 01/22/2014>, zfmask ! <OR CMOC 01/22/2014> add zfmask for iron limitation mask
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE :: zfmask                                      ! <OR CMOC 01/22/2014> add zfmask for iron limitation mask
      !
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > TYPE(FLD_N) ::   sn_dust, sn_riverdoc, sn_riverdic, sn_ndepo, sn_ironsed, sn_fmsk    ! <OR CMOC 01/21/2014> add sn_fmsk (iron limitation mask)        ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_fmsk, sn_riverdoc, sn_riverdic    ! <OR CMOC 01/21/2014> add sn_fmsk (iron limitation mask)        ! informations about the fields to be read
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > NAMELIST/nampissed/cn_dir, sn_dust, sn_riverdic, sn_riverdoc, sn_ndepo, sn_ironsed, &
      NAMELIST/nampissed/cn_dir, sn_riverdic, sn_riverdoc,                         &
        &                ln_river,                                                 & ! <CMOC OR 06/02/2014 Streamlining p4zsed > &                sn_fmsk, ln_dust, ln_river, ln_ndepo, ln_ironsed,         &       ! <OR CMOC 02/05/2014> sn_fmsk was missing from the namelist block
        &                sn_fmsk                ! <OR CMOC 02/05/2014> sn_fmsk was missing from the namelist block
        ! <CMOC OR 06/02/2014 Streamlining p4zsed > &                sedfeinput, dustsolub, wdust, nitrfix, diazolight, concfediaz 
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_sed_init')
      !
      !                                    ! number of seconds per year and per month
        ryyss    = nyear_len(1) * rday
! <CMOC OR 06/02/2014 Streamlining p4zsed >        rmtss    = ryyss / raamo
! <CMOC OR 06/02/2014 Streamlining p4zsed >        r1_rday  = 1. / rday
! <CMOC OR 06/02/2014 Streamlining p4zsed >        r1_ryyss = 1. / ryyss
      !                            !* set file information
      cn_dir  = './'            ! directory in which the model is executed
      ! ... default values (NB: frequency positive => hours, negative => months)
      !                  !   file       ! frequency !  variable   ! time intep !  clim   ! 'yearly' or ! weights  ! rotation   !
      !                  !   name       !  (hours)  !   name      !   (T/F)    !  (T/F)  !  'monthly'  ! filename ! pairs      !
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > sn_dust     = FLD_N( 'dust'       ,    -1     ,  'dust'     ,  .true.    , .true.  ,   'yearly'  , ''       , ''         )
      sn_riverdic = FLD_N( 'river'      ,   -12     ,  'riverdic' ,  .false.   , .true.  ,   'yearly'  , ''       , ''         )
      sn_riverdoc = FLD_N( 'river'      ,   -12     ,  'riverdoc' ,  .false.   , .true.  ,   'yearly'  , ''       , ''         )
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > sn_ndepo    = FLD_N( 'ndeposition',   -12     ,  'ndep'     ,  .false.   , .true.  ,   'yearly'  , ''       , ''         )
      ! <CMOC OR 06/02/2014 Streamlining p4zsed > sn_ironsed  = FLD_N( 'ironsed'    ,   -12     ,  'bathy'    ,  .false.   , .true.  ,   'yearly'  , ''       , ''         )
      sn_fmsk     = FLD_N( 'fermask'    ,   -12     ,  'femask'   ,  .false.   , .true.  ,   'yearly'  , ''       , ''         )   ! <OR CMOC 01/21/2014> add the iron limitation mask

      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampissed )

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' namelist : nampissed '
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~ '
         ! <CMOC OR 06/02/2014 Streamlining p4zsed > WRITE(numout,*) '    dust input from the atmosphere           ln_dust     = ', ln_dust
         WRITE(numout,*) '    river input of nutrients                 ln_river    = ', ln_river
         ! <CMOC OR 06/02/2014 Streamlining p4zsed > WRITE(numout,*) '    atmospheric deposition of n              ln_ndepo    = ', ln_ndepo
         ! <CMOC OR 06/02/2014 Streamlining p4zsed > WRITE(numout,*) '    fe input from sediments                  ln_sedinput = ', ln_ironsed
         ! <CMOC OR 06/02/2014 Streamlining p4zsed > WRITE(numout,*) '    coastal release of iron                  sedfeinput  = ', sedfeinput
         ! <CMOC OR 06/02/2014 Streamlining p4zsed > WRITE(numout,*) '    solubility of the dust                   dustsolub   = ', dustsolub
         ! <CMOC OR 06/02/2014 Streamlining p4zsed > WRITE(numout,*) '    sinking speed of the dust                wdust       = ', wdust
         ! <CMOC OR 06/02/2014 Streamlining p4zsed > WRITE(numout,*) '    nitrogen fixation rate                   nitrfix     = ', nitrfix
         ! <CMOC OR 06/02/2014 Streamlining p4zsed > WRITE(numout,*) '    nitrogen fixation sensitivty to light    diazolight  = ', diazolight
         ! <CMOC OR 06/02/2014 Streamlining p4zsed > WRITE(numout,*) '    fe half-saturation cste for diazotrophs  concfediaz  = ', concfediaz
       END IF

      ! <CMOC OR 06/02/2014 Streamlining p4zsed > IF( ln_dust .OR. ln_river .OR. ln_ndepo ) THEN
      IF( ln_river ) THEN
           ll_sbc = .TRUE.
      ELSE
           ll_sbc = .FALSE.
      ENDIF

      ! dust input from the atmosphere
      ! ------------------------------
! <CMOC OR 06/02/2014 Streamlining p4zsed >        IF( ln_dust ) THEN 
!         IF(lwp) WRITE(numout,*) '    initialize dust input from atmosphere '
!         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
!         !
!         ALLOCATE( sf_dust(1), STAT=ierr )           !* allocate and fill sf_sst (forcing structure) with sn_sst
!         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_init: unable to allocate sf_apr structure' )
!         !
!         CALL fld_fill( sf_dust, (/ sn_dust /), cn_dir, 'p4z_sed_init', 'Iron from sediment ', 'nampissed' )
!                                   ALLOCATE( sf_dust(1)%fnow(jpi,jpj,1)   )
!         IF( sn_dust%ln_tint )     ALLOCATE( sf_dust(1)%fdta(jpi,jpj,1,2) )
!         !
!         ! Get total input dust ; need to compute total atmospheric supply of Si in a year
!         CALL iom_open (  TRIM( sn_dust%clname ) , numdust )
!         CALL iom_gettime( numdust, zsteps, kntime=ntimes_dust)  ! get number of record in file
!         ALLOCATE( zdust(jpi,jpj,ntimes_dust) )
!         DO jm = 1, ntimes_dust
!            CALL iom_get( numdust, jpdom_data, TRIM( sn_dust%clvar ), zdust(:,:,jm), jm )
!         END DO
!         CALL iom_close( numdust )
!         sumdepsi = 0.e0
!         DO jm = 1, ntimes_dust
!            sumdepsi = sumdepsi + glob_sum( zdust(:,:,jm) * e1e2t(:,:) * tmask(:,:,1) ) 
!         ENDDO
!         sumdepsi = sumdepsi * r1_ryyss * 8.8 * 0.075 / 28.1 
!         DEALLOCATE( zdust)
!      ELSE
!         dust(:,:) = 0._wp
!         sumdepsi  = 0._wp
! <CMOC OR 06/02/2014 Streamlining p4zsed >        END IF

      ! iron mask/iron limitation for the ocean ! <CMOC OR 01/21/2014>
      ! ---------------------------------------
         IF(lwp) WRITE(numout,*) '    initialize CMOC iron limitation'
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         !
         ierr = 0 
         ALLOCATE( sf_fmsk(1), STAT=ierr )           !* allocate and fill sf_sst (forcing structure) with sn_sst
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_init: unable to allocate sf_fmsk structure' )
         ierr = 0
         CALL fld_fill( sf_fmsk, (/ sn_fmsk /), cn_dir, 'p4z_sed_init', 'Iron limitation mask', 'nampissed' )
                                   ALLOCATE( sf_fmsk(1)%fnow(jpi,jpj,1), STAT=ierr )
                                   IF( ierr > 0 ) THEN
                                            CALL ctl_stop( 'p4zsed: iron limitation mask, unable to allocate iron limitation array' )     ;    RETURN 
                                   ENDIF 
         ! IF( sn_fmsk%ln_tint )     ALLOCATE( sf_fmsk(1)%fdta(jpi,jpj,1,2) )
         !
         ! 
         CALL iom_open (  TRIM( sn_fmsk%clname ) , numfmsk )
         IF(lwp) WRITE(numout,*) '    iron mask file opened' 
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         CALL flush(numout)
         ALLOCATE(zfmask(jpi,jpj))
         CALL iom_get  ( numfmsk, jpdom_data, TRIM( sn_fmsk%clvar ), zfmask(:,:), 1 ) ! <CMOC OR 01/22/2014> get the data from the iron mask file without worrying about the time step
         !CALL fld_read( 1 , 1, sf_fmsk )
         IF(lwp) WRITE(numout,*) '    iron mask data read' 
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~'
         CALL flush(numout) 
         CALL iom_close( numfmsk )
           DO ji = 1, jpi
              DO jj = 1, jpj
                    xlimnfecmoc(ji,jj) = zfmask(ji,jj)
         !           IF(lwp) WRITE(numout,*) '    ji =', ji, ', jj=', jj, ', LFe=', xlimnfecmoc(ji,jj)
         !           IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         !           CALL flush(numout) 
              END DO
           END DO
         IF(lwp) WRITE(numout,*) '    iron mask data save in memory' 
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         CALL flush(numout)
         DEALLOCATE(zfmask) 
         IF(lwp) WRITE(numout,*) '    iron limitation initialization: done'
         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         CALL flush(numout)  
         !CALL fld_read( 1, 1, sf_fmsk )
         !xlimnfecmoc(:,:) = sf_fmsk(1)%fnow(:,:) ! No mask
         ! <CMOC OR 01/21/2014>

      ! nutrient input from rivers
      ! --------------------------
      IF( ln_river ) THEN
         ALLOCATE( sf_riverdic(1), STAT=ierr1 )           !* allocate and fill sf_sst (forcing structure) with sn_sst
         ALLOCATE( sf_riverdoc(1), STAT=ierr2 )           !* allocate and fill sf_sst (forcing structure) with sn_sst
         IF( ierr1 + ierr2 > 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_init: unable to allocate sf_apr structure' )
         !
         CALL fld_fill( sf_riverdic, (/ sn_riverdic /), cn_dir, 'p4z_sed_init', 'Input DOC from river ', 'nampissed' )
         CALL fld_fill( sf_riverdoc, (/ sn_riverdoc /), cn_dir, 'p4z_sed_init', 'Input DOC from river ', 'nampissed' )
                                   ALLOCATE( sf_riverdic(1)%fnow(jpi,jpj,1)   )
                                   ALLOCATE( sf_riverdoc(1)%fnow(jpi,jpj,1)   )
         IF( sn_riverdic%ln_tint ) ALLOCATE( sf_riverdic(1)%fdta(jpi,jpj,1,2) )
         IF( sn_riverdoc%ln_tint ) ALLOCATE( sf_riverdoc(1)%fdta(jpi,jpj,1,2) )
         ! Get total input rivers ; need to compute total river supply in a year
         CALL iom_open ( TRIM( sn_riverdic%clname ), numriv )
         CALL iom_gettime( numriv, zsteps, kntime=ntimes_riv)
         ALLOCATE( zriverdic(jpi,jpj,ntimes_riv) )   ;     ALLOCATE( zriverdoc(jpi,jpj,ntimes_riv) )
         DO jm = 1, ntimes_riv
            CALL iom_get( numriv, jpdom_data, TRIM( sn_riverdic%clvar ), zriverdic(:,:,jm), jm )
            CALL iom_get( numriv, jpdom_data, TRIM( sn_riverdoc%clvar ), zriverdoc(:,:,jm), jm )
         END DO
         CALL iom_close( numriv )
         ! N/P and Si releases due to coastal rivers
         ! -----------------------------------------
! <CMOC OR 06/02/2014 Streamlining p4zsed >  
!         rivpo4input = 0._wp 
!         rivalkinput = 0._wp 
!         DO jm = 1, ntimes_riv
!            rivpo4input = rivpo4input + glob_sum( ( zriverdic(:,:,jm) + zriverdoc(:,:,jm) ) * tmask(:,:,1) ) 
!            rivalkinput = rivalkinput + glob_sum(   zriverdic(:,:,jm)                       * tmask(:,:,1) ) 
!         END DO
!         rivpo4input = rivpo4input * 1E9 / 31.6_wp
!         rivalkinput = rivalkinput * 1E9 / 12._wp 
!         DEALLOCATE( zriverdic)   ;    DEALLOCATE( zriverdoc) 
! <CMOC OR 06/02/2014 Streamlining p4zsed >  
      ELSE
         rivinp(:,:) = 0._wp
         cotdep(:,:) = 0._wp
! <CMOC OR 06/02/2014 Streamlining p4zsed >  
!         rivpo4input = 0._wp
!         rivalkinput = 0._wp
! <CMOC OR 06/02/2014 Streamlining p4zsed >  
      END IF 

      ! nutrient input from dust
      ! ------------------------
! <CMOC OR 06/02/2014 Streamlining p4zsed >        IF( ln_ndepo ) THEN
!         IF(lwp) WRITE(numout,*) '    initialize the nutrient input by dust from ndeposition.orca.nc'
!         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!         ALLOCATE( sf_ndepo(1), STAT=ierr3 )           !* allocate and fill sf_sst (forcing structure) with sn_sst
!         IF( ierr3 > 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_init: unable to allocate sf_apr structure' )
!         !
!         CALL fld_fill( sf_ndepo, (/ sn_ndepo /), cn_dir, 'p4z_sed_init', 'Iron from sediment ', 'nampissed' )
!                                   ALLOCATE( sf_ndepo(1)%fnow(jpi,jpj,1)   )
!         IF( sn_ndepo%ln_tint )    ALLOCATE( sf_ndepo(1)%fdta(jpi,jpj,1,2) )
!         !
!         ! Get total input dust ; need to compute total atmospheric supply of N in a year
!         CALL iom_open ( TRIM( sn_ndepo%clname ), numdepo )
!         CALL iom_gettime( numdepo, zsteps, kntime=ntimes_ndep)
!         ALLOCATE( zndepo(jpi,jpj,ntimes_ndep) )
!         DO jm = 1, ntimes_ndep
!            CALL iom_get( numdepo, jpdom_data, TRIM( sn_ndepo%clvar ), zndepo(:,:,jm), jm )
!         END DO
!         CALL iom_close( numdepo )
!         nitdepinput = 0._wp
!         DO jm = 1, ntimes_ndep
!           nitdepinput = nitdepinput + glob_sum( zndepo(:,:,jm) * e1e2t(:,:) * tmask(:,:,1) ) 
!         ENDDO
!         nitdepinput = nitdepinput * 7.6 / 14E6 
!         DEALLOCATE( zndepo)
!      ELSE
!         nitdep(:,:) = 0._wp
!         nitdepinput = 0._wp
! <CMOC OR 06/02/2014 Streamlining p4zsed >        ENDIF

      ! coastal and island masks
      ! ------------------------
! <CMOC OR 06/02/2014 Streamlining p4zsed >        IF( ln_ironsed ) THEN     
!         IF(lwp) WRITE(numout,*) '    computation of an island mask to enhance coastal supply of iron'
!         IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!         CALL iom_open ( TRIM( sn_ironsed%clname ), numiron )
!         ALLOCATE( zcmask(jpi,jpj,jpk) )
!         CALL iom_get  ( numiron, jpdom_data, TRIM( sn_ironsed%clvar ), zcmask(:,:,:), 1 )
!         CALL iom_close( numiron )
!         !
!         DO jk = 1, 5
!            DO jj = 2, jpjm1
!               DO ji = fs_2, fs_jpim1
!                  IF( tmask(ji,jj,jk) /= 0. ) THEN
!                     zmaskt = tmask(ji+1,jj,jk) * tmask(ji-1,jj,jk) * tmask(ji,jj+1,jk)    &
!                        &                       * tmask(ji,jj-1,jk) * tmask(ji,jj,jk+1)
!                     IF( zmaskt == 0. )   zcmask(ji,jj,jk ) = MAX( 0.1, zcmask(ji,jj,jk) ) 
!                  END IF
!               END DO
!            END DO
!         END DO
!         CALL lbc_lnk( zcmask , 'T', 1. )      ! lateral boundary conditions on cmask   (sign unchanged)
!         DO jk = 1, jpk
!            DO jj = 1, jpj
!               DO ji = 1, jpi
!                  zexpide   = MIN( 8.,( fsdept(ji,jj,jk) / 500. )**(-1.5) )
!                  zdenitide = -0.9543 + 0.7662 * LOG( zexpide ) - 0.235 * LOG( zexpide )**2
!                  zcmask(ji,jj,jk) = zcmask(ji,jj,jk) * MIN( 1., EXP( zdenitide ) / 0.5 )
!               END DO
!            END DO
!         END DO
!         ! Coastal supply of iron
!         ! -------------------------
!         ironsed(:,:,jpk) = 0._wp
!         DO jk = 1, jpkm1
!            ironsed(:,:,jk) = sedfeinput * zcmask(:,:,jk) / ( fse3t(:,:,jk) * rday )
!         END DO
!         DEALLOCATE( zcmask)
!      ELSE
!         ironsed(:,:,:) = 0._wp
!      ENDIF
!      !
      IF( ll_sbc ) CALL p4z_sbc( nit000 ) 
!      !
!      IF(lwp) THEN 
!         WRITE(numout,*)
!         WRITE(numout,*) '    Total input of elements from river supply'
!         WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!         WRITE(numout,*) '    N Supply   : ', rivpo4input/7.6*1E3/1E12*14.,' TgN/yr'
!         WRITE(numout,*) '    Si Supply  : ', rivalkinput/6.*1E3/1E12*32.,' TgSi/yr'
!         WRITE(numout,*) '    Alk Supply : ', rivalkinput*1E3/1E12,' Teq/yr'
!         WRITE(numout,*) '    DIC Supply : ', rivpo4input*2.631*1E3*12./1E12,'TgC/yr'
!         WRITE(numout,*) 
!         WRITE(numout,*) '    Total input of elements from atmospheric supply'
!         WRITE(numout,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!         WRITE(numout,*) '    N Supply   : ', nitdepinput/7.6*1E3/1E12*14.,' TgN/yr'
!         WRITE(numout,*) 
! <CMOC OR 06/02/2014 Streamlining p4zsed >        ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_sed_init')
      !
   END SUBROUTINE p4z_sed_init

     INTEGER FUNCTION p4z_sed_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed_alloc  ***
      !!----------------------------------------------------------------------

      ! <CMOC OR 06/05/2014> Code Streamlining ! ALLOCATE( dust  (jpi,jpj), rivinp(jpi,jpj)     , cotdep(jpi,jpj),      &
      ! <CMOC OR 06/05/2014> &       nitdep(jpi,jpj), ironsed(jpi,jpj,jpk), STAT=p4z_sed_alloc )  
      ALLOCATE( rivinp(jpi,jpj)     , cotdep(jpi,jpj), STAT=p4z_sed_alloc )  

      IF( p4z_sed_alloc /= 0 ) CALL ctl_warn('p4z_sed_alloc : failed to allocate arrays.')

     END FUNCTION p4z_sed_alloc
#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_sed                         ! Empty routine
   END SUBROUTINE p4z_sed
#endif 

   !!======================================================================
END MODULE  p4zsed
