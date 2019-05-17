MODULE p4zlys
   !!======================================================================
   !!                         ***  MODULE p4zlys  ***
   !! TOP :   PISCES 
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr)  Calcon salinity dependence
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Improvment of calcite dissolution
   !!----------------------------------------------------------------------
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   'key_pisces'                                       PISCES bio-model
   !!----------------------------------------------------------------------
   !!   p4z_lys        :   Compute the CaCO3 dissolution 
   !!   p4z_lys_init   :   Read the namelist parameters
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_lys         ! called in trcsms_pisces.F90
   PUBLIC   p4z_lys_init    ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC :: kdca = 7.4E-3_wp   !: diss. rate constant calcite
   REAL(wp), PUBLIC :: nca  = 1.0_wp      !: order of reaction for calcite dissolution (not used)

   !! * Module variables
   REAL(wp) :: calcon = 1.03E-2           !: mean calcite concentration [Ca2+] in sea water [mole/kg solution]
   REAL(wp) :: r1_rday = 1.0_wp/86400._wp !: 1 / rday
 
   INTEGER  :: rmtss                      !: number of seconds per month 

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: p4zlys.F90 3321 2012-03-05 17:10:55Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p4z_lys( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_lys  ***
      !!
      !! ** Purpose :   CALCULATES DEGREE OF CACO3 SATURATION IN THE WATER
      !!                COLUMN, DISSOLUTION/PRECIPITATION OF CACO3 AND LOSS
      !!                OF CACO3 TO THE CACO3 SEDIMENT POOL.
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt ! ocean time step
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zalk, zdic, zph, zah2
      REAL(wp) ::   zdispot, zfact, zcalcon, zalka, zbot
      REAL(wp) ::   zph2, zph3, zpo4, zsi, zpd, zp0, zp1, zp3        ! coefficients added to account for P and Si contribution to TA
      REAL(wp) ::   zomegaca, zexcess, zexcess0
      REAL(wp) ::   ztmas, ztmas1
      REAL(wp) ::   zrfact
      CHARACTER (len=25) :: charout
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zco3, zcaldiss   
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('p4z_lys')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zco3, zcaldiss )
      !
      zco3    (:,:,:) = 0.
      zcaldiss(:,:,:) = 0.
      !     -------------------------------------------
      !     COMPUTE [CO3--] and [H+] CONCENTRATIONS
      !     -------------------------------------------
      
      DO jn = 1, 5                               !  BEGIN OF ITERATION
         !
!CDIR NOVERRCHK
         DO jk = 1, jpkm1
!CDIR NOVERRCHK
            DO jj = 1, jpj
!CDIR NOVERRCHK
               DO ji = 1, jpi

                  ztmas   = tmask_bgc_closea(ji,jj,jk)
                  ztmas1  = 1. - tmask_bgc_closea(ji,jj,jk)
                  zbot  = borat(ji,jj,1) * ztmas + 0.000416 * ztmas1 
                  zfact = rhop(ji,jj,1) / 1000. + rtrn
                  zdic  = trn(ji,jj,1,jpdic) / zfact * ztmas + 2. * ztmas1
                  zph   = MAX( hi(ji,jj,1), 1.e-10 ) / zfact * ztmas + 1.e-9 * ztmas1
                  zalka = trn(ji,jj,1,jptal) / zfact * ztmas + 2.2 * ztmas1
                  zph2 = zph*zph
                  zph3 = zph*zph2
                  zpo4 = (trn(ji,jj,1,jpno3)+trn(ji,jj,1,jpnh4)) / 16. *0.000001 / zfact
                  zsi = asi3(ji,jj,1) * 0.000001 / zfact                        ! silica is a static array based on initialization file, not a carried tracer

               ! CALCULATE P AND Si ION CONCENTRATIONS AS PER ORR ET AL (BPG EQUATIONS 43-47)
               ! zp3 = H3PO4, zp1 = HPO4(2-), zp0 = PO4(3-): denominator is the same for all 3 equations
                  zpd = 1./ ( zph3 + akp13(ji,jj,1)*zph2 + akp13(ji,jj,1)*akp23(ji,jj,1)*zph + akp13(ji,jj,1)*akp23(ji,jj,1)*akp33(ji,jj,1) )
                  zp3 = zph3*zpo4 * zpd
                  zp1 = zph*zpo4*akp13(ji,jj,1)*akp23(ji,jj,1) * zpd
                  zp0 = zpo4*akp13(ji,jj,1)*akp23(ji,jj,1)*akp33(ji,jj,1) * zpd
                  zsi = zsi / (1. + zph / aksi3(ji,jj,1))

               ! CALCULATE [ALK]([CO3--], [HCO3-])
                  zalk  = zalka - (  akw3(ji,jj,1) / zph - zph + zbot / ( 1.+ zph / akb3(ji,jj,1) ) + 2.*zp0 + zp1 - zp3 + zsi )

               ! CALCULATE [H+] AND [H2CO3]
                  zah2   = SQRT(  (zdic-zalk)*(zdic-zalk) + 4.* ( zalk * ak23(ji,jj,1)   &
                     &                                        / ak13(ji,jj,1) ) * ( 2.* zdic - zalk )  )
                  zah2   = 0.5 * ak13(ji,jj,1) / zalk * ( ( zdic - zalk ) + zah2 )
                  zco3(ji,jj,jk) = zalk / ( 2. + zah2 / ak23(ji,jj,jk) ) * zfact
                  hi(ji,jj,jk)   = zah2 * zfact

               END DO
            END DO
         END DO
         !
      END DO 

      !     ---------------------------------------------------------
      !        CALCULATE DEGREE OF CACO3 SATURATION AND CORRESPONDING
      !        DISSOLOUTION AND PRECIPITATION OF CACO3 (BE AWARE OF
      !        MGCO3)
      !     ---------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! DEVIATION OF [CO3--] FROM SATURATION VALUE
               ! Salinity dependance in zomegaca and divide by rhop/1000 to have good units
               zcalcon  = calcon * ( tsn(ji,jj,jk,jp_sal) / 35._wp )
               zfact    = rhop(ji,jj,jk) / 1000._wp
               zomegaca = ( zcalcon * zco3(ji,jj,jk) * zfact ) / aksp(ji,jj,jk) 

               ! SET DEGREE OF UNDER-/SUPERSATURATION
               excess(ji,jj,jk) = 1._wp - zomegaca
               !zexcess0 = MAX( 0., excess(ji,jj,jk) )
               !zexcess  = zexcess0

               ! AMOUNT CACO3 (12C) THAT RE-ENTERS SOLUTION
               !       (ACCORDING TO THIS FORMULATION ALSO SOME PARTICULATE
               !       CACO3 GETS DISSOLVED EVEN IN THE CASE OF OVERSATURATION)
               zdispot = kdca * r1_rday * trn(ji,jj,jk,jpcal)
              !  CHANGE OF [CO3--] , [ALK], PARTICULATE [CACO3],
              !       AND [SUM(CO2)] DUE TO CACO3 DISSOLUTION/PRECIPITATION
              zcaldiss(ji,jj,jk)  = zdispot                         ! calcite dissolution (first order, no saturation state dependence)
              zco3(ji,jj,jk)      = zco3(ji,jj,jk) + zcaldiss(ji,jj,jk) * rfact
              !
              tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) -         zcaldiss(ji,jj,jk)
              tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 2.E-6 * zcaldiss(ji,jj,jk)       
              tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + 1.E-6 * zcaldiss(ji,jj,jk) 
            END DO
         END DO
      END DO
      !
      IF( ln_diatrc )  THEN
         !
         IF( lk_iomput ) THEN
            CALL iom_put( "PH"    , -1. * LOG10( hi(:,:,:) )                * tmask_bgc_closea(:,:,:) )
            CALL iom_put( "CO3"   ,        zco3    (:,:,:) * 1e+3           * tmask_bgc_closea(:,:,:) )
            CALL iom_put( "CO3sat",        aksp    (:,:,:) * 1e+3 / calcon  * tmask_bgc_closea(:,:,:) )
            CALL iom_put( "DCAL"  ,        zcaldiss(:,:,:) * 1.e-3          * tmask_bgc_closea(:,:,:) )       ! conversion from umol/L/s to mol/m3/s
         ELSE
            trc3d(:,:,:,jp_pcs0_3d    ) = -1. * LOG10( hi(:,:,:) ) * tmask_bgc_closea(:,:,:)
            trc3d(:,:,:,jp_pcs0_3d + 1) = zco3(:,:,:)              * tmask_bgc_closea(:,:,:)
            trc3d(:,:,:,jp_pcs0_3d + 2) = aksp(:,:,:) / calcon     * tmask_bgc_closea(:,:,:)
         ENDIF
         !
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('lys ')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zco3, zcaldiss )
      !
      IF( nn_timing == 1 )  CALL timing_stop('p4z_lys')
      !
   END SUBROUTINE p4z_lys

   SUBROUTINE p4z_lys_init

      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_lys_init  ***
      !!
      !! ** Purpose :   Initialization of CaCO3 dissolution parameters
      !!
      !! ** Method  :   Read the nampiscal namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiscal
      !!
      !!----------------------------------------------------------------------

      NAMELIST/nampiscal/ kdca, nca

      REWIND( numnatp )                     ! read numnatp
      READ  ( numnatp, nampiscal )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for CaCO3 dissolution, nampiscal'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    diss. rate constant calcite (per month)   kdca      =', kdca
         WRITE(numout,*) '    order of reaction for calcite dissolution nca       =', nca
      ENDIF

      ! Number of seconds per month 
      rmtss =  nyear_len(1) * rday / raamo

   END SUBROUTINE p4z_lys_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_lys( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p4z_lys: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_lys
#endif 
   !!======================================================================
END MODULE  p4zlys
