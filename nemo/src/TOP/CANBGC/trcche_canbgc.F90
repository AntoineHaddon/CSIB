MODULE trcche_canbgc
   !!======================================================================
   !!                         ***  MODULE trcche  ***
   !! TOP :   Sea water chemistry computed following OMIP protocol
   !!======================================================================
   !! History :   OPA  !  1988     (E. Maier-Reimer)  Original code
   !!              -   !  1998     (O. Aumont)  addition
   !!              -   !  1999     (C. Le Quere)  modification
   !!   NEMO      1.0  !  2004     (O. Aumont)  modification
   !!              -   !  2006     (R. Gangsto)  modification
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J.Orr ) update O2 solubility constants
   !!                  !  2022-2023(J. Christian, O. Riche) NEMO4 integration
   !!----------------------------------------------------------------------
   !!   trcche      :  Sea water chemistry computed following OMIP protocol
   !!----------------------------------------------------------------------

   USE oce_trc           !  shared variables between ocean and passive tracers
   USE trc               !  passive tracers common variables
   USE sms_top_canbgc    !  TOP Source Minus Sink variables and miscellaneous.
   USE lib_mpp           !  MPP library

   USE in_out_manager    ! in_out_manager grants access to numout file ID
   !USE iom                       ! to access iom_put for diagnostics
   
   USE trc_closea_canbgc ! bgc-specific closea mask
   USE trcsrc_canbgc     ! external sources module

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_che_2D         !
   PUBLIC   trc_che_3D         !
   PUBLIC   trc_che_init_2D    !
   PUBLIC   trc_che_init_3D    !
   PUBLIC   trc_che_alloc      !

   ! O2 and CO2 arrays to be passed to trc_flx
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   K0CO2    
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   K0O2
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qh2co3
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qco3
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qomegac
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qaksp

! Constants and conversion factors
   REAL(wp), PUBLIC ::     xconv0   = 0.01_wp / 3600._wp !: coefficients for conversion 
   REAL(wp), PUBLIC ::     atcoxy   = 0.20946_wp         !: O2 fraction of air partial pressure (atm / atm)
   REAL(wp), PUBLIC :: no3_sf         ! scaling factor for NO3 (for estimation of PO4)

   REAL(wp) ::   salchl = 1. / 1.80655    ! conversion factor for salinity --> chlorinity (Wooster et al. 1969)
   REAL(wp) ::   o2atm  = 1. / ( 1000. * 0.20946 )  

   REAL(wp) ::   akcc1  = -171.9065       ! coeff. for apparent solubility equilibrium
   REAL(wp) ::   akcc2  =   -0.077993     ! Millero et al. 1995 from Mucci 1983
   REAL(wp) ::   akcc3  = 2839.319        
   REAL(wp) ::   akcc4  =   71.595        
   REAL(wp) ::   akcc5  =   -0.77712      
   REAL(wp) ::   akcc6  =    0.00284263   
   REAL(wp) ::   akcc7  =  178.34        
   REAL(wp) ::   akcc8  =   -0.07711     
   REAL(wp) ::   akcc9  =    0.0041249   

   REAL(wp) ::   rgas   = 83.143         ! universal gas constants
   REAL(wp) ::   oxyco  = 1. / 22.4144   ! converts from liters of an ideal gas to moles

   REAL(wp) ::   bor1   = 0.000232       ! constants for calculating borate concentration
   REAL(wp) ::   bor2   = 1. / 10.811

   REAL(wp) ::   calcium = 1.03e-2       ! calcium fraction of sea salt

   REAL(wp) ::   ca0    = -160.7333   ! WEISS & PRICE 1980, units mol/(kg atm)
   REAL(wp) ::   ca1    =  215.4152
   REAL(wp) ::   ca2    =  89.8920
   REAL(wp) ::   ca3    = - 1.47759
   REAL(wp) ::   ca4    =   0.029941
   REAL(wp) ::   ca5    = - 0.027455
   REAL(wp) ::   ca6    =   0.0053407

   REAL(wp) ::   c10    = -3633.86        ! Coeff. for 1. dissoc. of carbonic acid (Dickson et al., 2007)
   REAL(wp) ::   c11    =    61.2172    
   REAL(wp) ::   c12    =    -9.67770  
   REAL(wp) ::   c13    =     0.011555     
   REAL(wp) ::   c14    =    -0.0001152

   REAL(wp) ::   c20    =  -471.78       ! coeff. for 2. dissoc. of carbonic acid (Dickson et al., 2007)   
   REAL(wp) ::   c21    =   -25.9290   
   REAL(wp) ::   c22    =     3.16967   
   REAL(wp) ::   c23    =     0.01781
   REAL(wp) ::   c24    =    -0.0001122

   REAL(wp) ::   cb0    = -8966.90      ! Coeff. for 1. dissoc. of boric acid 
   REAL(wp) ::   cb1    = -2890.53      ! (Dickson and Goyet, 1994)
   REAL(wp) ::   cb2    =   -77.942
   REAL(wp) ::   cb3    =     1.728
   REAL(wp) ::   cb4    =    -0.0996
   REAL(wp) ::   cb5    =   148.0248
   REAL(wp) ::   cb6    =   137.1942
   REAL(wp) ::   cb7    =     1.62142
   REAL(wp) ::   cb8    =   -24.4344
   REAL(wp) ::   cb9    =   -25.085
   REAL(wp) ::   cb10   =    -0.2474 
   REAL(wp) ::   cb11   =     0.053105

   REAL(wp) ::   cw0    = -13847.26     ! Coeff. for dissoc. of water (Dickson and Riley, 1979 )
   REAL(wp) ::   cw1    =    148.9652  
   REAL(wp) ::   cw2    =    -23.6521
   REAL(wp) ::   cw3    =    118.67 
   REAL(wp) ::   cw4    =     -5.977 
   REAL(wp) ::   cw5    =      1.0495  
   REAL(wp) ::   cw6    =     -0.01615

   REAL(wp) ::   cp10    = -4576.752    ! Coeff. for dissoc. of H3PO4 (Dickson et al 2007)
   REAL(wp) ::   cp11    = 115.525
   REAL(wp) ::   cp12    = -18.453
   REAL(wp) ::   cp13    = -106.736
   REAL(wp) ::   cp14    = 0.69171
   REAL(wp) ::   cp15    = -0.65643
   REAL(wp) ::   cp16    = -0.01844
       
   REAL(wp) ::   cp20    = -8814.715    ! Coeff. for dissoc. of H2PO4- (Dickson et al 2007)
   REAL(wp) ::   cp21    = 172.0883
   REAL(wp) ::   cp22    = -27.927
   REAL(wp) ::   cp23    = -160.340
   REAL(wp) ::   cp24    = 1.3566
   REAL(wp) ::   cp25    = 0.37335
   REAL(wp) ::   cp26    = -0.05778
      
   REAL(wp) ::   cp30    = -3070.75     ! Coeff. for dissoc. of HPO4-- (Dickson et al 2007)
   REAL(wp) ::   cp31    = -18.141
   REAL(wp) ::   cp32    = 17.27039
   REAL(wp) ::   cp33    = 2.81197
   REAL(wp) ::   cp34    = -44.99486
   REAL(wp) ::   cp35    = -0.09984
   
   REAL(wp) ::   csi0    = -8904.2      ! Coeff. for dissoc. of Si(OH)4 (Dickson et al 2007)
   REAL(wp) ::   csi1    = 117.385
   REAL(wp) ::   csi2    = -19.334
   REAL(wp) ::   csi3    = -458.79
   REAL(wp) ::   csi4    = 3.5913
   REAL(wp) ::   csi5    = 188.74
   REAL(wp) ::   csi6    = -1.5998
   REAL(wp) ::   csi7    = -12.1652
   REAL(wp) ::   csi8    = 0.07871
   REAL(wp) ::   csi9    = -0.001005

   !                                    ! volumetric solubility constants for o2 in ml/L  
   REAL(wp) ::   ox0    =  2.00856      ! from Table 1 for Eq 8 of Garcia and Gordon, 1992.
   REAL(wp) ::   ox1    =  3.22400      ! corrects for moisture and fugacity, but not total atmospheric pressure
   REAL(wp) ::   ox2    =  3.99063      !      Original PISCES code noted this was a solubility, but 
   REAL(wp) ::   ox3    =  4.80299      ! was in fact a bunsen coefficient with units L-O2/(Lsw atm-O2)
   REAL(wp) ::   ox4    =  9.78188e-1   ! Hence, need to divide EXP( zoxy ) by 1000, ml-O2 => L-O2
   REAL(wp) ::   ox5    =  1.71069      ! and atcoxy = 0.20946 to add the 1/atm dimension.
   REAL(wp) ::   ox6    = -6.24097e-3   
   REAL(wp) ::   ox7    = -6.93498e-3 
   REAL(wp) ::   ox8    = -6.90358e-3
   REAL(wp) ::   ox9    = -4.29155e-3 
   REAL(wp) ::   ox10   = -3.11680e-7 

   REAL(wp) :: devk10  = -25.5          ! coefficients for pressure correction: this is done quite differently in NEMO4 PISCES vs NEMO3.4
   REAL(wp) :: devk11  = -15.82         ! 0, 1 = K1, K2
   REAL(wp) :: devk12  = -29.48         ! 2 = KB
   REAL(wp) :: devk13  = -20.02         ! 3 = Kw
   REAL(wp) :: devk14  = -18.03         ! 4 = K_SO4 (not used)
   REAL(wp) :: devk15  = -9.78          ! 5 = K_F (not used)
   REAL(wp) :: devk16  = -48.76         ! 6 = Ksp
   REAL(wp) :: devk17  = -14.51         ! 7 = K_P1
   REAL(wp) :: devk18  = -23.12         ! 8 = K_P2
   REAL(wp) :: devk19  = -26.57         ! 9 = K_P3
   REAL(wp) :: devk110  = -29.48        ! 10 = K_Si

   REAL(wp) :: devk20  = 0.1271
   REAL(wp) :: devk21  = -0.0219
   REAL(wp) :: devk22  = 0.1622
   REAL(wp) :: devk23  = 0.1119
   REAL(wp) :: devk24  = 0.0466
   REAL(wp) :: devk25  = -0.0090
   REAL(wp) :: devk26  = 0.5304
   REAL(wp) :: devk27  = 0.1211
   REAL(wp) :: devk28  = 0.1758
   REAL(wp) :: devk29  = 0.2020
   REAL(wp) :: devk210  = 0.1622

   REAL(wp) :: devk30  = 0.
   REAL(wp) :: devk31  = 0.
   REAL(wp) :: devk32  = 2.608E-3
   REAL(wp) :: devk33  = -1.409e-3
   REAL(wp) :: devk34  = 0.316e-3
   REAL(wp) :: devk35  = -0.942e-3
   REAL(wp) :: devk36  = 0.
   REAL(wp) :: devk37  = -0.321e-3
   REAL(wp) :: devk38  = -2.647e-3
   REAL(wp) :: devk39  = -3.042e-3
   REAL(wp) :: devk310  = -2.6080e-3

   REAL(wp) :: devk40  = -3.08E-3
   REAL(wp) :: devk41  = 1.13E-3
   REAL(wp) :: devk42  = -2.84E-3
   REAL(wp) :: devk43  = -5.13E-3
   REAL(wp) :: devk44  = -4.53e-3
   REAL(wp) :: devk45  = -3.91e-3
   REAL(wp) :: devk46  = -11.76e-3
   REAL(wp) :: devk47  = -2.67e-3
   REAL(wp) :: devk48  = -5.15e-3
   REAL(wp) :: devk49  = -4.08e-3
   REAL(wp) :: devk410  = -2.84e-3

   REAL(wp) :: devk50  = 0.0877E-3
   REAL(wp) :: devk51  = -0.1475E-3     
   REAL(wp) :: devk52  = 0.
   REAL(wp) :: devk53  = 0.0794E-3      
   REAL(wp) :: devk54  = 0.09e-3
   REAL(wp) :: devk55  = 0.054e-3
   REAL(wp) :: devk56  = 0.3692E-3
   REAL(wp) :: devk57  = 0.0427e-3
   REAL(wp) :: devk58  = 0.09e-3
   REAL(wp) :: devk59  = 0.0714e-3
   REAL(wp) :: devk510  = 0.0

   !!* Substitution
#  include "domzgr_substitute.h90"

CONTAINS

   SUBROUTINE trc_che_2D( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_che_2D  ***
      !!
      !! ** Purpose :   Sea water chemistry and CO2/O2 solubility computed following OMIP protocol (surface only)
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER, INTENT(in) ::   Kmm  ! time level indices
      INTEGER  ::   ji, jj, jk, jm
      REAL(wp) ::   ztkel, zt, zt2, zsal, zsal2
      REAL(wp) ::   ztgg, ztgg2, ztgg3, ztgg4, ztgg5
      REAL(wp) ::   zoxy, zlogt, zcek1
      REAL(wp) ::   zph, zah2, zbot, zdic, zcalk, ztalk, zfact
      REAL(wp) ::   zpo4, zsi
      REAL(wp) ::   zak1, zak2, zakb, zakw, zakp1, zakp2, zakp3, zaksi
      REAL(wp) ::   ztmas, ztmas1
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: hi
      !REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zph0
      REAL(wp), DIMENSION(2) :: hion_CA
      !!---------------------------------------------------------------------

      IF( ln_timing )  CALL timing_start('trc_che_2D')
      !
      ! ----------------------------------
      !
      ALLOCATE( hi(jpi, jpj) )
      !
      hi(:,:)=1.e-9
      CALL trc_src3d(kt,js3d_si)
      qasi3=src3d_dta(:,:,:,js3d_si)

! solubility of CO2 and O2 in seawater (note the names of these arrays are different from the PISCES-based one used in CanESM5)
      DO jj = 1, jpj
         DO ji = 1, jpi

            ztmas   = tmask_bgc_closea(ji,jj,1)
            ztmas1  = 1. - tmask_bgc_closea(ji,jj,1)
            !                             ! SET ABSOLUTE TEMPERATURE
            ztkel = ts(ji,jj,1,jp_tem,Kmm) + 273.15
            zt    = ztkel * 0.01
            zt2   = zt * zt
            !
            zsal  = ts(ji,jj,1,jp_sal,Kmm)*ztmas + ztmas1*35.
            zsal2 = zsal * zsal
            zlogt = LOG( zt )
            !                             ! LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1980)
            !                             !     AND FOR THE ATMOSPHERE FOR NON IDEAL GAS
            zcek1 = ca0 + ca1 / zt + ca2 * zlogt + ca3 * zt2 + zsal * ( ca4 + ca5 * zt + ca6 * zt2 )
            !                             ! LN(K0) OF SOLUBILITY OF O2 and N2 in ml/L (EQ. 8, GARCIA AND GORDON, 1992)
            ztgg  = LOG( ( 298.15 - ts(ji,jj,1,jp_tem,Kmm) ) / ztkel )  ! Set the GORDON & GARCIA scaled temperature
            ztgg2 = ztgg  * ztgg
            ztgg3 = ztgg2 * ztgg
            ztgg4 = ztgg3 * ztgg
            ztgg5 = ztgg4 * ztgg
            zoxy  = ox0 + ox1 * ztgg + ox2 * ztgg2 + ox3 * ztgg3 + ox4 * ztgg4 + ox5 * ztgg5   &
                   + zsal * ( ox6 + ox7 * ztgg + ox8 * ztgg2 + ox9 * ztgg3 ) +  ox10 * zsal2

            !                             ! SET SOLUBILITIES OF O2 AND CO2 
            K0CO2(ji,jj) = EXP( zcek1 ) * rhop(ji,jj,1) / 1000.         ! mol/(L atm)
            K0O2(ji,jj) = ( EXP( zoxy  ) * o2atm ) * oxyco              ! mol/(L atm)
            !
         END DO
      END DO

      ! SURFACE CHEMISTRY (PCO2 AND [H+] IN
      !     SURFACE LAYER); THE RESULT OF THIS CALCULATION
      !     IS USED TO COMPUTE AIR-SEA FLUX OF CO2

      DO jm = 1, 10
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! DUMMY VARIABLES FOR DIC, TA, AND BORATE
               ztmas = tmask_bgc_closea(ji,jj,1)
               ztmas1 = 1. - tmask_bgc_closea(ji,jj,1)
               zfact = rhop(ji,jj,1) / 1000. + rtrn
               zbot = qborat2(ji,jj) * ztmas + 0.000416 * ztmas1 
               zdic = tr(ji,jj,1,jqdic, Kmm) / zfact * ztmas + 0.002 * ztmas1
               ztalk = tr(ji,jj,1,jqtal, Kmm) / zfact * ztmas + 0.0024 * ztmas1
               ! add 400 uM "guardrail" to prevent extreme pCO2 in runoff-dominated environments
               zdic = MAX(zdic,0.0004)
               ztalk = MAX(ztalk,0.0004)

               zpo4 = tr(ji,jj,1,jqno3, Kmm) * no3_sf / 16. / zfact                        ! needs to include NH4 for CanOE when available
               zsi = qasi3(ji,jj,1) * 0.000001 / zfact                        ! silica is a static array based on initialization file, not a carried tracer

               ! initialize local scalar variables so that calculations are identical in 2D and 3D SRs
               zakp1=qakp12(ji,jj)
               zakp2=qakp22(ji,jj)
               zakp3=qakp32(ji,jj)
               zaksi=qaksi2(ji,jj)
               zakw=qakw2(ji,jj)
               zakb=qakb2(ji,jj)
               zak1=qak12(ji,jj)
               zak2=qak22(ji,jj)

               ! C chemistry calculations have been placed in a FUNCTION to prevent divergence between trc_che_2D and trc_che_3D
               ! this function returns the [H+] (1) and the carbonate alkalinity (2)
               zph = MAX( hi(ji,jj), 1.e-10 ) / zfact * ztmas + 1.e-9 * ztmas1
               hion_CA=hplus(zdic,ztalk,zph,zbot,zpo4,zsi,zak1,zak2,zakb,zakw,zakp1,zakp2,zakp3,zaksi)
               zah2=hion_CA(1)
               zcalk=hion_CA(2)
               ! record [H+] in L^-1 for next iteration
               hi(ji,jj) = zah2 * zfact
               ! calculate [CO2*] for export to gas flux SR
               qh2co3(ji,jj) = ( 2.* zdic - zcalk ) / ( 2.+ zak1 / zah2 ) * zfact
               qhi(ji,jj,1) = hi(ji,jj)    ! OR Jan 19th 2023
            END DO
         END DO
      END DO
      !
      ! OR Jan 24th 2023
      ! Moving pH diagnostics here
!      IF ( ln_cmoc ) THEN
!       ALLOCATE( zph0(jpi,jpj,jpk) )
!       zph0(:,:,:) = 0. 
!       DO jj= 1, jpj
!        DO ji= 1, jpi
!          zph0(ji,jj,1) = -1. * LOG10( MAX( qhi(ji,jj,1) + rtrn , rtrn ) ) 
!        END DO
!       END DO
!       CALL iom_put("pH", zph0(:,:,:) * tmask_bgc_closea(:,:,:))
!       DEALLOCATE( zph0 )
!      END IF
      ! 
      DEALLOCATE( hi )
      !
      IF( ln_timing )  CALL timing_stop('trc_che_2D')
      !
   END SUBROUTINE trc_che_2D

   SUBROUTINE trc_che_3D( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_che_3D  ***
      !!
      !! ** Purpose :   Sea water chemistry computed following OMIP protocol (all z levels)
      !!
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER, INTENT(in) ::   Kmm  ! time level indices
      INTEGER  ::   ji, jj, jk, jm
      REAL(wp) ::   zph, zah2, zbot, zdic, zcalk, ztalk, zfact, zcalcon
      REAL(wp) ::   zpo4, zsi
      REAL(wp) ::   zak1, zak2, zakb, zakw, zakp1, zakp2, zakp3, zaksi
      REAL(wp) ::   ztmas, ztmas1
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: hi
      !REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zph0
      REAL(wp), DIMENSION(2) :: hion_CA
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('trc_che_3D')
      !
      ALLOCATE( hi(jpi, jpj, jpk) )
      !
      !     -------------------------------------------
      !     COMPUTE [CO3--] and [H+] CONCENTRATIONS
      !     -------------------------------------------
      
      hi(:,:,:)=1.e-9
      qomegac(:,:,:) = 0.e0
      CALL trc_src3d(kt,js3d_si)
      qasi3=src3d_dta(:,:,:,js3d_si)

      DO jm = 1, 5                              !  BEGINNING OF ITERATION
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi

                  ztmas = tmask_bgc_closea(ji,jj,jk)
                  ztmas1 = 1. - tmask_bgc_closea(ji,jj,jk)
                  zfact = rhop(ji,jj,jk) / 1000. + rtrn
                  zbot = qborat3(ji,jj,jk) * ztmas + 0.000416 * ztmas1 
                  zdic = tr(ji,jj,jk,jqdic, Kmm) / zfact * ztmas + 0.002 * ztmas1
                  ztalk = tr(ji,jj,jk,jqtal, Kmm) / zfact * ztmas + 0.0024 * ztmas1
                  zdic = MAX(zdic,0.0004)
                  ztalk = MAX(ztalk,0.0004)
                  zpo4 = tr(ji,jj,jk,jqno3, Kmm) * no3_sf / 16. / zfact               ! needs to include NH4 for CanOE when available
                  zsi = qasi3(ji,jj,jk) * 0.000001 / zfact                        ! silica is a static array based on initialization file, not a carried tracer

               ! initialize local scalar variables so that calculations are identical in 2D and 3D SRs
                  zakp1=qakp13(ji,jj,jk)
                  zakp2=qakp23(ji,jj,jk)
                  zakp3=qakp33(ji,jj,jk)
                  zaksi=qaksi3(ji,jj,jk)
                  zakw=qakw3(ji,jj,jk)
                  zakb=qakb3(ji,jj,jk)
                  zak1=qak13(ji,jj,jk)
                  zak2=qak23(ji,jj,jk)

               ! C chemistry calculations have been placed in a FUNCTION to prevent divergence between trc_che_2D and trc_che_3D
               ! this function returns the [H+] (1) and the carbonate alkalinity (2)
                  zph = MAX( hi(ji,jj,jk), 1.e-10 ) / zfact * ztmas + 1.e-9 * ztmas1
                  hion_CA=hplus(zdic,ztalk,zph,zbot,zpo4,zsi,zak1,zak2,zakb,zakw,zakp1,zakp2,zakp3,zaksi)
                  zah2=hion_CA(1)
                  zcalk=hion_CA(2)
               ! record [H+] in L^-1 for next iteration
                  hi(ji,jj,jk) = zah2 * zfact
               ! calculate [CO3--] and [H+] for export to other SR's
                  qco3(ji,jj,jk) = zcalk / ( 2. + zah2 / zak2 )     ! no conversion to mol L^-1 as it is not applied to Ksp
                  zcalcon  = calcium * ( ts(ji,jj,jk,jp_sal,Kmm) / 35._wp )
                  zfact    = rhop(ji,jj,jk) / 1000._wp
                  qomegac(ji,jj,jk) = ( zcalcon * qco3(ji,jj,jk) * zfact ) / qaksp(ji,jj,jk)
                  WRITE(numout,*) "ji,jj,jk,qomegac(ji,jj,jk):",ji,jj,jk,qomegac(ji,jj,jk)
                  qhi(ji,jj,jk) = hi(ji,jj,jk)     ! OR Jan 19th 2023

               END DO
            END DO
         END DO
         !
      END DO 
      !
!      ALLOCATE( zph0(jpi,jpj,jpk) )
!      zph0(:,:,:) = rtrn 
!      DO jk= 1, jpk
!        DO jj= 1, jpj
!          DO ji= 1, jpi
!          zph0(ji,jj,jk) = -1. * LOG10( MAX( qhi(ji,jj,jk) + rtrn , rtrn ) ) 
!          END DO
!        END DO
!      END DO
!      CALL iom_put("pH", zph0(:,:,:) * tmask_bgc_closea(:,:,:))
!      DEALLOCATE( zph0 )      
      !
      DEALLOCATE( hi )
      !
      IF( ln_timing )  CALL timing_stop('trc_che_3D')

   END SUBROUTINE trc_che_3D

   SUBROUTINE trc_che_init_2D( Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_che_init_2D  ***
      !!
      !! ** Purpose :   Calculate surface values of dissociation constants
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kmm  ! time level indices
      INTEGER  ::   ji, jj
      REAL(wp) ::   ztkel, zsal, zsal2
      REAL(wp) ::   ztc, zcl
      REAL(wp) ::   zsqrt, ztr, zlogt
      REAL(wp) ::   zis, zis2, zsal15, zisqrt
      REAL(wp) ::   zak1, zak2, zakb, zakw, zakp1, zakp2, zakp3, zaksi
      REAL(wp) ::   zck1, zck2, zckw, zckb, zckp1, zckp2, zckp3, zcksi
      REAL(wp) ::   ztmas, ztmas1
      !!---------------------------------------------------------------------

      IF( ln_timing )  CALL timing_start('trc_che_init_2D')
      !
      ! ----------------------------------

          DO jj = 1, jpj
             DO ji = 1, jpi
 
                ztmas   = tmask_bgc_closea(ji,jj,1)
                ztmas1  = 1. - tmask_bgc_closea(ji,jj,1)
                ! SET ABSOLUTE TEMPERATURE
                ztkel   = ts(ji,jj,1,jp_tem,Kmm) + 273.15
                zsal  = ts(ji,jj,1,jp_sal,Kmm)*ztmas + ztmas1*35.
                zsqrt  = SQRT( zsal )
                zsal15  = zsqrt * zsal
                zlogt  = LOG( ztkel )
                ztr    = 1. / ztkel
                zis    = 19.924 * zsal / ( 1000.- 1.005 * zsal )
                zis2   = zis * zis
                zisqrt = SQRT( zis )
                ztc     = ts(ji,jj,1,jp_tem,Kmm)
                ! CHLORINITY (WOOSTER ET AL., 1969)
                zcl     = zsal * salchl
 
                ! DISSOCIATION CONSTANT FOR CARBONATE AND BORATE
                zckb    = ( cb0 + cb1 * zsqrt + cb2  * zsal + cb3 * zsal15 + cb4 * zsal * zsal ) * ztr   &
                   &    + ( cb5 + cb6 * zsqrt + cb7  * zsal )                                            &
                   &    + ( cb8 + cb9 * zsqrt + cb10 * zsal ) * zlogt + cb11 * zsqrt * ztkel 
 
                zck1    = c10 * ztr + c11 + c12 * zlogt + c13 * zsal + c14 * zsal * zsal
                zck2    = c20 * ztr + c21 + c22 * zlogt + c23 * zsal + c24 * zsal * zsal
 
                ! PKW (H2O) (DICKSON AND RILEY, 1979)
                zckw    = cw0 * ztr + cw1 + cw2 * zlogt + ( cw3 * ztr + cw4 + cw5 * zlogt ) * zsqrt + cw6 * zsal
 
                ! DISSOCIATION CONSTANTS FOR PHOSPHATE AND SILICATE
                zckp1    = cp10 * ztr + cp11 + cp12 * zlogt + (cp13 * ztr + cp14) * zsqrt + (cp15 * ztr + cp16) * zsal
                zckp2    = cp20 * ztr + cp21 + cp22 * zlogt + (cp23 * ztr + cp24) * zsqrt + (cp25 * ztr + cp26) * zsal
                zckp3    = cp30 * ztr + cp31 + (cp32 * ztr + cp33) * zsqrt + (cp34 * ztr + cp35) * zsal
 
                zcksi = csi0 * ztr + csi1 + csi2 * zlogt + (csi3 * ztr + csi4) * zisqrt + (csi5 * ztr + csi6) * zis &
                   & + (csi7 * ztr + csi8) * zis2 + LOG(1. + csi9 * zsal)
 
                zak1    = 10.**(zck1)
                zak2    = 10.**(zck2)
                zakb    = EXP( zckb  )
                zakw    = EXP( zckw )
                zakp1    = EXP( zckp1 )
                zakp2    = EXP( zckp2 )
                zakp3    = EXP( zckp3 )
                zaksi    = EXP( zcksi )

                ! 2D arrays of dissociation constants (no pressure correction)
                qak12(ji,jj) = zak1
                qak22(ji,jj) = zak2
                qakb2(ji,jj) = zakb
                qakw2(ji,jj) = zakw
                qakp12(ji,jj) = zakp1
                qakp22(ji,jj) = zakp2
                qakp32(ji,jj) = zakp3
                qaksi2(ji,jj) = zaksi
 
                ! TOTAL BORATE CONCENTR. (mol kg^-1)
                qborat2(ji,jj) = bor1 * zcl * bor2
 
             END DO
          END DO

      ! -------------------------------
      !
      IF( ln_timing )  CALL timing_stop('trc_che_init_2D')
      !
   END SUBROUTINE trc_che_init_2D

   SUBROUTINE trc_che_init_3D( Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_che_init_3D  ***
      !!
      !! ** Purpose :   Calculate values of dissociation constants at all depths
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kmm  ! time level indices
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   ztkel, zt, zt2, zsal, zbuf1, zbuf2
      REAL(wp) ::   zpres, ztc, zcl, zcpexp, zcpexp2, zfact
      REAL(wp) ::   zsqrt, ztr, zlogt
      REAL(wp) ::   zis, zis2, zsal15, zisqrt
      REAL(wp) ::   zak1, zak2, zakb, zakw, zakp1, zakp2, zakp3, zaksi
      REAL(wp) ::   zck1, zck2, zckw, zckb, zckp1, zckp2, zckp3, zcksi
      REAL(wp) ::   zaksp0, zaksp1
      REAL(wp) ::   ztmas, ztmas1
      !!---------------------------------------------------------------------

      IF( ln_timing )  CALL timing_start('trc_che_init_3D')
      !
      ! ----------------------------------

       DO jk = 1, jpk
          DO jj = 1, jpj
             DO ji = 1, jpi
 
                ztmas   = tmask_bgc_closea(ji,jj,jk)
                ztmas1  = 1. - tmask_bgc_closea(ji,jj,jk)
                ! PRESSURE in dbar
                zpres   = 1.025e-1 *gdept(ji,jj,jk,Kmm)

                ! ABSOLUTE TEMPERATURE
                ztkel   = ts(ji,jj,jk,jp_tem,Kmm) + 273.15
                zsal  = ts(ji,jj,1,jp_sal,Kmm)*ztmas + ztmas1*35.
                zsqrt  = SQRT( zsal )
                zsal15  = zsqrt * zsal
                zlogt  = LOG( ztkel )
                ztr    = 1. / ztkel
                zis    = 19.924 * zsal / ( 1000.- 1.005 * zsal )
                zis2   = zis * zis
                zisqrt = SQRT( zis )
                ztc     = ts(ji,jj,jk,jp_tem,Kmm)
                ! CHLORINITY (WOOSTER ET AL., 1969)
                zcl     = zsal * salchl
 
                ! DISSOCIATION CONSTANT FOR CARBONATE AND BORATE
                zckb    = ( cb0 + cb1 * zsqrt + cb2  * zsal + cb3 * zsal15 + cb4 * zsal * zsal ) * ztr   &
                   &    + ( cb5 + cb6 * zsqrt + cb7  * zsal )                                            &
                   &    + ( cb8 + cb9 * zsqrt + cb10 * zsal ) * zlogt + cb11 * zsqrt * ztkel 
 
                zck1    = c10 * ztr + c11 + c12 * zlogt + c13 * zsal + c14 * zsal * zsal
                zck2    = c20 * ztr + c21 + c22 * zlogt + c23 * zsal + c24 * zsal * zsal
 
                ! PKW (H2O) (DICKSON AND RILEY, 1979)
                zckw    = cw0 * ztr + cw1 + cw2 * zlogt + ( cw3 * ztr + cw4 + cw5 * zlogt ) * zsqrt + cw6 * zsal
 
                ! DISSOCIATION CONSTANTS FOR PHOSPHATE AND SILICATE
                zckp1    = cp10 * ztr + cp11 + cp12 * zlogt + (cp13 * ztr + cp14) * zsqrt + (cp15 * ztr + cp16) * zsal
                zckp2    = cp20 * ztr + cp21 + cp22 * zlogt + (cp23 * ztr + cp24) * zsqrt + (cp25 * ztr + cp26) * zsal
                zckp3    = cp30 * ztr + cp31 + (cp32 * ztr + cp33) * zsqrt + (cp34 * ztr + cp35) * zsal
 
                zcksi = csi0 * ztr + csi1 + csi2 * zlogt + (csi3 * ztr + csi4) * zisqrt + (csi5 * ztr + csi6) * zis &
                   & + (csi7 * ztr + csi8) * zis2 + LOG(1. + csi9 * zsal)

                ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN SEAWATER
                !       (S=27-43, T=2-25 DEG C) at pres =0 (atmos. pressure) (MUCCI 1983)
                zaksp0  = akcc1 + akcc2 * ztkel + akcc3 * ztr + akcc4 * LOG10( ztkel )   &
                   & + ( akcc5 + akcc6 * ztkel + akcc7 * ztr ) * zsqrt + akcc8 * zsal + akcc9 * zsal15
 
                zak1    = 10.**(zck1)
                zak2    = 10.**(zck2)
                zakb    = EXP( zckb  )
                zakw    = EXP( zckw )
                zakp1    = EXP( zckp1 )
                zakp2    = EXP( zckp2 )
                zakp3    = EXP( zckp3 )
                zaksi    = EXP( zcksi )
                zaksp1  = 10**(zaksp0)

                ! FORMULA FOR CPEXP AFTER EDMOND & GIESKES (1970)
                zcpexp  = zpres / (rgas*ztkel)
                zcpexp2 = zpres * zcpexp
                
                ! 3D arrays of dissociation constants (pressure corrected)
                zbuf1  = -     ( devk10 + devk20 * ztc + devk30 * ztc * ztc )
                zbuf2  = 0.5 * ( devk40 + devk50 * ztc )
                qak13(ji,jj,jk) = zak1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

                zbuf1  =     - ( devk11 + devk21 * ztc + devk31 * ztc * ztc )
                zbuf2  = 0.5 * ( devk41 + devk51 * ztc )
                qak23(ji,jj,jk) = zak2 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

                zbuf1  =     - ( devk12 + devk22 * ztc + devk32 * ztc * ztc )
                zbuf2  = 0.5 * ( devk42 + devk52 * ztc )
                qakb3(ji,jj,jk) = zakb * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

                zbuf1  =     - ( devk13 + devk23 * ztc + devk33 * ztc * ztc )
                zbuf2  = 0.5 * ( devk43 + devk53 * ztc )
                qakw3(ji,jj,jk) = zakw * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

                zbuf1  =     - ( devk17 + devk27 * ztc + devk37 * ztc * ztc )
                zbuf2  = 0.5 * ( devk47 + devk57 * ztc )
                qakp13(ji,jj,jk) = zakp1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

                zbuf1  =     - ( devk18 + devk28 * ztc + devk38 * ztc * ztc )
                zbuf2  = 0.5 * ( devk48 + devk58 * ztc )
                qakp23(ji,jj,jk) = zakp2 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

                zbuf1  =     - ( devk19 + devk29 * ztc + devk39 * ztc * ztc )
                zbuf2  = 0.5 * ( devk49 + devk59 * ztc )
                qakp33(ji,jj,jk) = zakp3 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

                zbuf1  =     - ( devk110 + devk210 * ztc + devk310 * ztc * ztc )
                zbuf2  = 0.5 * ( devk410 + devk510 * ztc )
                qaksi3(ji,jj,jk) = zaksi * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

                zbuf1  =     - ( devk16 + devk26 * ztc + devk36 * ztc * ztc )
                zbuf2  = 0.5 * ( devk46 + devk56 * ztc )
                qaksp(ji,jj,jk) = zaksp1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )
 
                ! TOTAL BORATE CONCENTR. (mol kg^-1)
                qborat3(ji,jj,jk) = bor1 * zcl * bor2
 
             END DO
          END DO
       END DO

      ! -------------------------------
      !
      IF( ln_timing )  CALL timing_stop('trc_che_init_3D')
      !
   END SUBROUTINE trc_che_init_3D

   PURE FUNCTION hplus(zdic,ztalk,zph,zbot,zpo4,zsi,zak1,zak2,zakb,zakw,zakp1,zakp2,zakp3,zaksi)

      ! calculates H+ conc at each iteration of the solver (localizing it this way prevents potential divergence between 2D and 3D versions)
      ! there are two outputs: [H+] and CA; these are both needed by the calling SR
      REAL(wp), INTENT(IN) :: zdic,ztalk,zph,zbot,zpo4,zsi,zak1,zak2,zakb,zakw,zakp1,zakp2,zakp3,zaksi
      REAL(wp) :: zph2, zph3, zpd, zp0, zp1, zp3, zsia, zcalk, zah2
      REAL(wp) :: hplus(2)
      zph2 = zph*zph
      zph3 = zph*zph2
      ! CALCULATE P AND Si ION CONCENTRATIONS AS PER ORR ET AL (BPG EQUATIONS 43-47)
      ! zp3 = H3PO4, zp1 = HPO4(2-), zp0 = PO4(3-): denominator is the same for all 3 equations
      zpd = 1./ ( zph3 + zakp1*zph2 + zakp1*zakp2*zph + zakp1*zakp2*zakp3 )
      zp3 = zph3*zpo4 * zpd
      zp1 = zph*zpo4*zakp1*zakp2 * zpd
      zp0 = zpo4*zakp1*zakp2*zakp3 * zpd
      zsia = zsi / (1. + zph / zaksi)

      ! CALCULATE [ALK]([CO3--], [HCO3-])
      zcalk = ztalk - ( zakw / zph - zph + zbot / ( 1.+ zph / zakb ) + 2.*zp0 + zp1 - zp3 + zsia )

      ! CALCULATE [H+] AND [H2CO3]
      zah2 = SQRT( (zdic-zcalk)*(zdic-zcalk) + 4.* ( zcalk * zak2 &
         &                                        / zak1 ) * ( 2.* zdic - zcalk ) )
      hplus(1) = 0.5 * zak1 / zcalk * ( ( zdic - zcalk ) + zah2 )
      hplus(2) = zcalk

   END FUNCTION hplus

   INTEGER FUNCTION trc_che_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_che_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER ::   ierr(6)        ! Local variables
      !!----------------------------------------------------------------------

      ierr(:)=0

      ALLOCATE( K0CO2 (jpi,jpj),    K0O2 (jpi,jpj),    & 
              & qh2co3(jpi,jpj),    qco3(jpi,jpj,jpk), & 
              & qaksp(jpi,jpj,jpk), qomegac(jpi,jpj,jpk),    STAT=ierr(1) )
      !
      trc_che_alloc = MAXVAL( ierr )

      IF( trc_che_alloc /= 0 )   CALL ctl_warn('trc_che_alloc : failed to allocate arrays.')
      !
   END FUNCTION trc_che_alloc

   !!======================================================================

END MODULE  trcche_canbgc

