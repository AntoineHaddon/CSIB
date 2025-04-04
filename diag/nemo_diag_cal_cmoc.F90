MODULE nemo_diag_cal_cmoc

   !!===============================================================
   !!              ***   MODULE nemo_diag_cal_cmoc   ***
   !!       nemo diagnostics calculations for biogeochemistry
   !!===============================================================
   !! 2018-04 (D. Yang): Original code
   !! 2019-01 (J. Christian): biogeochemistry (CMOC) version
   !!---------------------------------------------------------------
   !!
   !!---------------------------------------------------------------
   !! density            : seawater potential density
   !! cmip6_co3sat       : [CO3--] at calcite/aragonite saturation
   !! cmip6_cchem        : [CO3--], pH, Omega_X
   !! cmip6_o2sol        : oxygen saturation concentration
   !! saturation_depth   : depth of saturation horizon
   !!---------------------------------------------------------------
   USE nemo_diag_glovars_cmoc

   IMPLICIT NONE
   PRIVATE      
   PUBLIC cmip6_co3sat           ! Called by nemo_diag_cmoc.F90
   PUBLIC cmip6_cchem
   PUBLIC cmip6_o2sol
   PUBLIC saturation_depth
   PUBLIC cmip6_zo2min
   PUBLIC density

   REAL, PARAMETER :: salchl = 1. / 1.80655    ! conversion factor for salinity --> chlorinity (Wooster et al. 1969)
   REAL, PARAMETER :: Ca=0.010280              ! concentration of Ca at S=35 in mol kg^-1 (from Zeebe Table 1.1.6)

   REAL, PARAMETER :: o2atm  = 1. / ( 1000. * 0.20946 )
   REAL, PARAMETER :: rgas   = 83.143         ! universal gas constant

   REAL, PARAMETER :: akcc1  = -171.9065     ! coefficients for solubility product (calcite)
   REAL, PARAMETER :: akcc2  =   -0.077993 
   REAL, PARAMETER :: akcc3  = 2839.319
   REAL, PARAMETER :: akcc4  =   71.595
   REAL, PARAMETER :: akcc5  =   -0.77712
   REAL, PARAMETER :: akcc6  =    0.00284263
   REAL, PARAMETER :: akcc7  =  178.34
   REAL, PARAMETER :: akcc8  =   -0.07711
   REAL, PARAMETER :: akcc9  =    0.0041249

   REAL, PARAMETER :: bkcc1  = -171.945       ! coefficients for solubility product (aragonite)
   REAL, PARAMETER :: bkcc2  =-0.077993
   REAL, PARAMETER :: bkcc3  = 2903.293
   REAL, PARAMETER :: bkcc4  =   71.595
   REAL, PARAMETER :: bkcc5  =   -0.068393
   REAL, PARAMETER :: bkcc6  =    0.0017276
   REAL, PARAMETER :: bkcc7  =  88.135
   REAL, PARAMETER :: bkcc8  =   -0.10018
   REAL, PARAMETER :: bkcc9  =    0.0059415

   REAL, PARAMETER :: bor1   = 0.000232       ! borat constants
   REAL, PARAMETER :: bor2   = 1. / 10.811

   REAL, PARAMETER :: c10    = -3633.86        ! Coeff. for 1. dissoc. of carbonic acid (Dickson et al., 2007)
   REAL, PARAMETER :: c11    =    61.2172
   REAL, PARAMETER :: c12    =    -9.67770
   REAL, PARAMETER :: c13    =     0.011555
   REAL, PARAMETER :: c14    =    -0.0001152

   REAL, PARAMETER :: c20    =  -471.78       ! coeff. for 2. dissoc. of carbonic acid (Dickson et al., 2007)   
   REAL, PARAMETER :: c21    =   -25.9290
   REAL, PARAMETER :: c22    =     3.16967
   REAL, PARAMETER :: c23    =     0.01781
   REAL, PARAMETER :: c24    =    -0.0001122

   REAL, PARAMETER :: cb0    = -8966.90      ! Coeff. for 1. dissoc. of boric acid 
   REAL, PARAMETER :: cb1    = -2890.53      ! (Dickson and Goyet, 1994)
   REAL, PARAMETER :: cb2    =   -77.942
   REAL, PARAMETER :: cb3    =     1.728
   REAL, PARAMETER :: cb4    =    -0.0996
   REAL, PARAMETER :: cb5    =   148.0248
   REAL, PARAMETER :: cb6    =   137.1942
   REAL, PARAMETER :: cb7    =     1.62142
   REAL, PARAMETER :: cb8    =   -24.4344
   REAL, PARAMETER :: cb9    =   -25.085
   REAL, PARAMETER :: cb10   =    -0.2474
   REAL, PARAMETER :: cb11   =     0.053105

   REAL, PARAMETER :: cw0    = -13847.26     ! Coeff. for dissoc. of water (Dickson and Riley, 1979 )
   REAL, PARAMETER :: cw1    =    148.9652
   REAL, PARAMETER :: cw2    =    -23.6521
   REAL, PARAMETER :: cw3    =    118.67
   REAL, PARAMETER :: cw4    =     -5.977
   REAL, PARAMETER :: cw5    =      1.0495
   REAL, PARAMETER :: cw6    =     -0.01615

   REAL, PARAMETER :: cp10    = -4576.752    ! Coeff. for dissoc. of H3PO4 (Dickson et al 2007)
   REAL, PARAMETER :: cp11    = 115.525
   REAL, PARAMETER :: cp12    = -18.453
   REAL, PARAMETER :: cp13    = -106.736
   REAL, PARAMETER :: cp14    = 0.69171
   REAL, PARAMETER :: cp15    = -0.65643
   REAL, PARAMETER :: cp16    = -0.01844

   REAL, PARAMETER :: cp20    = -8814.715    ! Coeff. for dissoc. of H2PO4- (Dickson et al 2007)
   REAL, PARAMETER :: cp21    = 172.0883
   REAL, PARAMETER :: cp22    = -27.927
   REAL, PARAMETER :: cp23    = -160.340
   REAL, PARAMETER :: cp24    = 1.3566
   REAL, PARAMETER :: cp25    = 0.37335
   REAL, PARAMETER :: cp26    = -0.05778

   REAL, PARAMETER :: cp30    = -3070.75     ! Coeff. for dissoc. of HPO4-- (Dickson et al 2007)
   REAL, PARAMETER :: cp31    = -18.141
   REAL, PARAMETER :: cp32    = 17.27039
   REAL, PARAMETER :: cp33    = 2.81197
   REAL, PARAMETER :: cp34    = -44.99486
   REAL, PARAMETER :: cp35    = -0.09984

   REAL, PARAMETER :: csi0    = -8904.2      ! Coeff. for dissoc. of Si(OH)4 (Dickson et al 2007)
   REAL, PARAMETER :: csi1    = 117.385
   REAL, PARAMETER :: csi2    = -19.334
   REAL, PARAMETER :: csi3    = -458.79
   REAL, PARAMETER :: csi4    = 3.5913
   REAL, PARAMETER :: csi5    = 188.74
   REAL, PARAMETER :: csi6    = -1.5998
   REAL, PARAMETER :: csi7    = -12.1652
   REAL, PARAMETER :: csi8    = 0.07871
   REAL, PARAMETER :: csi9    = -0.001005

   !                                    ! volumetric solubility constants for o2 in ml/L  
   REAL, PARAMETER :: ox0    =  2.00856      ! from Table 1 for Eq 8 of Garcia and Gordon, 1992.
   REAL, PARAMETER :: ox1    =  3.22400      ! corrects for moisture and fugacity, but not total atmospheric pressure
   REAL, PARAMETER :: ox2    =  3.99063      !      Original PISCES code noted this was a solubility, but 
   REAL, PARAMETER :: ox3    =  4.80299      ! was in fact a bunsen coefficient with units L-O2/(Lsw atm-O2)
   REAL, PARAMETER :: ox4    =  9.78188e-1   ! Hence, need to divide EXP( zoxy ) by 1000, ml-O2 => L-O2
   REAL, PARAMETER :: ox5    =  1.71069      ! and atcox = 0.20946 to add the 1/atm dimension.
   REAL, PARAMETER :: ox6    = -6.24097e-3
   REAL, PARAMETER :: ox7    = -6.93498e-3
   REAL, PARAMETER :: ox8    = -6.90358e-3
   REAL, PARAMETER :: ox9    = -4.29155e-3
   REAL, PARAMETER :: ox10   = -3.11680e-7

   REAL, DIMENSION(6)  :: devk1, devk2, devk3, devk4, devk5   ! coeff. for seawater pressure correction 
   !                                                              ! (millero 95)
   DATA devk1 / -25.5    , -15.82    , -29.48  , -25.60     , -48.76   , -46.00    /
   DATA devk2 / 0.1271   , -0.0219   , 0.1622  , 0.2324     , 0.5304   , 0.5304    /
   DATA devk3 / 0.       , 0.        ,-2.608E-3,  -3.6246E-3, 0.       , 0.        /
   DATA devk4 / -3.08E-3 , 1.13E-3   , -2.84E-3, -5.13E-3   , -11.76E-3, -11.76E-3 /
   DATA devk5 / 0.0877E-3, -0.1475E-3,  0.     , 0.0794E-3  , 0.3692E-3, 0.3692E-3 /

CONTAINS

   SUBROUTINE cmip6_co3sat(Kspa, Kspc) 
      !!
      !!-------------------------------------------------------------
      !! Purpose: Compute [CO3--] at calcite/aragonite saturation  = 1
      !! Reference: Orr et al, 2017 
      !! Input fields:
      !!            3D: TT, SS, tmask
      !! Output:    carbonate ion concentration in mol m^-3
      !!-------------------------------------------------------------
      INTEGER                 :: i,j,k,l
      INTEGER, DIMENSION (1)  :: ierr                   ! local variable
      REAL :: ztkel, zsal, zsqrt, zsal15, zlogt, ztr, zis, zis2, zisqrt, ztc
      REAL :: zaksp0, zaksp1, zbuf1, zbuf2, zpres, zcpexp, zcpexp2,zca
      REAL, DIMENSION(imt,jmt,km,lm) :: Kspa, Kspc

       DO l=1,lm
        DO k=1,km
         DO j=1,jmt
          DO i=1,imt

               ! SET PRESSION
               zpres   = 1.025e-1 * deptht(k)

      ! SET ABSOLUTE TEMPERATURE
                ztkel   = TT(i,j,k,l) + 273.15
                zsal    = SS(i,j,k,l) + ( 1.-tmask(i,j,k) ) * 35.
                zsqrt  = SQRT( zsal )
                zsal15  = zsqrt * zsal
                zlogt  = LOG( ztkel )
                ztr    = 1. / ztkel
                zis    = 19.924 * zsal / ( 1000.- 1.005 * zsal )
                zis2   = zis * zis
                zisqrt = SQRT( zis )
                ztc     = TT(i,j,k,l) + ( 1.- tmask(i,j,k) ) * 20.
                zca     = Ca*zsal/35.                   ! [Ca++] in mol kg-1

      ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN SEAWATER
      !       (S=27-43, T=2-25 DEG C) at pres =0 (atmos. pressure) (MUCCI 1983)
                zaksp0  = akcc1 + akcc2 * ztkel + akcc3 * ztr + akcc4 * LOG10( ztkel )   &
              &   + ( akcc5 + akcc6 * ztkel + akcc7 * ztr ) * zsqrt + akcc8 * zsal + akcc9 * zsal15
                zaksp1  = 10.**(zaksp0)
               ! FORMULA FOR CPEXP AFTER EDMOND & GIESKES (1970)
               !        (REFERENCE TO CULBERSON & PYTKOQICZ (1968) AS MADE
               !        IN BROECKER ET AL. (1982) IS INCORRECT; HERE RGAS IS
               !        TAKEN TENFOLD TO CORRECT FOR THE NOTATION OF pres  IN
               !        DBAR INSTEAD OF BAR AND THE EXPRESSION FOR CPEXP IS
               !        MULTIPLIED BY LN(10.) TO ALLOW USE OF EXP-FUNCTION
               !        WITH BASIS E IN THE FORMULA FOR AKSPP (CF. EDMOND
               !        & GIESKES (1970), P. 1285-1286 (THE SMALL
               !        FORMULA ON P. 1286 IS RIGHT AND CONSISTENT WITH THE
               !        SIGN IN PARTIAL MOLAR VOLUME CHANGE AS SHOWN ON P. 1285))

               zcpexp  = zpres /(rgas*ztkel)
               zcpexp2 = zpres * zpres/(rgas*ztkel)

      ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE 
      !        AS FUNCTION OF PRESSURE FOLLOWING MILLERO
      !        (P. 1285) AND BERNER (1976)
                zbuf1  =     - ( devk1(5) + devk2(5) * ztc + devk3(5) * ztc * ztc )
                zbuf2  = 0.5 * ( devk4(5) + devk5(5) * ztc )
! [CO3--] in mol m^-3 
                Kspc(i,j,k,l) = zaksp1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )
                Kspc(i,j,k,l) = Kspc(i,j,k,l) + ( 1.- tmask(i,j,k) ) * 1.e-7
                co3_satc(i,j,k,l) = Kspc(i,j,k,l) / zca * prhop(i,j,k,l)
                co3_satc(i,j,k,l) = co3_satc(i,j,k,l) * tmask(i,j,k)

      ! APPARENT SOLUBILITY PRODUCT K'SP OF ARAGONITE IN SEAWATER
                zaksp0  = bkcc1 + bkcc2 * ztkel + bkcc3 * ztr + bkcc4 * LOG10( ztkel )   &
              &   + ( bkcc5 + bkcc6 * ztkel + bkcc7 * ztr ) * zsqrt + bkcc8 * zsal + bkcc9 * zsal15
                zaksp1  = 10.**(zaksp0)
                zbuf1  =     - ( devk1(6) + devk2(6) * ztc + devk3(6) * ztc * ztc )
                zbuf2  = 0.5 * ( devk4(6) + devk5(6) * ztc )
                Kspa(i,j,k,l) = Kspa(i,j,k,l) + ( 1.- tmask(i,j,k) ) * 1.e-7
                Kspa(i,j,k,l) = zaksp1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )
                co3_sata(i,j,k,l) = Kspa(i,j,k,l) / zca * prhop(i,j,k,l)
                co3_sata(i,j,k,l) = co3_sata(i,j,k,l) * tmask(i,j,k)

          ENDDO
         ENDDO
        ENDDO
       ENDDO

   END SUBROUTINE cmip6_co3sat

   SUBROUTINE cmip6_cchem(XDIC, XTA, Kspc, Kspa, CO3, pH, Om_C, Om_A)

      !!-------------------------------------------------------------
      !! Purpose: Solve carbon chemistry to generate [H+], which can then be used to calculate [CO3--] etc
      !! Reference: Orr et al, 2017 
      !! Input fields:
      !!            3D: TT, SS, DIC, TA, NO3, NH4, SI, tmask
      !! Output:    hydrogen ion concentration in mol kg^-1
      !!-------------------------------------------------------------
      INTEGER                 :: i,j,k,l,jm
      INTEGER, DIMENSION (1)  :: ierr                   ! local variable
      REAL :: ztkel, zsal, zsqrt, zsal15, zlogt, ztr, zis, zis2, zisqrt, ztc
      REAL :: zpres, zcl, zckb, zck1, zck2, zckw, zckp1, zckp2, zckp3, zcksi
      REAL :: zak1, zak2, zakb, zakw, zakp1, zakp2, zakp3, zaksi, zaksp1, zaksp2
      REAL :: zaksp0, zbuf1, zbuf2, zcpexp, zcpexp2, zbot, zfact, zdic, zph
      REAL :: zalka, zph2, zph3, zpo4, zsi, zpd, zp3, zp1, zp0, zalk, zah2, hion
      REAL :: zrhop, zr1, zr2, zr3, zr4, zt, zs, zsr, zca
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: hi, borat, ak13, ak23, akb3, akw3
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: akp13, akp23, akp33, aksi3
      REAL, DIMENSION(imt,jmt,km,lm) :: XDIC, XTA, CO3, pH, Om_A, Om_C, Kspc, Kspa

      !!----------------
      !! Allocate Arrays
      !!----------------
      ALLOCATE( hi(imt,jmt,km,lm),borat(imt,jmt,km,lm),ak13(imt,jmt,km,lm),ak23(imt,jmt,km,lm),     &
                akb3(imt,jmt,km,lm),akw3(imt,jmt,km,lm),akp13(imt,jmt,km,lm),akp23(imt,jmt,km,lm),  &
                akp33(imt,jmt,km,lm),aksi3(imt,jmt,km,lm), STAT=ierr(1) )

      IF (MAXVAL(ierr) /=0) THEN
         STOP 'Memory allocation error in cmip6_cchem'
      ENDIF

! this part is from p4zche.F90 (define the equilibrium constants, borate concentration etc as function of T, S, P)
       DO l=1,lm
        DO k=1,km
         DO j=1,jmt
          DO i=1,imt

               ! SET PRESSION
               zpres   = 1.025e-1 * deptht(k)

               ! SET ABSOLUTE TEMPERATURE
               ! Apply a data funnel to ensure that temperatures and salinities
               ! fall within a reasonable range
               ztkel   = MIN(MAX(TT(i,j,k,l),-2.),40.) + 273.15
               zsal    = MIN(MAX(SS(i,j,k,l),5.),44.) + ( 1.-tmask(i,j,k) ) * 35.
               zsqrt  = SQRT( zsal )
               zsal15  = zsqrt * zsal
               zlogt  = LOG( ztkel )
               ztr    = 1. / ztkel
               zis    = 19.924 * zsal / ( 1000.- 1.005 * zsal )
               zis2   = zis * zis
               zisqrt = SQRT( zis )
               ztc     = TT(i,j,k,l) + ( 1.- tmask(i,j,k) ) * 20.

               ! CHLORINITY (WOOSTER ET AL., 1969)
               zcl     = zsal * salchl

               ! TOTAL BORATE CONCENTR. [MOLES/L]
               borat(i,j,k,l) = bor1 * zcl * bor2

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

               ! K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW (H2O) (LIT.?)
               zak1    = 10.**(zck1)
               zak2    = 10.**(zck2)
               zakb    = EXP( zckb  )
               zakw    = EXP( zckw )
               zakp1    = EXP( zckp1 )
               zakp2    = EXP( zckp2 )
               zakp3    = EXP( zckp3 )
               zaksi    = EXP( zcksi )
               zaksp1  = 10.**(zaksp0)

               ! FORMULA FOR CPEXP AFTER EDMOND & GIESKES (1970)
               !        (REFERENCE TO CULBERSON & PYTKOQICZ (1968) AS MADE
               !        IN BROECKER ET AL. (1982) IS INCORRECT; HERE RGAS IS
               !        TAKEN TENFOLD TO CORRECT FOR THE NOTATION OF pres  IN
               !        DBAR INSTEAD OF BAR AND THE EXPRESSION FOR CPEXP IS
               !        MULTIPLIED BY LN(10.) TO ALLOW USE OF EXP-FUNCTION
               !        WITH BASIS E IN THE FORMULA FOR AKSPP (CF. EDMOND
               !        & GIESKES (1970), P. 1285-1286 (THE SMALL
               !        FORMULA ON P. 1286 IS RIGHT AND CONSISTENT WITH THE
               !        SIGN IN PARTIAL MOLAR VOLUME CHANGE AS SHOWN ON P. 1285))

               zcpexp  = zpres /(rgas*ztkel)
               zcpexp2 = zpres * zpres/(rgas*ztkel)

              ! KB OF BORIC ACID, K1,K2 OF CARBONIC ACID PRESSURE
               !        CORRECTION AFTER CULBERSON AND PYTKOWICZ (1968)
               !        (CF. BROECKER ET AL., 1982)

               zbuf1  =     - ( devk1(1) + devk2(1) * ztc + devk3(1) * ztc * ztc )
               zbuf2  = 0.5 * ( devk4(1) + devk5(1) * ztc )
               ak13(i,j,k,l) = zak1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

               zbuf1  =     - ( devk1(2) + devk2(2) * ztc + devk3(2) * ztc * ztc )
               zbuf2  = 0.5 * ( devk4(2) + devk5(2) * ztc )
               ak23(i,j,k,l) = zak2 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

               zbuf1  =     - ( devk1(3) + devk2(3) * ztc + devk3(3) * ztc * ztc )
               zbuf2  = 0.5 * ( devk4(3) + devk5(3) * ztc )
               akb3(i,j,k,l) = zakb * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

               zbuf1  =     - ( devk1(4) + devk2(4) * ztc + devk3(4) * ztc * ztc )
               zbuf2  = 0.5 * ( devk4(4) + devk5(4) * ztc )
               akw3(i,j,k,l) = zakw * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

               ! K_Px and K_Si (NO PRESSURE CORRECTION)
               akp13(i,j,k,l) = zakp1
               akp23(i,j,k,l) = zakp2
               akp33(i,j,k,l) = zakp3
               aksi3(i,j,k,l) = zaksi

          ENDDO
         ENDDO
        ENDDO
       ENDDO

! initialize [H+] as a constant 1.e-8
       DO l=1,lm
        DO k=1,km
         DO j=1,jmt
          DO i=1,imt
              hi(i,j,k,l)=1.e-8
          ENDDO
         ENDDO
        ENDDO
       ENDDO

! this part is from p4zflx.F90

       DO jm = 1, 10
        DO l=1,lm
         DO k=1,km
          DO j=1,jmt
           DO i=1,imt

               ! DUMMY VARIABLES FOR DIC, H+, AND BORATE
               zbot  = borat(i,j,k,l)
               zfact = prhop(i,j,k,l)*0.001 + (1.-tmask(i,j,k))
               zdic  = XDIC(i,j,k,l) * 0.000001 / zfact + (1.-tmask(i,j,k))*0.002
               zalka = XTA(i,j,k,l) * 0.000001 / zfact + (1.-tmask(i,j,k))*0.002
               zph   = MAX( hi(i,j,k,l), 1.e-10 ) / zfact
               zph2 = zph*zph
               zph3 = zph*zph2
               zpo4 = NO3(i,j,k,l) / 16. * 0.000001 / zfact
               zsi = asi3(i,j,k,l) * 0.000001 / zfact                        ! silica is a static array based on initialization file, not a carried tracer

               ! CALCULATE P AND Si ION CONCENTRATIONS AS PER ORR ET AL (BPG EQUATIONS 43-47)
               ! zp3 = H3PO4, zp1 = HPO4(2-), zp0 = PO4(3-): denominator is the same for all 3 equations
               zpd = 1./ ( zph3 + akp13(i,j,k,l)*zph2 + akp13(i,j,k,l)*akp23(i,j,k,l)*zph + &
                            akp13(i,j,k,l)*akp23(i,j,k,l)*akp33(i,j,k,l) )
               zp3 = zph3*zpo4 * zpd
               zp1 = zph*zpo4*akp13(i,j,k,l)*akp23(i,j,k,l) * zpd
               zp0 = zpo4*akp13(i,j,k,l)*akp23(i,j,k,l)*akp33(i,j,k,l) * zpd
               zsi = zsi / (1. + zph / aksi3(i,j,k,l))

               ! CALCULATE [ALK]([CO3--], [HCO3-])
               zalk  = zalka - (  akw3(i,j,k,l) / zph - zph + zbot / ( 1.+ zph / akb3(i,j,k,l) ) + 2.*zp0 + zp1 - zp3 + zsi )

               ! CALCULATE [H+] AND [H2CO3]
               zah2   = SQRT(  (zdic-zalk)*(zdic-zalk) + 4.* ( zalk * ak23(i,j,k,l)   &
                  &                                        / ak13(i,j,k,l) ) * ( 2.* zdic - zalk )  )
               zah2   = 0.5 * ak13(i,j,k,l) / zalk * ( ( zdic - zalk ) + zah2 )
               hi(i,j,k,l)   = zah2 * zfact + (1.-tmask(i,j,k))*1.e-7

           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       pH=0.
       CO3=0.
       Om_A=0.
       Om_C=0.

       DO l=1,lm
        DO k=1,km
         DO j=1,jmt
          DO i=1,imt
              if (tmask(i,j,k).eq.0. .or. hi(i,j,k,l).le.0.) cycle
              ! pH calculated from [H+] in mol L^-1
              pH(i,j,k,l)=ALOG10(hi(i,j,k,l))*(-1.)*tmask(i,j,k)
              ! convert [H+] to  mol kg^-1 (XDIC is in mmol m^-3; CO3 is in mol m^-3; hion and ak* are in mol kg^-1)
              zfact = prhop(i,j,k,l)*0.001 + (1.-tmask(i,j,k))
              hion=hi(i,j,k,l)/zfact
              CO3(i,j,k,l)=XDIC(i,j,k,l)*ak13(i,j,k,l)*ak23(i,j,k,l)/ &
                            (hion*hion + ak13(i,j,k,l)*hion + ak13(i,j,k,l)*ak23(i,j,k,l))*0.001
              CO3(i,j,k,l)=CO3(i,j,k,l)*tmask(i,j,k)
              zca     = Ca*SS(i,j,k,l)/35.                   ! [Ca++] in mol kg-1
              Om_A(i,j,k,l) = zca * (CO3(i,j,k,l) * 0.001 / zfact) / Kspa(i,j,k,l)      ! to calculate Omega, [CO3--] must be in mol kg^-1
              Om_C(i,j,k,l) = zca * (CO3(i,j,k,l) * 0.001 / zfact) / Kspc(i,j,k,l)
          ENDDO
         ENDDO
        ENDDO
       ENDDO

      DEALLOCATE( hi, borat, ak13, ak23, akb3, akw3, akp13, akp23, akp33, aksi3 )

   END SUBROUTINE cmip6_cchem

   SUBROUTINE cmip6_o2sol
      !!
      !!-------------------------------------------------------------
      !! Purpose: Compute solubility of oxygen in sea water
      !! Reference: Orr et al, 2017 
      !! Input fields:
      !!            3D: TT, SS, tmask
      !! Output:    solubility in mol m^-3
      !!-------------------------------------------------------------
      INTEGER                 :: i,j,k,l
      INTEGER, DIMENSION (1)  :: ierr                   ! local variable
      REAL :: ztkel, zsal, zt, zt2, zsal2, zlogt, zcek1, ztgg, ztgg2, ztgg3, ztgg4, ztgg5, zoxy 

       DO l=1,lm
        DO k=1,km
         DO j=1,jmt
          DO i=1,imt

      ! SET ABSOLUTE TEMPERATURE
                ztkel   = TT(i,j,k,l) + 273.15
                zt    = ztkel * 0.01
                zt2   = zt * zt
                zsal    = SS(i,j,k,l) + ( 1.-tmask(i,j,k) ) * 35.
                zsal2 = zsal * zsal
                zlogt = LOG( zt )
            !                             ! LN(K0) OF SOLUBILITY OF O2 and N2 in ml/L (EQ. 8, GARCIA AND GORDON, 1992)
                ztgg  = LOG( ( 298.15 - TT(i,j,k,l) ) / ztkel )  ! Set the GORDON & GARCIA scaled temperature
                ztgg2 = ztgg  * ztgg
                ztgg3 = ztgg2 * ztgg
                ztgg4 = ztgg3 * ztgg
                ztgg5 = ztgg4 * ztgg
                zoxy  = ox0 + ox1 * ztgg + ox2 * ztgg2 + ox3 * ztgg3 + ox4 * ztgg4 + ox5 * ztgg5   &
                   + zsal * ( ox6 + ox7 * ztgg + ox8 * ztgg2 + ox9 * ztgg3 ) +  ox10 * zsal2

                o2sol(i,j,k,l) =  EXP( zoxy  ) * 0.0445919       ! convert to mol m^-3
                o2sol(i,j,k,l) = o2sol(i,j,k,l) * tmask(i,j,k)

          ENDDO
         ENDDO
        ENDDO
       ENDDO

   END SUBROUTINE cmip6_o2sol

   SUBROUTINE saturation_depth(Om_C, Om_A)

      !!-------------------------------------------------------------
      !! Purpose: Calculate depth of calcite and aragonite saturation horizon based on criteria given in footnotes to Table 15 in Orr et al, 2017 
      !! If supersaturated water is never encountered then z_sat=0; if all levels are supersaturated then z_sat=MDV (99999.); otherwise, z_sat is depth of shallowest layer where Omega<1
      !! A value of 0 will be registered in the unlikely case that there is undersaturation at the surface and a subsurface supersaturation zone
      !! Input fields:
      !!            3D: Omega_A, Omega_C
      !! Output:    depth of saturation horizon
      !!-------------------------------------------------------------
      INTEGER                 :: i,j,k,l,isw1,isw2,kk
      REAL, DIMENSION(imt,jmt,km,lm) :: Om_C, Om_A

       DO l=1,lm
         DO j=1,jmt
          DO i=1,imt

           isw1=0
           isw2=0
           kk=0
           DO k=1,km
            IF (Om_C(i,j,k,l) .GT. 1. .AND. isw1 .EQ. 0) THEN
             kk=k+1
            ELSE
             isw1=1
             IF (k .EQ. 1) isw2=1
            ENDIF
           ENDDO
           zsat_c(i,j,l)=deptht(kk+1)       ! the RHS here will be garbage in case "all levels are supersaturated" because kk+1>km, but it will be overwritten in the next line
           IF (isw1 .EQ. 0) zsat_c(i,j,l)=99999.
           IF (isw1 .EQ. 1 .AND. isw2 .EQ. 1) zsat_c(i,j,l)=0.

           isw1=0
           isw2=0
           kk=0
           DO k=1,km
            IF (Om_A(i,j,k,l) .GT. 1. .AND. isw1 .EQ. 0) THEN
             kk=k+1
            ELSE
             isw1=1
             IF (k .EQ. 1) isw2=1
            ENDIF
           ENDDO
           zsat_a(i,j,l)=deptht(kk+1)
           IF (isw1 .EQ. 0) zsat_a(i,j,l)=99999.
           IF (isw1 .EQ. 1 .AND. isw2 .EQ. 1) zsat_a(i,j,l)=0.

          ENDDO
         ENDDO
        ENDDO

   END SUBROUTINE saturation_depth

   SUBROUTINE cmip6_zo2min

      !!-------------------------------------------------------------------------------------------------------
      !! Purpose: Calculate minimum oxygen concentration over the water column and the depth at which it occurs
      !! Input fields:
      !!            3D: O2
      !! Output:    minimum oxygen concentration, depth of minimum oxygen concentration
      !!-------------------------------------------------------------------------------------------------------
      INTEGER                 :: i,j,k,l

       DO l=1,lm
         DO j=1,jmt
          DO i=1,imt

           o2min(i,j,l)=O2(i,j,1,l)*1.e-3
           zo2min(i,j,l)=0.
           DO k=1,km
            IF (O2(i,j,k,l)*1.e-3 .LT. o2min(i,j,l) .AND. tmask(i,j,k) .NE. 0.) THEN
             o2min(i,j,l)=O2(i,j,k,l)*1.e-3
             zo2min(i,j,l)=deptht(k)
            ENDIF
           ENDDO

          ENDDO
         ENDDO
        ENDDO

   END SUBROUTINE cmip6_zo2min

   SUBROUTINE density
      !!-------------------------------------------------------------
      !! Purpose: Calculate seawater potential density (from eosbn2.F90)
      !! Input fields:
      !!            3D: T, S
      !! Output:    potential density in kg m^-3
      !!-------------------------------------------------------------
      INTEGER                 :: i,j,k,l
      REAL                    :: zt, zs, zsr, zr1, zr2, zr3, zr4, zrhop

        DO l=1,lm
         DO k=1,km
          DO j=1,jmt
           DO i=1,imt

                zt   = TT(i,j,k,l)
                zs    = SS(i,j,k,l) + ( 1.-tmask(i,j,k) ) * 35.
                zsr  = SQRT( zs )
                  !
                  ! compute volumic mass pure water at atm pressure
                  zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4 )*zt   &
                     &                          -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
                  ! seawater volumic mass atm pressure
                  zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
                     &                                         -4.0899e-3 ) *zt+0.824493
                  zr3= ( -1.6546e-6*zt+1.0227e-4 )    *zt-5.72466e-3
                  zr4= 4.8314e-4
                  !
                  ! potential volumic mass (reference to the surface)
                  zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
                  !
                  ! save potential volumic mass
                  prhop(i,j,k,l) = zrhop * tmask(i,j,k)
                  !
          ENDDO
         ENDDO
        ENDDO
       ENDDO

   END SUBROUTINE density

END MODULE nemo_diag_cal_cmoc

