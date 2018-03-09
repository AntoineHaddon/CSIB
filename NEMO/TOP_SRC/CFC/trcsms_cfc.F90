

MODULE trcsms_cfc
   !!======================================================================
   !! *** MODULE trcsms_cfc ***
   !! TOP : CFC main model
   !!======================================================================
   !! History : OPA ! 1999-10 (JC. Dutay) original code
   !! NEMO 1.0 ! 2004-03 (C. Ethe) free form + modularity
   !! 2.0 ! 2007-12 (C. Ethe, G. Madec) reorganisation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !! 'key_cfc' CFC tracers
   !!----------------------------------------------------------------------
   !! trc_sms_cfc : compute and add CFC suface forcing to CFC trends
   !! trc_cfc_cst : sets constants for CFC surface forcing computation
   !!----------------------------------------------------------------------
   USE oce_trc ! Ocean variables
   USE par_trc ! TOP parameters
   USE sbc_oce, only : fr_i ! Fraction of surface covered by ice
   USE trc ! TOP variables
   USE trdmod_oce
   USE trdmod_trc
   USE iom ! I/O library
   USE par_cfc, only : jp_cfc
   USE dom_oce, only : nsec_year, nyear_len, nyear

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_sms_cfc ! called in ???
   PUBLIC trc_sms_cfc_alloc ! called in trcini_cfc.F90

!   #include "domzgr_substitute.h90"

   INTEGER , PUBLIC, PARAMETER :: jphem = 2 ! parameter for the 2 hemispheres
   INTEGER , PUBLIC :: jpyear ! Number of years read in CFC1112 file
   INTEGER , PUBLIC :: cfc_year_offset ! Offset from model year.

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) :: p_cfc_year ! Year associated with the atmospheric partial pressure
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: p_cfc ! partial hemispheric pressure for CFC
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: xphem ! spatial interpolation factor for patm
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: qtr_cfc ! flux at surface
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: qint_cfc ! cumulative flux
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: patm ! atmospheric function

   REAL(wp), DIMENSION(4,jp_cfc) :: soa ! coefficient for solubility of CFC [mol/l/atm]
   REAL(wp), DIMENSION(3,jp_cfc) :: sob ! "               "
   REAL(wp), DIMENSION(5,jp_cfc) :: sca ! coefficients for schmidt number in degre Celcius

   ! ! coefficients for conversion
   REAL(wp), parameter :: xconv1 = 1.0 ! conversion from to
   REAL(wp), parameter :: xconv2 = 0.01/3600. ! conversion from cm/h to m/s:
   REAL(wp), parameter :: xconv3 = 1.0e+3 ! conversion from mol/l/atm to mol/m3/atm
   REAL(wp), parameter :: xconv4 = 1.0e-12 ! conversion from mol/m3/atm to mol/m3/pptv
   REAL(wp), parameter :: kw_scale = 0.251 !< Scale factor used when calculating Schmidt number
                                            !! 0.251 is the updated Wanninkhof 2014 value
                                            !! 0.39 is the Wanninkhof 1992 values

CONTAINS
   SUBROUTINE trc_sms_cfc( kt )
      !!----------------------------------------------------------------------
      !! *** ROUTINE trc_sms_cfc ***
      !!
      !! ** Purpose : Compute the surface boundary contition on CFC 11
      !! passive tracer associated with air-mer fluxes and add it
      !! to the general trend of tracers equations.
      !!
      !! ** Method : - get the atmospheric partial pressure - given in pico -
      !! - computation of solubility ( in 1.e-12 mol/l then in 1.e-9 mol/m3)
      !! - computation of transfert speed ( given in cm/hour ----> cm/s )
      !! - the input function is given by :
      !! speed * ( concentration at equilibrium - concentration at surface )
      !! - the input function is in pico-mol/m3/s and the
      !! CFC concentration in pico-mol/m3
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: kt ! ocean time-step index
      !
      INTEGER :: ji, jj, jn, jl, jm, js
      INTEGER :: iyear_beg, iyear_end
      INTEGER :: ierr
      INTEGER :: nyears, cfc_year
      REAL(wp) :: yearfrac, wt1, wt2
      REAL(wp) :: ztap, zdtap, iztap
      REAL(wp) :: zt1, zt2, zt3, zv2
      REAL(wp) :: zsol ! solubility
      REAL(wp) :: zsch ! schmidt number
      REAL(wp) :: zca_cfc ! concentration at equilibrium
      REAL(wp) :: zak_cfc ! transfert coefficients
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zpp_cfc ! atmospheric partial pressure of CFC
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zpatm ! atmospheric function
      !!----------------------------------------------------------------------
      !
      !
      IF( nn_timing == 1 ) CALL timing_start('trc_sms_cfc')
      !
      ALLOCATE( zpatm(jphem,jp_cfc), STAT=ierr )
      IF( ierr > 0 ) THEN
         CALL ctl_stop( 'trc_sms_cfc: unable to allocate zpatm array' ) ; RETURN
      ENDIF
      ALLOCATE( zpp_cfc(jpi,jpj,jp_cfc), STAT=ierr )
      IF( ierr > 0 ) THEN
         CALL ctl_stop( 'trc_sms_cfc: unable to allocate zpp_cfc array' ) ; RETURN
      ENDIF
      CALL trc_cfc_cst
      ! Temporal interpolation
      ! ----------------------
      nyears = SIZE(p_cfc_year(:))
      ! Calculate the given year that the CFC module sees
      cfc_year = nyear + cfc_year_offset
      yearfrac = ( nsec_year ) / ( 86400. * nyear_len(1) )
      ! Check to make sure that the current 'cfc_year' is within the observational range
      ! and set time interpolation factors keeping in mind that atmospheric values are
      ! annual averages with time at the middle of the year. If earlier than the record,
      ! atmospheric concentratio assumed to be equal to the first entry. If later than
      ! the record, then its assumed to be equal to the last entry
      IF ( cfc_year < p_cfc_year(1) ) THEN
        iyear_beg = 1
        wt1 = 1.
      ELSEIF ( cfc_year > p_cfc_year(nyears) ) THEN
        iyear_beg = nyears
        wt1 = 1.
      ELSEIF ( yearfrac > 0.5 ) THEN
        iyear_beg = cfc_year - FLOOR(p_cfc_year(1)) + 1
        wt1 = 1.5 - yearfrac
      ELSEIF ( yearfrac < 0.5 ) THEN
        iyear_beg = cfc_year-FLOOR(p_cfc_year(1))
        wt1 = 0.5 - yearfrac
      ELSEIF ( yearfrac == 0.5 ) THEN
        iyear_beg = cfc_year - FLOOR(p_cfc_year(1)) + 1
        wt1 = 1.
      ENDIF
      wt2 = 1. - wt1
      iyear_end = iyear_beg + 1
      ! !------------!
      DO jl = 1, jp_cfc ! CFC loop !
         ! !------------!
         ! Index of tracer in big tracer array
         jn = jp_cfc0 + jl - 1
         ! time interpolation at time kt
         DO jm = 1, jphem
            zpatm(jm,jl) = p_cfc(iyear_beg, jm, jl)*wt1 + p_cfc(iyear_end, jm, jl)*wt2
         END DO
         ! !------------!
         DO jj = 1, jpj ! i-j loop !
            DO ji = 1, jpi !------------!
               ! space interpolation
               zpp_cfc(ji,jj,jl) = xphem(ji,jj) * zpatm(1,jl) &
                  & + ( 1.- xphem(ji,jj) ) * zpatm(2,jl)
               ! Computation of concentration at equilibrium : in picomol/l
               ! coefficient for solubility for CFC-11/12 in mol/l/atm
               IF( tmask(ji,jj,1) .GE. 0.5 ) THEN
                  ztap = ( tsn(ji,jj,1,jp_tem ) + 273.15 ) * 0.01
                  iztap = 100. / ( tsn(ji,jj,1,jp_tem ) + 273.15 )
                  ! Temperature terms first
                  zdtap = soa(1,jl) + soa(2,jl)*iztap + soa(3,jl)*LOG(ztap) + soa(4,jl)*(ztap*ztap)
                  zdtap = zdtap + tsn(ji,jj,1,jp_sal)*( sob(1,jl) + sob(2,jl)*ztap + sob(3,jl)*(ztap*ztap))
                  zsol = EXP( zdtap )
               ELSE
                  zsol = 0.e0
               ENDIF
               ! conversion from mol/l/atm to mol/m3/atm and from mol/m3/atm to mol/m3/pptv
               zsol = xconv4 * xconv3 * zsol * tmask(ji,jj,1)
               ! concentration at equilibrium
               zca_cfc = xconv1 * zpp_cfc(ji,jj,jl) * zsol * tmask(ji,jj,1)
               ! Computation of speed transfert
               ! Schmidt number
               zsch = calc_schmidt_number(5, sca(:,jl), tsn(ji, jj, 1, jp_tem))
               ! speed transfert : formulae of wanninkhof 1992
               zv2 = wndm(ji,jj) * wndm(ji,jj)
               zsch = zsch / 660.
               zak_cfc = (( kw_scale * xconv2 * zv2 / SQRT(zsch) )*(1.-fr_i(ji,jj))) * tmask(ji,jj,1)
               ! Input function : speed *( conc. at equil - concen at surface )
               ! trn in pico-mol/l idem qtr; ak in en m/a
               qtr_cfc(ji,jj,jl) = -zak_cfc * ( trb(ji,jj,1,jn) - zca_cfc ) &
                  & * tmask(ji,jj,1) * ( 1. - fr_i(ji,jj) )
               ! Add the surface flux to the trend
               tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + qtr_cfc(ji,jj,jl) / e3t(ji,jj,1)
               ! cumulation of surface flux at each time step
               qint_cfc(ji,jj,jl) = qint_cfc(ji,jj,jl) + qtr_cfc(ji,jj,jl) * rdt
               ! !----------------!
            END DO ! end i-j loop !
         END DO !----------------!
         ! !----------------!
      END DO ! end CFC loop !
      ! !----------------!
      IF( ln_diatrc ) THEN
        !
        IF( lk_iomput ) THEN
           CALL iom_put( "cfc11qtr" , qtr_cfc (:,:,1) )
           CALL iom_put( "cfc11qint" , qint_cfc(:,:,1) )
           CALL iom_put( "cfc11patm" , zpp_cfc(:,:,1) )
           CALL iom_put( "cfc12qtr" , qtr_cfc (:,:,2) )
           CALL iom_put( "cfc12qint" , qint_cfc(:,:,2) )
           CALL iom_put( "cfc12patm" , zpp_cfc(:,:,2) )
           CALL iom_put( "sf6qtr" , qtr_cfc (:,:,3) )
           CALL iom_put( "sf6qint" , qint_cfc(:,:,3) )
           CALL iom_put( "sf6patm" , zpp_cfc(:,:,3) )
        ELSE
           DO jl = 1, jp_cfc
             trc2d(:,:,jp_cfc0_2d + 2*jl-2 ) = qtr_cfc (:,:,jl)
             trc2d(:,:,jp_cfc0_2d + 2*jl-1 ) = qint_cfc(:,:,jl)
           ENDDO
        END IF
        !
      END IF
      IF( l_trdtrc ) THEN
          DO jn = jp_cfc0, jp_cfc1
            CALL trd_mod_trc( tra(:,:,:,jn), jn, jptra_trd_sms, kt ) ! save trends
          END DO
      END IF
      !
      IF( nn_timing == 1 ) CALL timing_stop('trc_sms_cfc')
      !
   END SUBROUTINE trc_sms_cfc
   SUBROUTINE trc_cfc_cst
      !!---------------------------------------------------------------------
      !! *** trc_cfc_cst ***
      !!
      !! ** Purpose : sets constants for CFC model
      !!---------------------------------------------------------------------
      ! coefficient for CFC11
      !----------------------
      ! Solubility
      soa(1,1) = -229.9261
      soa(2,1) = 319.6552
      soa(3,1) = 119.4471
      soa(4,1) = -1.39165
      sob(1,1) = -0.142382
      sob(2,1) = 0.091459
      sob(3,1) = -0.0157274
      ! Schmidt number
      ! Values from Wanninkhof [2014]
      sca(1,1) = 3579.2
      sca(2,1) = -222.63
      sca(3,1) = 7.5749
      sca(4,1) = -0.14595
      sca(5,1) = 0.0011874
      ! coefficient for CFC12
      !----------------------
      ! Solubility
      soa(1,2) = -218.0971
      soa(2,2) = 298.9702
      soa(3,2) = 113.8049
      soa(4,2) = -1.39165
      sob(1,2) = -0.143566
      sob(2,2) = 0.091015
      sob(3,2) = -0.0153924
      ! schmidt number
      ! Values from Wanninkhof [2014]
      sca(1,2) = 3828.1
      sca(2,2) = -249.86
      sca(3,2) = 8.7603
      sca(4,2) = -0.1716
      sca(5,2) = 0.001408
      ! coefficient for SF6
      !----------------------
      ! Table 3, column 2 of Bullister, Wisegarver, and Menzia [Deep-Sea Research, 2001]
      ! Solubility
      soa(1,3) = -80.0343
      soa(2,3) = 117.232
      soa(3,3) = 29.5817
      soa(4,3) = 0. ! Unlike CFCs Solubility function for SF6 has no fourth temperature term
      sob(1,3) = 0.0335183
      sob(2,3) = -0.0373942
      sob(3,3) = 0.00774862
      ! schmidt number
      ! Values from Wanninkhof [2014]
      sca(1,3) = 3177.5
      sca(2,3) = -200.57
      sca(3,3) = 6.8865
      sca(4,3) = -0.13335
      sca(5,3) = 0.0010877
   END SUBROUTINE trc_cfc_cst
   INTEGER FUNCTION trc_sms_cfc_alloc()
      !!----------------------------------------------------------------------
      !! *** ROUTINE trc_sms_cfc_alloc ***
      !!----------------------------------------------------------------------
      ALLOCATE( xphem (jpi,jpj) , &
         & qtr_cfc (jpi,jpj,jp_cfc) , &
         & qint_cfc(jpi,jpj,jp_cfc) , STAT=trc_sms_cfc_alloc )
         !
      IF( trc_sms_cfc_alloc /= 0 ) CALL ctl_warn('trc_sms_cfc_alloc : failed to allocate arrays.')
      !
   END FUNCTION trc_sms_cfc_alloc
   !> Calculates the Schmidt numbr as a function of temperature via Horner's method
   FUNCTION calc_schmidt_number(n, coeffs, temp) RESULT(sc)
     INTEGER :: n !< Number of polynomial coefficients
     REAL, DIMENSION(n) :: coeffs !< Polynomial coefficieints in increasing order of degree
     REAL :: temp !< Temperature (in degC)
     REAL :: sc !< Schmidt number
     INTEGER k
     sc = coeffs(n)
     do k=n-1,1,-1
      sc = sc*temp + coeffs(k)
     enddo
   END FUNCTION calc_schmidt_number
   !!======================================================================
END MODULE trcsms_cfc
