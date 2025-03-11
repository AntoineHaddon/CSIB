MODULE trcsms_canoe
   !!======================================================================
   !!                         ***  MODULE trcsms_canoe  ***
   !! TOP :   Main module of the CANOE tracers
   !!======================================================================
   !! History :      !  2007  (C. Ethe, G. Madec)  Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!                !  2022  (O. Riche) NEMO4 integration  
   !!----------------------------------------------------------------------
   !! trc_sms_canoe       : CANOE model main routine
   !! trc_sms_canoe_alloc : allocate arrays specific to CANOE sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc
   !
   USE lbclnk             ! exchange fields over tile boundaries
   !USE lib_fortran
   !
   USE trcopt_canbgc      ! PAR attenuation
   USE trcche_canbgc      ! carbon chemistry eq. constants
   USE trcflx_canbgc      ! air-flux gas exch.
   USE trcsink_canbgc     ! particule sinking
   !
   USE sms_top_canbgc     ! basic shared TOP variables, also contains ext. src array declarations
   USE sms_canoe
   !USE sms_canoe, ONLY    : ln_canoenegtr
   !
   USE canoetemp          ! CanOE temperature dependencies module   
   USE canoeprod          ! CanOE PP module
   USE canoenzd           ! CanOE grazing/mortality/remineralization module

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_canoe       ! called by trcsms.F90 module
   PUBLIC   trc_sms_canoe_alloc ! called by trcini_canoe.F90 module
   !PUBLIC   trc_sms_canoe_dmp

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: rnegtr2     ! Array used to indicate negative tracer values 
   LOGICAL , PUBLIC ::   ll_sbc  ! trigger for external sources (ln_dust0, ln_river0, and ln_ndepo0)
 
   ! Defined HERE the arrays specific to CANOE sms and ALLOCATE them in trc_sms_canoe_alloc

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms_canoe.F90 12841 2020-05-01 10:52:40Z cetlod $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_canoe( kt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_canoe  ***
      !!
      !! ** Purpose :   main routine of CANOE model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      USE par_canoe
      USE trcsrc_canbgc             ! loading external files/sources
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      INTEGER  ::  jnt              ! time (-step) splitting index
      INTEGER  ::  jn, ji, jj, jk   ! dummy loop indices
      INTEGER  ::  zrfact           ! working variable      
      INTEGER  ::  jp_tot           ! jp_bgc+jp_cmoc
      REAL(wp) ::  ztra

      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrmyt
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: rtrbbio
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_canoe')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_canoe:  CANOE model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      !
      ! Sum of all the tracers shared TOP + CANOE
      jp_tot = jp_bgc + jp_canoe
      !
      ALLOCATE(rnegtr2(jpi,jpj,jpk))      
      ALLOCATE(rtrbbio(jpi,jpj,jpk,jp_tot))
      !
      IF( ln_dust0 .OR. ln_river0 .OR. ln_ndepo0 ) THEN   ;   ll_sbc = .TRUE.
      ELSE                                                ;   ll_sbc = .FALSE.
      ENDIF

      !
      ! IF ( kt == nit000) THEN
        ! ! Calling external sources
        ! IF(lwp) WRITE(numout,*)
        ! IF(lwp) WRITE(numout,*) '   number of external sources to be read:', nb_src2d+nb_src3d
        ! IF(lwp) WRITE(numout,*) '   3D fields: ', nb_src3d
        ! IF(lwp) WRITE(numout,*) '   2D fields: ', nb_src2d
        ! IF(lwp) WRITE(numout,*)
        ! IF(lwp) WRITE(numout,*) '   ext. source array indices:'
        ! IF(lwp) WRITE(numout,*) '   js3d_no3 /NO3    : ', js3d_no3
        ! IF(lwp) WRITE(numout,*) '   js3d_si  /Si     : ', js3d_si
        ! IF(lwp) WRITE(numout,*) '   js3d_po4 /PO4    : ', js3d_po4
        ! IF(lwp) WRITE(numout,*) '   js3d_doc /DOC    : ', js3d_doc
        ! IF(lwp) WRITE(numout,*) '   js3d_fe  /Fer    : ', js3d_fe
        ! IF(lwp) WRITE(numout,*) '   js3d_hyfe/Fe bott: ', js3d_hyfe
        ! IF(lwp) WRITE(numout,*) '   js2d_chla/CHLA   : ', js2d_chla
        ! IF(lwp) WRITE(numout,*) '   js2d_dust /dust  : ', js2d_dust
        ! IF(lwp) WRITE(numout,*) '   js2d_par  /fr_par: ', js2d_par
        ! IF(lwp) WRITE(numout,*) '   js2d_femask      : ', js2d_femask
        ! IF(lwp) WRITE(numout,*) '   js2d_ndep        : ', js2d_ndep
        ! IF(lwp) WRITE(numout,*) '   js2d_rdic        : ', js2d_rdic
        ! IF(lwp) WRITE(numout,*) '   js2d_rdoc        : ', js2d_rdoc
        ! IF(lwp) WRITE(numout,*) '   js2d_rpoc        : ', js2d_rpoc
        ! IF(lwp) WRITE(numout,*) '   js2d_fsol1       : ', js2d_fsol1
        ! IF(lwp) WRITE(numout,*) '   js2d_fsol2       : ', js2d_fsol2
        ! CALL FLUSH(numout)
      ! ENDIF

      IF( kt == nittrc000 ) THEN       
        !
        IF( .NOT. ln_rsttr ) THEN
          !
          qndayflxtr = nday_year
          !
          IF(lwp) write(numout,*)
          IF(lwp) write(numout,*) ' New chemical constants and various rates for biogeochemistry at new day : ', nday_year
          IF(lwp) write(numout,*) '~~~~~~'
          !
          CALL trc_che_2D( kt, Kmm )   ! initialization of carbon chemistry constants
        ELSE
            WRITE(numout,*)
            WRITE(numout,*) 'Should something be done for the restart mode here? Nothing coded here yet, some code exists in TOP/trcini.F90 to take care of this though.'
            WRITE(numout,*)
        ENDIF
         !
      ENDIF
      ! 
      IF( qndayflxtr /= nday_year ) THEN      ! New days
        !
        qndayflxtr = nday_year
  
        IF(lwp) write(numout,*)
        IF(lwp) write(numout,*) ' New chemical constants and various rates for biogeochemistry at new day : ', nday_year
        IF(lwp) write(numout,*) '~~~~~~'
  
        CALL trc_che_2D( kt,Kmm )           ! computation of carbon chemistry constants
        CALL trc_che_3D( kt,Kmm )           ! computation of carbon chemistry constants
            !
      ENDIF  
      !
      ! Update temperature dependencies
      CALL canoe_temp( Kmm )
      !
      DO jn = 1, jp_tot                    !   Store the tracer concentrations before entering CMOC
        rtrbbio(:,:,:,jn) = tr(:,:,:,jn, Kbb)
      END DO
      !  
      DO jnt = 1, qnrdttrc             ! Potential time splitting if requested
        !
        !!!!!!! Start of "p4zbio" block !!!!!!! 
        ! This is the equivalent of p4z_bio call
        ! in trcsms_pisces.F90/CanESM5/CANOE
        !
        WRITE(numout,*)
        WRITE(numout,*) 'time step            #', kt
        WRITE(numout,*) 'split loop iteration #', jnt
        WRITE(numout,*) 'narea                #', narea
        WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        CALL FLUSH(numout)
       !
        CALL canoe_sink( kt , jnt, Kbb, Kmm, Krhs )     ! particle sinking 
       !
        CALL trc_opt( kt, jnt, Kbb, Kmm )       ! 3-band PAR attenuation
       !
       ! call CanOE production S/R
        CALL canoe_prod( kt, jnt, Kbb, Kmm, Krhs  )
        CALL canoe_meso( kt, jnt, Kbb, Kmm, Krhs )
        CALL canoe_mzoo( kt, jnt, Kbb, Kmm, Krhs )
        CALL canoe_mort1( kt, jnt, Kbb, Kmm, Krhs  )
        CALL canoe_mort2( kt, jnt, Kbb, Kmm, Krhs  )
        CALL canoe_rem( kt, jnt, Kbb, Kmm, Krhs  )
        CALL trc_src_fedep( kt, Kbb, Kmm, Krhs )
        CALL trc_src_fesed( Kbb, Kmm, Krhs )
        CALL trc_n2fx_canoe( kt, jnt, Kbb, Kmm, Krhs )
        CALL trc_sms_canoe_dmp( kt, Kbb, Kmm, Krhs )
        !
        ! Initialize rnegtr2, if no call to trc_xnegtr tra used w/o correction
        rnegtr2(:,:,:) = 1._wp
        !
        IF( ln_canoenegtr )  CALL trc_xnegtr( 1, jp_tot, Kbb, Kmm, Krhs, rnegtr2 )   !!! O Riche Nov 8th 2022 ! reside in sms_top_canbgc.F90
        DO jn = 1, jp_tot
          tr(:,:,:,jn, Kbb) = tr(:,:,:,jn, Kbb) + rnegtr2(:,:,:) * tr(:,:,:,jn, Krhs)        
          tr(:,:,:,jn, Krhs) = 0._wp
        END DO
        !  
        CALL trc_flx( kt, Kbb, Kmm, Krhs )     ! compute air-sea gas exchange 
      !
      END DO
      !
      DO jn = 1, jp_tot
        tr(:,:,:,jn, Krhs) = ( tr(:,:,:,jn, Kbb) - rtrbbio(:,:,:,jn) ) * qfactr
        tr(:,:,:,jn, Kbb) = rtrbbio(:,:,:,jn)
        rtrbbio(:,:,:,jn) = 0._wp
      END DO
      !
      DEALLOCATE(rnegtr2)
      DEALLOCATE(rtrbbio)
      ! 
      !CALL total_element(totfe,totn)
      !WRITE(numout,*) totfe, totn
      !
      IF( l_trdtrc )  ALLOCATE( ztrmyt(jpi,jpj,jpk) )

      ! Save the trends in the mixed layer
      IF( l_trdtrc ) THEN
          DO jn = 1, jp_tot
            ztrmyt(:,:,:) = tr(:,:,:,jn, Krhs)
            CALL trd_trc( ztrmyt, jn, jptra_sms, kt, Kmm )   ! save trends
          END DO
          DEALLOCATE( ztrmyt )
      END IF
  
      !
      IF( ln_timing )   CALL timing_stop('trc_sms_canoe')
      !
   END SUBROUTINE trc_sms_canoe
   
   SUBROUTINE trc_sms_canoe_dmp( kt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                    ***  trc_sms_canoe_dmp  ***
      !!
      !! adjust global total N and TA while leaving spatial distribution intact
      !! enforces global conservation in case that DNF is larger or smaller than denitrification
      !!----------------------------------------------------------------------
      !
      USE lib_fortran,   ONLY: glob_sum

      INTEGER, INTENT( in )  ::     kt, Kbb, Kmm, Krhs
      !
      REAL(wp) :: zdntrsum, zdnfsum, ztau, znsum, zalksum  !   , zalksum0
      !!---------------------------------------------------------------------

         ! adjust NO3 according to difference between global total rates of denitrification and N2 fixation
         ! in CanESM5 this adjustment was made multiplicative rather than additive to prevent negative concentrations
         ! however that approach can not work in NEMO4 as we are avoiding any direct alteration to the tracer concentration arrays
         ! it might work if we apply it only to a subset of the model layers
         ! at present denitr is in mmol/m^3/tstep and zn2fix is in mmol/m^3/d

         znsum = glob_sum('totalNO3', tr(:,:,35:74,jqno3,Kbb)*cvol(:,:,35:74) )               ! global total N in mmol^-1
         zdnfsum = glob_sum('totalDNF',  zn2fix(:,:,:) * cvol(:,:,:)  )*xstepb                ! global total in mmolN tstep^-1
         zdntrsum = glob_sum('totalDenitr',  denitr(:,:,:) * cvol(:,:,:)  )                   ! global total in mmolN tstep^-1
         ztau = (zdntrsum-zdnfsum)/(znsum*1.0008)                                             ! 1.0008 is an approximate correction for non-NO3 N
         tr(:,:,35:74,jqno3,Krhs) = tr(:,:,35:74,jqno3,Krhs) + tr(:,:,35:74,jqno3,Kbb) * ztau           ! correction to NO3
! add/remove alkalinity to compensate for processes bypassed in nonphysical addition/subtraction of NO3 (- 1 mol / mol N)
         zalksum = glob_sum('totalTA',  tr(:,:,35:74,jqtal,Kbb)*cvol(:,:,35:74) )
         tr(:,:,35:74,jqtal,Krhs) = tr(:,:,35:74,jqtal,Krhs) - tr(:,:,35:74,jqtal,Kbb) * ztau * znsum/zalksum * 1.e-6
         !tr(:,:,35:74,jqtal,Krhs) = tr(:,:,35:74,jqtal,Krhs) - tr(:,:,35:74,jqno3,Kbb) * ztau * 1.e-6

   END SUBROUTINE trc_sms_canoe_dmp

   SUBROUTINE total_element(totfe,totn, Kmm)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE total_element  ***
      !!
      !! ** Purpose :   calculate model total N and Fe   2015/03/23 JRC
      !!---------------------------------------------------------------------
      !
      USE lib_fortran,   ONLY: glob_sum
      USE sms_canoe,     ONLY: rr_c2n, rr_fe2c
      !
      INTEGER, INTENT(in) ::    Kmm  ! time level indices
      REAL(wp) :: totfe, totn
      !!---------------------------------------------------------------------

      totfe = 0._wp
      totn  = 0._wp

      totn = glob_sum( 'total_element',                                   &
                         (   (tr(:,:,:,jqno3, Kmm)   + tr(:,:,:,jrnh4, Kmm)       &
      &                     + tr(:,:,:,jrnn, Kmm)  + tr(:,:,:,jrdn, Kmm))*rr_c2n  &
      &                     + tr(:,:,:,jrzoo, Kmm) + tr(:,:,:,jrmes, Kmm)         &
      &                     + tr(:,:,:,jrpoc, Kmm) + tr(:,:,:,jrgoc, Kmm)  ) * cvol(:,:,:)  )
      totn = totn/rr_c2n

      ! zoo and poc concs are in C units

      totfe = glob_sum( 'total_element',                                               &
                          (   tr(:,:,:,jrfer, Kmm) + tr(:,:,:,jrdfe, Kmm) + tr(:,:,:,jrnfe, Kmm)   &
      &                     + tr(:,:,:,jrzoo, Kmm)*rr_fe2c + tr(:,:,:,jrmes, Kmm)*rr_fe2c      &
      &                     + tr(:,:,:,jrpoc, Kmm)*rr_fe2c + tr(:,:,:,jrgoc, Kmm)*rr_fe2c  ) * cvol(:,:,:)  )

      !
   END SUBROUTINE total_element

   INTEGER FUNCTION trc_sms_canoe_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_canoe_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to CANOE
      ! ALLOCATE( tab(...) , STAT=trc_sms_canoe_alloc )
      trc_sms_canoe_alloc = 0      ! set to zero if no array to be allocated
      !
      IF( trc_sms_canoe_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_canoe_alloc : failed to allocate arrays' )
      !

   END FUNCTION trc_sms_canoe_alloc

   !!======================================================================
END MODULE trcsms_canoe


