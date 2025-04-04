MODULE trcsms_cmoc
   !!======================================================================
   !!                         ***  MODULE trcsms_cmoc  ***
   !! TOP :   Main module of the CMOC tracers
   !!======================================================================
   !! History :      !  2007  (C. Ethe, G. Madec)  Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!                !  2022  (O. Riche) NEMO4 integration  
   !!----------------------------------------------------------------------
   !! trc_sms_cmoc       : CMOC model main routine
   !! trc_sms_cmoc_alloc : allocate arrays specific to CMOC sms
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc
   !
   USE lbclnk             ! exchange fields over tile boundaries
   !
   USE trcopt_canbgc      ! PAR attenuation
   USE trcche_canbgc      ! carbon chemistry eq. constants
   USE trcflx_canbgc      ! air-flux gas exch.
   USE trcsink_canbgc     ! particule sinking
   !
   USE sms_top_canbgc     ! basic shared TOP variables, also contains ext. src array declarations
   USE sms_cmoc, ONLY     : ln_cmocnegtr
   !
   USE cmocprod           ! CMOC PP module
   ! USE cmocmort            ! CMOC phyto mortality module
   ! USE cmocrem             ! CMOC carbon remineralization
   ! USE cmoczoo             ! CMOC zooplankton grazing
   USE cmocnzd             ! consolidated module containing remineralization to (N)itrate
                           ! (Z)ooplankton grazing, and (D)etritus for phytoplankton mortality
   !
   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_sms_cmoc       ! called by trcsms.F90 module
   PUBLIC trc_sms_cmoc_alloc ! called by trcini_cmoc.F90 module 
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: qnegtr2     ! Array used to indicate negative tracer values 
   LOGICAL , PUBLIC ::   ll_sbc  ! trigger for external sources (ln_dust0, ln_river0, and ln_ndepo0)
   
   !
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms_cmoc.F90 12841 2020-05-01 10:52:40Z cetlod $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_cmoc( kt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_cmoc  ***
      !!
      !! ** Purpose :   main routine of CMOC model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      USE par_cmoc
      USE trcsrc_canbgc             ! loading external files/sources
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      INTEGER  ::  jnt			        ! time (-step) splitting index
      INTEGER  ::  jn, ji, jj, jk   ! dummy loop indices
      INTEGER  ::  zrfact           ! working variable
      INTEGER  ::  jp_tot           ! jp_bgc+jp_cmoc
      REAL(wp) ::  ztra
      
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrmyt
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: qtrbbio
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_cmoc')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_cmoc:  CMOC model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      !
      ! Sum of all the tracers shared TOP + CMOC
      jp_tot = jp_bgc + jp_cmoc
      !
      ALLOCATE(qnegtr2(jpi,jpj,jpk))
      ALLOCATE(qtrbbio(jpi,jpj,jpk,jp_tot))
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
          CALL trc_che_2D( kt, Kmm )   ! initialize the carbon chemistry constants
          !
        !  
        ELSE
            WRITE(numout,*)
            WRITE(numout,*) 'Should something be done for the restart mode here? Nothing coded here yet, some code exists in TOP/trcini.F90 to take care of this though.'
            WRITE(numout,*)
        ENDIF
          !
      ENDIF
      ! 
      ! Do we need this or is this covered at least partly by all the new
      ! external sources subroutines, e.g. trcsrc.F90 modules.
      ! anything else to add?
      ! ?IF( ll_sbc ) CALL p4z_sbc( kt )   ! external sources of nutrients
      ! Do we need a CMOC- and CanOE-specific *_sbc.F90 file?
      !
      IF( qndayflxtr /= nday_year ) THEN      ! New days
        !
        qndayflxtr = nday_year

        IF(lwp) write(numout,*)
        IF(lwp) write(numout,*) ' New chemical constants and various rates for biogeochemistry at new day : ', nday_year
        IF(lwp) write(numout,*) '~~~~~~'

        CALL trc_che_2D( kt, Kmm )   ! computation of carbon chemistry constants
        !
      ENDIF
      !
      DO jn = 1, jp_tot                    !   Store the tracer concentrations before entering CMOC
        qtrbbio(:,:,:,jn) = tr(:,:,:,jn, Kbb)
      END DO
      !
      ! O Riche Sept 14th 2022
      ! Move here before cmoc_prod as issue with PAR being set to 0s
      ! at initialization (current state)
      ! also need to add time splitting loop 1=> qnrdttrc
      ! and so trc_opt_1band and trc_opt (CanOE)
      ! needs jnt index/input arg along with kt see below
      ! for cmoc_prod.
      !
      DO jnt = 1, qnrdttrc             ! Potential time splitting if requested
        !
        !!!!!!! Start of "p4zbio" block !!!!!!! 
        ! This is the equivalent of p4z_bio call
        ! in trcsms_pisces.F90/CanESM5/CMOC
        !
        ! trcsink calls go here according to p4z_bio
        ! O Riche Jan 18th 2023
        ! will have to make sure the time stepping scheme
        ! properly work with NEMO4/PISCES/TRP infrastructure
        ! before activating.
        !
        ! Test print narea
        WRITE(numout,*)
        WRITE(numout,*) 'time step            #', kt
        WRITE(numout,*) 'split loop iteration #', jnt
        WRITE(numout,*) 'narea                #', narea
        WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        CALL FLUSH(numout)
        !
        CALL cmoc_sink( kt , jnt, Kbb, Kmm, Krhs  )     ! particule sinking 
        !
        CALL trc_opt_1band( kt, jnt, Kmm )  ! 1-band PAR attenuation
        !
        CALL cmoc_prod( kt, jnt, Kbb, Kmm, Krhs )      ! Primary production
        !
        CALL cmoc_rem( kt, jnt, Kbb, Kmm, Krhs )       ! NO3 remineralization
        !
        CALL cmoc_mort( kt, jnt, Kbb, Kmm, Krhs )      ! Phyto mortality
        !
        CALL cmoc_zoo( kt, jnt, Kbb, Kmm, Krhs )       ! Zooplankton grazing and mortality
        !
        ! Here CMOC would call the new subroutines that
        ! compute the various sources that were scattered
        ! within CanESM5/CMOC p4zsed.F90 code, e.g.
        ! river sources
        ! Formely p4z_sbc in p4zsed.F90
        !! OR Jan19 23 ! IF ( jnt == 1 .AND. ll_sbc ) CALL trc_src_criver( kt, Krhs )
        ! POC bottom instant. rem
        CALL trc_bott_cmoc( Kmm, Krhs)
        ! n2 fixation/denitrification
        !! OR Jan19 23 ! CALL cmoc_rem_denit
        CALL cmoc_rem_denit( Kmm )
        CALL trc_n2fx_denit_cmoc( par_1band, kt, jnt, Kbb, Kmm, Krhs )
        !
        ! Move here to be consistent with NEMO4 and sidestepping from CanESM5 CMOC NEMO
        CALL trc_flx(kt, Kbb, Kmm, Krhs)               ! compute air-sea gas exchange
        ! IF the radioactive tracer was added there would be also a call to p4z_dcy( kt ) equivalent (trc_dcy?) here. 
        !       
        !!!!!! O Riche Nov 8th 2022
        !!!!!! replace this by a call to trc_xnegtr subroutine
        !!!!!! sitting higher in CANBGC
        ! Enforce conservation and positive values of tracers
        ! by adjusting the time step using tra trend
        !
        ! Initialize qnegtr2, if no call to trc_xnegtr tra used w/o correction
        qnegtr2(:,:,:) = 1._wp
        !
        IF( ln_cmocnegtr )  CALL trc_xnegtr( 1, jp_tot, Kbb, Kmm, Krhs, qnegtr2 )   !!! O Riche Nov 8th 2022 ! reside in sms_top_canbgc.F90
        DO jn = 1, jp_tot
          tr(:,:,:,jn, Kbb) = tr(:,:,:,jn, Kbb) + qnegtr2(:,:,:) * tr(:,:,:,jn, Krhs)        
          tr(:,:,:,jn, Krhs) = 0._wp
        END DO
        !  
        !!!!!!! End   of "p4zbio" block !!!!!!!        
        !
        !!!!!!! Start of "p4zsed" block !!!!!!!
        !
        !!!!!!! End   of "p4zsed" block !!!!!!!
        !
      END DO
      !
      DO jn = 1, jp_tot
        tr(:,:,:,jn, Krhs) = ( tr(:,:,:,jn, Kbb) - qtrbbio(:,:,:,jn) ) * qfactr
        tr(:,:,:,jn, Kbb) = qtrbbio(:,:,:,jn)
        qtrbbio(:,:,:,jn) = 0._wp
      END DO
      !
      DEALLOCATE(qnegtr2)
      DEALLOCATE(qtrbbio)
      !
      !! OR Jan 24th 2023
      !! markers to track values of the tracers
      
      !! End of OR Jan 24th 2023
      !
      ! Save the trends in the mixed layer
      IF( l_trdtrc ) THEN
          ALLOCATE( ztrmyt(jpi,jpj,jpk) )
          DO jn = 1, jp_tot
            ztrmyt(:,:,:) = tr(:,:,:,jn, Krhs)
            CALL trd_trc( ztrmyt, jn, jptra_sms, kt, Kmm )   ! save trends
          END DO
          DEALLOCATE( ztrmyt )
      END IF
      !
      ! ! O Riche DBG Oct 21st 2022
      ! IF( lwp .AND. kt == nittrc000 ) THEN
        ! WRITE(numout,*)
        ! WRITE(numout,*) 'Checking trn index attribution:'
        ! WRITE(numout,*) 'jqdic = ', jqdic
        ! WRITE(numout,*) 'jqtal = ', jqtal
        ! WRITE(numout,*) 'jqoxy = ', jqoxy
        ! WRITE(numout,*) 'jqno3 = ', jqno3
        ! WRITE(numout,*) 'jqpoc = ', jqpoc
        ! WRITE(numout,*) 'jqphy = ', jqphy
        ! WRITE(numout,*) 'jqnch = ', jqnch
        ! WRITE(numout,*) 'jqzoo = ', jqzoo
        ! WRITE(numout,*) 'jp_age =', jp_age
      ! CALL FLUSH(numout)
      ! ENDIF
      ! !
      ! ! O Riche Oct 25th 2022
      ! ! test value of jp_tot to see if jp_age is involved
      ! IF( lwp .AND. kt == nittrc000 ) THEN
        ! WRITE(numout,*) 'trc_sms_cmoc: jp_age and jp_tot check'
        ! WRITE(numout,*) 'jp_age = ', jp_age
        ! WRITE(numout,*) 'jp_tot = ', jp_tot
      ! ENDIF
      ! !
      IF( ln_timing )   CALL timing_stop('trc_sms_cmoc')
      !
   END SUBROUTINE trc_sms_cmoc


   INTEGER FUNCTION trc_sms_cmoc_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_cmoc_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to CMOC
      ! ALLOCATE( tab(...) , STAT=trc_sms_cmoc_alloc )
      trc_sms_cmoc_alloc = 0      ! set to zero if no array to be allocated
      !
      IF( trc_sms_cmoc_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_cmoc_alloc : failed to allocate arrays' )
      !
   END FUNCTION trc_sms_cmoc_alloc
    
END MODULE trcsms_cmoc
