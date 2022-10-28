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

   USE trcopt_canbgc      ! PAR attenuation
   USE trcche_canbgc      ! carbon chemistry eq. constants
   USE trcflx_canbgc      ! air-flux gas exch.
   USE sms_top_canbgc     ! basic shared TOP variables, also contains ext. src array declarations

   USE cmocprod           ! CMOC PP module
   USE cmocmort           ! CMOC phyto mortality module
   USE cmocrem            ! CMOC carbon remineralization

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_cmoc       ! called by trcsms.F90 module
   PUBLIC   trc_sms_cmoc_alloc ! called by trcini_cmoc.F90 module 
   
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: qnegtr     ! Array used to indicate negative tracer values 

   ! Defined HERE the arrays specific to CMOC sms and ALLOCATE them in trc_sms_cmoc_alloc

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms_cmoc.F90 12841 2020-05-01 10:52:40Z cetlod $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_cmoc( kt )
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
      INTEGER  ::  jnt			        ! time (-step) splitting index
      INTEGER  ::  jn, ji, jj, jk   ! dummy loop indices
      INTEGER  ::  zrfact           ! working variable
      INTEGER  ::  jp_tot           ! jp_bgc+jp_cmoc
      REAL(wp) ::  ztra
      
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrmyt
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
      IF ( kt == nit000) THEN
        ! Calling external sources
        IF(lwp) WRITE(numout,*)
        IF(lwp) WRITE(numout,*) '   number of external sources to be read:', nb_src2d+nb_src3d
        IF(lwp) WRITE(numout,*) '   3D fields: ', nb_src3d
        IF(lwp) WRITE(numout,*) '   2D fields: ', nb_src2d
        IF(lwp) WRITE(numout,*)
        IF(lwp) WRITE(numout,*) '   ext. source array indices:'
        IF(lwp) WRITE(numout,*) '   js3d_no3 /NO3    : ', js3d_no3
        IF(lwp) WRITE(numout,*) '   js3d_si  /Si     : ', js3d_si
        IF(lwp) WRITE(numout,*) '   js3d_po4 /PO4    : ', js3d_po4
        IF(lwp) WRITE(numout,*) '   js3d_doc /DOC    : ', js3d_doc
        IF(lwp) WRITE(numout,*) '   js3d_fe  /Fer    : ', js3d_fe
        IF(lwp) WRITE(numout,*) '   js3d_hyfe/Fe bott: ', js3d_hyfe
        IF(lwp) WRITE(numout,*) '   js2d_chla/CHLA   : ', js2d_chla
        IF(lwp) WRITE(numout,*) '   js2d_dust /dust  : ', js2d_dust
        IF(lwp) WRITE(numout,*) '   js2d_par  /fr_par: ', js2d_par
        IF(lwp) WRITE(numout,*) '   js2d_femask      : ', js2d_femask
        IF(lwp) WRITE(numout,*) '   js2d_ndep        : ', js2d_ndep
        IF(lwp) WRITE(numout,*) '   js2d_rdic        : ', js2d_rdic
        IF(lwp) WRITE(numout,*) '   js2d_rdoc        : ', js2d_rdoc
        IF(lwp) WRITE(numout,*) '   js2d_rpoc        : ', js2d_rpoc
        IF(lwp) WRITE(numout,*) '   js2d_fsol1       : ', js2d_fsol1
        IF(lwp) WRITE(numout,*) '   js2d_fsol2       : ', js2d_fsol2
        CALL FLUSH(numout)
      ENDIF

      IF( kt == nittrc000 ) THEN
        !
        ALLOCATE( qnegtr(jpi,jpj,jpk) )
        !
        IF( .NOT. ln_rsttr ) THEN
          !
          qndayflxtr = nday_year
          !
          IF(lwp) write(numout,*)
          IF(lwp) write(numout,*) ' New chemical constants and various rates for biogeochemistry at new day : ', nday_year
          IF(lwp) write(numout,*) '~~~~~~'
          !
          CALL trc_che           ! computation of carbon chemistry constants
          ! initialize the chemical constants
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

        CALL trc_che           ! computation of carbon chemistry constants
      !
      ENDIF                            ! initialize the chemical constants

      CALL trc_flx( kt )     ! compute air-sea gas exchange    
      !
      ! O Riche Sept 14th 2022
      ! Move here before cmoc_prod as issue with PAR being set to 0s
      ! at initialization (current state)
      ! also need to add time splitting loop 1=> qnrdttrc
      ! and so trc_opt_1band and trc_opt (CanOE)
      ! needs jnt index/input arg along with kt see below
      ! for cmoc_prod.
      !
      CALL trc_opt_1band( kt )        ! 1-band PAR attenuation
      !
      DO jnt = 1, qnrdttrc             ! Potential time splitting if requested
        CALL cmoc_prod( kt, jnt )      ! PP subroutine
      END DO
      !
      CALL cmoc_mort( kt )
      !
      CALL cmoc_rem( kt )
      !
      ! Is this below necessary? (NEMO3.4.1 code)
      ! DO jn = jp_bgc+1, jp_bgc+jp_cmoc
        ! CALL lbc_lnk( trn(:,:,:,jn), 'T', 1. )
        ! CALL lbc_lnk( trb(:,:,:,jn), 'T', 1. )
        ! CALL lbc_lnk( tra(:,:,:,jn), 'T', 1. )
      ! END DO
      ! !
      ! IF( l_trdtrc )  ALLOCATE( ztrmyt(jpi,jpj,jpk) )
      ! !
      ! ! Save the trends in the mixed layer
      IF( l_trdtrc ) THEN
          DO jn = 1, jp_tot
            ztrmyt(:,:,:) = tra(:,:,:,jn)
            CALL trd_trc( ztrmyt, jn, jptra_sms, kt )   ! save trends
          END DO
          DEALLOCATE( ztrmyt )
      END IF
      !
      ! O Riche DBG Oct 21st 2022
      IF( lwp .AND. kt == nittrc000 ) THEN
        WRITE(numout,*)
        WRITE(numout,*) 'Checking trn index attribution:'
        WRITE(numout,*) 'jqdic = ', jqdic
        WRITE(numout,*) 'jqtal = ', jqtal
        WRITE(numout,*) 'jqoxy = ', jqoxy
        WRITE(numout,*) 'jqno3 = ', jqno3
        WRITE(numout,*) 'jqpoc = ', jqpoc
        WRITE(numout,*) 'jqphy = ', jqphy
        WRITE(numout,*) 'jqnch = ', jqnch
        WRITE(numout,*) 'jqzoo = ', jqzoo
        WRITE(numout,*) 'jp_age =', jp_age
      CALL FLUSH(numout)
      ENDIF
      !
      qnegtr(:,:,:) = 1.e0
      ! O Riche Oct 25th 2022
      ! test value of jp_tot to see if jp_age is involved
      IF( lwp .AND. kt == nittrc000 ) THEN
        WRITE(numout,*) 'trc_sms_cmoc: jp_age and jp_tot check'
        WRITE(numout,*) 'jp_age = ', jp_age
        WRITE(numout,*) 'jp_tot = ', jp_tot
      ENDIF
      
      DO jn = 1, jp_tot
        DO jk = 1, jpk
           DO jj = 1, jpj
              DO ji = 1, jpi
                 IF( ( trb(ji,jj,jk,jn) + tra(ji,jj,jk,jn) ) < 0.e0 ) THEN
                    ztra             = ABS( trb(ji,jj,jk,jn) ) / ( ABS( tra(ji,jj,jk,jn) ) + rtrn )
                    qnegtr(ji,jj,jk) = MIN( qnegtr(ji,jj,jk),  ztra )
                 ENDIF
             END DO
           END DO
        END DO
      END DO
      !                                ! where at least 1 tracer concentration becomes negative
      !                                ! and by tracer we mean only the CMOC or shared BGC tracer.
      DO jn = 1, jp_tot 
        trb(:,:,:,jn) = trb(:,:,:,jn) + qnegtr(:,:,:) * tra(:,:,:,jn)
        tra(:,:,:,jn) = 0._wp
      END DO
      !
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
