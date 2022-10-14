MODULE trcsms_canoe
   !!======================================================================
   !!                         ***  MODULE trcsms_canoe  ***
   !! TOP :   Main module of the CANOE tracers
   !!======================================================================
   !! History :      !  2007  (C. Ethe, G. Madec)  Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_sms_canoe       : CANOE model main routine
   !! trc_sms_canoe_alloc : allocate arrays specific to CANOE sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc

   USE trcopt_canbgc      ! PAR attenuation
   USE trcche_canbgc      ! carbon chemistry eq. constants
   USE trcflx_canbgc      ! air-flux gas exch.
   USE sms_top_canbgc     ! basic shared TOP variables, also contains ext. src array declarations

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_canoe       ! called by trcsms.F90 module
   PUBLIC   trc_sms_canoe_alloc ! called by trcini_canoe.F90 module

   ! Defined HERE the arrays specific to CANOE sms and ALLOCATE them in trc_sms_canoe_alloc

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms_canoe.F90 12841 2020-05-01 10:52:40Z cetlod $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_canoe( kt )
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
      INTEGER ::   jn				        ! dummy loop index
      INTEGER ::   zrfact           ! working variable
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrmyt
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_canoe')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_canoe:  CANOE model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~'

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

      ! ! O Riche Aug 16th 2022
      ! ! for now this is a placeholder
      ! ! surface chlorophyll (and phosphate will in trcflx)
      ! ! are the only fields used so far
      ! zrfact = 86400 / rdttrc * 14 !!!! time steps per 2 weeks
      ! IF ( MOD(kt,zrfact) == 0 ) THEN
        ! CALL trc_src3d( kt, js3d_si    )    
        ! CALL trc_src3d( kt, js3d_no3   )   
        ! CALL trc_src3d( kt, js3d_po4   )   
        ! CALL trc_src3d( kt, js3d_doc   )   
        ! CALL trc_src3d( kt, js3d_fe    )    
        ! CALL trc_src3d( kt, js3d_hyfe  )  
        ! !
        ! CALL trc_src2d( kt, js2d_chla  )      
        ! CALL trc_src2d( kt, js2d_par   )      
        ! CALL trc_src2d( kt, js2d_dust  )     
        ! CALL trc_src2d( kt, js2d_femask)    
        ! CALL trc_src2d( kt, js2d_ndep  )     
        ! CALL trc_src2d( kt, js2d_rdoc  )      
        ! CALL trc_src2d( kt, js2d_rdic  )      
        ! CALL trc_src2d( kt, js2d_rpoc  )      
        ! CALL trc_src2d( kt, js2d_fsol1 )      
        ! CALL trc_src2d( kt, js2d_fsol2 )      
      ! END IF

      IF( qndayflxtr /= nday_year ) THEN      ! New days
        !
        qndayflxtr = nday_year
  
        IF(lwp) write(numout,*)
        IF(lwp) write(numout,*) ' New chemical constants and various rates for biogeochemistry at new day : ', nday_year
        IF(lwp) write(numout,*) '~~~~~~'
  
        CALL trc_che           ! computation of carbon chemistry constants
            !
      ENDIF
  
      CALL trc_flx( kt )     ! compute air-sea gas exchange    
      !
      CALL trc_opt( kt )     ! test PAR attenuation

      ! ! O Riche Aug 26th 2022
      ! ! test trc_src_fedep and trc_src_fesed
      zrfact = 86400 / rdttrc * 14 !!!! time steps per 2 weeks
      CALL trc_src_fedep( kt )
      ! IF ( kt == nit000 .OR. MOD(kt,zrfact) == 0 )  CALL trc_src_fedep( kt )
      CALL trc_src_fesed        ! This source does not vary with time
      ! IF ( kt == nit000 .OR. MOD(kt,zrfact) == 0 )  CALL trc_src_fesed        ! This source does not vary with time
      
      
      !
      IF( l_trdtrc )  ALLOCATE( ztrmyt(jpi,jpj,jpk) )

      ! Save the trends in the mixed layer
      IF( l_trdtrc ) THEN
          DO jn = 1, jp_canoe
            ztrmyt(:,:,:) = tra(:,:,:,jn)
            CALL trd_trc( ztrmyt, jn, jptra_sms, kt )   ! save trends
          END DO
          DEALLOCATE( ztrmyt )
      END IF
	  
      !
      IF( ln_timing )   CALL timing_stop('trc_sms_canoe')
      !
   END SUBROUTINE trc_sms_canoe


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
