MODULE trc_closea_canbgc

	!!----------------------------------------------------------------------
	!!                     ***  trc_closeabgc.F90  ***  
	!! TOP :   create a close sea mask to be used by BGCMs
	!!----------------------------------------------------------------------
	!! History :   0.0  !  2022-08 (O. Riche) 
  !!                  !  code reads domain_cfg.nc that should contain
  !!                  !  a variable named 'closea_bgc_mask'
  !!                  !  using info in make_closea_masks.py 
  !!                  !  location ./CanNEMO_tmp_src/nemo/tools/DOMAINcfg/
	!!----------------------------------------------------------------------

    USE par_oce         !: access jq* indices declaration
    USE dom_oce         !: access tmask declaration
    USE iom            !  I/O manager

	IMPLICIT NONE
    PRIVATE
  
    REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tmask_bgc_closea !: make the mask real(wp)

    PUBLIC trc_closea_init

CONTAINS
     
    SUBROUTINE trc_closea_init(read_var_flag)
    !
    USE lib_mpp , ONLY: ctl_stop
    !
    LOGICAL, OPTIONAL, INTENT(in) :: read_var_flag   ! 
    LOGICAL                       :: read_var_flag0  ! IF .FALSE. do not read closea_bgc_mask
    !
		INTEGER ::   ierr    ! Local variables  
    INTEGER ::   inum    ! input file identifier
    INTEGER ::   varid   ! variable identier
    INTEGER ::   jk
    REAL(wp), DIMENSION(jpi,jpj) :: zdata_in  ! temporary real array for input

    ! O Riche Aug 29th 2022
    ! create a solution to be
    ! able to switch on and off closea
    ! mask for BGCM.
    ! O Riche Aug 29th 2022
    ! Add the same failsafe than in the ocean physics
    ! in case the subroutine is called but no mask is
    ! available in domain_cfg.nc.
    ! O Riche Aug 22nd 2022
    ! Move closea_mask_bgc code from closea to
    ! to here.
    ! O Riche Aug 15th 2022
    ! This can be moved elsewhere but it is necessary
    ! as closea_mask_bgc/closea_mask are only allocated space
    ! when closea_mask exists in the domain_cfg.nc symlinked
    ! file.
		! initialize tmask_bgc_closea
		!
    IF ( .NOT. PRESENT(read_var_flag) ) THEN
      read_var_flag0 = .false.
    ELSE
      read_var_flag0 = read_var_flag
    ENDIF
    !
    ALLOCATE( tmask_bgc_closea(jpi,jpj,jpk) , STAT=ierr )
    !
    IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'trc_closea_init: failed to allocate tmask_bgc_closea array')
    tmask_bgc_closea(:,:,:) = tmask(:,:,:)
    !
    IF ( lwp )  WRITE(numout,*) 'trc_closea_init: tmask_bgc_closea initialized as tmask'
    IF ( lwp )  WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    IF ( lwp )  WRITE(numout,*)
    !
    IF ( read_var_flag0 ) THEN
      IF ( lwp )  WRITE(numout,*) 'trc_closea_init: prepping closed sea mask for BGCMs'
      IF ( lwp )  WRITE(numout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      IF ( lwp )  WRITE(numout,*)
      CALL iom_open(cn_domcfg, inum)
      varid = iom_varid(inum, 'closea_bgc_mask', ldstop = .false.)
      IF( varid > 0 ) THEN  ! the mask exists
        CALL iom_get(inum, jpdom_global,'closea_bgc_mask',zdata_in(:,:))
        CALL iom_close(inum)
        DO jk = 1, jpk
          tmask_bgc_closea(:,:,jk) = zdata_in(:,:)*tmask(:,:,jk)
        END DO
      ENDIF
    ENDIF
    !

    END SUBROUTINE trc_closea_init
     
END MODULE trc_closea_canbgc
