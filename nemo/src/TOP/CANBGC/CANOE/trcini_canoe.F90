MODULE trcini_canoe
   !!======================================================================
   !!                         ***  MODULE trcini_canoe  ***
   !! TOP :   initialisation of the CANOE tracers
   !!======================================================================
   !! History :        !  2007  (C. Ethe, G. Madec) Original code
   !!                  !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_ini_canoe   : CANOE model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE par_canoe
   USE trcnam_canoe     ! CANOE SMS namelist
   USE trcsms_canoe

   USE trc_closea_canbgc ! bgc closea mask
   USE trcflx_canbgc     ! air-sea gas exch.
   USE trcche_canbgc     ! carbon chemistry 
   USE trcsrc_canbgc     ! external sources/other data 
   USE sms_top_canbgc    ! access ext. source arrays declaration

   USE trcopt_canbgc     ! PAR attenuation
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_canoe   ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcini_canoe.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_canoe
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_canoe  ***  
      !!
      !! ** Purpose :   initialization for CANOE model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
      !
      INTEGER  :: jn
      CHARACTER(len = 20)  ::  cltra
      REAL(wp), SAVE ::   sco2   =  2.312e-3_wp
      REAL(wp), SAVE ::   alka0  =  2.426e-3_wp
      REAL(wp), SAVE ::   oxyg0  =  177.6e-6_wp
      REAL(wp), SAVE ::   no30   =    5.0e-6_wp
      ! !
      CALL trc_nam_canoe   
	    !                       ! Allocate CANOE arrays
      IF( trc_sms_canoe_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_canoe: unable to allocate CANOE arrays' )

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_canoe: passive tracer unit vector'
      IF(lwp) WRITE(numout,*) ' To check conservation : '
      IF(lwp) WRITE(numout,*) '   1 - No sea-ice model '
      IF(lwp) WRITE(numout,*) '   2 - No runoff ' 
      IF(lwp) WRITE(numout,*) '   3 - precipitation and evaporation equal to 1 : E=P=1 ' 
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      ! assign an index in trc array for each prognostic variable
      DO jn = 1,jp_canoe
       write(numout,*) ctrcnm(jn)
       cltra = ctrcnm(jn)
       IF( cltra == 'DIC'      )   jqdic = jn      !: dissolved inorganic carbon concentration
       IF( cltra == 'Alkalini' )   jqtal = jn      !: total alkalinity
       IF( cltra == 'O2'       )   jqoxy = jn      !: oxygen concentration
       IF( cltra == 'NO3'      )   jqno3 = jn      !: NO3 concentration
      END DO
      !
      IF( .NOT. ln_rsttr ) THEN
        trn(:,:,:,jqdic) = sco2
        trn(:,:,:,jqtal) = alka0 
        trn(:,:,:,jqoxy) = oxyg0
        trn(:,:,:,jqno3) = no30
      ENDIF
      !
      ! closea mask for BGCM
      CALL trc_closea_bgc(read_var_flag=.true.)
      !
      ! Test allocation of space for CanOE arrays before initialization
      CALL canoe_alloc ! allocate arrays space, see end of this module
      !
      ! O Riche Aug 4th 2022
      ! Initialise external sources reading
      ! Check namelist_top_cfg for &trcsrc_dta section
      ! call trc_src_init
      ! prep reading external sources
      ! open the files
      CALL trc_src_init
      !
      ! call all the BGC initialization subroutines in TOP tier
      CALL trc_flx_init
      !
      ! O Riche Aug 16th 2022
      CALL trc_opt_init
      !
      ! !
   END SUBROUTINE trc_ini_canoe

   SUBROUTINE canoe_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE pisces_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of CANOE 
      !!----------------------------------------------------------------------
      !
      USE lib_mpp, only: mpp_sum         ! ierr sum over all processors      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      !ierr =        sms_canoe_alloc()          ! Start of CANOE-related alloc routines...
      ierr =        trc_opt_alloc()
      ierr = ierr + sms_top_alloc()
      ierr = ierr + trc_che_alloc()
      ierr = ierr + trc_flx_alloc()
      !
      IF( lk_mpp    )   CALL mpp_sum( 'canoe_alloc', ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'canoe_alloc: unable to allocate CANOE arrays' )
      !
   END SUBROUTINE canoe_alloc

   !!======================================================================
END MODULE trcini_canoe
