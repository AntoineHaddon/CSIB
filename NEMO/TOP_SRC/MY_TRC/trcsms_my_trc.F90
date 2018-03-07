MODULE trcsms_my_trc
   !!======================================================================
   !!                         ***  MODULE trcsms_my_trc  ***
   !! TOP :   Main module of the MY_TRC tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_my_trc
   !!----------------------------------------------------------------------
   !!   'key_my_trc'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_sms_my_trc       : MY_TRC model main routine
   !! trc_sms_my_trc_alloc : allocate arrays specific to MY_TRC sms
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trdmod_oce
   USE trdmod_trc
   USE phycst
   USE fldread         ! read input fields

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_my_trc       ! called by trcsms.F90 module
   PUBLIC   trc_sms_my_trc_alloc ! called by trcini_my_trc.F90 module

   INTEGER , PARAMETER ::   jpfld = 2   ! maximum number of files to read
   INTEGER , PARAMETER ::   jpo1  = 1   ! index of tracer o1
   INTEGER , PARAMETER ::   jpo2  = 2   ! index of tracer o2
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf    ! structure of input fields (file informations, fields read)

   !! * Substitution
#  include "domzgr_substitute.h90"

   ! Defined HERE the arrays specific to MY_TRC sms and ALLOCATE them in trc_sms_my_trc_alloc

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_my_trc.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_my_trc( kt )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_my_trc  ***
      !!
      !! ** Purpose :   main routine of MY_TRC model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER ::   jn   ! dummy loop index
      REAL(wp) :: zfact = 1._wp             ! temporary scalar
      REAL(wp) :: dtyrs
      CHARACTER(len=100) ::  cn_dir             ! Root directory for location of tracer files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i  ! array of namelist information structures
      TYPE(FLD_N) :: sn_o1, sn_o2               ! informations about the fields to be read
      NAMELIST/namsbc_olo2/ cn_dir, ln_o1, ln_o2, sn_o1, sn_o2
!!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_my_trc')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_my_trc:  MY_TRC model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
! IF PISCES is used, time-stepping is Euler, so use a factor of 2, relative
! to leapfrog stepping when PISCES is not used. (stupid, but beyond CCCma control)
#if defined key_pisces 
      dtyrs = 2.0_wp / (3600._wp * 24. * 365.) ! fraction of a year per time step 
#else
      dtyrs = 1.0_wp / (3600._wp * 24. * 365.) ! fraction of a year per time step 
#endif
      tra(:,:,:,jpage) = tra(:,:,:, jpage) + dtyrs ! Add the time to the tendancy.
      tra(:,:,1,jpage) = 0._wp  ! Hard restoring to counter E-P & river dilution, equivalent to relaxation time=0
      trn(:,:,1,jpage) = 0._wp  ! Hard restoring to counter E-P & river dilution
      ! WRITE(numout,*) 'Max surface, ocean age', maxval(trn(:,:,1,jpage)), maxval(trn(:,:,:,jpage))
      ! WRITE(numout,*) 'Max surface, ocean age tra:', maxval(tra(:,:,1,jpage)), maxval(tra(:,:,:,jpage))

! Special tracers o1 & o2
      IF( kt == nit000 ) THEN                  ! First call kt=nit000
      ! set file information
         cn_dir = './'                         ! directory in which the model is executed
         ! ... default values (NB: frequency positive => hours, negative => months)
         !            !  file  ! frequency !  variable      ! time intep !  clim   ! 'yearly' or ! weights  ! rotation  !
         !            !  name  ! (hours)   !   name         !   (T/F)    !  (T/F)  ! 'monthly'   ! filename ! pairs     !
         sn_o1 = FLD_N(  'o1' ,    -1     , 'heat_flux_anom',   .false.  , .true. ,   'yearly'  ,  ''       , ''        )
         sn_o2 = FLD_N(  'o2' ,    -1     , 'heat_flux_anom',   .false.  , .true. ,   'yearly'  ,  ''       , ''        )
         !
         REWIND ( numnam )                     ! read in namlist namsbc_o1o2
         READ   ( numnam, namsbc_o1o2 )           
         !                                     ! Control print
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,'trc_sms_my_trc : o1 o2 tracers'
            WRITE(numout,*) '~~~~~~~ '
            WRITE(numout,*) '   Namelist namsbc_o1o2'
            WRITE(numout,*) 'o1 flux in a file to be read  ln_o1 = ', ln_o1
            WRITE(numout,*) 'o2 flux in a file to be read  ln_o2 = ', ln_o2
         ENDIF
         !
         !                                     ! store namelist information in an array
         slf_i(jpo1) = sn_o1   ;   slf_i(jpo2) = sn_o2
         !
         ALLOCATE( sf(jpfld), STAT=ierror )    ! set sf structure
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'trc_sms_my_trc: unable to allocate sf structure' )   ;   RETURN
         ENDIF
         DO ji= 1, jpfld
            ALLOCATE( sf(ji)%fnow(jpi,jpj,1) )
            IF( slf_i(ji)%ln_tint ) ALLOCATE( sf(ji)%fdta(jpi,jpj,1,2) )
         END DO
         !
         CALL fld_fill( sf, slf_i, cn_dir, 'trc_sms_my_trc', 'o1 o2 tracers', 'namsbc_o1o2' )
         !
      ENDIF
      IF( ln_o1 .or. ln_o2 ) THEN
         CALL fld_read( kt, nn_fsbc, sf )      ! input fields provided at the current time-step
      !
         IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN  ! update o1 o2 fluxes at each SBC frequency
            tra(:,:,1,jpo1) = zfact * sf(jpo1)%fnow(:,:,1) * ro0cpr / fse3t(:,:,1)
            tra(:,:,1,jpo2) = zfact * sf(jpo2)%fnow(:,:,1) * ro0cpr / fse3t(:,:,1)
         ENDIF
      ELSE
      ! Apply uniform fluxes of 1 W/m2 and 2 W/m2 to tracers.
         tra(:,:,1,jpo1) = 1._wp * ro0cpr / fse3t(:,:,1)
         tra(:,:,1,jpo2) = 2._wp * ro0cpr / fse3t(:,:,1) 
      ENDIF
  
      ! To get the heat flux anomaly field from netcdf, do something like for runoff:
      ! CALL fld_read ( kt, nn_fsbc, sf_rnf   )    ! Read Runoffs data and provide it at kt
      ! rnf(:,:) =  sf_rnf(1)%fnow(:,:,1) 

      ! To increment temperature with the heat flux anomaly, this is what we want to do:
      ! zfact = 0.5e0
      ! sbc_tsc_b(:,:,:) = sbc_tsc(:,:,:)
      ! sbc_tsc(ji,jj,jp_tem) = ro0cpr * qns(ji,jj) ! where qns is our heat flux anomaly
      ! z1_e3t = zfact / fse3t(ji,jj,1)
      ! tsa(ji,jj,1,jn) = tsa(ji,jj,1,jn) + ( sbc_tsc_b(ji,jj,jn) + sbc_tsc(ji,jj,jn) ) * z1_e3t

      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_my_trc')
      !
   END SUBROUTINE trc_sms_my_trc


   INTEGER FUNCTION trc_sms_my_trc_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_my_trc_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to MY_TRC
      ! ALLOCATE( tab(...) , STAT=trc_sms_my_trc_alloc )
      trc_sms_my_trc_alloc = 0      ! set to zero if no array to be allocated
      !
      IF( trc_sms_my_trc_alloc /= 0 ) CALL ctl_warn('trc_sms_my_trc_alloc : failed to allocate arrays')
      !
   END FUNCTION trc_sms_my_trc_alloc


#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No MY_TRC model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_my_trc( kt )             ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_my_trc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_my_trc
#endif

   !!======================================================================
END MODULE trcsms_my_trc
