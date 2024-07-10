MODULE trcsms_csib
   !!======================================================================
   !!                         ***  MODULE trcsms_csib  ***
   !! TOP :   Main module of the CSIB tracers
   !!======================================================================
   !! History :      !  2007  (C. Ethe, G. Madec)  Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_sms_csib       : CSIB model main routine
   !! trc_sms_csib_alloc : allocate arrays specific to CSIB sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc

   USE ice              ! ice variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_csib       ! called by trcsms.F90 module
   PUBLIC   trc_sms_csib_alloc ! called by trcini_csib.F90 module

   ! Defined HERE the arrays specific to CSIB sms and ALLOCATE them in trc_sms_csib_alloc
   REAL(wp), PUBLIC, PARAMETER ::   epsi10 = 1.e-10_wp  !: small number

   !! Sea ice tracers are like ice model variables and have 2 equivalent variables: 
   !!    - one extenisve for dynamics
   !!    - one intensive for thermodynanics and biogeochemistry
   !!
   !! **********************************************************************|
   !! ***         Category dependent state variables (prognostic)        ***|
   !! **********************************************************************|
   !!                                                                       |
   !! ** Global variables                                                   |
   !!-------------|-------------|---------------------------------|---------|
   !! icedia_gca  |      -      |    Ice algae grid cell average  | mmol/m3 |
   !!                                                                       |
   !!-------------|-------------|---------------------------------|---------|
   !!                                                                       |
   !! ** Equivalent variables                                               |
   !!-------------|-------------|---------------------------------|---------|
   !! icedia      | -           |    Ice algae per ice area       | mmol/m3 |



   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: icedia            !  Ice algae per ice area
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: icedia_gca        !  Ice algae grid cell average
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: icediagca_2d      !  Ice algae grid cell average, 2d version for ice model

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsms_my_trc.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_sms_csib( kt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_csib  ***
      !!
      !! ** Purpose :   main routine of CSIB model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      INTEGER ::   jn   ! dummy loop index
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_sms_csib')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_csib:  CSIB model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'



      ! Conversion from global to equivalent variables
      WHERE( a_i(:,:,:) < epsi10 .OR. v_i(:,:,:) < epsi10 .OR. h_i(:,:,:) < epsi10 )
         icedia(:,:,:)=0._wp
      ELSEWHERE
         icedia(:,:,:) = icedia_gca(:,:,:) / a_i(:,:,:)
      END WHERE


      

      ! Conversion from equivalent to global variables
      icedia_gca(:,:,:) = icedia(:,:,:) * a_i(:,:,:)


      IF( ln_timing )   CALL timing_stop('trc_sms_csib')
      !
   END SUBROUTINE trc_sms_csib


   INTEGER FUNCTION trc_sms_csib_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_csib_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to CSIB
      ! ALLOCATE( tab(...) , STAT=trc_sms_csib_alloc )
      trc_sms_csib_alloc = 0      ! set to zero if no array to be allocated
      
      ALLOCATE(icedia(jpi,jpj,jpl), icedia_gca(jpi,jpj,jpl), icediagca_2d(jpij,jpl), &
         &     STAT=trc_sms_csib_alloc)

      IF( trc_sms_csib_alloc /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_csib_alloc : failed to allocate arrays' )
      !
   END FUNCTION trc_sms_csib_alloc

   !!======================================================================
END MODULE trcsms_csib
