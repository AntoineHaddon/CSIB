MODULE cmocsink
   !!!!!! O Riche Nov 8th 2022
   !!!!!! move this into CANBGC as trcsink_canbgc.F90
   !!======================================================================
   !!                         ***  MODULE cmocsink  ***
   !! TOP :  CANBGC  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!           CMOC1  !  2013-15  (O. Riche) POC sinking, Oct 28th 2015 simplified sink scheme according to Jim's own modifications for CMOC2
   !!           CMOC1  !  2016-02  (N. Swart) Adds calcite sinking.
   !!           CMOC2  !  2022-23  NEMO4 implementation
   !!----------------------------------------------------------------------
   !!   cmoc_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   cmoc_sink_init  :  Unitialisation of sinking speed parameters
   !!   cmoc_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   USE in_out_manager  ! I/O manager
   USE dom_oce         ! ocean space and time domain 
   USE timing          ! Timing
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   USE sms_top_canbgc
   USE sms_cmoc
   USE trc_closea_canbgc ! tmask_bgc_closea


   IMPLICIT NONE
   PRIVATE

   PUBLIC   cmoc_sink         ! called in trcsms_cmoc.F90
   PUBLIC   cmoc_sink_init    ! called in trcsms_cmoc.F90
   ! PUBLIC   cmoc_sink_alloc

   ! REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed 
   ! REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking  !: POC sinking fluxes

   ! INTEGER  :: iksed  = 10

   !!* Substitution
!#  include "top_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: cmocsink.F90 3685 2012-11-27 15:39:02Z cetlod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE cmoc_sink_init  ***
      !!----------------------------------------------------------------------
 
  SUBROUTINE cmoc_sink( kt , jnt )

      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - NEED TO BE DETAILED FOR ONCE
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, jnt  
  
  END SUBROUTINE cmoc_sink


  SUBROUTINE cmoc_sink_init
      !
      INTEGER :: ji, jj, ikt, ios     !: working integers for loops and I/O
      !
      NAMELIST/namcmoccal/ rmcico_cmoc, trcico_cmoc, aci_cmoc, dci_cmoc
      !
      REWIND( numnatp_refb )              ! Namelist namcmoccal in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcmoccal, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcmoccal in reference namelist_cmoc' )
      REWIND( numnatp_cfgb )              ! Namelist namcmoccal in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcmoccal, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcmoccal in configuration namelist_cmoc' )

      IF(lwm) WRITE( numonpb, namcmoccal )

      !
      IF( lwp ) THEN
        WRITE(numout,*) ' Namelist parameters for calcite export  , namcmoccal'
        WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        WRITE(numout,*) '    Maximum rain ratio                      rmcico_cmoc =', rmcico_cmoc
        WRITE(numout,*) '    Rain ratio half-point temperature       trcico_cmoc =', trcico_cmoc
        WRITE(numout,*) '    Rain ratio scaling factor                  aci_cmoc =',    aci_cmoc
        WRITE(numout,*) '    CaCO3 redissolution depth scale            dci_cmoc =',    dci_cmoc
        WRITE(numout,*) ' '      
      END IF
      !
      ! Define an open ocean mask, based on where mbkt > nk_bal_cmoc == 28 (or as define in namelist_cmoc_ref)
      ! Move to here instead of p4z_sink2 in the original CanESM5/CMOC code, since this only needs to be defined
      ! once
      !
      ALLOCATE( oomask(jpi, jpj) )
      oomask(:,:) = 0.0_wp
      !
      DO jj = 1, jpj
         DO ji = 1,jpi
            ikt = mbkt(ji,jj)
            IF ( ikt > nk_bal_cmoc ) THEN
               oomask(ji, jj) = 1.0_wp
            ENDIF
         ENDDO                                                               
      ENDDO

      !
  END SUBROUTINE cmoc_sink_init
  

END MODULE cmocsink