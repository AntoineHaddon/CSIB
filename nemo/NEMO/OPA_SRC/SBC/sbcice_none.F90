MODULE sbcice_none
   !!======================================================================
   !!                       ***  MODULE  sbcice  ***
   !! Surface module :  Keep SST above freezing point
   !!======================================================================
   !! History :   3.4.1! 2013-08  (D. Yang)  Original code 
   !!
   !!----------------------------------------------------------------------
   !!   sbc_ice_none   : Prevent SST from dropping below freezing point
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE eosbn2          ! equation of state
   USE sbc_oce         ! surface boundary condition: ocean fields
   USE sbccpl
   USE fldread         ! read input field
   USE iom             ! I/O manager library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ice_none    ! routine called in sbcmod

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: sbcice_none.F90 ???? 2013-08-19 19:15:05Z D.YANG $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_ice_none( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_ice_none  ***
      !!
      !! ** Purpose :   Keep SST above freezing point when sea-ice model are not used
      !!
      !! ** Method  : - blah blah blah, ...
      !!
      !! ** Action  :   utau, vtau : remain unchanged
      !!                taum, wndm : remain unchanged
      !!                qns, qsr   : update heat flux below sea-ice
      !!                emp, emps  : update freshwater flux below sea-ice
      !!                fr_i       : update the ice fraction
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      !INTEGER  ::   ierror     ! return error code
      REAL(wp) ::   Tocnfrz
      !!
      !!---------------------------------------------------------------------
      !IF( MOD( kt-1, nn_fsbc) == 0 ) THEN
         !
         Tocnfrz=-1.96_wp
      !  WRITE(numout,*) 'Tocnfrz',Tocnfrz
         CALL FLUSH(numout)

         DO jj = 1, jpj
            DO ji = 1, jpi
               !
               tsn(ji,jj,1,jp_tem) = MAX( tsn(ji,jj,1,jp_tem), Tocnfrz )     ! avoid over-freezing point temperature
            END DO
         END DO
         CALL lbc_lnk( tsn (:,:,1,jp_tem), 'T', -1. )
         !
      !ENDIF
      !
   END SUBROUTINE sbc_ice_none

   !!======================================================================
END MODULE sbcice_none
