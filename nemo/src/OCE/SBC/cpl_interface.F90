!> Contains the interfaces to the CanCPL or OASIS couplers
MODULE CPL_INTERFACE

   use sbc_oce,  only : lk_oasis, lk_cancpl
   use par_kind, only : wp
   use lib_mpp,  only : ctl_stop

   use cpl_oasis3, only : cpl_oasis3_init  => cpl_init
   use cpl_oasis3, only : cpl_define
   use cpl_oasis3, only : cpl_snd
   use cpl_oasis3, only : cpl_rcv
   use cpl_oasis3, only : cpl_freq
   use cpl_oasis3, only : cpl_finalize
   use cpl_oasis3, only : OASIS_idle, OASIS_Rcv

   use cpl_cancpl, only : cpl_cancpl_init
   use cpl_cancpl, only : cpl_cancpl_define
   use cpl_cancpl, only : cpl_cancpl_snd
   use cpl_cancpl, only : cpl_cancpl_rcv
   use cpl_cancpl, only : cpl_cancpl_freq
   use cpl_cancpl, only : cpl_cancpl_finalize

   use cpl_types, only : COUPLER_idle, COUPLER_Rcv

   implicit none; private

   !> Abstract interfaces to allow for either the OASIS or CanCPL routines to be used
   ABSTRACT INTERFACE
      SUBROUTINE cpl_define_sub( krcv, ksnd, kcplmodel )
         INTEGER, INTENT(in) :: krcv ! Number of coupling fields to receive
         INTEGER, INTENT(in) :: ksnd ! Number of coupling fields to send
         INTEGER, INTENT(in) :: kcplmodel      ! Maximum number of models to/from which NEMO is potentialy sending/receiving data
      END SUBROUTINE
   END INTERFACE

   ABSTRACT INTERFACE
      SUBROUTINE cpl_rcv_sub( kid, kstep, pdata, pmask, kinfo )
         import wp
         INTEGER                   , INTENT(in   ) ::   kid       ! variable index in the array
         INTEGER                   , INTENT(in   ) ::   kstep     ! ocean time-step in seconds
         REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pdata     ! IN to keep the value if nothing is done
         REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pmask     ! coupling mask
         INTEGER                   , INTENT(  out) ::   kinfo     ! OASIS3 info argument
       END SUBROUTINE
   END INTERFACE

   ABSTRACT INTERFACE
      SUBROUTINE cpl_snd_sub( kid, kstep, pdata, kinfo )
         import wp
         INTEGER                   , INTENT(in   ) ::   kid       ! variable index in the array
         INTEGER                   , INTENT(  out) ::   kinfo     ! OASIS3 info argument
         INTEGER                   , INTENT(in   ) ::   kstep     ! ocean time-step in seconds
         REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdata
      END SUBROUTINE cpl_snd_sub
   END INTERFACE

   ABSTRACT INTERFACE
      INTEGER FUNCTION cpl_freq_sub( cdfieldname )
         CHARACTER(len = *), INTENT(in) ::   cdfieldname    ! field name as set in namcouple file
      END FUNCTION cpl_freq_sub
   END INTERFACE

   ABSTRACT INTERFACE
      SUBROUTINE cpl_finalize_sub( )
      END SUBROUTINE cpl_finalize_sub
   END INTERFACE

   PUBLIC   cpl_init
   PROCEDURE(cpl_snd_sub),      POINTER, PUBLIC, SAVE :: cpl_snd => NULL()
   PROCEDURE(cpl_rcv_sub),      POINTER, PUBLIC, SAVE :: cpl_rcv => NULL()
   PROCEDURE(cpl_freq_sub),     POINTER, PUBLIC, SAVE :: cpl_freq => NULL()
   PROCEDURE(cpl_define_sub),   POINTER, PUBLIC, SAVE :: cpl_define => NULL()
   PROCEDURE(cpl_finalize_sub), POINTER, PUBLIC, SAVE :: cpl_finalize => NULL()

   CONTAINS

   !> Set the correct procedure pointers for the requested coupler
   SUBROUTINE cpl_init( cd_modname, kl_comm )
      CHARACTER(len = *), INTENT(in   ) ::   cd_modname   ! model name as set in namcouple file
      INTEGER           , INTENT(  out) ::   kl_comm      ! local communicator of the model

      if (lk_cancpl) then
        call cpl_cancpl_init( kl_comm )
        cpl_snd => cpl_cancpl_snd
        cpl_rcv => cpl_cancpl_rcv
        cpl_freq => cpl_cancpl_freq
        cpl_define => cpl_cancpl_define
        cpl_finalize => cpl_cancpl_finalize
        COUPLER_idle = OASIS_idle
        COUPLER_Rcv  = OASIS_Rcv
      elseif (lk_oasis) then
        call cpl_oasis3_init( cd_modname, kl_comm )
        COUPLER_idle = OASIS_idle
        COUPLER_Rcv  = OASIS_Rcv
      else
        CALL ctl_stop( 'STOP', 'cpl_interface_cpl_init: No coupler selected, set lk_cancpl or lk_oasis3' )
      endif

   end subroutine cpl_init


END MODULE CPL_INTERFACE
