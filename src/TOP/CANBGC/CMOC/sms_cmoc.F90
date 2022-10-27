MODULE sms_cmoc   
   !!----------------------------------------------------------------------
   !!                     ***  sms_pisces.F90  ***  
   !! TOP :   PISCES Source Minus Sink variables
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-02 (O. Aumont) original code
   !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
   !!            CMOC1 !  2013-15 (O. Riche) add CMOC global variables and namelist parameters
   !!----------------------------------------------------------------------

   USE par_kind   !: access wp kind

   IMPLICIT NONE
   PUBLIC
   
   ! PP module
   REAL(wp), SAVE :: achl_cmoc  
   REAL(wp), SAVE :: thm_cmoc   
   REAL(wp), SAVE :: tau_cmoc   
   REAL(wp), SAVE :: itau_cmoc
   REAL(wp), SAVE :: ep_cmoc    
   REAL(wp), SAVE :: tvm_cmoc   
   REAL(wp), SAVE :: vm_cmoc    
   REAL(wp), SAVE :: kn_cmoc 

   !   Redfield ratio and euphotic zone
   REAL(wp), SAVE :: cnrr_cmoc  
   REAL(wp), SAVE :: ncrr_cmoc

   ! Mortality module
   REAL(wp), SAVE :: mpd_cmoc
   REAL(wp), SAVE :: mpd2_cmoc

   ! Remineralization module
   REAL(wp), SAVE :: ed_cmoc
   REAL(wp), SAVE :: reref_cmoc
 
END MODULE sms_cmoc