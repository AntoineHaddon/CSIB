MODULE sms_canoe   
   !!----------------------------------------------------------------------
   !!                     ***  sms_pisces.F90  ***  
   !! TOP :   PISCES Source Minus Sink variables
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-02 (O. Aumont) original code
   !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
   !!            CANOE2!  2022-23 (O. Riche & J Christian ) add CANOE global variables and namelist parameters
   !!----------------------------------------------------------------------

   USE par_kind   !: access wp kind

   IMPLICIT NONE
   PUBLIC
   
   ! Particles sinking speed
   REAL(wp), SAVE :: ws_canoe
   REAL(wp), SAVE :: wsbio, wsbio2, wsbioc

   ! qnegtr block skipping switch
   LOGICAL, SAVE :: ln_canoenegtr

    
   
   
END MODULE sms_canoe