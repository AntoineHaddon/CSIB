MODULE par_fafmip
   !!======================================================================
   !!                        ***  par_fafmip  ***
   !! TOP :   set the fafmip parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: par_fafmip.F90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   USE par_lobster, ONLY : jp_lobster      !: number of tracers in LOBSTER
   USE par_lobster, ONLY : jp_lobster_2d   !: number of 2D diag in LOBSTER
   USE par_lobster, ONLY : jp_lobster_3d   !: number of 3D diag in LOBSTER
   USE par_lobster, ONLY : jp_lobster_trd  !: number of biological diag in LOBSTER

   USE par_pisces , ONLY : jp_pisces       !: number of tracers in PISCES
   USE par_pisces , ONLY : jp_pisces_2d    !: number of 2D diag in PISCES
   USE par_pisces , ONLY : jp_pisces_3d    !: number of 3D diag in PISCES
   USE par_pisces , ONLY : jp_pisces_trd   !: number of biological diag in PISCES

   USE par_cfc    , ONLY : jp_cfc          !: number of tracers in CFC
   USE par_cfc    , ONLY : jp_cfc_2d       !: number of tracers in CFC
   USE par_cfc    , ONLY : jp_cfc_3d       !: number of tracers in CFC
   USE par_cfc    , ONLY : jp_cfc_trd      !: number of tracers in CFC

   USE par_c14b   , ONLY : jp_c14b         !: number of tracers in C14
   USE par_c14b   , ONLY : jp_c14b_2d      !: number of tracers in C14
   USE par_c14b   , ONLY : jp_c14b_3d      !: number of tracers in C14
   USE par_c14b   , ONLY : jp_c14b_trd     !: number of tracers in C14

   IMPLICIT NONE

   INTEGER, PARAMETER ::   jp_ln      = jp_lobster     + jp_pisces     + jp_cfc     + jp_c14b     + jp_my_trc     !: 
   INTEGER, PARAMETER ::   jp_ln_2d   = jp_lobster_2d  + jp_pisces_2d  + jp_cfc_2d  + jp_c14b_2d  + jp_my_trc_2d  !:
   INTEGER, PARAMETER ::   jp_ln_3d   = jp_lobster_3d  + jp_pisces_3d  + jp_cfc_3d  + jp_c14b_3d  + jp_my_trc_3d  !:
   INTEGER, PARAMETER ::   jp_ln_trd  = jp_lobster_trd + jp_pisces_trd + jp_cfc_trd + jp_c14b_trd + jp_my_trc_trd !:

#if defined key_fafmip
   !!---------------------------------------------------------------------
   !!   'key_fafmip'                     user defined tracers (fafmip)
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_fafmip     = .TRUE.   !: PTS flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_fafmip     =  2       !: number of PTS tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_fafmip_2d  =  0       !: additional 2d output arrays ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_fafmip_3d  =  0       !: additional 3d output arrays ('key_trc_diaadd')
   INTEGER, PUBLIC, PARAMETER ::   jp_fafmip_trd =  0       !: number of sms trends for fafmip

   ! assign an index in trc arrays for each PTS prognostic variables
   INTEGER, PUBLIC, PARAMETER ::   jpTr  = jp_ln + 1     !: FAFMIP "redistributed heat" tracer
   INTEGER, PUBLIC, PARAMETER ::   jpTa  = jp_ln + 2     !: FAFMIP "added heat" tracer

#else
   !!---------------------------------------------------------------------
   !!   Default                           No user defined tracers (fafmip)
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_fafmip     = .FALSE.  !: fafmip flag 
   INTEGER, PUBLIC, PARAMETER ::   jp_fafmip     =  0       !: No fafmip tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_fafmip_2d  =  0       !: No fafmip additional 2d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_fafmip_3d  =  0       !: No fafmip additional 3d output arrays 
   INTEGER, PUBLIC, PARAMETER ::   jp_fafmip_trd =  0       !: number of sms trends for fafmip
#endif

   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC, PARAMETER ::   jp_faf0     = jp_ln     + 1              !: First index of fafmip passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_faf1     = jp_ln     + jp_fafmip      !: Last  index of fafmip passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_faf0_2d  = jp_ln_2d  + 1              !: First index of fafmip passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_faf1_2d  = jp_ln_2d  + jp_fafmip_2d   !: Last  index of fafmip passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_faf0_3d  = jp_ln_3d  + 1              !: First index of fafmip passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_faf1_3d  = jp_ln_3d  + jp_fafmip_3d   !: Last  index of fafmip passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_faf0_trd = jp_ln_trd + 1              !: First index of fafmip passive tracers
   INTEGER, PUBLIC, PARAMETER ::   jp_faf1_trd = jp_ln_trd + jp_fafmip_trd  !: Last  index of fafmip passive tracers

   !!======================================================================
END MODULE par_fafmip
