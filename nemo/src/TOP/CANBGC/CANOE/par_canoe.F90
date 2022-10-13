MODULE par_canoe
   !!======================================================================
   !!                        ***  par_canoe  ***
   !! TOP :   set the CANOE parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: par_canoe.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

   IMPLICIT NONE

   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   ! INTEGER, PUBLIC ::   jp_knu0             !: First index of CANOE passive tracers
   ! INTEGER, PUBLIC ::   jp_knu1             !: Last  index of CANOE passive tracers
   ! ! INTEGER, PUBLIC ::   jqdic_can           !: DIC
   ! INTEGER, PUBLIC ::   jpalk_can           !: TA
   ! INTEGER, PUBLIC ::   jqoxy_can           !: O2 
   !!======================================================================
END MODULE par_canoe