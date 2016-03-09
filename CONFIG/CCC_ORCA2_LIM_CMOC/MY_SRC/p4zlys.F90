  MODULE p4zlys
   !!======================================================================
   !!                         ***  MODULE p4zlys  ***
   !! TOP :   PISCES 
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr)  Calcon salinity dependence
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Improvment of calcite dissolution
   !!            CMOC1 !  2013-15  (O. Riche) no more PISCES CaCO3 dissolution diagnostics; keep CO3 and pH diagnostics
   !!            CMOC1 !  2013-15  (N. Swart)) Removed, surface ph diagnostic moved to pz4flx
   !!           
   !!----------------------------------------------------------------------
#if defined key_pisces
#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_lys( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p4z_lys: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_lys
#endif 
   !!======================================================================
END MODULE  p4zlys
