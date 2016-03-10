MODULE p4zopt
   !!======================================================================
   !!                         ***  MODULE p4zopt  ***
   !! TOP - PISCES : Compute the light availability in the water column
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.2  !  2009-04  (C. Ethe, G. Madec)  optimisation
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Improve light availability of nano & diat
   !!           CMOC1  !  2013-15  (O. Riche) change total chlorophyll (l.#84), bug fix euph. zone bottom (l.#172), NEMO list
   !!----------------------------------------------------------------------
#if defined  key_pisces
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No PISCES bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE p4z_opt                   ! Empty routine
   END SUBROUTINE p4z_opt
#endif 

   !!======================================================================
END MODULE  p4zopt
