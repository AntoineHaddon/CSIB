MODULE par_csib
   !!======================================================================
   !!                        ***  par_csib  ***
   !! TOP :   set the CSIB parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: par_my_trc.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

   USE par_kind   !: access wp kind

   IMPLICIT NONE
   PUBLIC
   
   REAL(wp), SAVE :: z_ia=0.03_wp               ! height of skeletal layer

   !!======================================================================
END MODULE par_csib