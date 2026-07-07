MODULE par_dms
   !!======================================================================
   !!                        ***  par_dms  ***
   !! TOP :   set the DMS parameters
   !!======================================================================
   !! History :      ! 2026 (T. Sou, A. Haddon) Ocean DMS 
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

   USE par_kind   !: access wp kind
   IMPLICIT NONE

   ! Indices of tracers
   INTEGER, PUBLIC        ::   jrdmspd                    ! ocean DMSP
   INTEGER, PUBLIC        ::   jrdms                      ! ocean DMS

END MODULE par_dms
