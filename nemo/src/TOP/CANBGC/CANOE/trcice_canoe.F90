MODULE trcice_canoe
   !!======================================================================
   !!                         ***  MODULE trcice_canoe  ***
   !!----------------------------------------------------------------------
   !! trc_ice_canoe       : CANOE model seaice coupling routine
   !!----------------------------------------------------------------------
   !! History :        !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ice_ini_canoe       ! called by trcice.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcice_canoe.F90 10069 2018-08-28 14:12:24Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ice_ini_canoe
      !!----------------------------------------------------------------------
      !!                     ***  trc_ice_canoe  ***
      !!
      !!----------------------------------------------------------------------
      !
      !
   END SUBROUTINE trc_ice_ini_canoe

   !!======================================================================
END MODULE trcice_canoe
