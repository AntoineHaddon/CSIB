MODULE par_canoe
   !!======================================================================
   !!                        ***  par_canoe  ***
   !! TOP :   set the CANOE parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!                  !  2022-23  (J. Christian, O. Riche) NEMO4 integration
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: par_canoe.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

   IMPLICIT NONE

   ! INTEGER, PUBLIC :: jrnch !    
   ! INTEGER, PUBLIC :: jrdch ! 
   ! INTEGER, PUBLIC :: jpphy   ! small sized phyto C biomass
   ! INTEGER, PUBLIC :: jpnch   ! small sized phyto chl-a
   ! INTEGER, PUBLIC :: jpzoo   ! small sized zoo C biomass
   ! INTEGER, PUBLIC :: jppoc   ! small sized POC
   INTEGER, PUBLIC :: jpdia   ! large sized phyto C by
   INTEGER, PUBLIC :: jpdch   ! large sized phyto chl-a
   INTEGER, PUBLIC :: jpmes   ! large sized zoo C biomass
   INTEGER, PUBLIC :: jpgoc   ! large sized POC
   INTEGER, PUBLIC :: jpcal   ! calcite particules

   !!======================================================================

END MODULE par_canoe
