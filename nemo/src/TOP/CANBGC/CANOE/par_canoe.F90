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

   ! first four are declared in par_trc
   !INTEGER, PUBLIC :: jrpoc  ! small sized POC
   !INTEGER, PUBLIC :: jrphy  ! small sized phyto C biomass
   !INTEGER, PUBLIC :: jrnch  ! small sized phyto chl-a
   !INTEGER, PUBLIC :: jrzoo  ! small sized zoo C biomass
   INTEGER, PUBLIC :: jrnh4  ! Ammonium Concentration
   INTEGER, PUBLIC :: jrfer  ! Dissolved Iron Concentration
   INTEGER, PUBLIC :: jrcal  ! calcite particules
   INTEGER, PUBLIC :: jrgoc  ! large sized POC
   INTEGER, PUBLIC :: jrnn   ! Nanophytoplankton N concentration 
   INTEGER, PUBLIC :: jrnfe  ! Nanophytoplankton Fe concentration 
   INTEGER, PUBLIC :: jrdia  ! large sized phyto C biomass
   INTEGER, PUBLIC :: jrdch  ! large sized phyto chl-a
   INTEGER, PUBLIC :: jrdn   ! Diatoms N concentration 
   INTEGER, PUBLIC :: jrdfe  ! Diatoms Fe concentration 
   INTEGER, PUBLIC :: jrmes  ! large sized zoo C biomass

   !!======================================================================

END MODULE par_canoe
