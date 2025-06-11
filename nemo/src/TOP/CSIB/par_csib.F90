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
   
   INTEGER, PUBLIC      ::    jp_csib =5                 ! number of ice tracers in CSIB model
   
   ! Indices of tracers
   INTEGER, PUBLIC      ::    jridiac=1                  ! Ice diatoms C biomass
   INTEGER, PUBLIC      ::    jridian=2                  ! Ice diatoms N biomass
   INTEGER, PUBLIC      ::    jridiach=3                 ! Ice diatoms Chl biomass
   INTEGER, PUBLIC      ::    jrino3=4                   ! Ice nitrate
   INTEGER, PUBLIC      ::    jrinh4=5                   ! Ice ammonium
   
   CHARACTER(len=lca), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   icetrcnm   !: tracer names (used for restart input and output. can't init here because if restart then ice model executed first and needs to read restart files. So init is in function trcini_csib.f90/trc_ini_csibnames.) 
   

   !!======================================================================
END MODULE par_csib