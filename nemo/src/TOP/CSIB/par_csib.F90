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
   
   LOGICAL , PUBLIC, SAVE ::   ln_csib              ! CSIB flag 
   LOGICAL , PUBLIC, SAVE ::   ln_ibgcspinup        ! flag to start ice BGC variables from ocean surface values, for spin-up in case of restart from a run without ice BGC


   INTEGER, PUBLIC      ::    jp_csib =5                 ! number of ice tracers in CSIB model
   
   ! Indices of tracers
   INTEGER, PUBLIC      ::    jridiac=1                  ! Ice diatoms C biomass
   INTEGER, PUBLIC      ::    jridian=2                  ! Ice diatoms N biomass
   INTEGER, PUBLIC      ::    jridiach=3                 ! Ice diatoms Chl biomass
   INTEGER, PUBLIC      ::    jrino3=4                   ! Ice nitrate
   INTEGER, PUBLIC      ::    jrinh4=5                   ! Ice ammonium
   
   ! tracer names (used for restart input and output). 
   CHARACTER(len=lca), PUBLIC, SAVE, DIMENSION(5) ::   icetrcnm  = (/ 'icediac','icedian','icediach','iceno3','icenh4'  /)
   

   !!======================================================================
END MODULE par_csib