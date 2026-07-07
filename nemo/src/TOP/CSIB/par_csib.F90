MODULE par_csib
   !!======================================================================
   !!                        ***  par_csib  ***
   !! TOP :   set the CSIB parameters
   !!======================================================================
   !! History :   !  2025  (A. Haddon)  Original code
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: par_my_trc.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

   USE par_kind   !: access wp kind

   IMPLICIT NONE
   
   LOGICAL , PUBLIC, SAVE ::   ln_csib=.false.        ! CSIB flag 
   LOGICAL , PUBLIC, SAVE ::   ln_dmsice=.false.      ! sea ice DMS flag
   LOGICAL , PUBLIC, SAVE ::   ln_ibgcspinup=.false.  ! flag to start ice BGC variables from ocean surface values, for spin-up in case of restart from a run without ice BGC

   INTEGER, PUBLIC      ::    jp_csib                    ! number of ice tracers in CSIB model
   
   ! Indices of tracers
   INTEGER, PUBLIC      ::    jridiac=1                  ! Ice algae C biomass
   INTEGER, PUBLIC      ::    jridian=2                  ! Ice algae N biomass
   INTEGER, PUBLIC      ::    jridiach=3                 ! Ice algae Chl biomass
   INTEGER, PUBLIC      ::    jrino3=4                   ! Ice nitrate
   INTEGER, PUBLIC      ::    jrinh4=5                   ! Ice ammonium
   INTEGER, PUBLIC      ::    jridmspd=6                 ! Ice DMSpd
   INTEGER, PUBLIC      ::    jridms=7                   ! Ice DMS

   ! tracer names (used for restart input and output). 
   CHARACTER(len=lca), PUBLIC, SAVE, DIMENSION(7) :: icetrcnm = (/ 'icediac','icedian','icediach','iceno3','icenh4','icedmspd','icedms'  /)

   

   !!======================================================================
END MODULE par_csib