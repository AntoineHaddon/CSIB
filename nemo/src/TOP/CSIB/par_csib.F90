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
   
   INTEGER, PUBLIC      ::    jp_csib =3                 ! number of ice tracers in CSIB model
   
   ! Indices of tracers
   INTEGER, PUBLIC      ::    jridia=1                   ! Ice diatoms
   INTEGER, PUBLIC      ::    jrino3=2                   ! Ice nitrate
   INTEGER, PUBLIC      ::    jrinh4=3                   ! Ice ammonium
   
   CHARACTER(len=lca), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   icetrcnm   !: tracer name 
   

   !!======================================================================
END MODULE par_csib