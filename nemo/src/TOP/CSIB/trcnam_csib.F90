MODULE trcnam_csib
   !!======================================================================
   !!                      ***  MODULE trcnam_csib  ***
   !! TOP :   initialisation of some run parameters for CSIB bio-model
   !!======================================================================
   !! History :      !  2007  (C. Ethe, G. Madec) Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_nam_csib      : CSIB model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_csib   ! called by trcini_csib

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcnam_my_trc.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_csib
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_csib  ***  
      !!
      !! ** Purpose :   read CSIB namelist
      !!
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_csib : read CSIB namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~'
      !
   END SUBROUTINE trc_nam_csib
   
   !!======================================================================
END MODULE trcnam_csib
