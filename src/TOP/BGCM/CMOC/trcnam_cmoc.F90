MODULE trcnam_cmoc
   !!======================================================================
   !!                      ***  MODULE trcnam_cmoc  ***
   !! TOP :   initialisation of some run parameters for CMOC bio-model
   !!======================================================================
   !! History :      !  2007  (C. Ethe, G. Madec) Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_nam_cmoc      : CMOC model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables

   USE iom             ! IO manager
   USE sms_top         ! shared arrays across BGCM code

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_cmoc   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcnam_cmoc.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_cmoc
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_cmoc  ***  
      !!
      !! ** Purpose :   read CMOC namelist
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER :: jl, jn
      INTEGER :: ios, ioptio         ! Local integer
      CHARACTER(LEN=20)::   clname
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      clname = 'namelist_cmoc'

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_cmoc : read CMOC namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~'
      !
      CALL ctl_opn( numnatp_refb, TRIM( clname )//'_ref', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatp_cfgb, TRIM( clname )//'_cfg', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) CALL ctl_opn( numonpb     , 'output.namelist.cmoc' , 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      !
      !
   END SUBROUTINE trc_nam_cmoc
   
   !!======================================================================
END MODULE trcnam_cmoc
