MODULE trcnam_cmoc
   !!======================================================================
   !!                      ***  MODULE trcnam_cmoc  ***
   !! TOP :   initialisation of some run parameters for CMOC bio-model
   !!======================================================================
   !! History :      !  2007  (C. Ethe, G. Madec) Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!                !  2022  (O. Riche) NEMO4 integration 
   !!----------------------------------------------------------------------
   !! trc_nam_cmoc      : CMOC model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables

   USE iom             ! IO manager
   USE sms_top_canbgc  ! shared arrays across BGCM code
   USE sms_cmoc, ONLY  : ws_cmoc, ln_cmocnegtr

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
      CHARACTER(LEN=20)::   clname
      !!----------------------------------------------------------------------
      INTEGER ::   ios       ! Local integer
      !!----------------------------------------------------------------------
      NAMELIST/namcmocws/ ws_cmoc
      NAMELIST/namcmocnegtr/ ln_cmocnegtr

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
      ! Reading particles sinking speed here as it is used in both cmocsink.F90 and trcsrc_canbgc.F90/trc_bott_cmoc subroutines
      REWIND( numnatp_refb )              ! Namelist namcmocws in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcmocws, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcmocws in reference namelist_cmoc' )
      REWIND( numnatp_cfgb )              ! Namelist namcmocws in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcmocws, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcmocws in configuration namelist_cmoc' )

      IF(lwm) WRITE( numonpb, namcmocws )    
      !
      ! Reading ln_cmocnegtr, if .false. the qnegtr block in trcsms_cmoc.F90 is skipped
      REWIND( numnatp_refb )              ! Namelist namcmocnegtr in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcmocnegtr, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcmocnegtr in reference namelist_cmoc' )
      REWIND( numnatp_cfgb )              ! Namelist namcmocnegtr in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcmocnegtr, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcmocnegtr in configuration namelist_cmoc' )
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_cmoc : namcmocws and namcmocnegtr'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) '    POC sinking speed                           ws_cmoc  =', ws_cmoc
      IF(lwp) WRITE(numout,*) '    xnegtrc block switch                    ln_cmocnegtr =', ln_cmocnegtr
      !
      IF(lwm) WRITE( numonpb, namcmocnegtr )    
      !
   END SUBROUTINE trc_nam_cmoc
   
   !!======================================================================
END MODULE trcnam_cmoc
