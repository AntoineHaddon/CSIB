MODULE trcnam_canoe
   !!======================================================================
   !!                      ***  MODULE trcnam_canoe  ***
   !! TOP :   initialisation of some run parameters for CANOE bio-model
   !!======================================================================
   !! History :      !  2007  (C. Ethe, G. Madec) Original code
   !!                !  2016  (C. Ethe, T. Lovato) Revised architecture
   !!----------------------------------------------------------------------
   !! trc_nam_canoe      : CANOE model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   !
   USE iom             ! IO manager
   !
   USE sms_top_canbgc  ! shared arrays across BGCM code
   USE sms_canoe
   !
   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_canoe   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcnam_canoe.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_canoe
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_canoe  ***  
      !!
      !! ** Purpose :   read CANOE namelist
      !!
      !!----------------------------------------------------------------------
      !
      CHARACTER(LEN=20)::   clname
      !!----------------------------------------------------------------------
      INTEGER ::   ios       ! Local integer
      !!----------------------------------------------------------------------
      NAMELIST/namcanbio/ ws_canoe, ws_canoe2, ws_canoec
      NAMELIST/namcanoenegtr/ ln_canoenegtr
      
      IF(lwp) WRITE(numout,*)
      clname = 'namelist_canoe'

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_canoe : read CANOE namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~'
      !
      CALL ctl_opn( numnatp_refb, TRIM( clname )//'_ref', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatp_cfgb, TRIM( clname )//'_cfg', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) CALL ctl_opn( numonpb     , 'output.namelist.can' , 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      !
      !
      ! Reading particles sinking speed here as it is used in both canoesink.F90
      REWIND( numnatp_refb )              ! Namelist namcanbio in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcanbio, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcanbio in reference namelist_canoe' )
      REWIND( numnatp_cfgb )              ! Namelist namcanbio in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcanbio, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcanbio in configuration namelist_canoe' )

      IF(lwm) WRITE( numonpb, namcanbio )    
      !
      ! Reading ln_canoenegtr, if .false. the qnegtr block in trcsms_canoe.F90 is skipped
      REWIND( numnatp_refb )              ! Namelist namcanoenegtr in reference namelist : Passive tracer variables
      READ  ( numnatp_refb, namcanoenegtr, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcanoenegtr in reference namelist_canoe' )
      REWIND( numnatp_cfgb )              ! Namelist namcanoenegtr in configuration namelist : Passive tracer variables
      READ  ( numnatp_cfgb, namcanoenegtr, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcanoenegtr in configuration namelist_canoe' )

      IF(lwm) WRITE( numonpb, namcanoenegtr )    
      !
      END SUBROUTINE trc_nam_canoe
      !
   !!======================================================================
END MODULE trcnam_canoe
