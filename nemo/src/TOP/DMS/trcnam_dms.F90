MODULE trcnam_dms
   !!======================================================================
   !!                      ***  MODULE trcnam_dms  ***
   !! TOP :   initialisation of some run parameters for DMS bio-model
   !!======================================================================
   !! History :      ! 2026 (T. Sou, A. Haddon) Ocean DMS 
   !!----------------------------------------------------------------------
   !! trc_nam_dms      : DMS model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables

   USE sms_top_canbgc  ! canoe namelist reference

   USE par_dms
   USE trcsms_dms


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_dms   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcnam_dms.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_dms
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_dms  ***  
      !!
      !! ** Purpose :   read DMS namelist
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER ::   ios       ! Local integer
   

      NAMELIST/namdmsoce/ q_p1,q_p2,f_z1,f_z2,f_e1,f_e2,f_yield,k_dmspd,k_dms,k_free,k_photo                    !t ocedms

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_dms : read ocean DMS namelist'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~'


      ! read namelists
      REWIND( numnatp_refb )              
      READ  ( numnatp_refb, namdmsoce, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdmsoce in reference namelist_canoe' )
      REWIND( numnatp_cfgb )              
      READ  ( numnatp_cfgb, namdmsoce, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdmsoce in configuration namelist_canoe' )

      IF(lwm) WRITE( numonpb, namdmsoce )  ! output parameters in output.namelist.can 


      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for ocean DMS, namdmsoce'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) ' q_p1    =',q_p1
         WRITE(numout,*) ' q_p2    =',q_p2
         WRITE(numout,*) ' f_z1    =',f_z1
         WRITE(numout,*) ' f_z2    =',f_z2
         WRITE(numout,*) ' f_e1    =',f_e1
         WRITE(numout,*) ' f_e2    =',f_e2
         WRITE(numout,*) ' f_yield =',f_yield
         WRITE(numout,*) ' k_dmspd =',k_dmspd
         WRITE(numout,*) ' k_dms   =',k_dms
         WRITE(numout,*) ' k_free  =',k_free
         WRITE(numout,*) ' k_photo =',k_photo
         
      ENDIF

      ! conversions of namelist constants from 1/day to 1/sec
      k_dmspd = k_dmspd / 86400._wp
      k_dms   = k_dms   / 86400._wp
      k_free  = k_free  / 86400._wp
      k_photo = k_photo / 86400._wp



!
   END SUBROUTINE trc_nam_dms
   
   !!======================================================================
END MODULE trcnam_dms
