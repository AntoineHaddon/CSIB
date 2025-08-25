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

   USE trcsms_csib

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
      CHARACTER(LEN=20)::   clname
      INTEGER ::   ios       ! Local integer
      INTEGER ::   numnatp_refcsib = -1           !! Logical unit for the ref namelist for the parameters of the CSIB model
      INTEGER ::   numnatp_cfgcsib = -1           !! Logical unit for the cfg namelist for the parameters of the CSIB model
      INTEGER ::   numonpbcsib      = -1          !! Logical unit for the above ref/cfg namelists output

      !!----------------------------------------------------------------------
      NAMELIST/namicedia/ z_ia, qnidiamin, qnidiamax, qchidiaref, alrefidia, pcrefidia, betaidia, vnref, knh4, kno3, etares, ch2nmax, min_icedia, t_ia, r_m1, r_m2, f_p2, f_flsh, f_slgh, dt_mo, t_mo, d_mo
      NAMELIST/namicenit/ f_rm, r_ni, c_di, c_nu
      NAMELIST/namicedic/ sicpump, icedicref, icetalref, f_dicsw, f_dicsw_melt
      
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_csib : read CSIB namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~'
     
      ! open namelists
      clname = 'namelist_csib'
      CALL ctl_opn( numnatp_refcsib, TRIM( clname )//'_ref', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatp_cfgcsib, TRIM( clname )//'_cfg', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      
      ! read ref namelists

      ! Namelist namicedia in reference namelist
      REWIND( numnatp_refcsib )              
      READ  ( numnatp_refcsib, namicedia, IOSTAT = ios, ERR = 901)
      901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namicedia in reference namelist_csib' )
      ! Namelist namicenit in reference namelist
      REWIND( numnatp_refcsib )              
      READ  ( numnatp_refcsib, namicenit, IOSTAT = ios, ERR = 902)
      902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namicenit in reference namelist_csib' )
      ! Namelist namicedic in reference namelist
      REWIND( numnatp_refcsib )              
      READ  ( numnatp_refcsib, namicedic, IOSTAT = ios, ERR = 903)
      903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namicedic in reference namelist_csib' )


      ! read cfg namelists

      ! Namelist namicedia in configuration
      REWIND( numnatp_cfgcsib )
      READ  ( numnatp_cfgcsib, namicedia, IOSTAT = ios, ERR = 904 )
      904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namicedia in configuration namelist_csib' )
       ! Namelist namicenit in configuration
      REWIND( numnatp_cfgcsib )
      READ  ( numnatp_cfgcsib, namicenit, IOSTAT = ios, ERR = 905 )
      905   IF( ios >  0 )   CALL ctl_nam ( ios , 'namicenit in configuration namelist_csib' )
       ! Namelist namicedic in configuration
      REWIND( numnatp_cfgcsib )
      READ  ( numnatp_cfgcsib, namicedic, IOSTAT = ios, ERR = 906 )
      906   IF( ios >  0 )   CALL ctl_nam ( ios , 'namicedic in configuration namelist_csib' )
      
      

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Parameters for CSIB'
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Bottom diatoms'
         WRITE(numout,*) ' z_ia =',z_ia 
         WRITE(numout,*) ' qnidiamin =',qnidiamin 
         WRITE(numout,*) ' qnidiamax =',qnidiamax 
         WRITE(numout,*) ' qchidiaref =',qchidiaref 
         WRITE(numout,*) ' alrefidia =',alrefidia
         WRITE(numout,*) ' pcrefidia =',pcrefidia
         WRITE(numout,*) ' betaidia =',betaidia
         WRITE(numout,*) ' vnref =',vnref
         WRITE(numout,*) ' knh4 =',knh4
         WRITE(numout,*) ' kno3 =',kno3
         WRITE(numout,*) ' etares =',etares
         WRITE(numout,*) ' ch2nmax =',ch2nmax
         WRITE(numout,*) ' min_icedia =',min_icedia
         WRITE(numout,*) ' t_ia =',t_ia
         WRITE(numout,*) ' r_m1 =',r_m1
         WRITE(numout,*) ' r_m2 =',r_m2
         WRITE(numout,*) ' f_p2 =',f_p2
         WRITE(numout,*) ' f_flsh =',f_flsh
         WRITE(numout,*) ' f_slgh =',f_slgh
         WRITE(numout,*) ' dt_mo =',dt_mo
         WRITE(numout,*) ' t_mo =',t_mo
         WRITE(numout,*) ' d_mo =',d_mo
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Ice Nitrogen'
         WRITE(numout,*) ' f_rm =',f_rm
         WRITE(numout,*) ' r_ni =',r_ni
         WRITE(numout,*) ' c_di =',c_di
         WRITE(numout,*) ' c_nu =',c_nu
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Ice DIC'
         WRITE(numout,*) ' sicpump =',sicpump
         WRITE(numout,*) ' icedicref =',icedicref
         WRITE(numout,*) ' icetalref =',icetalref
         WRITE(numout,*) ' f_dicsw =',f_dicsw
         WRITE(numout,*) ' f_dicsw_melt =',f_dicsw_melt


      ENDIF

      ! convert time unit from /day to /sec
      alrefidia = alrefidia / 86400._wp 
      pcrefidia = pcrefidia / 86400._wp
      betaidia = betaidia / 86400._wp
      vnref = vnref / 86400._wp
      r_m1 = r_m1 / 86400._wp
      r_m2 = r_m2 / 86400._wp
      r_ni = r_ni  / 86400._wp
      dt_mo = dt_mo / 86400._wp
      ! convert temperature units
      t_mo = t_mo + 273.15_wp

      ! output namelists
      IF(lwm) CALL ctl_opn( numonpbcsib     , 'output.namelist.csib' , 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) WRITE( numonpbcsib, namicedia )
      IF(lwm) WRITE( numonpbcsib, namicenit )
      IF(lwm) WRITE( numonpbcsib, namicedic )

      
   END SUBROUTINE trc_nam_csib
   
   !!======================================================================
END MODULE trcnam_csib
