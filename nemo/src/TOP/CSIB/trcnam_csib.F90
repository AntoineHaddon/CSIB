MODULE trcnam_csib
   !!======================================================================
   !!                      ***  MODULE trcnam_csib  ***
   !! TOP :   initialisation of some run parameters for CSIB bio-model
   !!======================================================================
   !! History :      !  2025 (A. Haddon) Original code
   !!----------------------------------------------------------------------
   !! trc_nam_csib      : CSIB model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables

   USE trcsms_csib
   USE par_csib

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
      INTEGER ::   ios       ! Local integer

      !!----------------------------------------------------------------------
      NAMELIST/namicetra/ ln_csib, &
         & z_ia, qnidiamin, cn_fct, cn_pow, qchidiaref, alrefidia, pcrefidia, betaidia, cigr, vnref, knh4, kno3, etares, ch2nmax, min_icedia, t_ia, r_m1, r_m2, f_p2, f_flsh, f_slgh, dt_mo, t_mo, d_mo, &
         & f_rm, r_ni, c_di, c_nu, &
         & sicpump, icedicref, icetalref, f_dicsw, f_dicsw_melt, &
         & ln_dmsice, q_pi, f_zi, f_ei, f_yieldi, k_dmspdi, k_dmsi, k_freei, k_photoi, h_ni
      
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '  trc_nam_csib : read icetra namelist in ice namelist file for CSIB2'
      IF(lwp) WRITE(numout,*) 
     
      READ  ( numnam_ice_ref, namicetra, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namicetra in reference namelist' )
      READ  ( numnam_ice_cfg, namicetra, IOSTAT = ios, ERR = 902 )
902   IF( ios > 0 )   CALL ctl_nam ( ios , 'namicetra in configuration namelist' )
      
      IF(.NOT. ln_csib) THEN 
         IF(lwm) WRITE(numout,*) 'CSIB not active'
         RETURN
      END IF

      IF(lwm) WRITE( numoni, namicetra )

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) 'Parameters for CSIB'
         WRITE(numout,*) ' '
         WRITE(numout,*) 'General parameters'
         WRITE(numout,*) '  ln_csib = ', ln_csib
         WRITE(numout,*) '  ln_dmsice = ', ln_dmsice
         WRITE(numout,*) '  '
         WRITE(numout,*) 'Ice algae'
         WRITE(numout,*) '  z_ia =',z_ia 
         WRITE(numout,*) '  qnidiamin =',qnidiamin 
         WRITE(numout,*) '  cn_fct =', cn_fct 
         WRITE(numout,*) '  cn_pow =', cn_pow 
         WRITE(numout,*) '  qchidiaref =',qchidiaref 
         WRITE(numout,*) '  alrefidia =',alrefidia
         WRITE(numout,*) '  pcrefidia =',pcrefidia
         WRITE(numout,*) '  betaidia =',betaidia
         WRITE(numout,*) '  cigr =',cigr
         WRITE(numout,*) '  vnref =',vnref
         WRITE(numout,*) '  knh4 =',knh4
         WRITE(numout,*) '  kno3 =',kno3
         WRITE(numout,*) '  etares =',etares
         WRITE(numout,*) '  ch2nmax =',ch2nmax
         WRITE(numout,*) '  min_icedia =',min_icedia
         WRITE(numout,*) '  t_ia =',t_ia
         WRITE(numout,*) '  r_m1 =',r_m1
         WRITE(numout,*) '  r_m2 =',r_m2
         WRITE(numout,*) '  f_p2 =',f_p2
         WRITE(numout,*) '  f_flsh =',f_flsh
         WRITE(numout,*) '  f_slgh =',f_slgh
         WRITE(numout,*) '  dt_mo =',dt_mo
         WRITE(numout,*) '  t_mo =',t_mo
         WRITE(numout,*) '  d_mo =',d_mo
         WRITE(numout,*) '  '
         WRITE(numout,*) 'Ice Nitrogen'
         WRITE(numout,*) '  f_rm =',f_rm
         WRITE(numout,*) '  r_ni =',r_ni
         WRITE(numout,*) '  c_di =',c_di
         WRITE(numout,*) '  c_nu =',c_nu
         WRITE(numout,*) '  '
         WRITE(numout,*) 'Ice DIC'
         WRITE(numout,*) '  sicpump =',sicpump
         WRITE(numout,*) '  icedicref =',icedicref
         WRITE(numout,*) '  icetalref =',icetalref
         WRITE(numout,*) '  f_dicsw =',f_dicsw
         WRITE(numout,*) '  f_dicsw_melt =',f_dicsw_melt
      ENDIF

      ! convert time unit from /day to /sec
      alrefidia = alrefidia / 86400._wp 
      pcrefidia = pcrefidia / 86400._wp
      betaidia = betaidia / 86400._wp
      cigr = cigr / 86400._wp
      vnref = vnref / 86400._wp
      r_m1 = r_m1 / 86400._wp
      r_m2 = r_m2 / 86400._wp
      r_ni = r_ni  / 86400._wp
      dt_mo = dt_mo / 86400._wp
      ! convert temperature units
      t_mo = t_mo + 273.15_wp


      IF (ln_dmsice .AND. lwp ) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) 'Ice DMS'
         WRITE(numout,*) '  q_pi     =',q_pi
         WRITE(numout,*) '  f_zi     =',f_zi
         WRITE(numout,*) '  f_ei     =',f_ei
         WRITE(numout,*) '  f_yieldi =',f_yieldi
         WRITE(numout,*) '  k_dmspdi =',k_dmspdi
         WRITE(numout,*) '  k_dmsi   =',k_dmsi
         WRITE(numout,*) '  k_freei  =',k_freei
         WRITE(numout,*) '  k_photoi =',k_photoi
         
         ! dmsice conversions of namelist constants from 1/day to 1/sec
         k_dmspdi = k_dmspdi /  86400._wp
         k_dmsi   = k_dmsi   /  86400._wp
         k_freei  = k_freei  /  86400._wp
         k_photoi = k_photoi /  86400._wp
      ENDIF

      

      
   END SUBROUTINE trc_nam_csib
   
   !!======================================================================
END MODULE trcnam_csib
