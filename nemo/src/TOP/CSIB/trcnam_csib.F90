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
      NAMELIST/namicedia/ z_ia, mu_max, t_ia, r_pp, h_ni, vnh4, C2N_dia, C2CH_dia, b_ia, r_m1, r_m2, f_p2, f_rm, r_ni, c_di, c_nu
      
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_csib : read CSIB namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~'
     
      ! open namelist
      clname = 'namelist_csib'
      CALL ctl_opn( numnatp_refcsib, TRIM( clname )//'_ref', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      CALL ctl_opn( numnatp_cfgcsib, TRIM( clname )//'_cfg', 'OLD'    , 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      
      ! read ref and cfg namelist
      REWIND( numnatp_refcsib )              ! Namelist namicedia in reference namelist
      READ  ( numnatp_refcsib, namicedia, IOSTAT = ios, ERR = 901)
      901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namicedia in reference namelist_csib' )
      REWIND( numnatp_cfgcsib )              ! Namelist namicedia in configuration namelist
      READ  ( numnatp_cfgcsib, namicedia, IOSTAT = ios, ERR = 902 )
      902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namicedia in configuration namelist_csib' )
      
      

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Parameters for CSIB'
         WRITE(numout,*) ' '
         WRITE(numout,*) ' z_ia =',z_ia 
         WRITE(numout,*) ' mu_max =',mu_max
         WRITE(numout,*) ' t_ia =',t_ia
         WRITE(numout,*) ' r_pp =',r_pp
         WRITE(numout,*) ' h_ni =',h_ni
         WRITE(numout,*) ' vnh4 =',vnh4
         WRITE(numout,*) ' C2N_dia =',C2N_dia
         WRITE(numout,*) ' C2CH_dia =',C2CH_dia
         WRITE(numout,*) ' b_ia =',b_ia
         WRITE(numout,*) ' r_m1 =',r_m1
         WRITE(numout,*) ' r_m2 =',r_m2
         WRITE(numout,*) ' f_p2 =',f_p2
         WRITE(numout,*) ' f_rm =',f_rm
         WRITE(numout,*) ' r_ni =',r_ni
         WRITE(numout,*) ' c_nu =',c_nu
         WRITE(numout,*) ' c_di =',c_di
         WRITE(numout,*) ' '
      ENDIF


      N2C_dia = 1._wp / C2N_dia ! N to C ratio
      CH2C_dia = 1._wp / C2CH_dia ! CH to C ratio

      ! convert time unit from /day to /sec
      mu_max = mu_max / 86400._wp 
      r_m1 = r_m1 / 86400._wp
      r_m2 = r_m2 / 86400._wp
      r_ni = r_ni  / 86400._wp

      ! output namelist
      IF(lwm) CALL ctl_opn( numonpbcsib     , 'output.namelist.csib' , 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
      IF(lwm) WRITE( numonpbcsib, namicedia )

      
   END SUBROUTINE trc_nam_csib
   
   !!======================================================================
END MODULE trcnam_csib
