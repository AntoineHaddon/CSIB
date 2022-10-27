MODULE cmocrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP/CANBGC :   CMOC Compute remineralization/scavenging of organic compounds
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!          CMOC 1  !  2013-2015(O. Riche) remineralization, PIC export and nitrogen fixation, based on Zahariev et al 2008
   !!          CMOC 1  !  2016-02  (N. Swart) Bugfixes and moves calcite flux to p4zsink; DNF to p4zsed.
   !!          CMOC 2  !  2022-23  (O. Riche) NEMO4 remineralization   
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   cmoc_rem       :  Compute remineralization/scavenging of organic compounds
   !!   cmoc_rem_init  :  Initialisation of parameters for remineralisation
   !!----------------------------------------------------------------------   
CONTAINS

  SUBROUTINE cmoc_rem
  
  END SUBROUTINE cmoc_rem
  
  
  SUBROUTINE cmoc_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the namcmocpoc namelist and check the parameters
      !!                called at the first timestep
      !!
      !! ** input   :   Namelist namcmocpoc
      !!
      !!----------------------------------------------------------------------

      ! <CMOC code OR 10/15/2015> CMOC namelist


      ! control print
      IF(lwp) THEN
         WRITE(numout,*) ' Namelist parameters for remineralization, namcmocpoc'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    Remineralisation rate of POC              reref_cmoc=', reref_cmoc
         WRITE(numout,*) '    Activation energy for remineralization    ed_cmoc   =', ed_cmoc   
         WRITE(numout,*) ' '
      ENDIF
      ! 
  END SUBROUTINE cmoc_rem_init

END MODULE cmocrem