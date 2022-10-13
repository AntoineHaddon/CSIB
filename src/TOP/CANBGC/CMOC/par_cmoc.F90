MODULE par_cmoc
   !!======================================================================
   !!                        ***  par_cmoc  ***
   !! TOP :   set the CMOC parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!                  !  2022     (O. Riche) NEMO4 integration
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: par_cmoc.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

   IMPLICIT NONE

   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   ! INTEGER, PUBLIC ::   jp_knu0             !: First index of CMOC passive tracers
   ! INTEGER, PUBLIC ::   jp_knu1             !: Last  index of CMOC passive tracers
   ! INTEGER, PUBLIC ::   jqdic_can           !: DIC
   ! INTEGER, PUBLIC ::   jpalk_can           !: TA
   ! INTEGER, PUBLIC ::   jqoxy_can           !: O2 
   INTEGER, PUBLIC :: jqpoc !
   INTEGER, PUBLIC :: jqphy !
   INTEGER, PUBLIC :: jqzoo !
   INTEGER, PUBLIC :: jqnch !
   ! INTEGER, PUBLIC :: jpdab !: abiotic DIC 
   ! INTEGER, PUBLIC :: jpaab !: abiotic Alkalinity
   ! INTEGER, PUBLIC :: jpoab !: abiotic oxygen
   ! INTEGER, PUBLIC :: jpdnt !: natural DIC
   ! INTEGER, PUBLIC :: jpdrc !: abiotic DI14C
   !!======================================================================
END MODULE par_cmoc
