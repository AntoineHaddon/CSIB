MODULE nemo_diag_glovars_cmoc

   !!========================================================
   !!          ***   MODULE nemo_diag_glovars   ***          
   !!          declarations of all global variables         
   !!========================================================
   !! 2018-04 (D. Yang): Original code
   !! 2019-01 (J. Christian): biogeochemistry (CMOC) version
   !!--------------------------------------------------------

      INTEGER            :: imt, jmt, km, lm, ly            ! i, j, k and time dimensions
                                                            ! for monthly and yearly
      REAL, DIMENSION(:,:), ALLOCATABLE       :: e2u, e1v
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE     :: e3u, e3v, e3t
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: umask, vmask, tmask, fsdept
      REAL, DIMENSION(:), ALLOCATABLE         :: deptht
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: TT, SS
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: CC, AA, CAB, CNT
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: asi3, NO3, O2
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: co3_sata, co3_satc
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: CO3full, pHfull
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: CO3abio, pHabio, CO3nat, pHnat
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: K_sp_cal, K_sp_arag, Omega_C, Omega_A
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: Omega_C_abio, Omega_A_abio
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: Omega_C_nat, Omega_A_nat
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: prhop, o2sol
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: zsat_a, zsat_c
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: o2min, zo2min
 
END MODULE nemo_diag_glovars_cmoc

