MODULE nemo_diag_glovars_canoe

   !!========================================================
   !!          ***   MODULE nemo_diag_glovars   ***          
   !!          declarations of all global variables         
   !!========================================================
   !! 2018-04 (D. Yang): Original code
   !! 2019-01 (J. Christian): biogeochemistry (CanOE) version
   !!--------------------------------------------------------

      INTEGER            :: imt, jmt, km, lm, ly            ! i, j, k and time dimensions
                                                            ! for monthly and yearly
      REAL, DIMENSION(:,:), ALLOCATABLE       :: e2u, e1v
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE     :: e3u, e3v, e3t
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: umask, vmask, tmask, fsdept
      REAL, DIMENSION(:), ALLOCATABLE         :: deptht
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: TT, SS, CC, AA
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: asi3, NO3, NH4, O2
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: co3_sata, co3_satc
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: CO3, pH
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: K_sp_cal, K_sp_arag, Omega_C, Omega_A
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: prhop, o2sol
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: zsat_a, zsat_c
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: o2min, zo2min

END MODULE nemo_diag_glovars_canoe

