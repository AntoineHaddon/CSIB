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
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: e3u, e3v, e3t
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: umask, vmask, tmask, fsdept
      REAL, DIMENSION(:), ALLOCATABLE         :: deptht
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: TT, SS
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: CC, AA, CAB, CNT
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: asi3, NO3
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: borat, ak13, ak23, akb3, akw3
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: akp13, akp23, akp33, aksi3
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: co3_sata, co3_satc
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: CO3full, pHfull
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: CO3abio, pHabio, CO3nat, pHnat
!      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: XDIC, XTA, CO3, pH
!      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: hi
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: o2sol
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: prhop

END MODULE nemo_diag_glovars_cmoc

