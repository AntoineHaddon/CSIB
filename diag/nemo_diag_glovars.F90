MODULE nemo_diag_glovars

   !!========================================================
   !!          ***   MODULE nemo_diag_glovars   ***          
   !!          declarations of all global variables         
   !!========================================================
   !! 2018-04 (D. Yang): Original code
   !!--------------------------------------------------------

      INTEGER            :: imt, jmt, km, lm, ly            ! i, j, k and time dimensions
                                                            ! for monthly and yearly
      INTEGER, PARAMETER :: nline = 15                      ! number of sections
      INTEGER, PARAMETER :: nline_ice = 4                   ! number of ice-sections
      REAL, DIMENSION(:,:), ALLOCATABLE       :: e2u, e1v
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE     :: e3u, e3v, e3t
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: umask, vmask, tmask
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: u, v
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: ssh       ! sossheig (m)
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: tnc, tnp, snc, snp  ! tn & sn at previous 
                                                                     ! and current time steps
      REAL, DIMENSION(:,:), ALLOCATABLE       :: mfo       ! mass transport through transects
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: msftbarot ! quasi-barotropic streamfuction
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE   :: opottemptend, osalttend ! T & S tendencies
      ! ice transport 
      !REAL, DIMENSION(:,:,:), ALLOCATABLE     :: xmtrpice,xmtrpsnw,xatrp ! i-direction inputs (not deeded in our sections)
      REAL, DIMENSION(:,:,:), ALLOCATABLE     :: ymtrpice,ymtrpsnw,yatrp ! j-direction inputs
      REAL, DIMENSION(:,:), ALLOCATABLE       :: simassacrossline, snmassacrossline, siareaacrossline  ! transport through ice-transects

END MODULE nemo_diag_glovars





