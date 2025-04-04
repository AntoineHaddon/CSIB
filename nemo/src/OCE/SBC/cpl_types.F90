module cpl_types

use par_kind, only : wp

implicit none; private

INTEGER, PUBLIC, PARAMETER ::   nmaxfld=80   ! Maximum number of coupling fields
INTEGER, PUBLIC, PARAMETER ::   nmaxcat=5    ! Maximum number of coupling fields
INTEGER, PUBLIC, PARAMETER ::   nmaxcpl=5    ! Maximum number of coupling fields
TYPE, PUBLIC ::   FLD_C                     !
   CHARACTER(len = 32) ::   cldes      ! desciption of the coupling strategy
   CHARACTER(len = 32) ::   clcat      ! multiple ice categories strategy
   CHARACTER(len = 32) ::   clvref     ! reference of vector ('spherical' or 'cartesian')
   CHARACTER(len = 32) ::   clvor      ! orientation of vector fields ('eastward-northward' or 'local grid')
   CHARACTER(len = 32) ::   clvgrd     ! grids on which is located the vector fields
END TYPE FLD_C

TYPE, PUBLIC ::   FLD_CPL               !: Type for coupling field information
   LOGICAL               ::   laction   ! To be coupled or not
   CHARACTER(len = 8)    ::   clname    ! Name of the coupling field
   CHARACTER(len = 1)    ::   clgrid    ! Grid type
   REAL(wp)              ::   nsgn      ! Control of the sign change
   INTEGER, DIMENSION(nmaxcat,nmaxcpl) ::   nid   ! Id of the field (no more than 9 categories and 9 extrena models)
   INTEGER               ::   nct       ! Number of categories in field
   INTEGER               ::   ncplmodel ! Maximum number of models to/from which this variable may be sent/received
   CHARACTER(len = 32)   ::   cldes      ! desciption of the coupling strategy
END TYPE FLD_CPL

TYPE(FLD_CPL), DIMENSION(nmaxfld), PUBLIC, TARGET, SAVE ::   srcv, ssnd   !: Coupling fields

INTEGER, PUBLIC, SAVE :: COUPLER_Rcv = 1
INTEGER, PUBLIC, SAVE :: COUPLER_Snd = 2
INTEGER, PUBLIC, SAVE :: COUPLER_idle = 0




end module cpl_types
