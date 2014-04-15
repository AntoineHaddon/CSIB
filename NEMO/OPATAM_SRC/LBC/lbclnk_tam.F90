MODULE lbclnk_tam
   !!======================================================================
   !!                       ***  MODULE  lbclnk_tam  ***
   !! Ocean        : TAM of lateral boundary conditions
   !!=====================================================================
   !! History :  OPA  ! 1997-06  (G. Madec)     Original code
   !!   NEMO     1.0  ! 2002-09  (G. Madec)     F90: Free form and module
   !!            3.2  ! 2009-03  (R. Benshila)  External north fold treatment
   !! History of TAM : 3.2 ! ???
   !!                  3.4 ! 2012-03 (P.-A. Bouttier) Update
   !!----------------------------------------------------------------------
#if   defined key_mpp_mpi
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'             MPI massively parallel processing library
   !!----------------------------------------------------------------------
   !!   lbc_lnk_adj      : generic interface for mpp_lnk_3d_adj and mpp_lnk_2d_adj routines defined in lib_mpp
   !!   lbc_lnk_e_adj    : generic interface for mpp_lnk_2d_e_adj routine defined in lib_mpp
   !!----------------------------------------------------------------------
   USE lib_mpp_tam          ! distributed memory computing library

   INTERFACE lbc_lnk_adj
      MODULE PROCEDURE mpp_lnk_3d_gather_adj, mpp_lnk_3d_adj, mpp_lnk_2d_adj
   END INTERFACE

   INTERFACE lbc_lnk_e_adj
      MODULE PROCEDURE mpp_lnk_2d_e_adj
   END INTERFACE

   PUBLIC lbc_lnk_adj       ! ocean lateral boundary conditions
   PUBLIC lbc_lnk_e_adj

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

#else
   !!----------------------------------------------------------------------
   !!   Default option                              shared memory computing
   !!----------------------------------------------------------------------
   !!   lbc_lnk_adj      : generic interface for lbc_lnk_3d and lbc_lnk_2d
   !!   lbc_lnk_3d_adj   : set the lateral boundary condition on a 3D variable on ocean mesh
   !!   lbc_lnk_2d_adj   : set the lateral boundary condition on a 2D variable on ocean mesh
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE lbcnfd_tam          ! north fold

   IMPLICIT NONE
   PRIVATE

   INTERFACE lbc_lnk_adj
      MODULE PROCEDURE lbc_lnk_3d_gather_adj, lbc_lnk_3d_adj, lbc_lnk_2d_adj
   END INTERFACE

   INTERFACE lbc_lnk_e_adj
      MODULE PROCEDURE lbc_lnk_2d_adj
   END INTERFACE

   PUBLIC   lbc_lnk_adj       ! ocean/ice  lateral boundary conditions
   PUBLIC   lbc_lnk_e_adj

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lbc_lnk_3d_gather_adj( pt3d1, cd_type1, pt3d2, cd_type2, psgn )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE lbc_lnk_3d_gather_adj  ***
      !!
      !! ** Purpose :   Adjoint of set lateral boundary conditions on two 3D arrays (non mpp case)
      !!
      !! ** Method  :   psign = -1 :    change the sign across the north fold
      !!                      =  1 : no change of the sign across the north fold
      !!                      =  0 : no change of the sign across the north fold and
      !!                             strict positivity preserved: use inner row/column
      !!                             for closed boundaries.
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                , INTENT(in   ) ::   cd_type1, cd_type2   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pt3d1   , pt3d2      ! 3D array on which the lbc is applied
      REAL(wp)                        , INTENT(in   ) ::   psgn                 ! control of the sign
      !!----------------------------------------------------------------------
      !
      CALL lbc_lnk_3d_adj( pt3d2, cd_type2, psgn)
      CALL lbc_lnk_3d_adj( pt3d1, cd_type1, psgn)
      !
   END SUBROUTINE lbc_lnk_3d_gather_adj


   SUBROUTINE lbc_lnk_3d_adj( pt3d, cd_type, psgn, cd_mpp, pval )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE lbc_lnk_3d_adj  ***
      !!
      !! ** Purpose :   Adjoint of set lateral boundary conditions on a 3D array (non mpp case)
      !!
      !! ** Method  :   psign = -1 :    change the sign across the north fold
      !!                      =  1 : no change of the sign across the north fold
      !!                      =  0 : no change of the sign across the north fold and
      !!                             strict positivity preserved: use inner row/column
      !!                             for closed boundaries.
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout)           ::   pt3d      ! 3D array on which the lbc is applied
      REAL(wp)                        , INTENT(in   )           ::   psgn      ! control of the sign
      CHARACTER(len=3)                , INTENT(in   ), OPTIONAL ::   cd_mpp    ! MPP only (here do nothing)
      REAL(wp)                        , INTENT(in   ), OPTIONAL ::   pval      ! background value (for closed boundaries)
      !!
      REAL(wp) ::   zland
      !!----------------------------------------------------------------------

      IF( PRESENT( pval ) ) THEN   ;   zland = pval      ! set land value (zero by default)
      ELSE                         ;   zland = 0.e0
      ENDIF


      IF( PRESENT( cd_mpp ) ) THEN
         ! only fill the overlap area and extra allows
         ! this is in mpp case. In this module, just do nothing
      ELSE
         !
         !                                     ! North-South boundaries
         !                                     ! ======================
         SELECT CASE ( nperio )
         !
         CASE ( 2 )                               !**  South symmetric  --  North closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'W' )                   ! T-, U-, W-points
               pt3d(:,jpj,:) = 0.0_wp
               pt3d(:, 3 ,:) = pt3d(:,3,:) + pt3d(:,1,:)
               pt3d(:, 1 ,:) = 0.0_wp
            CASE ( 'V' , 'F' )                         ! V-, F-points
               pt3d(:,jpj,:) = 0.0_wp
               pt3d(:, 2 ,:) = pt3d(:,2,:) + psgn * pt3d(:,1,:)
               pt3d(:, 1 ,:) = 0.0_wp
            END SELECT
            !
         CASE ( 3 , 4 , 5 , 6 )                   !**  North fold  T or F-point pivot  --  South closed
            CALL lbc_nfd_adj( pt3d(:,:,:), cd_type, psgn  )  ! North fold
            pt3d(jpi,jpj,:) = 0.0_wp
            pt3d( 1 ,jpj,:) = 0.0_wp
            !
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' , 'I' )             ! all points except F-point
               pt3d(:, 1 ,:) = 0.0_wp
            END SELECT
            !
         CASE DEFAULT                             !**  North closed  --  South closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
               pt3d(:,jpj,:) = 0.0_wp
               pt3d(:, 1 ,:) = 0.0_wp
            CASE ( 'F' )                               ! F-point
               pt3d(:,jpj,:) = 0.0_wp
            END SELECT
            !
         END SELECT
         !
         !                                     !  East-West boundaries
         !                                     ! ======================
         SELECT CASE ( nperio )
         !
         CASE ( 1 , 4 , 6 )                       !**  cyclic east-west
            pt3d(  2  ,:,:) = pt3d(  2  ,:,:) + pt3d(jpi,:,:)
            pt3d(jpi,:,:) = 0.0_wp
            pt3d(jpim1,:,:) = pt3d(jpim1,:,:) + pt3d( 1 ,:,:)
            pt3d( 1 ,:,:) = 0.0_wp
            !
         CASE DEFAULT                             !**  East closed  --  West closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
               pt3d( 1 ,:,:) = 0.0_wp
               pt3d(jpi,:,:) = 0.0_wp
            CASE ( 'F' )                               ! F-point
               pt3d(jpi,:,:) = 0.0_wp
            END SELECT
            !
         END SELECT
         !
      ENDIF
      !
   END SUBROUTINE lbc_lnk_3d_adj


   SUBROUTINE lbc_lnk_2d_adj( pt2d, cd_type, psgn, cd_mpp, pval )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE lbc_lnk_2d_adj  ***
      !!
      !! ** Purpose :   Adjoint of set lateral boundary conditions on a 2D array (non mpp case)
      !!
      !! ** Method  :   psign = -1 :    change the sign across the north fold
      !!                      =  1 : no change of the sign across the north fold
      !!                      =  0 : no change of the sign across the north fold and
      !!                             strict positivity preserved: use inner row/column
      !!                             for closed boundaries.
      !!----------------------------------------------------------------------
      CHARACTER(len=1)            , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout)           ::   pt2d      ! 2D array on which the lbc is applied
      REAL(wp)                    , INTENT(in   )           ::   psgn      ! control of the sign
      CHARACTER(len=3)            , INTENT(in   ), OPTIONAL ::   cd_mpp    ! MPP only (here do nothing)
      REAL(wp)                    , INTENT(in   ), OPTIONAL ::   pval      ! background value (for closed boundaries)
      !!
      REAL(wp) ::   zland
      !!----------------------------------------------------------------------
      IF( PRESENT( pval ) ) THEN   ;   zland = pval      ! set land value (zero by default)
      ELSE                         ;   zland = 0.e0
      ENDIF
      IF (PRESENT(cd_mpp)) THEN
         ! only fill the overlap area and extra allows
         ! this is in mpp case. In this module, just do nothing
      ELSE
         !
         !                                     ! North-South boundaries
         !                                     ! ======================
         SELECT CASE ( nperio )
         !
         CASE ( 2 )                               !**  South symmetric  --  North closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'W' )                   ! T-, U-, W-points
               pt2d(:,jpj) = 0.0_wp
               pt2d(:, 3 ) = pt2d(:,3) + pt2d(:,1)
               pt2d(:, 1 ) = 0.0_wp
            CASE ( 'V' , 'F' )                         ! V-, F-points
               pt2d(:,jpj) = 0.e0
               pt2d(:, 2 ) = pt2d(:,2) + psgn * pt2d(:,1)
               pt2d(:, 1 ) = 0.e0
            END SELECT
            !
         CASE ( 3 , 4 , 5 , 6 )                   !**  North fold  T or F-point pivot  --  South closed
            CALL lbc_nfd_adj( pt2d(:,:), cd_type, psgn )
            pt2d(jpi,jpj) = 0.0_wp
            pt2d( 1 ,jpj) = 0.0_wp
            pt2d( 1 , 1 ) = 0.0_wp
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' ,  'I')             ! all points except F-point
               pt2d(:, 1 ) = 0.0_wp
            END SELECT
            !
         CASE DEFAULT                             !**  North closed  --  South closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
               pt2d(:, 1 ) = 0.0_wp
               pt2d(:,jpj) = 0.0_wp
            CASE ( 'F' )                               ! F-point
               pt2d(:,jpj) = 0.0_wp
            END SELECT
            !
         END SELECT
         !
         !                                     ! East-West boundaries
         !                                     ! ====================
         SELECT CASE ( nperio )
         !
         CASE ( 1 , 4 , 6 )                       !** cyclic east-west
            pt2d(  2  ,:) = pt2d(  2  ,:) + pt2d( jpi ,:)
            pt2d( jpi ,:) = 0.0_wp
            pt2d(jpim1,:) = pt2d(jpim1,:) + pt2d( 1 ,:)
            pt2d(  1  ,:) = 0.0_wp
            !
         CASE DEFAULT                             !** East closed  --  West closed
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'U' , 'V' , 'W' )            ! T-, U-, V-, W-points
               pt2d( 1 ,:) = 0.0_wp
               pt2d(jpi,:) = 0.0_wp
            CASE ( 'F' )                              ! F-point
               pt2d(jpi,:) = 0.0_wp
            END SELECT
            !
         END SELECT
         !
      ENDIF
      !
   END SUBROUTINE lbc_lnk_2d_adj

#endif

   !!======================================================================
END MODULE lbclnk_tam
