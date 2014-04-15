MODULE dynspg_exp_tam
   !!======================================================================
   !!                   ***  MODULE  dynspg_exp_tam  TANGENT/ADJOINT OF MODULE dynspg_exp***
   !! Ocean dynamics:  surface pressure gradient trend
   !!======================================================================
   !! History of the direct module:
   !!            2.0  !  2005-11  (V. Garnier, G. Madec, L. Bessieres) Original code
   !!            3.2  !  2009-06  (G. Madec, M. Leclair, R. Benshila) introduce sshwzv module
   !! History of the tam module:
   !!            3.2  !  2010-06  (A. Vidard) tam of the 2009-06 version
   !!----------------------------------------------------------------------
#if defined key_dynspg_exp   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_dynspg_exp'                              explicit free surface
   !!----------------------------------------------------------------------
   !!   dyn_spg_exp  : update the momentum trend with the surface
   !!                      pressure gradient in the free surface constant
   !!                      volume case with vector optimization
   !!----------------------------------------------------------------------
   USE par_kind
   USE phycst
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE gridrandom
   USE dotprodfld
   USE paresp
   USE in_out_manager
   USE tstool_tam
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_spg_exp_tan       ! routine called by step.F90
   PUBLIC   dyn_spg_exp_adj       ! routine called by step.F90
   PUBLIC   dyn_spg_exp_adj_tst   ! routine called by tamtst.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LOCEAN-IPSL (2009)
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_spg_exp_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_exp_tan  ***
      !!
      !! ** Purpose :   Compute the now trend due to the surface pressure
      !!              gradient in case of explicit free surface formulation and
      !!              add it to the general trend of momentum equation.
      !!
      !! ** Method  :   Explicit free surface formulation. Add to the general
      !!              momentum trend the surface pressure gradient :
      !!                      (ua,va) = (ua,va) + (spgu,spgv)
      !!              where spgu = -1/rau0 d/dx(ps) = -g/e1u di( sshn )
      !!                    spgv = -1/rau0 d/dy(ps) = -g/e2v dj( sshn )
      !!
      !! ** Action :   (ua,va)   trend of horizontal velocity increased by
      !!                         the surf. pressure gradient trend
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kt         ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk               ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg_exp_tan')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_exp_tan : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   (explicit free surface)'
         !
         spgu_tl(:,:) = 0._wp  ;   spgv_tl(:,:) = 0._wp
         !
         IF( lk_vvl .AND. lwp ) WRITE(numout,*) '              lk_vvl=T : spg is included in dynhpg'
      ENDIF
      IF( .NOT. lk_vvl ) THEN          !* fixed volume : add the surface pressure gradient trend
         !
         DO jj = 2, jpjm1                    ! now surface pressure gradient
            DO ji = fs_2, fs_jpim1   ! vector opt.
               spgu_tl(ji,jj) = - grav * ( sshn_tl(ji+1,jj) - sshn_tl(ji,jj) ) / e1u(ji,jj)
               spgv_tl(ji,jj) = - grav * ( sshn_tl(ji,jj+1) - sshn_tl(ji,jj) ) / e2v(ji,jj)
            END DO
         END DO
         DO jk = 1, jpkm1                    ! Add it to the general trend
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ua_tl(ji,jj,jk) = ua_tl(ji,jj,jk) + spgu_tl(ji,jj)
                  va_tl(ji,jj,jk) = va_tl(ji,jj,jk) + spgv_tl(ji,jj)
               END DO
            END DO
         END DO
         !
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_spg_exp_tan')
      !
   END SUBROUTINE dyn_spg_exp_tan
   SUBROUTINE dyn_spg_exp_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_exp_adj  ***
      !!
      !! ** Purpose :   Compute the now trend due to the surface pressure
      !!              gradient in case of explicit free surface formulation and
      !!              add it to the general trend of momentum equation.
      !!
      !! ** Method  :   Explicit free surface formulation. Add to the general
      !!              momentum trend the surface pressure gradient :
      !!                      (ua,va) = (ua,va) + (spgu,spgv)
      !!              where spgu = -1/rau0 d/dx(ps) = -g/e1u di( sshn )
      !!                    spgv = -1/rau0 d/dy(ps) = -g/e2v dj( sshn )
      !!
      !! ** Action :   (ua,va)   trend of horizontal velocity increased by
      !!                         the surf. pressure gradient trend
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kt         ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk               ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg_exp_adj')
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_exp_adj : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   (explicit free surface)'
      END IF

      spgu_ad(:,:) = 0._wp   ;   spgv_ad(:,:) = 0._wp

      IF( .NOT. lk_vvl ) THEN          !* fixed volume : add the surface pressure gradient trend
         !
         DO jk = 1, jpkm1                    ! Add it to the general trend
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  spgu_ad(ji,jj) = spgu_ad(ji,jj) + ua_ad(ji,jj,jk)
                  spgv_ad(ji,jj) = spgv_ad(ji,jj) + va_ad(ji,jj,jk)
               END DO
            END DO
         END DO
         DO jj = jpjm1, 2, -1                    ! now surface pressure gradient
            DO ji = fs_jpim1, fs_2, -1           ! vector opt.

               spgu_ad(ji,jj) = - grav * spgu_ad(ji,jj) / e1u(ji,jj)
               spgv_ad(ji,jj) = - grav * spgv_ad(ji,jj) / e2v(ji,jj)
               sshn_ad(ji+1,jj) = sshn_ad(ji+1,jj) + spgu_ad(ji,jj)
               sshn_ad(ji,jj+1) = sshn_ad(ji,jj+1) + spgv_ad(ji,jj)
               sshn_ad(ji,jj)   = sshn_ad(ji,jj) - spgu_ad(ji,jj) - spgv_ad(ji,jj)
            END DO
         END DO
         !
      ENDIF
      !
   !
   IF( nn_timing == 1 )  CALL timing_stop('dyn_spg_exp_adj')
   !
   END SUBROUTINE dyn_spg_exp_adj
   SUBROUTINE dyn_spg_exp_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE dyn_spg_exp_adj_tst ***
      !!
      !! ** Purpose : Test the adjoint routine.
      !!
      !! ** Method  : Verify the scalar product
      !!
      !!                 ( L dx )^T W dy  =  dx^T L^T W dy
      !!
      !!              where  L   = tangent routine
      !!                     L^T = adjoint routine
      !!                     W   = diagonal matrix of scale factors
      !!                     dx  = input perturbation (random field)
      !!                     dy  = L dx
      !!
      !!
      !! History :
      !!        ! 2010-06 (A. Vidard)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      !! * Local declarations
      INTEGER ::  &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d        ! 2D seed for the random number generator
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zua_tlin ,     & ! Tangent input
         & zva_tlin ,     & ! Tangent input
         & zua_tlout,     & ! Tangent output
         & zva_tlout,     & ! Tangent output
         & zua_adin ,     & ! Adjoint input
         & zva_adin ,     & ! Adjoint input
         & zua_adout,     & ! Adjoint output
         & zva_adout,     & ! Adjoint output
         & zr3d           ! 3D random field
     REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
        & zsshn_tlin,    &
        & zsshn_adout,   &
        & zr2d
       CHARACTER(LEN=14) :: &
         & cl_name
      ! Allocate memory

      ALLOCATE( &
         & zua_tlin(   jpi,jpj,jpk),     &
         & zva_tlin(   jpi,jpj,jpk),     &
         & zsshn_tlin( jpi,jpj    ),     &
         & zua_tlout(  jpi,jpj,jpk),     &
         & zva_tlout(  jpi,jpj,jpk),     &
         & zua_adin(   jpi,jpj,jpk),     &
         & zva_adin(   jpi,jpj,jpk),     &
         & zua_adout(  jpi,jpj,jpk),     &
         & zva_adout(  jpi,jpj,jpk),     &
         & zsshn_adout(jpi,jpj    ),     &
         & zr3d(       jpi,jpj,jpk),     &
         & zr2d(       jpi,jpj    )      &
         & )
      !==================================================================
      ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
      !    dy = ( hdivb_tl, hdivn_tl )
      !==================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------

      ua_ad( :,:,:) = 0.0_wp
      va_ad( :,:,:) = 0.0_wp
      sshn_ad( :,:) = 0.0_wp
      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random( zr3d, 'U', 0.0_wp, stdu )
      zua_tlin(:,:,:) = zr3d(:,:,:)
      CALL grid_random( zr3d, 'V', 0.0_wp, stdv )
      zva_tlin(:,:,:) = zr3d(:,:,:)
      CALL grid_random( zr2d, 'T', 0.0_wp, stdssh )
      zsshn_tlin(:,:) = zr2d(:,:)


      ua_tl = zua_tlin
      va_tl = zva_tlin
      sshn_tl = zsshn_tlin
      CALL dyn_spg_exp_tan( nit000 )
      zua_tlout = ua_tl
      zva_tlout = va_tl
      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zua_adin(ji,jj,jk) = zua_tlout(ji,jj,jk) &
                 &               * e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk) &
                 &               * umask(ji,jj,jk)
              zva_adin(ji,jj,jk) = zva_tlout(ji,jj,jk) &
                 &               * e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk) &
                 &               * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1 = DOT_PRODUCT( zua_tlout, zua_adin ) &
         & + DOT_PRODUCT( zva_tlout, zva_adin )

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      ua_ad = zua_adin
      va_ad = zva_adin

      CALL dyn_spg_exp_adj( nit000 )

      zua_adout   = ua_ad
      zva_adout   = va_ad
      zsshn_adout = sshn_ad

      zsp2 = DOT_PRODUCT( zua_tlin  , zua_adout   ) &
         & + DOT_PRODUCT( zva_tlin  , zva_adout   ) &
         & + DOT_PRODUCT( zsshn_tlin, zsshn_adout )

      ! 14 char:'12345678901234'
      cl_name = 'dyn_spg_exp   '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE(         &
         & zua_tlin,      &
         & zva_tlin,      &
         & zsshn_tlin,    &
         & zua_tlout,     &
         & zva_tlout,     &
         & zua_adin,      &
         & zva_adin,      &
         & zua_adout,     &
         & zva_adout,     &
         & zsshn_adout,   &
         & zr3d,          &
         & zr2d           &
         & )

   END SUBROUTINE dyn_spg_exp_adj_tst

#else
   !!----------------------------------------------------------------------
   !!   Default case :   Empty module   No standart explicit free surface
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE dyn_spg_exp_tan( kt )       ! Empty routine
      WRITE(*,*) 'dyn_spg_exp: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_exp_tan
   SUBROUTINE dyn_spg_exp_adj( kt )       ! Empty routine
      WRITE(*,*) 'dyn_spg_exp: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_exp_adj
   SUBROUTINE dyn_spg_exp_adj_tst( kt )       ! Empty routine
      WRITE(*,*) 'dyn_spg_exp: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_exp_adj_tst
#endif

   !!======================================================================
END MODULE dynspg_exp_tam
