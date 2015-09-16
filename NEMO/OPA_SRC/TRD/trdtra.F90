MODULE trdtra
   !!======================================================================
   !!                       ***  MODULE  trdtra  ***
   !! Ocean diagnostics:  ocean tracers trends
   !!=====================================================================
   !! History :  1.0  !  2004-08  (C. Talandier) Original code
   !!            2.0  !  2005-04  (C. Deltel)    Add Asselin trend in the ML budget
   !!            3.3  !  2010-06  (C. Ethe) merge TRA-TRC 
   !!            3.4.1!  2015-08  (D. Yang) output 3D tracer trends
   !!            3.4.1!  2015-09  (D. Yang) added diagnostics for the "PURE" Kz trend
   !!                                       in case of iso-neutral diffusion
   !!----------------------------------------------------------------------
#if  defined key_trdtra || defined key_trdtrc || defined key_trdmld || defined key_trdmld_trc 
   !!----------------------------------------------------------------------
   !!   trd_tra      : Call the trend to be computed
   !!   trd_tra_adv   : transform a div(U.T) trend into a U.grad(T) trend
   !!   trd_tra_mng   : tracer trend manager: dispatch to the diagnostic modules
   !!   trd_tra_iom   : output 3D tracer trends using IOM
   !!----------------------------------------------------------------------
   USE oce              ! ocean dynamics and tracers variables
   USE dom_oce          ! ocean domain 
   USE zdf_oce          ! ocean vertical physics
   USE trdmod_oce       ! ocean active mixed layer tracers trends 
   USE trdmod           ! ocean active mixed layer tracers trends 
   USE trdmod_trc       ! ocean passive mixed layer tracers trends 
   USE zdfddm           ! vertical physics: double diffusion
   USE in_out_manager   ! I/O manager
   USE iom              ! I/O manager library
   USE lib_mpp          ! MPP library
   USE wrk_nemo         ! Memory allocation


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trd_tra          ! called by all  traXX modules
 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: trdtx, trdty, trdt  !:

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: trdtra.F90 3425 2012-07-05 11:33:22Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trd_tra_alloc()
      !!----------------------------------------------------------------------------
      !!                  ***  FUNCTION trd_tra_alloc  ***
      !!----------------------------------------------------------------------------
      ALLOCATE( trdtx(jpi,jpj,jpk) , trdty(jpi,jpj,jpk) , trdt(jpi,jpj,jpk) , STAT= trd_tra_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum ( trd_tra_alloc )
      IF( trd_tra_alloc /= 0 )   CALL ctl_warn('trd_tra_alloc: failed to allocate arrays')
   END FUNCTION trd_tra_alloc


   SUBROUTINE trd_tra( kt, ctype, ktra, ktrd, ptrd, pun, ptra )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra  ***
      !! 
      !! ** Purpose : Dispatch all trends computation, e.g. vorticity, mld or 
      !!              integral constraints
      !!
      !! ** Method/usage : For the mixed-layer trend, the control surface can be either
      !!       a mixed layer depth (time varying) or a fixed surface (jk level or bowl). 
      !!      Choose control surface with nn_ctls in namelist NAMTRD :
      !!        nn_ctls = 0  : use mixed layer with density criterion 
      !!        nn_ctls = 1  : read index from file 'ctlsurf_idx'
      !!        nn_ctls > 1  : use fixed level surface jk = nn_ctls
      !!----------------------------------------------------------------------
      !
      INTEGER                         , INTENT(in)           ::  kt      ! time step
      CHARACTER(len=3)                , INTENT(in)           ::  ctype   ! tracers trends type 'TRA'/'TRC'
      INTEGER                         , INTENT(in)           ::  ktra    ! tracer index
      INTEGER                         , INTENT(in)           ::  ktrd    ! tracer trend index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)           ::  ptrd    ! tracer trend  or flux
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::  pun     ! velocity 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::  ptra    ! Tracer variablea
      !
      REAL(wp), POINTER, DIMENSION(:,:,:)  :: zwt, zws, ztrdt, ztrds     ! 3D workspace
      !!----------------------------------------------------------------------

      CALL wrk_alloc( jpi, jpj, jpk, ztrds )

      IF( .NOT. ALLOCATED( trdtx ) ) THEN       ! allocate trdtra arrays
         IF( trd_tra_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trd_tra : unable to allocate arrays' )
      ENDIF
      
      ! Control of optional arguments
      IF( ctype == 'TRA' .AND. ktra == jp_tem ) THEN   !==  Temperature trend  ==!
      !   IF( PRESENT( ptra ) ) THEN    
            SELECT CASE( ktrd )            ! shift depending on the direction
            CASE( jptra_trd_xad )  ;  CALL trd_tra_adv( ptrd, pun, ptra, 'X', trdtx ) 
            CASE( jptra_trd_yad )  ;  CALL trd_tra_adv( ptrd, pun, ptra, 'Y', trdty ) 
            CASE( jptra_trd_zad )  ;  CALL trd_tra_adv( ptrd, pun, ptra, 'Z', trdt  ) 
            CASE( jptra_trd_bbc,   &       ! qsr, bbc: on temperature only, send to trd_tra_mng
               &  jptra_trd_qsr )  ;  trdt(:,:,:) = ptrd(:,:,:) * tmask(:,:,:)
                                      ztrds(:,:,:) = 0.
                                      CALL trd_tra_mng( trdt, ztrds, ktrd, kt   )
                                      CALL trd_mod( trdt, ztrds, ktrd, ctype, kt )
            CASE DEFAULT
               trdt(:,:,:) = ptrd(:,:,:) * tmask(:,:,:)
            END SELECT
      !   ELSE
      !   END IF
      END IF

      IF( ctype == 'TRA' .AND. ktra == jp_sal ) THEN 
      !   IF( PRESENT( ptra ) ) THEN    
            SELECT CASE( ktrd )       ! shift depending on the direction
            !                         ! advection: transform the advective flux into a trend
            !                         !            and send T & S trends to trd_tra_mng
            CASE( jptra_trd_xad )  
                                CALL trd_tra_adv( ptrd, pun, ptra, 'X', ztrds ) 
                                CALL trd_tra_mng( trdtx, ztrds, ktrd, kt   )
                                CALL trd_mod( trdtx, ztrds, ktrd, ctype, kt   )
            CASE( jptra_trd_yad )  
                                CALL trd_tra_adv( ptrd, pun, ptra, 'Y', ztrds ) 
                                CALL trd_tra_mng( trdty, ztrds, ktrd, kt   )
                                CALL trd_mod( trdty, ztrds, ktrd, ctype, kt   )
            CASE( jptra_trd_zad )    
                                CALL trd_tra_adv( ptrd, pun, ptra, 'Z', ztrds ) 
                                CALL trd_tra_mng( trdt , ztrds, ktrd, kt   )
                                CALL trd_mod( trdt , ztrds, ktrd, ctype, kt   )
            CASE( jptra_trd_zdfp )       ! diagnose the "PURE" Kz trend (here: just before the swap)
               !                         ! iso-neutral diffusion case otherwise jptra_zdf is "PURE"
               CALL wrk_alloc( jpi, jpj, jpk, zwt, zws, ztrdt )
               !
               zwt(:,:, 1 ) = 0._wp   ;   zws(:,:, 1 ) = 0._wp            ! vertical diffusive fluxes
               zwt(:,:,jpk) = 0._wp   ;   zws(:,:,jpk) = 0._wp
               DO jk = 2, jpk
                  zwt(:,:,jk) =   avt(:,:,jk) * ( tsa(:,:,jk-1,jp_tem) - tsa(:,:,jk,jp_tem) ) / fse3w(:,:,jk) * tmask(:,:,jk)
                  zws(:,:,jk) = fsavs(:,:,jk) * ( tsa(:,:,jk-1,jp_sal) - tsa(:,:,jk,jp_sal) ) / fse3w(:,:,jk) * tmask(:,:,jk)
               END DO
               !
               ztrdt(:,:,jpk) = 0._wp   ;   ztrds(:,:,jpk) = 0._wp
               DO jk = 1, jpkm1
                  ztrdt(:,:,jk) = ( zwt(:,:,jk) - zwt(:,:,jk+1) ) / fse3t(:,:,jk)
                  ztrds(:,:,jk) = ( zws(:,:,jk) - zws(:,:,jk+1) ) / fse3t(:,:,jk) 
               END DO
               CALL trd_tra_mng( ztrdt, ztrds, jptra_trd_zdfp, kt )  
               !
               CALL wrk_dealloc( jpi, jpj, jpk, zwt, zws, ztrdt )
               !
            CASE DEFAULT                 ! other trends: mask and send T & S trends to trd_tra_mng
       !  ELSE
            ztrds(:,:,:) = ptrd(:,:,:) * tmask(:,:,:)
            !IF( ktrd == jptra_trd_ldf .OR. ktrd == jptra_trd_zdf ) THEN
              CALL trd_tra_mng( trdt , ztrds, ktrd, kt   )
            !END IF
            CALL trd_mod( trdt, ztrds, ktrd, ctype, kt )  
            END SELECT
       !  END IF
      END IF

      IF( ctype == 'TRC' ) THEN
         !
         IF( PRESENT( ptra ) ) THEN  
            SELECT CASE( ktrd )            ! shift depending on the direction
            CASE( jptra_trd_xad )  
                                CALL trd_tra_adv( ptrd, pun, ptra, 'X', ztrds ) 
                                CALL trd_mod_trc( ztrds, ktra, ktrd, kt       )
            CASE( jptra_trd_yad )  
                                CALL trd_tra_adv( ptrd, pun, ptra, 'Y', ztrds ) 
                                CALL trd_mod_trc( ztrds, ktra, ktrd, kt       )
            CASE( jptra_trd_zad )    
                                CALL trd_tra_adv( ptrd, pun, ptra, 'Z', ztrds ) 
                                CALL trd_mod_trc( ztrds, ktra, ktrd, kt       )
            END SELECT
         ELSE
            ztrds(:,:,:) = ptrd(:,:,:)
            CALL trd_mod_trc( ztrds, ktra, ktrd, kt       )  
         END IF
         !
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, jpk, ztrds )
      !
   END SUBROUTINE trd_tra


   SUBROUTINE trd_tra_adv( pf, pun, ptn, cdir, ptrd )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra_adv  ***
      !! 
      !! ** Purpose :   transformed the i-, j- or k-advective flux into thes
      !!              i-, j- or k-advective trends, resp.
      !! ** Method  :   i-advective trends = -un. di-1[T] = -( di-1[fi] - tn di-1[un] )
      !!                k-advective trends = -un. di-1[T] = -( dj-1[fi] - tn dj-1[un] )
      !!                k-advective trends = -un. di+1[T] = -( dk+1[fi] - tn dk+1[un] )
      !!----------------------------------------------------------------------
      REAL(wp)        , INTENT(in ), DIMENSION(jpi,jpj,jpk) ::   pf      ! advective flux in one direction
      REAL(wp)        , INTENT(in ), DIMENSION(jpi,jpj,jpk) ::   pun     ! now velocity  in one direction
      REAL(wp)        , INTENT(in ), DIMENSION(jpi,jpj,jpk) ::   ptn     ! now or before tracer 
      CHARACTER(len=1), INTENT(in )                         ::   cdir    ! X/Y/Z direction
      REAL(wp)        , INTENT(out), DIMENSION(jpi,jpj,jpk) ::   ptrd    ! advective trend in one direction
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ii, ij, ik   ! index shift function of the direction
      REAL(wp) ::   zbtr         ! local scalar
      !!----------------------------------------------------------------------

      SELECT CASE( cdir )            ! shift depending on the direction
      CASE( 'X' )   ;   ii = 1   ; ij = 0   ;   ik = 0      ! i-advective trend
      CASE( 'Y' )   ;   ii = 0   ; ij = 1   ;   ik = 0      ! j-advective trend
      CASE( 'Z' )   ;   ii = 0   ; ij = 0   ;   ik =-1      ! k-advective trend
      END SELECT

      !                              ! set to zero uncomputed values
      ptrd(jpi,:,:) = 0.e0   ;   ptrd(1,:,:) = 0.e0
      ptrd(:,jpj,:) = 0.e0   ;   ptrd(:,1,:) = 0.e0
      ptrd(:,:,jpk) = 0.e0
      !
      !
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbtr    = 1.e0/ ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
               ptrd(ji,jj,jk) = - zbtr * (      pf (ji,jj,jk) - pf (ji-ii,jj-ij,jk-ik)                    &
                 &                          - ( pun(ji,jj,jk) - pun(ji-ii,jj-ij,jk-ik) ) * ptn(ji,jj,jk)  )&
                 &                     * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE trd_tra_adv

   SUBROUTINE trd_tra_mng( ptrdx, ptrdy, ktrd, kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra_mng  ***
      !!
      !! ** Purpose :   Dispatch all tracer trends computation, e.g. 3D output,
      !!                integral constraints, potential energy, and/or
      !!                mixed layer budget.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdx   ! Temperature or U trend
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdy   ! Salinity    or V trend
      INTEGER                   , INTENT(in   ) ::   ktrd    ! tracer trend index
      INTEGER                   , INTENT(in   ) ::   kt      ! time step
      !!----------------------------------------------------------------------

      !Duo IF( neuler == 0 .AND. kt == nit000    ) THEN   ;   r2dt =      rdt      ! = rdtra (restart with Euler time stepping)
      !Duo ELSEIF(               kt <= nit000 + 1) THEN   ;   r2dt = 2. * rdt      ! = 2 rdttra (leapfrog)
      !Duo ENDIF

      !                   ! 3D output of tracers trends using IOM interface
      IF( l_trdtra )          CALL trd_tra_iom ( ptrdx, ptrdy, ktrd, kt )

      !                   ! Integral Constraints Properties for tracers trends                                       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !Duo IF( ln_glo_trd )   CALL trd_glo( ptrdx, ptrdy, ktrd, 'TRA', kt )

      !                   ! Potential ENergy trends
      !Duo IF( ln_PE_trd  )   CALL trd_pen( ptrdx, ptrdy, ktrd, kt, r2dt )

      !                   ! Mixed layer trends for active tracers
      !Duo IF( ln_tra_mxl )   THEN
         !-----------------------------------------------------------------------------------------------
         ! W.A.R.N.I.N.G :
         ! jptra_ldf : called by traldf.F90
         !                 at this stage we store:
         !                  - the lateral geopotential diffusion (here, lateral = horizontal)
         !                  - and the iso-neutral diffusion if activated
         ! jptra_zdf : called by trazdf.F90
         !                 * in case of iso-neutral diffusion we store the vertical diffusion component in the
         !                   lateral trend including the K_z contrib, which will be removed later (see trd_mxl)
         !-----------------------------------------------------------------------------------------------

      !Duo   SELECT CASE ( ktrd )
      !Duo   CASE ( jptra_xad )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_xad, '3D' )   ! zonal    advection
      !Duo   CASE ( jptra_yad )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_yad, '3D' )   ! merid.   advection
      !Duo   CASE ( jptra_zad )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_zad, '3D' )   ! vertical advection
      !Duo   CASE ( jptra_ldf )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_ldf, '3D' )   ! lateral  diffusion
      !Duo   CASE ( jptra_bbl )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_bbl, '3D' )   ! bottom boundary layer
      !Duo   CASE ( jptra_zdf )
      !Duo      IF( ln_traldf_iso ) THEN ; CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_ldf, '3D' )   ! lateral  diffusion (K_z)
      !Duo      ELSE                   ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_zdf, '3D' )   ! vertical diffusion (K_z)
      !Duo      ENDIF
      !Duo   CASE ( jptra_dmp )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_dmp, '3D' )   ! internal 3D restoring (tradmp)
      !Duo   CASE ( jptra_qsr )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_for, '3D' )   ! air-sea : penetrative sol radiat
      !Duo   CASE ( jptra_nsr )        ;   ptrdx(:,:,2:jpk) = 0._wp   ;   ptrdy(:,:,2:jpk) = 0._wp
      !Duo                                 CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_for, '2D' )   ! air-sea : non penetr sol radiation
      !Duo   CASE ( jptra_bbc )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_bbc, '3D' )   ! bottom bound cond (geoth flux)
      !Duo   CASE ( jptra_npc )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_npc, '3D' )   ! non penetr convect adjustment
      !Duo   CASE ( jptra_atf )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_atf, '3D' )   ! asselin time filter (last trend)
                                   !
      !Duo                                 CALL trd_mxl( kt, r2dt )                             ! trends: Mixed-layer (output)
      !Duo   END SELECT
         !
      !Duo ENDIF
      !
   END SUBROUTINE trd_tra_mng
      

   SUBROUTINE trd_tra_iom( ptrdx, ptrdy, ktrd, kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra_iom  ***
      !!
      !! ** Purpose :   output 3D tracer trends using IOM
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdx   ! Temperature or U trend
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdy   ! Salinity    or V trend
      INTEGER                   , INTENT(in   ) ::   ktrd    ! tracer trend index
      INTEGER                   , INTENT(in   ) ::   kt      ! time step
      !!
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      INTEGER ::   ikbu, ikbv   ! local integers
      REAL(wp), POINTER, DIMENSION(:,:)   ::   z2dx, z2dy   ! 2D workspace
      !!----------------------------------------------------------------------
      !
!!gm Rq: mask the trends already masked in trd_tra, but lbc_lnk should probably be added
      !
      SELECT CASE( ktrd )
      CASE( jptra_trd_xad  )   ;   CALL iom_put( "ttrd_xad" , ptrdx )        ! x- horizontal advection
                                   CALL iom_put( "strd_xad" , ptrdy )
      CASE( jptra_trd_yad  )   ;   CALL iom_put( "ttrd_yad" , ptrdx )        ! y- horizontal advection
                                   CALL iom_put( "strd_yad" , ptrdy )
      CASE( jptra_trd_zad  )   ;   CALL iom_put( "ttrd_zad" , ptrdx )        ! z- vertical   advection
                                   CALL iom_put( "strd_zad" , ptrdy )
                               IF( .NOT. lk_vvl ) THEN                   ! cst volume : adv flux through z=0 surface
                                 CALL wrk_alloc( jpi, jpj, z2dx, z2dy )
                                  z2dx(:,:) = wn(:,:,1) * tsn(:,:,1,jp_tem) / fse3t(:,:,1)
                                  z2dy(:,:) = wn(:,:,1) * tsn(:,:,1,jp_sal) / fse3t(:,:,1)
                                  CALL iom_put( "ttrd_sad", z2dx )
                                  CALL iom_put( "strd_sad", z2dy )
                                  CALL wrk_dealloc( jpi, jpj, z2dx, z2dy )
                               ENDIF
      CASE( jptra_trd_ldf  )   ;   CALL iom_put( "ttrd_ldf" , ptrdx )        ! lateral diffusion
                                   CALL iom_put( "strd_ldf" , ptrdy )
      CASE( jptra_trd_zdf  )   ;   CALL iom_put( "ttrd_zdf" , ptrdx )        ! vertical diffusion (including Kz contribution)
                                   CALL iom_put( "strd_zdf" , ptrdy )
      !CASE( jptra_trd_zdfp )   ;   CALL iom_put( "ttrd_zdfp", ptrdx )        ! PURE vertical diffusion (no isoneutral contribution)
      !                             CALL iom_put( "strd_zdfp", ptrdy )
      CASE( jptra_trd_dmp  )   ;   CALL iom_put( "ttrd_dmp" , ptrdx )        ! internal restoring (damping)
                                   CALL iom_put( "strd_dmp" , ptrdy )
      CASE( jptra_trd_bbl  )   ;   CALL iom_put( "ttrd_bbl" , ptrdx )        ! bottom boundary layer
                                   CALL iom_put( "strd_bbl" , ptrdy )
      CASE( jptra_trd_npc  )   ;   CALL iom_put( "ttrd_npc" , ptrdx )        ! static instability mixing
                                   CALL iom_put( "strd_npc" , ptrdy )
      CASE( jptra_trd_nsr  )   ;   CALL iom_put( "ttrd_qns" , ptrdx )        ! surface forcing + runoff (ln_rnf=T)
                                   CALL iom_put( "strd_cdt" , ptrdy )
      CASE( jptra_trd_qsr  )   ;   CALL iom_put( "ttrd_qsr" , ptrdx )        ! penetrative solar radiat. (only on temperature)
      CASE( jptra_trd_bbc  )   ;   CALL iom_put( "ttrd_bbc" , ptrdx )        ! geothermal heating   (only on temperature)
      CASE( jptra_trd_atf  )   ;   CALL iom_put( "ttrd_atf" , ptrdx )        ! asselin time Filter
                                   CALL iom_put( "strd_atf" , ptrdy )
      END SELECT
      !
   END SUBROUTINE trd_tra_iom

#   else
   !!----------------------------------------------------------------------
   !!   Default case :          Dummy module           No trend diagnostics
   !!----------------------------------------------------------------------
   USE par_oce      ! ocean variables trends
CONTAINS
   SUBROUTINE trd_tra( kt, ctype, ktra, ktrd, ptrd, pu, ptra )
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in)           ::  kt      ! time step
      CHARACTER(len=3)                , INTENT(in)           ::  ctype   ! tracers trends type 'TRA'/'TRC'
      INTEGER                         , INTENT(in)           ::  ktra    ! tracer index
      INTEGER                         , INTENT(in)           ::  ktrd    ! tracer trend index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)           ::  ptrd    ! tracer trend 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::  pu      ! velocity 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::  ptra    ! Tracer variable 
      WRITE(*,*) 'trd_3d: You should not have seen this print! error ?', ptrd(1,1,1), ptra(1,1,1), pu(1,1,1),   &
         &                                                               ktrd, ktra, ctype, kt
   END SUBROUTINE trd_tra
#   endif
   !!======================================================================
END MODULE trdtra
