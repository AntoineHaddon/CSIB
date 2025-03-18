MODULE sbcssm
   !!======================================================================
   !!                       ***  MODULE  sbcssm  ***
   !! Off-line : interpolation of the physical fields
   !!======================================================================
   !! History :  3.4  ! 2012-03 (S. Alderson)  original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_ssm_init  : initialization, namelist read, and SAVEs control
   !!   sbc_ssm       : Interpolation of the fields
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE c1d            ! 1D configuration: ln_c1d
   USE dom_oce        ! ocean domain: variables
   USE zdf_oce        ! ocean vertical physics: variables
   USE sbc_oce        ! surface module: variables
   USE phycst         ! physical constants
   USE eosbn2         ! equation of state - Brunt Vaisala frequency
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE zpshde         ! z-coord. with partial steps: horizontal derivatives
   USE closea         ! for ln_closea
   USE icb_oce        ! for icebergs
#if defined key_si3
   USE ice             ! sea-ice: variables
   USE icevar          ! sea-ice: operations
   USE icecor          ! sea-ice: corrections
#endif
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lib_mpp        ! distributed memory computing library
   USE prtctl         ! print control
   USE fldread        ! read input fields
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ssm_init       ! called by sbc_init
   PUBLIC   sbc_ssm            ! called by sbc
   PUBLIC   sbc_ssm_ice_init   ! called by sbc_init
   PUBLIC   sbc_ssm_ice        ! called by sbc

   CHARACTER(len=100) ::   cn_dir        ! Root directory for location of ssm files
   LOGICAL            ::   ln_3d_uve     ! specify whether input velocity data is 3D
   LOGICAL            ::   ln_read_frq   ! specify whether we must read frq or not

   LOGICAL            ::   l_sasread     ! Ice intilisation: =T read a file ; =F anaytical initilaistion
   LOGICAL            ::   l_initdone = .false.
   INTEGER     ::   nfld_3d
   INTEGER     ::   nfld_2d
   INTEGER     ::   nfld_ice

   INTEGER     ::   jf_tem         ! index of temperature
   INTEGER     ::   jf_sal         ! index of salinity
   INTEGER     ::   jf_usp         ! index of u velocity component
   INTEGER     ::   jf_vsp         ! index of v velocity component
   INTEGER     ::   jf_ssh         ! index of sea surface height
   INTEGER     ::   jf_e3t         ! index of first T level thickness
   INTEGER     ::   jf_frq         ! index of fraction of qsr absorbed in the 1st T level

   INTEGER     ::   jf_ifr         ! index of ice fraction
   INTEGER     ::   jf_ims         ! index of ice mass
   INTEGER     ::   jf_tic         ! index of ice surface temperature
   INTEGER     ::   jf_ial         ! index of sea-ice albedo (not implemented yet)

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_ssm_3d  ! structure of input fields (file information, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_ssm_2d  ! structure of input fields (file information, fields read)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_ssm_ice  ! structure of input fields (file information, fields read)

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/SAS 4.0 , NEMO Consortium (2018)
   !! $Id: sbcssm.F90 15023 2021-06-18 14:35:25Z gsamson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_ssm( kt, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ssm  ***
      !!
      !! ** Purpose :  Prepares dynamics and physics fields from a NEMO run
      !!               for an off-line simulation using surface processes only
      !!
      !! ** Method : calculates the position of data
      !!             - interpolates data if needed
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm   ! ocean time level indices
      ! (not needed for SAS but needed to keep a consistent interface in sbcmod.F90)
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) ::   ztinta     ! ratio applied to after  records when doing time interpolation
      REAL(wp) ::   ztintb     ! ratio applied to before records when doing time interpolation
      REAL(wp), DIMENSION(jpi,jpj)     ::  sstfrz
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start( 'sbc_ssm')

      IF ( l_sasread ) THEN
         IF( nfld_3d > 0 ) CALL fld_read( kt, 1, sf_ssm_3d )      !==   read data at kt time step   ==!
         IF( nfld_2d > 0 ) CALL fld_read( kt, 1, sf_ssm_2d )      !==   read data at kt time step   ==!
         !
         IF( ln_3d_uve ) THEN
            IF( .NOT. ln_linssh ) THEN
               e3t_m(:,:) = sf_ssm_3d(jf_e3t)%fnow(:,:,1) * tmask(:,:,1) ! vertical scale factor
            ELSE
               e3t_m(:,:) = e3t_0(:,:,1)                                 ! vertical scale factor
            ENDIF
            ssu_m(:,:) = sf_ssm_3d(jf_usp)%fnow(:,:,1) * umask(:,:,1)    ! u-velocity
            ssv_m(:,:) = sf_ssm_3d(jf_vsp)%fnow(:,:,1) * vmask(:,:,1)    ! v-velocity
         ELSE
            IF( .NOT. ln_linssh ) THEN
               e3t_m(:,:) = sf_ssm_2d(jf_e3t)%fnow(:,:,1) * tmask(:,:,1) ! vertical scale factor
            ELSE
               e3t_m(:,:) = e3t_0(:,:,1)                                 ! vertical scale factor
            ENDIF
            IF( TRIM(sf_ssm_2d(jf_usp)%clrootname) == 'NOT USED' ) &
               &     sf_ssm_2d(jf_usp)%fnow(:,:,1) = 0._wp
            IF( TRIM(sf_ssm_2d(jf_vsp)%clrootname) == 'NOT USED' ) &
               &     sf_ssm_2d(jf_vsp)%fnow(:,:,1) = 0._wp
            ssu_m(:,:) = sf_ssm_2d(jf_usp)%fnow(:,:,1) * umask(:,:,1)    ! u-velocity
            ssv_m(:,:) = sf_ssm_2d(jf_vsp)%fnow(:,:,1) * vmask(:,:,1)    ! v-velocity
         ENDIF
         !
         IF( TRIM(sf_ssm_2d(jf_sal)%clrootname) == 'NOT USED' ) &
            &     sf_ssm_2d(jf_sal)%fnow(:,:,1) = 33.252_wp
         IF( TRIM(sf_ssm_2d(jf_tem)%clrootname) == 'NOT USED' ) &
            &     CALL eos_fzp( sf_ssm_2d(jf_sal)%fnow(:,:,1), sf_ssm_2d(jf_tem)%fnow(:,:,1) )
         IF( TRIM(sf_ssm_2d(jf_ssh)%clrootname) == 'NOT USED' ) &
            &     sf_ssm_2d(jf_ssh)%fnow(:,:,1) = 0._wp
         sst_m(:,:) = sf_ssm_2d(jf_tem)%fnow(:,:,1) * tmask(:,:,1)    ! temperature
         sss_m(:,:) = sf_ssm_2d(jf_sal)%fnow(:,:,1) * tmask(:,:,1)    ! salinity
         ssh_m(:,:) = sf_ssm_2d(jf_ssh)%fnow(:,:,1) * tmask(:,:,1)    ! sea surface height
         IF( ln_read_frq ) THEN
            frq_m(:,:) = sf_ssm_2d(jf_frq)%fnow(:,:,1) * tmask(:,:,1) ! solar penetration
         ELSE
            frq_m(:,:) = 1._wp
         ENDIF
      ELSE
         sss_m(:,:) = 33.252_wp                          ! =33.252 to obtain a physical value for the freezing point of 271.2K (teos10)
         CALL eos_fzp( sss_m(:,:), sst_m(:,:) )          ! sst_m is set at the freezing point
         ssu_m(:,:) = 0._wp
         ssv_m(:,:) = 0._wp
         ssh_m(:,:) = 0._wp
         e3t_m(:,:) = e3t_0(:,:,1)                       !clem: necessary at least for sas2D
         frq_m(:,:) = 1._wp                              !              - -
         ssh  (:,:,Kmm) = 0._wp                          !              - -
      ENDIF
      CALL eos_fzp( sss_m(:,:), sstfrz(:,:) )          ! set min sst_m to the freezing point
      WHERE(sst_m(:,:).le. sstfrz(:,:) ) sst_m(:,:)=sstfrz(:,:)

      IF ( nn_ice == 1.or.ln_cpl ) THEN
         ts(:,:,1,jp_tem,Kmm) = sst_m(:,:)
         ts(:,:,1,jp_sal,Kmm) = sss_m(:,:)
         ts(:,:,1,jp_tem,Kbb) = sst_m(:,:)
         ts(:,:,1,jp_sal,Kbb) = sss_m(:,:)
      ENDIF
      uu (:,:,1,Kbb) = ssu_m(:,:)
      vv (:,:,1,Kbb) = ssv_m(:,:)

      IF(sn_cfctl%l_prtctl) THEN            ! print control
         CALL prt_ctl(tab2d_1=sst_m, clinfo1=' sst_m   - : ', mask1=tmask   )
         CALL prt_ctl(tab2d_1=sss_m, clinfo1=' sss_m   - : ', mask1=tmask   )
         CALL prt_ctl(tab2d_1=ssu_m, clinfo1=' ssu_m   - : ', mask1=umask   )
         CALL prt_ctl(tab2d_1=ssv_m, clinfo1=' ssv_m   - : ', mask1=vmask   )
         CALL prt_ctl(tab2d_1=ssh_m, clinfo1=' ssh_m   - : ', mask1=tmask   )
         IF( .NOT.ln_linssh )   CALL prt_ctl(tab2d_1=ssh_m, clinfo1=' e3t_m   - : ', mask1=tmask   )
         IF( ln_read_frq    )   CALL prt_ctl(tab2d_1=frq_m, clinfo1=' frq_m   - : ', mask1=tmask   )
      ENDIF
      !
      IF( l_initdone ) THEN          !   Mean value at each nn_fsbc time-step   !
         CALL iom_put( 'ssu_m', ssu_m )
         CALL iom_put( 'ssv_m', ssv_m )
         CALL iom_put( 'sst_m', sst_m )
         CALL iom_put( 'sss_m', sss_m )
         CALL iom_put( 'ssh_m', ssh_m )
         CALL iom_put( 'e3t_m', e3t_m )
         IF( ln_read_frq    )   CALL iom_put( 'frq_m', frq_m )
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop( 'sbc_ssm')
      !
   END SUBROUTINE sbc_ssm

   SUBROUTINE sbc_ssm_ice( kt, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ssm_ice  ***
      !!
      !! ** Purpose :  Prepares ice fields from a NEMO run
      !!               for an off-line simulation using surface processes only
      !!
      !! ** Method : calculates the position of data
      !!             - interpolates data if needed
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm   ! ocean time level indices
      ! (not needed for SAS but needed to keep a consistent interface in sbcmod.F90)
      !
      INTEGER  ::   ji, jj, jk, jl     ! dummy loop indices
      REAL(wp) ::   ztinta     ! ratio applied to after  records when doing time interpolation
      REAL(wp) ::   ztintb     ! ratio applied to before records when doing time interpolation
      REAL(wp) ::   ztmelts
      REAL(wp), DIMENSION(jpi,jpj)     ::  at_i_read, vt_i_read
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start( 'sbc_ssm_ice')

      IF ( l_sasread ) THEN
         IF( nfld_ice > 0 ) CALL fld_read( kt, 1, sf_ssm_ice )      !==   read data at kt time step   ==!
         !

                                     ! 1. -- Change the category concentrations according to the input (will be rebin later)
         IF( TRIM(sf_ssm_ice(jf_tic)%clrootname) /= 'NOT USED' ) THEN
            t_su(:,:,:) = 0.
            t_su (:,:,1) = sf_ssm_ice(jf_tic)%fnow(:,:,1)
         ENDIF
         IF( TRIM(sf_ssm_ice(jf_ifr)%clrootname) /= 'NOT USED' ) THEN 
             at_i_read(:,:) = sf_ssm_ice(jf_ifr)%fnow(:,:,1)
         ELSE ! if ice volume not read, the volume does not change
             at_i_read = at_i
         ENDIF
         IF( TRIM(sf_ssm_ice(jf_ims)%clrootname) /= 'NOT USED' ) THEN
            vt_i_read = sf_ssm_ice(jf_ims)%fnow(:,:,1) / rhoi
         ELSE ! if ice volume not read, the volume does not change
            vt_i_read = vt_i
         ENDIF

         ! limit the input sea-ice concentration to the rn_max
         WHERE(at_i_read(:,:).gt. rn_amax_2d(:,:) ) at_i_read(:,:)=rn_amax_2d(:,:)
         ! defined the sea ice variables for different cases. 
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
             IF (at_i(ji,jj).gt.epsi20) THEN ! change in concentration
                 a_i(ji,jj,:)=a_i(ji,jj,:)*at_i_read(ji,jj)/at_i(ji,jj)
                 IF (vt_i(ji,jj).gt.epsi20) THEN !
                     v_i(ji,jj,:)=v_i(ji,jj,:)*vt_i_read(ji,jj)/vt_i(ji,jj)
                 ELSE
                     v_i(ji,jj  ,:) = 0._wp ; v_i(ji,jj  ,1) = max(vt_i_read(ji,jj),ht_i_new(ji,jj)*at_i_read(ji,jj)  )
                 ENDIF
             ELSEIF (at_i_read(ji,jj).gt.epsi20) THEN !new ice (put in first category) 
                 a_i(ji,jj  ,:) = 0._wp ; a_i(ji,jj  ,1) = at_i_read(ji,jj)
                 v_i(ji,jj  ,:) = 0._wp ; v_i(ji,jj  ,1) = max(vt_i_read(ji,jj) ,ht_i_new(ji,jj)*at_i_read(ji,jj) )
                 h_s(ji,jj,:)   = 0._wp
                 t_s(ji,jj,:,:) = rt0 
                 t_i(ji,jj,:,:) = rt0 
                 t_su(ji,jj ,:) = rt0 
                 s_i (ji,jj ,:) = rn_simin 
                 o_i (ji,jj ,:) = 0._wp
             ELSE
                 a_i(ji,jj  ,:) = 0._wp 
                 v_i(ji,jj  ,:) = 0._wp
                 h_s(ji,jj,:)   = 0._wp
                 o_i (ji,jj ,:) = 0._wp
             ENDIF
             v_s (ji,jj ,:) = h_s(ji,jj,:) * a_i(ji,jj,:)
             sv_i(ji,jj ,:) = s_i(ji,jj,:) * v_i(ji,jj,:)
         END_2D

         CALL ice_cor( kt , 0 )      ! 2. -- Check for thickness <rn_himin  and >rn_amax
                                     ! 3. -- Rebin categories with thickness out of bounds     
                                     ! 4. -- Check for salinity in bounds [Simin,Simax] 
         DO jl  = 1, jpl             ! 5. -- Re-calculate the enthalpy (snow & ice)
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                DO jk = 1, nlay_s           
                    t_s(ji,jj,jk,jl) = MIN( t_s(ji,jj,jk,jl), -0.15_wp + rt0 )           ! Force t_s to be lower than -0.15deg (arbitrary) => likely conservation issue
                    !                                                                    !       otherwise instant melting can occur
                    e_s(ji,jj,jk,jl) = rhos * ( rcpi * ( rt0 - t_s(ji,jj,jk,jl) ) + rLfus )   ! enthalpy in J/m3
                    e_s(ji,jj,jk,jl) = e_s(ji,jj,jk,jl) * v_s(ji,jj,jl) * r1_nlay_s           ! enthalpy in J/m2
                END DO               
            END_2D
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               t_su(ji,jj,jl) = MIN( t_su(ji,jj,jl), -0.15_wp + rt0 )                  ! Force t_su to be lower than -0.15deg (arbitrary)
               DO jk = 1, nlay_i
                    ztmelts          = - rTmlt  * sz_i(ji,jj,jk,jl)             ! Melting temperature in C
                    t_i(ji,jj,jk,jl) = MIN( t_i(ji,jj,jk,jl), (ztmelts-0.15_wp) + rt0 )  ! Force t_i to be lower than melting point (-0.15) => likely conservation issue
                    !
                    e_i(ji,jj,jk,jl) = rhoi * ( rcpi  * ( ztmelts - ( t_i(ji,jj,jk,jl) - rt0 ) )           &   ! enthalpy in J/m3
                       &                      + rLfus * ( 1._wp - ztmelts / ( t_i(ji,jj,jk,jl) - rt0 ) )   &
                       &                      - rcp   *   ztmelts )                  
                    e_i(ji,jj,jk,jl) = e_i(ji,jj,jk,jl) * v_i(ji,jj,jl) * r1_nlay_i                            ! enthalpy in J/m2
               END DO
            END_2D
         END DO               
         CALL ice_var_agg(1)         ! 6. -- integrate variables over layers and categories post inputs
      ENDIF


      IF(sn_cfctl%l_prtctl) THEN            ! print control
         CALL prt_ctl(tab3d_1=a_i , clinfo1=' a_i     - : ', kdim=jpl      )
         CALL prt_ctl(tab3d_1=t_su, clinfo1=' t_su    - : ', kdim=jpl      )
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop( 'sbc_ssm_ice')
      !
   END SUBROUTINE sbc_ssm_ice

   SUBROUTINE sbc_ssm_init( Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ssm_init  ***
      !!
      !! ** Purpose :   Initialisation of sea surface mean data
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Kmm   ! ocean time level indices
      ! (not needed for SAS but needed to keep a consistent interface in sbcmod.F90)
      INTEGER  :: ierr, ierr0, ierr1, ierr2, ierr3   ! return error code
      INTEGER  :: ifpr                               ! dummy loop indice
      INTEGER  :: inum, idv, idimv, jpm              ! local integer
      INTEGER  ::   ios                              ! Local integer output status for namelist read
      !!
      CHARACTER(len=100)                     ::  cn_dir       ! Root directory for location of core files
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::  slf_3d       ! array of namelist information on the fields to read
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::  slf_2d       ! array of namelist information on the fields to read
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::  slf_ice      ! array of namelist information on the fields to read
      TYPE(FLD_N) ::   sn_tem, sn_sal                     ! information about the fields to be read
      TYPE(FLD_N) ::   sn_usp, sn_vsp
      TYPE(FLD_N) ::   sn_ssh, sn_e3t, sn_frq
      !!
      TYPE(FLD_N) ::   sn_ifr, sn_ims, sn_tic, sn_ial
      !!
      NAMELIST/namsbc_sas/ l_sasread, cn_dir, ln_3d_uve, ln_read_frq,   &
         &                 sn_tem, sn_sal, sn_usp, sn_vsp, sn_ssh, sn_e3t, sn_frq, &
         &                 sn_ifr, sn_ims, sn_tic, sn_ial
      !!----------------------------------------------------------------------
      !
      IF( ln_rstart .AND. nn_components == jp_iam_sas )   RETURN
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_ssm_init : sea surface mean data initialisation '
         WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF
      !
      READ  ( numnam_ref, namsbc_sas, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_sas in reference namelist' )
      READ  ( numnam_cfg, namsbc_sas, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_sas in configuration namelist' )
      IF(lwm) WRITE ( numond, namsbc_sas )
      !
      IF(lwp) THEN                              ! Control print
         WRITE(numout,*) '   Namelist namsbc_sas'
         WRITE(numout,*) '      Initialisation using an input file                                 l_sasread   = ', l_sasread
         WRITE(numout,*) '      Are we supplying a 3D u,v and e3 field                             ln_3d_uve   = ', ln_3d_uve
         WRITE(numout,*) '      Are we reading frq (fraction of qsr absorbed in the 1st T level)   ln_read_frq = ', ln_read_frq
      ENDIF
      !
      !! switch off stuff that isn't sensible with a standalone module
      !! note that we need sbc_ssm called first in sbc
      !
      IF( ln_apr_dyn ) THEN
         IF( lwp ) WRITE(numout,*) '         ==>>>   No atmospheric gradient needed with StandAlone Surface scheme'
         ln_apr_dyn = .FALSE.
      ENDIF
      IF( ln_rnf ) THEN
         IF( lwp ) WRITE(numout,*) '         ==>>>   No runoff needed with StandAlone Surface scheme'
         ln_rnf = .FALSE.
      ENDIF
      IF( ln_ssr ) THEN
         IF( lwp ) WRITE(numout,*) '         ==>>>   No surface relaxation needed with StandAlone Surface scheme'
         ln_ssr = .FALSE.
      ENDIF
      IF( nn_fwb > 0 ) THEN
         IF( lwp ) WRITE(numout,*) '         ==>>>   No freshwater budget adjustment needed with StandAlone Surface scheme'
         nn_fwb = 0
      ENDIF
      IF( ln_closea ) THEN
         IF( lwp ) WRITE(numout,*) '         ==>>>   No closed seas adjustment needed with StandAlone Surface scheme'
         ln_closea = .false.
      ENDIF
      !
      IF( l_sasread ) THEN                       ! store namelist information in an array
         !
         !! following code is a bit messy, but distinguishes between when u,v are 3d arrays and
         !! when we have other 3d arrays that we need to read in
         !! so if a new field is added i.e. jf_new, just give it the next integer in sequence
         !! for the corresponding dimension (currently if ln_3d_uve is true, 4 for 2d and 3 for 3d,
         !! alternatively if ln_3d_uve is false, 6 for 2d and 1 for 3d), reset nfld_3d, nfld_2d,
         !! and the rest of the logic should still work
         !
         jf_tem = 1   ;   jf_ssh = 3   ! default 2D fields index
         jf_sal = 2   ;   jf_frq = 4   !
         
         !
         IF( ln_3d_uve ) THEN
            jf_usp = 1   ;   jf_vsp = 2   ;   jf_e3t = 3     ! define 3D fields index
            nfld_3d  = 2 + COUNT( (/.NOT.ln_linssh/) )       ! number of 3D fields to read
            nfld_2d  = 3 + COUNT( (/ln_read_frq/) )          ! number of 2D fields to read
         ELSE
            jf_usp = 4   ;   jf_e3t = 6                      ! update 2D fields index
            jf_vsp = 5   ;   jf_frq = 6 + COUNT( (/.NOT.ln_linssh/) )
            !
            nfld_3d  = 0                                     ! no 3D fields to read
            nfld_2d  = 5 + COUNT( (/.NOT.ln_linssh/) ) + COUNT( (/ln_read_frq/) )    ! number of 2D fields to read
         ENDIF
#if defined key_si3
         jf_ifr =  1   ;   jf_ims = 2 ; jf_tic =  3   ;  jf_ial =  4  ! Sea-Ice 2D fields
         nfld_ice = 4
#else
         jf_ifr = -1 ;   jf_ims = -1;   jf_tic = -1 ;  jf_ial = -1 ! Sea-Ice 2D fields (dummy value to avoid bad matching)
         nfld_ice = 0
#endif
         !
         IF( nfld_3d > 0 ) THEN
            ALLOCATE( slf_3d(nfld_3d), STAT=ierr )         ! set slf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'sbc_ssm_init: unable to allocate slf 3d structure' )   ;   RETURN
            ENDIF
            slf_3d(jf_usp) = sn_usp
            slf_3d(jf_vsp) = sn_vsp
            IF( .NOT.ln_linssh )   slf_3d(jf_e3t) = sn_e3t
         ENDIF
         !
         IF( nfld_2d > 0 ) THEN
            ALLOCATE( slf_2d(nfld_2d), STAT=ierr )         ! set slf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'sbc_ssm_init: unable to allocate slf 2d structure' )   ;   RETURN
            ENDIF
            slf_2d(jf_tem) = sn_tem   ;   slf_2d(jf_sal) = sn_sal   ;   slf_2d(jf_ssh) = sn_ssh
            IF( ln_read_frq )   slf_2d(jf_frq) = sn_frq
            IF( .NOT. ln_3d_uve ) THEN
               slf_2d(jf_usp) = sn_usp ; slf_2d(jf_vsp) = sn_vsp
               IF( .NOT.ln_linssh )   slf_2d(jf_e3t) = sn_e3t
            ENDIF
         ENDIF
         !
#if defined key_si3
         IF( nfld_ice > 0 ) THEN ! sf_ssm_ice filled here, but allocation done in sbc_ssm_ice_init
            ALLOCATE( slf_ice(nfld_ice), STAT=ierr )         ! set slf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'sbc_ssm_init: unable to allocate slf 2d structure' )   ;   RETURN
            ENDIF
            slf_ice(jf_ifr) = sn_ifr   ; slf_ice(jf_ims) = sn_ims   ;   slf_ice(jf_tic) = sn_tic   ;   slf_ice(jf_ial) = sn_ial
            ALLOCATE( sf_ssm_ice(nfld_ice), STAT=ierr )         ! set sf structure
            CALL fld_fill( sf_ssm_ice, slf_ice, cn_dir, 'sbc_ssm_init', 'Ice Data in file', 'namsbc_ssm' )
         ENDIF
#endif
         !
         ierr1 = 0    ! default definition if slf_?d(ifpr)%ln_tint = .false.
         IF( nfld_3d > 0 ) THEN
            ALLOCATE( sf_ssm_3d(nfld_3d), STAT=ierr )         ! set sf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'sbc_ssm_init: unable to allocate sf structure' )   ;   RETURN
            ENDIF
            DO ifpr = 1, nfld_3d
               ALLOCATE( sf_ssm_3d(ifpr)%fnow(jpi,jpj,jpk)    , STAT=ierr0 )
               IF( slf_3d(ifpr)%ln_tint )   ALLOCATE( sf_ssm_3d(ifpr)%fdta(jpi,jpj,jpk,2)  , STAT=ierr1 )
               IF( ierr0 + ierr1 > 0 ) THEN
                  CALL ctl_stop( 'sbc_ssm_init : unable to allocate sf_ssm_3d array structure' )   ;   RETURN
               ENDIF
            END DO
            !                                         ! fill sf with slf_i and control print
            CALL fld_fill( sf_ssm_3d, slf_3d, cn_dir, 'sbc_ssm_init', '3D Data in file', 'namsbc_ssm' )
            sf_ssm_3d(jf_usp)%cltype = 'U'   ;   sf_ssm_3d(jf_usp)%zsgn = -1._wp
            sf_ssm_3d(jf_vsp)%cltype = 'V'   ;   sf_ssm_3d(jf_vsp)%zsgn = -1._wp
         ENDIF
         !
         IF( nfld_2d > 0 ) THEN
            ALLOCATE( sf_ssm_2d(nfld_2d), STAT=ierr )         ! set sf structure
            IF( ierr > 0 ) THEN
               CALL ctl_stop( 'sbc_ssm_init: unable to allocate sf 2d structure' )   ;   RETURN
            ENDIF
            DO ifpr = 1, nfld_2d
               ALLOCATE( sf_ssm_2d(ifpr)%fnow(jpi,jpj,1)    , STAT=ierr0 )
               IF( slf_2d(ifpr)%ln_tint )   ALLOCATE( sf_ssm_2d(ifpr)%fdta(jpi,jpj,1,2)  , STAT=ierr1 )
               IF( ierr0 + ierr1 > 0 ) THEN
                  CALL ctl_stop( 'sbc_ssm_init : unable to allocate sf_ssm_2d array structure' )   ;   RETURN
               ENDIF
            END DO
            !
            CALL fld_fill( sf_ssm_2d, slf_2d, cn_dir, 'sbc_ssm_init', '2D Data in file', 'namsbc_ssm' )
            IF( .NOT. ln_3d_uve ) THEN
               sf_ssm_2d(jf_usp)%cltype = 'U'   ;   sf_ssm_2d(jf_usp)%zsgn = -1._wp
               sf_ssm_2d(jf_vsp)%cltype = 'V'   ;   sf_ssm_2d(jf_vsp)%zsgn = -1._wp
            ENDIF
         ENDIF
         !
         IF( nfld_3d > 0 )   DEALLOCATE( slf_3d, STAT=ierr )
         IF( nfld_2d > 0 )   DEALLOCATE( slf_2d, STAT=ierr )
         IF( nfld_2d > 0 )   DEALLOCATE( slf_ice, STAT=ierr )
         !
      ENDIF
      !
      CALL sbc_ssm( nit000, Kbb, Kmm )   ! need to define ss?_m arrays used in iceistate
      l_initdone = .TRUE.
      !
   END SUBROUTINE sbc_ssm_init

   SUBROUTINE sbc_ssm_ice_init( Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ssm_init  ***
      !!
      !! ** Purpose :   Initialisation of sea surface mean ice data
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Kmm   ! ocean time level indices
      ! (not needed for SAS but needed to keep a consistent interface in sbcmod.F90)
      INTEGER  :: ierr, ierr0, ierr1, ierr2, ierr3   ! return error code
      INTEGER  :: ifpr                               ! dummy loop indice
      INTEGER  :: inum, idv, idimv, jpm              ! local integer
      INTEGER  ::   ios                              ! Local integer output status for namelist read


      IF( ln_rstart .AND. nn_components == jp_iam_oce )   RETURN
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_ssm_ice_init : sea surface mean ice data initialisation '
         WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF
      !
      IF( l_sasread ) THEN                       ! store namelist information in an array

         IF( nfld_ice > 0 ) THEN 
            DO ifpr = 1, nfld_ice
               ALLOCATE( sf_ssm_ice(ifpr)%fnow(jpi,jpj,1)    , STAT=ierr0 )
               IF( sf_ssm_ice(ifpr)%ln_tint )   ALLOCATE( sf_ssm_ice(ifpr)%fdta(jpi,jpj,1,2)  , STAT=ierr1 )
               IF( ierr0 + ierr1 > 0 ) THEN
                  CALL ctl_stop( 'sbc_ssm_ice_init : unable to allocate sf_ssm_ice array structure' )   ;   RETURN
               ENDIF
            END DO
         ENDIF


      ENDIF
      !
      CALL sbc_ssm_ice( nit000, Kbb, Kmm )   ! need to define ss?_m arrays used in iceistate
      !
   END SUBROUTINE sbc_ssm_ice_init
   !!======================================================================
END MODULE sbcssm
