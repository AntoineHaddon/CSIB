MODULE trcini_csib
   !!======================================================================
   !!                         ***  MODULE trcini_csib  ***
   !! TOP :   initialisation of the CSIB tracers
   !!======================================================================
   !! History :      !  2025 (A. Haddon) Original code
   !!----------------------------------------------------------------------
   !! trc_ini_csib   : CSIB model initialisation
   !!----------------------------------------------------------------------
   USE par_kind   !: access wp kind
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE par_csib
   USE trcnam_csib     ! csib SMS namelist
   USE trcsms_csib

   USE dom_oce, ONLY: glamt, gphit               ! latitude/longitude for funky initiation
   USE ice , ONLY: a_i, jpl, t_i, nlay_i
   
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_csib         ! called by trcini.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcini_my_trc.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_csib
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_csib  ***  
      !!
      !! ** Purpose :   initialization for CSIB model
      !!
      !! ** Method  : - Read the csib namelist and check the parameter values
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl,jn          ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z3d   ! 3D workspace
      CHARACTER(len=25) ::   znam

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_csib:'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
            
      CALL trc_nam_csib ! read namelist 
      
      ! Set number of tracers
      IF (ln_dmsice) THEN
         jp_csib=7
      ELSE
         jp_csib=5
      END IF
      
      ! Allocate sms_CSIB arrays
      IF( trc_sms_csib_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_csib: unable to allocate CSIB arrays' )

      ! Init of sea ice BGC variables
      icetra(:,:,:,:) = 0._wp

      IF( ln_rstart ) THEN ! if restart 

         IF( iom_varid( numrir, icetrcnm(1), ldstop = .FALSE. ) > 0 ) THEN ! check that ice restart file has CSIB variables
            IF(lwp) WRITE(numout,*) '  read sea ice tracer variables (CSIB) from ice restart file'
            DO jn = 1,jp_csib
               znam = icetrcnm(jn)
               CALL iom_get( numrir, jpdom_auto, znam , z3d ) ; icetra(:,:,:,jn) = z3d(:,:,:)
            ENDDO
         ELSE ! restart but no CSIB variables in ice restart
            IF(lwp) WRITE(numout,*) '  no sea ice tracers (CSIB) in ice restart file - set equal to ocean surface '
            ln_ibgcspinup=.true. ! Init is done in trcsms_csib because CanOE init is done after CSIB init
         ENDIF
         
      ELSE ! Full model spin up
         
         IF(lwp) WRITE(numout,*) '  Full model spin up - set equal to ocean surface '
         ln_ibgcspinup=.true. ! Init is done in trcsms_csib because CanOE init is done after CSIB init
         
         ! init with constant value where latitude > ...
         ! WHERE( gphit(:,:) > 85._wp )   ;   icetra(:,:,3,jridia)=1._wp
         ! END WHERE
         
         ! WHERE( gphit(:,:)>75._wp .AND. gphit(:,:)<80._wp .AND. glamt(:,:)>100._wp .AND. glamt(:,:)<150._wp )  
         !    iceno3(:,:,2)=1._wp
         ! END WHERE
         ! WHERE( gphit(:,:)>70._wp .AND. gphit(:,:)<75._wp .AND. glamt(:,:)>-160._wp .AND. glamt(:,:)<-130._wp )  
         !    icenh4(:,:,2)=1._wp
         ! END WHERE
            
      ENDIF
      
      DO jn = 1,jp_csib
         icetra_gca(:,:,:,jn) = icetra(:,:,:,jn) * a_i(:,:,:)
      ENDDO
            
      ! initialize ratios, fluxes and process rates
      qnidia(:,:,:) = 0._wp
      qchidia(:,:,:) = 0._wp
      qnidiamax(:,:,:) = 0._wp

      flushrate(:,:,:) = 0._wp
      bogup(:,:,:) = 0._wp
      lagup(:,:,:) = 0._wp

      flush_dia(:,:,:) = 0._wp
      slough_dia(:,:,:) = 0._wp
      lamloss_dia(:,:,:) = 0._wp
      t_i_b(:,:,:) = 0._wp
      dt_i(:,:,:) = 0._wp
      meltoff_dia(:,:,:) = 0._wp
      bogup_dia(:,:,:) = 0._wp
      lagup_dia(:,:,:) = 0._wp
      nxsicedia(:,:,:) = 0._wp
      cxsicedia(:,:,:) = 0._wp
      
      flush_no3(:,:,:) = 0._wp
      slough_no3(:,:,:) = 0._wp
      lamloss_no3(:,:,:) = 0._wp
      moldif_no3(:,:,:) = 0._wp
      lagup_no3(:,:,:) = 0._wp
      bogup_no3(:,:,:) = 0._wp
      
      flush_nh4(:,:,:) = 0._wp
      slough_nh4(:,:,:) = 0._wp
      lamloss_nh4(:,:,:) = 0._wp
      moldif_nh4(:,:,:) = 0._wp
      lagup_nh4(:,:,:) = 0._wp
      bogup_nh4(:,:,:) = 0._wp
      
      fric_vel(:,:) = 0._wp
      
      phot_dia(:,:,:) = 0._wp
      lim_PAR(:,:,:) = 0._wp
      lim_nut(:,:,:) = 0._wp
      lim_ice(:,:,:) = 0._wp
      diaupn(:,:,:) = 0._wp
      chlsyn(:,:,:) = 0._wp
      mortlin_dia(:,:,:) = 0._wp
      mortquad_dia(:,:,:) = 0._wp
      remin_dia(:,:,:) = 0._wp
      nitri(:,:,:) = 0._wp

      ! init sea ice temp at previous time step with current temp: used for sea ice temp change, init at 0 can cause large derivative
       DO jj = 1, jpj
         DO ji = 1, jpi
            DO jl = 1, jpl 
               t_i_b(ji,jj,jl) = SUM(t_i(ji,jj,:,jl)) / nlay_i  
            ENDDO
         ENDDO
      ENDDO

      IF (ln_dmsice) THEN
         flush_dmspd(:,:,:) = 0._wp
         slough_dmspd(:,:,:) = 0._wp
         lamloss_dmspd(:,:,:) = 0._wp
         lagup_dmspd(:,:,:) = 0._wp
         bogup_dmspd(:,:,:) = 0._wp
         bogup_dmspd(:,:,:) = 0._wp
         lagup_dmspd(:,:,:) = 0._wp
         
         flush_dms(:,:,:)   =0._wp
         slough_dms(:,:,:) = 0._wp
         lamloss_dms(:,:,:) = 0._wp
         lagup_dms(:,:,:) = 0._wp
         bogup_dms(:,:,:) = 0._wp
         bogup_dms(:,:,:) = 0._wp
         lagup_dms(:,:,:) = 0._wp

         dmsp_exud(:,:,:) = 0._wp
         dmsp_lysis(:,:,:) = 0._wp
         dms_phot(:,:,:) = 0._wp
      ENDIF

   END SUBROUTINE trc_ini_csib


   !!======================================================================
END MODULE trcini_csib
