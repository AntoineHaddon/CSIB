MODULE trcini_csib
   !!======================================================================
   !!                         ***  MODULE trcini_csib  ***
   !! TOP :   initialisation of the CSIB tracers
   !!======================================================================
   !! History :        !  2007  (C. Ethe, G. Madec) Original code
   !!                  !  2016  (C. Ethe, T. Lovato) Revised architecture
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
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_csib         ! called by trcini.F90 module
   PUBLIC   trc_ini_csibnames   

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcini_my_trc.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_ini_csib( Kmm )
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_csib  ***  
      !!
      !! ** Purpose :   initialization for CSIB model
      !!
      !! ** Method  : - Read the csib namelist and check the parameter values
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kmm         ! time level indices
      INTEGER  ::   ji, jj, jl,jn          ! dummy loop indices
      
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_csib:'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      
      ! Allocate sms_CSIB arrays
      IF( trc_sms_csib_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_csib: unable to allocate CSIB arrays' )
      
      ! read namelist
      CALL trc_nam_csib

      
      IF( .NOT. ln_rsttr ) THEN
         CALL trc_ini_csibnames() ! when there is a restart this is done by the sea ice model
         
         icetra(:,:,:,:) = 0._wp
         
         ! init with constant value where latitude > ...
         ! WHERE( gphit(:,:) > 85._wp )   ;   icetra(:,:,3,jridia)=1._wp
         ! END WHERE
         
         ! WHERE( gphit(:,:)>75._wp .AND. gphit(:,:)<80._wp .AND. glamt(:,:)>100._wp .AND. glamt(:,:)<150._wp )  
         !    iceno3(:,:,2)=1._wp
         ! END WHERE
         ! WHERE( gphit(:,:)>70._wp .AND. gphit(:,:)<75._wp .AND. glamt(:,:)>-160._wp .AND. glamt(:,:)<-130._wp )  
         !    icenh4(:,:,2)=1._wp
         ! END WHERE
         
         
         DO jn = 1,jp_csib
            icetra_gca(:,:,:,jn) = icetra(:,:,:,jn) * a_i(:,:,:)
         ENDDO
         
      ENDIF
      
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

   END SUBROUTINE trc_ini_csib


   SUBROUTINE trc_ini_csibnames()
      ! allocate and initiate array of ice tracer variables
      ! called if no restart by trc_ini_csib 
      ! or if restart called by sea ice model in icedyn_adv_pra/adv_pra_rst when reading advection moments

      INTEGER :: trc_sms_csib_allocnames = 0
      ALLOCATE( icetrcnm(jp_csib) ,    STAT=trc_sms_csib_allocnames)
      IF( trc_sms_csib_allocnames /= 0 ) CALL ctl_stop( 'STOP', 'trc_sms_csib_allocnames : failed to allocate arrays' )

      icetrcnm(jridiac) = 'icediac'
      icetrcnm(jridian) = 'icedian'
      icetrcnm(jridiach) = 'icediach'
      icetrcnm(jrino3) = 'iceno3'
      icetrcnm(jrinh4) = 'icenh4'

   END SUBROUTINE trc_ini_csibnames

   !!======================================================================
END MODULE trcini_csib
