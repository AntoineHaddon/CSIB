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
   USE ice , ONLY: a_i, jpl
   USE par_canoe        ! indices of CanOE model variables, e.g. jrdia: diatoms

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_ini_csib   ! called by trcini.F90 module

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
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kmm     ! time level indices
      INTEGER  ::   ji, jj , jl        ! dummy loop indices
      
      REAL(wp) :: zmax !for diagnostics/debug
      !
      CALL trc_nam_csib
      !
      !                       ! Allocate CSIB arrays
      IF( trc_sms_csib_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_csib: unable to allocate CSIB arrays' )

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_csib: passive tracer unit vector'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
      
      IF( .NOT. ln_rsttr ) THEN
         icedia(:,:,:)=0._wp
         iceno3(:,:,:)=0._wp
         icenh4(:,:,:)=0._wp
         
         ! init from ocean surface concentration
         ! icedia(:,:,:) = tr(:,:,1,jrdia,Kmm)
         DO jl = 1, jpl
            WHERE( a_i(:,:,jl) > 1e-4 )
               iceno3(:,:,jl) = tr(:,:,1,jqno3,Kmm)
               icenh4(:,:,jl) = tr(:,:,1,jrnh4,Kmm)
            END WHERE
         ENDDO

         ! init with constant value where latitude > ...
         ! WHERE( gphit(:,:) > 85._wp )   ;   icedia(:,:,3)=1._wp
         ! ELSEWHERE                     ;   icedia(:,:,3)=0._wp
         ! END WHERE

         ! WHERE( gphit(:,:)>75._wp .AND. gphit(:,:)<80._wp .AND. glamt(:,:)>100._wp .AND. glamt(:,:)<150._wp )  
         !    iceno3(:,:,2)=1._wp
         ! END WHERE
         ! WHERE( gphit(:,:)>70._wp .AND. gphit(:,:)<75._wp .AND. glamt(:,:)>-160._wp .AND. glamt(:,:)<-130._wp )  
         !    icenh4(:,:,2)=1._wp
         ! END WHERE

         icedia_gca(:,:,:) = icedia(:,:,:) * a_i(:,:,:)
         iceno3_gca(:,:,:) = iceno3(:,:,:) * a_i(:,:,:)
         icenh4_gca(:,:,:) = icenh4(:,:,:) * a_i(:,:,:)
      ENDIF

      ! initialize fluxes and process rates
      flushrate(:,:,:) = 0._wp
      bogup(:,:,:) = 0._wp
      lagup(:,:,:) = 0._wp

      flush_dia(:,:,:) = 0._wp
      lamloss_dia(:,:,:) = 0._wp
      bogup_dia(:,:,:) = 0._wp
      lagup_dia(:,:,:) = 0._wp
      
      growth_dia(:,:,:) = 0._wp
      lim_lig(:,:,:) = 0._wp

      flush_no3(:,:,:) = 0._wp
      lamloss_no3(:,:,:) = 0._wp
      moldif_no3(:,:,:) = 0._wp
      lagup_no3(:,:,:) = 0._wp

      flush_nh4(:,:,:) = 0._wp
      lamloss_nh4(:,:,:) = 0._wp
      moldif_nh4(:,:,:) = 0._wp
      lagup_nh4(:,:,:) = 0._wp
      
      fric_vel(:,:) = 0._wp


       ! For debug/diagnostics: print max 
      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) 'init, max N hemisphere : '
      DO jl = 1, jpl
         zmax = MAXVAL( iceno3(:,:,jl), MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,jl) > 1e-4 )
         CALL mpp_max( "trc_sms_csib", zmax )
         IF(lwp) WRITE(numout,*) 'ice category ', jl , ' no3_i : ' , zmax

         zmax = MAXVAL( tr(:,:,1,jqno3,Kmm) - iceno3(:,:,jl), MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,jl) > 1e-4 )
         CALL mpp_max( "trc_sms_csib", zmax )
         IF(lwp) WRITE(numout,*) 'no3_o(Kmm) - no3_i : ' , zmax
      ENDDO
      zmax = MAXVAL( tr(:,:,1,jqno3,Kmm) , MASK= gphit(:,:) > 0._wp .AND. a_i(:,:,jl) > 1e-4 )
      CALL mpp_max( "trc_sms_csib", zmax )
      IF(lwp) WRITE(numout,*) 'no3_o(Kmm): ' , zmax

      !
   END SUBROUTINE trc_ini_csib

   !!======================================================================
END MODULE trcini_csib
