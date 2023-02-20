MODULE nemo_diag_cal

   !!===============================================================
   !!                ***   MODULE nemo_diag_cal   ***
   !!                 nemo diagnostics calculations
   !!===============================================================
   !! 2018-08 (D. Yang):   Original code 
   !!---------------------------------------------------------------
   !!
   !!---------------------------------------------------------------
   !! cmip6_mfo      : mass transports through sections (mfo)   
   !! cmip6_msftbarot: quasi-barotropic streamfunction (msftbarot)
   !! cmip6_tstend   : T & S tendencies (opottemptend & osalttend)
   !!---------------------------------------------------------------
   USE nemo_diag_glovars

   IMPLICIT NONE
   PRIVATE      
   PUBLIC cmip6_mfo                   ! Called by nemo_diag.F90
   PUBLIC cmip6_mfo_ice               ! Called by nemo_diag.F90
   PUBLIC cmip6_msftbarot             ! Called by nemo_diag.F90
   PUBLIC cmip6_tstend                ! Called by nemo_diag.F90

   REAL, PARAMETER :: rho0 = 1035.    ! LIM (1026. if CICE)

CONTAINS

   SUBROUTINE cmip6_mfo
      !!
      !!-------------------------------------------------------------
      !! Purpose:   Compute mass transports through sections (mfo)
      !! Reference: Griffies et al, 2016: OMIP contribution to CMIP6: 
      !!            experimental and diagnostic protocol for the 
      !!            physical component of the Ocean Model Intercomparison 
      !!            Project, p3269-3270.
      !! Section mask file (mfo_line_mask): Provided by O. Saenko
      !!            mfo(1) - Barents Opening            (pseudo-north) 
      !!            mfo(2) - Bering Strait              (pseudo-north)
      !!            mfo(3) - Canadian Archipelago       (pseudo-south)
      !!            mfo(4) - Caribbean Windward Passage (pseudo-north)
      !!            mfo(5) - Denmark Strait             (pseudo-south)  
      !!            mfo(6) - Drake Passage              (east)  
      !!            mfo(7) - English Channel            (pseudo-north)
      !!            mfo(8) - Faroe-Scotland Channel     (pseudo-east)
      !!            mfo(9) - Florida-Bahamas Strait     (pseudo-east)  
      !!            mfo(10)- Fram Strait                (pseudo-south)  
      !!            mfo(11)- Iceland Faroe Channel      (pseudo-north) 
      !!            mfo(12)- Indonesian Throughflow     (south)   
      !!            mfo(13)- Mozambique Channel         (south)   
      !!            mfo(14)- Pacific Eq. Undercurrent   (east)   
      !!            mfo(15)- Taiwan-Luzon Straits       (*) 
      !!            (*) negative means net inflow to South China Sea 
      !! Input fields: 
      !!            2D: e1v, e2u
      !!            3D: u, v, umask, vmask, e3u, e3v
      !! Output:    mfo(15) in kg/s 
      !!-------------------------------------------------------------
      INTEGER                       :: i, j, k, l, i0   ! dummy loop indices
      INTEGER, DIMENSION (1)        :: ierr             ! local variable
      CHARACTER, DIMENSION(imt,jmt) :: secmask          ! section mask for transports
      !!----------------
      !! Allocate Arrays
      !!----------------
      ALLOCATE( mfo(nline,lm), STAT=ierr(1) )
      mfo(:,:) =0. 

      IF (MAXVAL(ierr) /=0) THEN
         STOP 'Memory allocation error in cmip6_mfo'
      ENDIF
      !!-------------------------------------------------------------
      !! Read in ASCII mfo_line_mask file which defines the sections  
      !!-------------------------------------------------------------
      OPEN (10, file='mfo_line_mask', status='unknown')
      if (jmt.eq.292) then
        i0=1
      elseif (jmt.eq.331) then
        i0=41
      else
        Print*, 'WARNING: mfo_line_mask made for ORCA1 and eORCA1 only for now. mfo is empty.'
        return
      endif
      secmask='#'
      DO j=jmt,41,-1 
         READ(10,'(362A1)')(secmask(i,j),i=1,imt)
      ENDDO
      CLOSE(10)
      !!-----------------------------------------
      !! Compute mass transports through sections
      !!-----------------------------------------
      WRITE(*,*) 'mfo', lm
      DO l=1,lm
         DO j=1,jmt
            DO i=1,imt
               DO k=1,km
                  IF (secmask(i,j).eq.'1') THEN !  Barents opening 
                     mfo(1,l) = mfo(1,l)+v(i,j,k,l)*e1v(i,j)*e3v(i,j,k,l)*vmask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'2') THEN !  Bering Strait  
                     mfo(2,l) = mfo(2,l)+v(i,j,k,l)*e1v(i,j)*e3v(i,j,k,l)*vmask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'3') THEN !  Canadian Archipelago  
                     mfo(3,l) = mfo(3,l)+v(i,j,k,l)*e1v(i,j)*e3v(i,j,k,l)*vmask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'4') THEN ! Caribbean Windward Passage 
                     mfo(4,l) = mfo(4,l)+v(i,j,k,l)*e1v(i,j)*e3v(i,j,k,l)*vmask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'5') THEN ! Denmark Strait  
                     mfo(5,l) = mfo(5,l)+v(i,j,k,l)*e1v(i,j)*e3v(i,j,k,l)*vmask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'6') THEN ! Drake Passage   
                     mfo(6,l) = mfo(6,l)+u(i,j,k,l)*e2u(i,j)*e3u(i,j,k,l)*umask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'7') THEN ! English Channel   
                     mfo(7,l) = mfo(7,l)+v(i,j,k,l)*e1v(i,j)*e3v(i,j,k,l)*vmask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'8') THEN ! Faroe-Scotland Channel    
                     mfo(8,l) = mfo(8,l)+u(i,j,k,l)*e2u(i,j)*e3u(i,j,k,l)*umask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'9') THEN !  Florida-Bahamas Strait   
                     mfo(9,l) = mfo(9,l)+u(i,j,k,l)*e2u(i,j)*e3u(i,j,k,l)*umask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'a') THEN ! Fram Strait  
                     mfo(10,l)= mfo(10,l)+v(i,j,k,l)*e1v(i,j)*e3v(i,j,k,l)*vmask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'b') THEN ! Iceland Faroe Channel    
                     mfo(11,l)= mfo(11,l)+v(i,j,k,l)*e1v(i,j)*e3v(i,j,k,l)*vmask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'c') THEN ! Indonesian Throughflow  
                     mfo(12,l)= mfo(12,l)+v(i,j,k,l)*e1v(i,j)*e3v(i,j,k,l)*vmask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'d') THEN ! Mozambique Channel  
                     mfo(13,l)= mfo(13,l)+v(i,j,k,l)*e1v(i,j)*e3v(i,j,k,l)*vmask(i,j,k)*rho0
                  ENDIF 
                  IF (secmask(i,j).eq.'e'.and.k.le.16) THEN ! Pacific Eq. Undercurrent   
                     IF (u(i,j,k,l).gt.0.) THEN ! only eastward transport 
                        mfo(14,l)= mfo(14,l)+u(i,j,k,l)*e2u(i,j)*e3u(i,j,k,l)*umask(i,j,k)*rho0
                     ENDIF 
                  ENDIF 
                  IF (secmask(i,j).eq.'f'.or.secmask(i,j).eq.'g') THEN 
                     IF (secmask(i,j).eq.'f') THEN ! Taiwan Strait  
                        mfo(15,l)= mfo(15,l)+v(i,j,k,l)*e1v(i,j)*e3v(i,j,k,l)*vmask(i,j,k)*rho0
                     ELSE                      ! Luzon Strait  
                        mfo(15,l)=mfo(15,l)+u(i,j,k,l)*e2u(i,j)*e3u(i,j,k,l)*umask(i,j,k)*rho0
                     ENDIF 
                  ENDIF       
               ENDDO
            ENDDO
         ENDDO
         WRITE(*,*) '-----',l
         WRITE(*,'(15f10.2)')(mfo(i,l)*1.e-9, i=1,15)   !  Sv 
      ENDDO
      ! 
   END SUBROUTINE cmip6_mfo

   SUBROUTINE cmip6_mfo_ice
      !!
      !!-------------------------------------------------------------
      !! Purpose:   Compute sea ice area (siareaacrossline), sea ice mass (simassacrossline) and snow mass () transport throug sections (mfo_ice)
      !! Reference: Notz et al. 2016: The CMIP6 Sea-Ice Model Intercomparison 
      !!            Project (SIMIP): understanding sea ice through climate-model 
      !!            simulations.  doi.org/10.5194/gmd-9-3427-2016
      !! Section mask file (mfo_line_mask): Provided by O. Saenko
      !!            mfo(10)- Fram Strait                (pseudo-south)  
      !!            mfo(3) - Canadian Archipelago       (pseudo-south)
      !!            mfo(1) - Barents Opening            (pseudo-north) 
      !!            mfo(2) - Bering Strait              (pseudo-north)
      !!        ** used the file mfo_line_mask for now, but it needs to
      !!        ** be separate from the physics in the future
      !! Input fields: 
      !!            2D: ymtrpice,ymtrpsnw,yatrp
      !! Output:    simassacrossline,snmassacrossline,siareaacrossline
      !!-------------------------------------------------------------
      INTEGER                       :: i, j, k, l, i0   ! dummy loop indices
      INTEGER, DIMENSION (1)        :: ierr             ! local variable
      CHARACTER, DIMENSION(imt,jmt) :: secmask          ! section mask for transports
      !!----------------
      !! Allocate Arrays
      !!----------------
      ALLOCATE( simassacrossline(nline_ice,lm),snmassacrossline(nline_ice,lm),siareaacrossline(nline_ice,lm), STAT=ierr(1) )

      siareaacrossline(:,:) =0.
      simassacrossline(:,:) =0.
      snmassacrossline(:,:) =0.

      IF (MAXVAL(ierr) /=0) THEN
         STOP 'Memory allocation error in cmip6_mfo'
      ENDIF
      !!-------------------------------------------------------------
      !! Read in ASCII mfo_line_mask file which defines the sections  
      !!-------------------------------------------------------------
      OPEN (10, file='mfo_line_mask', status='unknown')
      if (jmt.eq.292) then
        i0=1
      elseif (jmt.eq.331) then
        i0=41
      else
        Print*, 'WARNING: mfo_line_mask made for ORCA1 and eORCA1 only for now. mfo is empty.'
        return
      endif
      secmask='#'
      DO j=jmt,41,-1 
         READ(10,'(362A1)')(secmask(i,j),i=1,imt)
      ENDDO
      CLOSE(10)
      !!-----------------------------------------
      !! Compute mass transports through sections
      !!-----------------------------------------

      DO l=1,lm
         DO j=1,jmt
            DO i=1,imt
                  IF (secmask(i,j).eq.'1') THEN !  Barents opening 
                     siareaacrossline(3,l) = siareaacrossline(3,l)+yatrp(i,j,l)
                     simassacrossline(3,l) = simassacrossline(3,l)+ymtrpice(i,j,l)
                     snmassacrossline(3,l) = snmassacrossline(3,l)+ymtrpsnw(i,j,l)
                  ENDIF 
                  IF (secmask(i,j).eq.'2') THEN !  Bering Strait  
                     siareaacrossline(4,l) = siareaacrossline(4,l)+yatrp(i,j,l)
                     simassacrossline(4,l) = simassacrossline(4,l)+ymtrpice(i,j,l)
                     snmassacrossline(4,l) = snmassacrossline(4,l)+ymtrpsnw(i,j,l)
                  ENDIF 
                  IF (secmask(i,j).eq.'3') THEN !  Canadian Archipelago  
                     siareaacrossline(2,l) = siareaacrossline(2,l)+yatrp(i,j,l)
                     simassacrossline(2,l) = simassacrossline(2,l)+ymtrpice(i,j,l)
                     snmassacrossline(2,l) = snmassacrossline(2,l)+ymtrpsnw(i,j,l)
                  ENDIF 
                  IF (secmask(i,j).eq.'a') THEN ! Fram Strait  
                     siareaacrossline(1,l) = siareaacrossline(1,l)+yatrp(i,j,l)
                     simassacrossline(1,l) = simassacrossline(1,l)+ymtrpice(i,j,l)
                     snmassacrossline(1,l) = snmassacrossline(1,l)+ymtrpsnw(i,j,l)
                  ENDIF 
            ENDDO
         ENDDO
      ENDDO
      ! 
   END SUBROUTINE cmip6_mfo_ice

   SUBROUTINE cmip6_msftbarot
      !!
      !!-------------------------------------------------------------
      !! Purpose:   Compute quasi-barotropic (mass) streamfunction
      !! Reference: Griffies et al, 2016: OMIP contribution to CMIP6: 
      !!            experimental and diagnostic protocol for the 
      !!            physical component of the Ocean Model Intercomparison 
      !!            Project, p3263-3264, equ.(H46).
      !! Input fields:
      !!            2D: e2u, ssh
      !!            3D: u, e3u, umask
      !! Output:    msftbarot in kg/s
      !! Comment:   Compute on u grid
      !!-------------------------------------------------------------
      INTEGER                       :: i, j, k, l       ! dummy loop indices
      INTEGER, DIMENSION (1)        :: ierr             ! local variable
      REAL                          :: ztmp             ! local scalar
      REAL, DIMENSION(imt,jmt)      :: uzint            ! local array
      !!----------------
      !! Allocate arrays
      !!----------------
      ALLOCATE( msftbarot(imt,jmt,lm), STAT=ierr(1) )

      IF (MAXVAL(ierr) /=0) THEN
         STOP 'Memory allocation error in msftbarot_cmip6'
      ENDIF
      !!---------------------------------------------------------
      !! uzint = vertical integal of u from free srface to bottom
      !!--------------------------------------------------------- 
      DO l=1,lm
         DO j=1,jmt
            DO i=1,imt
               uzint(i,j) = 0. 
               msftbarot(i,j,l) = 0.
               ztmp = 0. 
               DO k=1,km
                  ztmp = ztmp+u(i,j,k,l)*e3u(i,j,k,l)*umask(i,j,k)
               ENDDO
               uzint(i,j) = ztmp + ssh(i,j,l)*u(i,j,1,l)*umask(i,j,1)
            ENDDO
         ENDDO
      !!-----------------------------------------
      !! meridional integal of uzint * e2u * rho0
      !!-----------------------------------------
         DO i=1,imt
            DO j=2,jmt 
               msftbarot(i,j,l) = (msftbarot(i,j-1,l)-uzint(i,j)*e2u(i,j)*rho0)
            ENDDO
         ENDDO 
      ENDDO
      ! msftbarot(144,45,1) = -6.18654e+09 [ kg/s ]
      ! msftbarot(144,45,4) = -3.07254e+09 [ kg/s ]
      ! WRITE(*,*) 'msftbarot'
      ! WRITE(*,*) msftbarot(144,45,1), msftbarot(144,45,4), msftbarot(225,135,1)

   END SUBROUTINE cmip6_msftbarot

   SUBROUTINE cmip6_tstend ( l )
      !!
      !!-------------------------------------------------------------
      !! Purpose: Compute T & S tendencies (opottemptend & osalttend)
      !! Reference: Griffies et al, 2016: OMIP contribution to CMIP6: 
      !!            experimental and diagnostic protocol for the 
      !!            physical component of the Ocean Model Intercomparison 
      !!            Project, p3282.
      !! Input fields:
      !!            3D: tnp, tnc, snp, snc, e3t, tmask
      !! Output:    opottemptend in W m-2 & osalttend in kg m-2 s-1
      !!-------------------------------------------------------------
      INTEGER                 :: l                      ! time index
      INTEGER, DIMENSION (1)  :: ierr                   ! local variable
      REAL                    :: coeft, coefs           ! coefficients
      !!----------------
      !! Allocate Arrays
      !!----------------
      ALLOCATE( opottemptend(imt,jmt,km,ly), osalttend(imt,jmt,km,ly), STAT=ierr(1) )

      IF (MAXVAL(ierr) /=0) THEN
         STOP 'Memory allocation error in cmip6_tstend'
      ENDIF
      !!-------------------------
      !! Compute T & S tendencies
      !!-------------------------
      ! rhozero=1035 kg/m3; cpocean=4000 J/(kgC)
      ! coeft=1035*4000/(365*24*3600)=0.13127854
      coeft = 0.13127854
      ! coefs=1035/(1000*365*24*3600)=3.3*e-8
      coefs = 3.3e-8
      !WRITE(*,*) 'tnp & tnc & snp & snc & e3t, tmask'
      !WRITE(*,*) tnp(144,45,1), tnc(144,45,1), snp(144,45,1), snc(144,45,1), e3t(144,45,1), tmask(144,45,1)
      opottemptend(:,:,:,l) = ( tnc(:,:,:) - tnp(:,:,:) ) * e3t(:,:,:,l*12) * tmask(:,:,:) * coeft
      osalttend(:,:,:,l) = ( snc(:,:,:) - snp(:,:,:) ) * e3t(:,:,:,l*12) * tmask(:,:,:) * coefs
      !WRITE(*,*) 'opottemptend'
      !opottemptend(144,45,1,l)=1.756493546 [ Wm-2 ]
      !WRITE(*,*) opottemptend(144,45,1,l)
      !WRITE(*,*) 'osalttend'
      !osalttend(144,45,1,l)=1.8e-8 [ kg m-2 s-1 ]
      !WRITE(*,*) osalttend(144,45,1,l)

   END SUBROUTINE cmip6_tstend

END MODULE nemo_diag_cal
