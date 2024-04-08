module ccc_nemo_rtd_utils
      implicit none
      private
      public area_ave, area_ave_flx, moc, noleap_days

contains

      !=========================================================
      ! Area averaging of 3d field over the selected regions 
      !=========================================================
      SUBROUTINE area_ave(e1,e2,e3, mask, a, imt,jmt,km,a_mean,ss,kk)
            implicit none
            integer imt, jmt, km, i, j, kk
            real e1(imt,jmt),e2(imt,jmt), e3(imt,jmt,km)
            real a(imt,jmt), mask(imt,jmt)
            real a_mean, ss, s1, vol

                s1=0.
                ss=0.
                do i=1,imt
                    do j=1,jmt 
                        if (mask(i,j).gt.0.5) then  ! mask the region of interst
                            vol = e1(i,j)*e2(i,j)*e3(i,j,kk)
                            ss=ss+vol
                            s1=s1+a(i,j)*vol
                        endif
                    enddo
                enddo

                a_mean = 0.
                if (ss.ne.0.) then
                    a_mean =s1/ss
                endif

            return
      END SUBROUTINE area_ave

      !=========================================================
      ! Area averaging of 2d field over the selected regions 
      !=========================================================
      SUBROUTINE area_ave_flx(e1,e2, mask, a, imt, jmt,a_mean,ss)
            implicit none
            integer imt, jmt, i, j
            real e1(imt,jmt),e2(imt,jmt)
            real a(imt,jmt), mask(imt,jmt)
            real a_mean, ss, s1, arc

                s1=0.
                ss=0.
                do i=1,imt
                    do j=1,jmt
                            arc = e1(i,j)*e2(i,j)*mask(i,j)
                            ss=ss+arc
                            s1=s1+a(i,j)*arc
                    enddo
                enddo

                a_mean = 0.
                if (ss.ne.0.) then
                    a_mean =s1/ss
                endif

            return
      END SUBROUTINE area_ave_flx

      !=========================================================
      ! Global meridional overturning (Sv) (valid only south of 20N)  
      !=========================================================
      SUBROUTINE moc(e1v, e3v, mask, v, imt, jmt, km, over_psi)
            implicit none
            integer imt, jmt, km, i, j, k
            REAL, DIMENSION(imt, jmt) :: e1v
            REAL, DIMENSION(imt, jmt, km) :: e3v, v, mask
            REAL, DIMENSION(jmt, km) :: over_tran, over_psi
            REAL s


            do j = 1, jmt
                do k = km, 1, -1
                    s=0.
                    do i = 1, imt
                      if (mask(i,j,k).gt.0.5) then  ! mask the region of interst
                        s = s + v(i, j, k)*e1v(i, j)*e3v(i, j, k) 
                      endif  
                    enddo
                    over_tran(j, k) = s
                    over_psi(j, k)  = 0.
                enddo
            enddo

            do j = 1, jmt
                do k = km, 1, -1
                    if (k.eq.km) then
                        over_psi(j,k) = -over_tran(j,k)
                    else
                        over_psi(j,k) = over_psi(j, k + 1) - over_tran(j,k)
                    endif
                enddo
            enddo

            do j = 1, jmt
                do k = 1, km
                    over_psi(j, k)  =  over_psi(j, k)*1.e-6 ! to Sv                      
                enddo
            enddo

            return
      END SUBROUTINE moc

      SUBROUTINE noleap_days(year, mon, day, days_elapsed)
      !    Given a year, mon, day, returns the number of days elapsed
      !    since 01-01-0000 (using a noleap/365_day calendar)
          IMPLICIT NONE
          INTEGER, INTENT(IN)  :: year, mon, day
          INTEGER, INTENT(OUT) :: days_elapsed
          INTEGER              :: doy, dpy, mmon, myear

          ! Adjust if mon > 12. This is a potential here. No correction
          ! for days being off. Im assuming day will mostly be 1 anyway.

          IF (mon > 12) THEN
              myear = year + mon/12
              mmon = MOD(mon, 12)
          ELSE
             myear = year
             mmon = mon
          ENDIF

          ! In a given (noleap) year, compute the day of year    
          ! http://www.davidgsimpson.com/software/greg2doy_f90.txt
          !
          doy = ((275*mmon)/9) - ((mmon+9)/6) + day - 30

          ! Compute the number of days in the preceeding years
          dpy = 365 * myear

          ! Tally for the final result
          days_elapsed = dpy + doy
          return
      END SUBROUTINE noleap_days


end module ccc_nemo_rtd_utils
