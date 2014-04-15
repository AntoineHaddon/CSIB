MODULE mt19937ar
   !!======================================================================
   !!                       ***  MODULE  mt19937ar ***
   !! Mersen Twister random generator
   !!======================================================================
   !!                 
   !! ** Purpose :
   !! ** Method  :  
   !! References : This a a Fortran90 translation of mt19937ar.c:
   !! A C-program for MT19937, with initialization improved 2002/1/26.
   !! Coded by Takuji Nishimura and Makoto Matsumoto.
   !!
   !! Before using, initialize the state by using init_genrand(seed)  
   !! or init_by_array(init_key, key_length).
   !!
   !! Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   !! All rights reserved.                          
   !!
   !! Redistribution and use in source and binary forms, with or without
   !! modification, are permitted provided that the following conditions
   !! are met:
   !!
   !!   1. Redistributions of source code must retain the above copyright
   !!      notice, this list of conditions and the following disclaimer.
   !!
   !!   2. Redistributions in binary form must reproduce the above copyright
   !!      notice, this list of conditions and the following disclaimer in the
   !!      documentation and/or other materials provided with the distribution.
   !!
   !!   3. The names of its contributors may not be used to endorse or promote 
   !!      products derived from this software without specific prior written 
   !!      permission.
   !!
   !! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   !! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   !! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   !! A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   !! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   !! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   !! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   !! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   !! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   !! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   !! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   !!
   !!
   !! Any feedback is very welcome.
   !! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   !! email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
   !!
   !! Until 2001/4/6, MT had been distributed under GNU Public License, 
   !! but after 2001/4/6, we decided to let MT be used for any purpose, 
   !! including commercial use. 2002-versions mt19937ar.c, mt19937ar-cok.c 
   !! are considered to be usable freely. 
   !!      
   !! History: 
   !!        !  97 - 02 (Makoto Matsumoto/Takuji Nishimura) original c version
   !!        !  10-03 (A. Vidard) F90 translation and NEMOTAM adaptation
   !-----------------------------------------------------------------------!
   USE par_kind, ONLY : wp
   PRIVATE
   INTEGER, PARAMETER :: jpn = 624           ! period parameters
   INTEGER, PARAMETER :: jpm = 397

   INTEGER*8, PARAMETER :: jpa   = Z'9908b0df' ! constant vector A
   INTEGER*8, PARAMETER :: jpmsb = Z'80000000' ! most significant w-r bit mask
   INTEGER*8, PARAMETER :: jplsb = Z'7fffffff' ! least significant r bit mask
!
   INTEGER*8, PARAMETER :: jpall = Z'ffffffff' ! all bits mask

   INTEGER*8, PARAMETER :: jptmp1= Z'9d2c5680' ! Tampering paramters
   INTEGER*8, PARAMETER :: jptmp2= Z'efc60000'

   INTEGER :: mti
   INTEGER*8, DIMENSION(0:jpn-1) :: mt
   LOGICAL :: l_mtinit = .FALSE.

   PUBLIC &
      & l_mtinit,   &
      & mtrand_int31, &
      & mtrand_int32, &
      & mtrand_real1, &
      & mtrand_real2, &
      & mtrand_real3, &
      & mtrand_real53,&
      & init_by_array,&
      & init_mtrand
CONTAINS
   SUBROUTINE init_mtrand(ks)
!-----------------------------------------------------------------------
!!    initialize mt(0:N-1) with a seed
!-----------------------------------------------------------------------
      INTEGER :: ks
      !
      mt(0) = IAND(INT(ks,kind=8),jpall)
      DO mti = 1, jpn-1
         mt(mti) = 1812433253 * IEOR(mt(mti-1),ISHFT(mt(mti-1),-30)) + mti
         ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. 
         ! In the previous versions, MSBs of the seed affect   
         ! only MSBs of the array mt[].                        
         ! modified by Makoto Matsumoto 
         mt(mti) = IAND(mt(mti),jpall)
         ! for >32 bit machines
      END DO
      l_mtinit = .TRUE.
!
   END SUBROUTINE init_mtrand
   SUBROUTINE init_by_array(kinit_key,key_length)
!-----------------------------------------------------------------------
!     initialize by an array with array-length
!     init_key is the array for initializing keys
!     key_length is its length
!-----------------------------------------------------------------------
      INTEGER :: kinit_key(0:*)
      INTEGER :: key_length
      INTEGER :: ji,jj,jk
!
      CALL init_mtrand(19650218)
      ji = 1
      jj = 0
      DO jk = MAX(jpn,key_length),1,-1
         mt(ji) = IEOR(mt(ji),IEOR(mt(ji-1),ISHFT(mt(ji-1),-30)) * 1664525) &
            &   + kinit_key(jj) + jj
         mt(ji) = IAND(mt(ji),jpall)
         ji = ji + 1
         jj = jj + 1
         IF (ji .GE. jpn) THEN
            mt(0) = mt(jpn-1)
            ji = 1
         ENDIF
         IF (jj .GE. key_length) THEN
            jj = 0
         ENDIF
      END DO
      DO  jk = jpn-1,1,-1
         mt(ji) = IEOR(mt(ji),IEOR(mt(ji-1),ISHFT(mt(ji-1),-30)) * 1566083941) - ji ! non linear
         mt(ji) = IAND(mt(ji),jpall) ! for WORDSIZE > 32 machines
         ji = ji + 1
         IF (ji .GE. jpn) THEN
            mt(0) = mt(jpn-1)
            ji = 1
         ENDIF
      END DO
      mt(0) = jpmsb ! MSB is 1; assuring non-zero initial array
!
   END SUBROUTINE init_by_array
   FUNCTION mtrand_int32()
!-----------------------------------------------------------------------
!     generates a random number on [0,0xffffffff]-interval
!-----------------------------------------------------------------------
      INTEGER :: mtrand_int32
      INTEGER :: jk
      INTEGER*8 :: iy
      INTEGER*8, DIMENSION(0:1) :: mag01
!
      mag01(0)    = 0
      mag01(1)    = jpa          ! mag01[x] = x * vector A  for x=0,1
!
      IF(.NOT. l_mtinit)THEN     ! if init_mtrand() has not been called
         CALL init_mtrand(5489)  ! a default initial seed is used
      ENDIF
!
      IF (mti .GE. jpn) THEN
         DO jk = 0,jpn-jpm-1
            iy     = IOR(IAND(mt(jk),jpmsb),IAND(mt(jk+1),jplsb))
            mt(jk) = IEOR(IEOR(mt(jk+jpm),ISHFT(iy,-1)),mag01(IAND(iy,INT(1,kind=8))))
         END DO
         DO jk = jpn-jpm,jpn-1-1
            iy     = IOR(IAND(mt(jk),jpmsb),IAND(mt(jk+1),jplsb))
            mt(jk) = IEOR(IEOR(mt(jk+(jpm-jpn)),ISHFT(iy,-1)),mag01(IAND(iy,INT(1,kind=8))))
         END DO
         iy     = IOR(IAND(mt(jpn-1),jpmsb),IAND(mt(0),jplsb))
         mt(jk) = IEOR(IEOR(mt(jpm-1),ISHFT(iy,-1)),mag01(IAND(iy,INT(1,kind=8))))
         mti    = 0
      ENDIF
!
      iy  = mt(mti)
      mti = mti + 1
!
      iy = IEOR(iy,ISHFT(iy,-11))
      iy = IEOR(iy,IAND(ISHFT(iy,7),jptmp1))
      iy = IEOR(iy,IAND(ISHFT(iy,15),jptmp2))
      iy = IEOR(iy,ISHFT(iy,-18))
      !
      mtrand_int32 = iy
   END FUNCTION mtrand_int32
   FUNCTION mtrand_int31()
!-----------------------------------------------------------------------
!     generates a random number on [0,0x7fffffff]-interval
!-----------------------------------------------------------------------
      INTEGER :: mtrand_int31
      mtrand_int31 = INT(ISHFT(mtrand_int32(),-1))
   END FUNCTION mtrand_int31
   FUNCTION mtrand_real1()
!-----------------------------------------------------------------------
!     generates a random number on [0,1]-real-interval
!-----------------------------------------------------------------------
      REAL(wp) :: mtrand_real1, zr
      zr = REAL(mtrand_int32(),wp)
      IF (zr.LT.0._wp) zr = zr + 2._wp**32
      mtrand_real1 = zr / 4294967295._wp
   END FUNCTION mtrand_real1
   FUNCTION mtrand_real2()
!-----------------------------------------------------------------------
!     generates a random number on [0,1)-real-interval
!-----------------------------------------------------------------------
      REAL(wp) :: mtrand_real2,zr
      zr = REAL(mtrand_int32(),wp)
      IF (zr.LT.0._wp) zr = zr + 2._wp**32
      mtrand_real2=zr / 4294967296._wp
   END FUNCTION mtrand_real2
   FUNCTION mtrand_real3()
!-----------------------------------------------------------------------
!     generates a random number on (0,1)-real-interval
!-----------------------------------------------------------------------
      REAL(wp) :: mtrand_real3, zr
      zr = REAL(mtrand_int32(),wp)
      IF (zr.LT.0._wp) zr = zr + 2._wp**32
      mtrand_real3 = ( zr + 0.5_wp ) / 4294967296._wp
   END FUNCTION mtrand_real3
   FUNCTION mtrand_res53()
!-----------------------------------------------------------------------
!     generates a random number on [0,1) with 53-bit resolution
!-----------------------------------------------------------------------
      REAL(wp) :: mtrand_res53
      REAL(wp) :: za,zb
      za = ISHFT(mtrand_int32(),-5)
      zb = ISHFT(mtrand_int32(),-6)
      IF (za.LT.0._wp) za = za + 2._wp**32
      IF (zb.LT.0._wp) zb = zb + 2._wp**32
      mtrand_res53 = ( za * 67108864._wp + zb ) / 9007199254740992._wp
   END FUNCTION mtrand_res53
END MODULE mt19937ar
