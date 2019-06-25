      subroutine setupmprod
c
c
c  @(#) setupmprod.f  Barth  Jun-1999
c  This routine sets up mass production
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
c
c
c  Include global constants and variables.
c
      include 'globaer.h'
c
c
c  Local delclarations
c
c  <pprofile> are mass production profile values
c
      dimension pprofile(60),q(60)
c
c
    1 format(i4,2(3x,1pe13.6))

c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupmprod'
c
c
c  Total mass production rate [g/cm2/s]
c
      tot_rate = 1.225d-14
c
c  Mass production function
c
      z0 = 250.0   ! altitude of maximum mass production [km]
      zs =  25.0   ! scale height of mass production [km]

      z0 = 400.0 

      totalprod = 0.
      do k = 1,NZ
        if(zl3(k) .gt. 350.d5 .and. zl3(k) .lt. 450.d5) then
       !if(zl3(k) .gt. 110.d5 .and. zl3(k) .lt. 390.d5) then
          q(k)  = exp(-0.5*(( zl3(k)/1.d5 - z0 )/zs)**2)
        else
          q(k) = 0.
        endif
        totalprod = totalprod + q(k)
      enddo
c
c  Mass production profile
c
      do k = 1,NZ
       partprod = 0.
       do ki = k,NZ
         partprod = partprod + q(ki)
       enddo
       pprofile(k) = partprod/totalprod
      enddo
 
      do k=1, NZ
       if(k .lt. NZ) pprofile(k) = abs(pprofile(k)-pprofile(k+1))

c  Use this when mass is added in <psolve>
        rprod(k) = tot_rate * pprofile(k) /
     $              (4./3.*PI*r(1,1)**3 *rhoelem(1) *10.d5)

c 		rprod(k) = (3./4./PI/r(1,1)**3) / (zs*1.0d5) 
c    $			   / rhoelem(1) * dtime * ratemp * pprofile(k)

         if(zl3(NZ) .lt. 110.d5) rprod(k) = 0.
       enddo

c  Look at mass production
c      write(*,*) 'Tholin mass production values:',tot_rate,r(1,1),
c    $              rhoelem(1)
c      do k=1,nxyz
c       write(*,1) k,rprod(k),pprofile(k)
c      enddo
c
c  Return to caller with mass production values.
c
      return
      end
