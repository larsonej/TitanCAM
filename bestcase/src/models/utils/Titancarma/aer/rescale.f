       subroutine rescale
c
c
c  @(#) rescale.f  Ackerman  Aug-1997
c
c  This routine rescales vertical winds and diffusion coefficients,
c  and fall velocities.
c
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter rescale'
c
c-------------------------------------------------------------------------------
c
c
c  Rescale fields defined at boundaries of vertical layers using metric factors.
c   Linearly interpolate/extrapolate nearest 2 vertical midpoint metric factors.
c
      do k=1,NZP1

        if( k .eq. 1 )then
         k1 = 1
         k2 = min( NZ, 2 )
        else if( k .eq. NZP1 )then
         k1 = NZ
         k2 = NZ
        else
         k1 = min( NZ, max( 1, k - 1 ) )
	 k2 = k
        endif

        do ixy=1,NXY

          a_mid_k1 = zc2(ixy,k1)
          a_mid_k2 = zc2(ixy,k2)
 
          if(  a_mid_k2 .ne. a_mid_k1 )then
            frac = ( zl2(ixy,k) - a_mid_k1 ) / ( a_mid_k2 - a_mid_k1 )
          else
            frac = 0.
          endif

          zmet_k1 = zmetl2(ixy,k1)
          zmet_k2 = zmetl2(ixy,k2)
          zm_k_old = zmet_k1 + frac * ( zmet_k2 - zmet_k1 )

          zmet_k1 = zmet2(ixy,k1)
          zmet_k2 = zmet2(ixy,k2)
          zm_k_new = zmet_k1 + frac * ( zmet_k2 - zmet_k1 )
   
          scalef = zm_k_old / zm_k_new

          dkz2(ixy,k) = dkz2(ixy,k) * scalef**2
          w2(ixy,k) = w2(ixy,k) * scalef

        enddo
c
c
c  No horizontal variation for <vf>: use column with ixy = NXY
c
        do ig = 1,NGROUP
          do ibin = 1,NBIN
            vf(k,ibin,ig) = vf(k,ibin,ig) * scalef
          enddo
        enddo

      enddo
c
c
c  Return to caller with winds and diffusion coefficients scaled.
c
      return
      end
