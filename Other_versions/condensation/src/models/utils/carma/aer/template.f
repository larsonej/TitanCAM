       subroutine metric
c
c
c  @(#) metric.f  Ackerman  Jul-1997
c
c    <pvap> is vapor pressure in units of [dyne/cm^2]
c
c  Uses temperature <t> as input.
c
c  Argument list input:
c
c    igas: gas identifier
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter metric'
c
c-------------------------------------------------------------------------------
c
c
c  Scale 3-D fields defined at midpoint of vertical layers using metric factors.
c   Rule of thumb:  Any quantity with unit of distance factor ds{x,y,z} should
c   have a factor of 1./{x,y,z}met.
c
      do ixyz=1,NXYZ
        u3(ixyz) = u3(ixyz) / xmet3(ixyz)
        v3(ixyz) = v3(ixyz) / ymet3(ixyz)
        rhoa3(ixyz) = rhoa3(ixyz) * xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz)
        dkx3(ixyz) = dkx3(ixyz) / ( xmet3(ixyz)**2 )
        dky3(ixyz) = dky3(ixyz) / ( ymet3(ixyz)**2 )
      enddo
c
c
c  Scale 3-D fields defined at boundaries of vertical layers using metric factors.
c   Linearly interpolate/extrapolate nearest 2 vertical midpoint metric factors.
c
      do k=1,NZP1
        if( k .eq. 1 )then
         k1 = 1
         k2 = 2
        else if( k .eq. NZP1 )then
         k1 = NZ - 1
         k2 = NZ
        else
         k1 = k - 1
         k2 = k
        endif
        do ixy=1,NXY
          a_mid_k1 = zc2(ixy,k1) 
          a_mid_k2 = zc2(ixy,k2) 
          frac = ( zl2(ixy,k) - a_mid_k1 ) / ( a_mid_k2 - a_mid_k1 ) 

          xmet_k1 = xmet2(ixy,k1)
          xmet_k2 = xmet2(ixy,k2)
          xmet_k = xmet_k1 + frac * ( xmet_k2 - xmet_k1 )
          dkx2(ixy,k) = dkx2(ixy,k) / ( xmet_k**2 )

          ymet_k1 = ymet2(ixy,k1)
          ymet_k2 = ymet2(ixy,k2)
          ymet_k = ymet_k1 + frac * ( ymet_k2 - ymet_k1 )
          dky2(ixy,k) = dky2(ixy,k) / ( ymet_k**2 )

          zmet_k1 = zmet2(ixy,k1)
          zmet_k2 = zmet2(ixy,k2)
          zmet_k = zmet_k1 + frac * ( zmet_k2 - zmet_k1 )
          dkz2(ixy,k) = dkz2(ixy,k) / ( zmet_k**2 )
          w2(ixy,k) = w2(ixy,k) / zmet_k
        enddo
      enddo
c
c
c  Return to caller with 3-D fields scaled.
c
      return
      end
