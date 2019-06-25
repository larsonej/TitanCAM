       subroutine vertadv
c
c
c  @(#) vertadv.f  Jensen  Dec-1996
c  This routine calculates vertrical advection rates using
c  Piecewise Polynomial Method [Colela and Woodard, J. Comp. Phys.,
c  54, 174-201, 1984]
c
c  <vertadvu(k)> is upward vertical transport rate into level k
c                from level k-1 [cm/s].
c  <vertadvd(k)> is downward vertical transport rate into level k-1 from level k.
c
c  Modified  Sep-1997  (McKie)
c  Remove <ixy> from arg list
c  <ixy> now available as a global var in common block.
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
c  Local declarations
c
      dimension dela(NZ), delma(NZ), aju(NZ),
     $  ar(NZ), al(NZ), a6(NZ), cold(NZ)
c
c
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter vertadv'
c
c
c  Store old values of cvert
c
      do k = 1,NZ
        cold(k) = cvert(k)
      enddo
c
c
c  Initialize fluxes to zero
c
      do k = 1,NZ+1
        vertadvu(k) = 0.
        vertadvd(k) = 0.
      enddo
c
c
c  Set some constants
c
      nzm1 = max( 1, NZ-1 )
      nzm2 = max( 1, NZ-2 )
      itwo = min( 2, NZ   )
c
c
c  First, use cubic fits to estimate concentration values at layer
c  boundaries
c
      do k = 2,NZ-1

        dpc = cvert(k) / dz2(ixy,k)
        dpc1 = cvert(k+1) / dz2(ixy,k+1)
        dpcm1 = cvert(k-1) / dz2(ixy,k-1)
        ratt1 = dz2(ixy,k) /
     $    ( dz2(ixy,k-1) + dz2(ixy,k) + dz2(ixy,k+1) )
        ratt2 = ( 2.*dz2(ixy,k-1) + dz2(ixy,k) ) /
     $          ( dz2(ixy,k+1) + dz2(ixy,k) )
        ratt3 = ( 2.*dz2(ixy,k+1) + dz2(ixy,k) ) /
     $          ( dz2(ixy,k-1) + dz2(ixy,k) )
        dela(k) = ratt1 *
     $           ( ratt2*(dpc1-dpc) + ratt3*(dpc-dpcm1) )

        if( (dpc1-dpc)*(dpc-dpcm1) .gt. 0. .and. dela(k) .ne. 0. ) then
          delma(k) = min( abs(dela(k)), 2.*abs(dpc-dpc1),
     $         2.*abs(dpc-dpcm1) ) * abs(dela(k))/dela(k)
        else
          delma(k) = 0.
        endif

      enddo     ! k = 2,NZ-2

      do k = 2,NZ-2

        dpc = cvert(k) / dz2(ixy,k)
        dpc1 = cvert(k+1) / dz2(ixy,k+1)
        dpcm1 = cvert(k-1) / dz2(ixy,k-1)
        rat1 = dz2(ixy,k) /
     $          ( dz2(ixy,k) + dz2(ixy,k+1) )
        rat2 = 2. * dz2(ixy,k+1) * dz2(ixy,k) /
     $         ( dz2(ixy,k) + dz2(ixy,k+1) )
        rat3 = ( dz2(ixy,k-1) + dz2(ixy,k) ) /
     $         ( 2.*dz2(ixy,k) + dz2(ixy,k+1) )
        rat4 = ( dz2(ixy,k+2) + dz2(ixy,k+1) ) /
     $         ( 2.*dz2(ixy,k+1) + dz2(ixy,k) )
        den1 = dz2(ixy,k-1) + dz2(ixy,k) +
     $         dz2(ixy,k+1) + dz2(ixy,k+2)
c
c  <aju(k)> is the estimate for concentration (dn/dz) at layer
c  boundary <k>+1/2.
c
        aju(k) = dpc + rat1*(dpc1-dpc) + 1./den1 *
     $           ( rat2*(rat3-rat4)*(dpc1-dpc) -
     $           dz2(ixy,k)*rat3*delma(k+1) +
     $           dz2(ixy,k+1)*rat4*delma(k) )

      enddo     ! k = 2,NZ-2
c
c
c  Now construct polynomial functions in each layer
c
      do k = 3,NZ-2

        al(k) = aju(k-1)
        ar(k) = aju(k)

      enddo
c
c  Use linear functions in first two and last two layers
c

      ar(itwo) = aju(itwo)
      al(itwo) = cvert(1)/dz2(ixy,1) +
     $        (zl2(ixy,itwo)-zc2(ixy,1)) /
     $        (zc2(ixy,itwo)-zc2(ixy,1)) *
     $        (cvert(itwo)/dz2(ixy,itwo)-
     $        cvert(1)/dz2(ixy,1))
      ar(1) = al(itwo)
      al(1) = cvert(1)/dz2(ixy,1) -
     $        (zc2(ixy,1)-zl2(ixy,1)) /
     $        (zc2(ixy,itwo)-zc2(ixy,1)) *
     $        (cvert(itwo)/dz2(ixy,itwo)-
     $        cvert(1)/dz2(ixy,1))

      al(nzm1) = aju(nzm2)
      ar(nzm1) = cvert(nzm1)/dz2(ixy,nzm1) +
     $        (zl2(ixy,NZ)-zc2(ixy,nzm1))
     $        / (zc2(ixy,NZ)-zc2(ixy,nzm1)) *
     $        (cvert(NZ)/dz2(ixy,NZ)-
     $        cvert(nzm1)/dz2(ixy,nzm1))
      al(NZ) = ar(nzm1)
      ar(NZ) = cvert(nzm1)/dz2(ixy,nzm1) +
     $        (zl2(ixy,NZ+1)-zc2(ixy,nzm1))
     $        / (zc2(ixy,NZ)-zc2(ixy,nzm1)) *
     $        (cvert(NZ)/dz2(ixy,NZ)-
     $        cvert(nzm1)/dz2(ixy,nzm1))
c
c  Ensure that boundary values are not negative
c
      al(1) = max( al(1), 0.*ONE )
      ar(NZ) = max( ar(NZ), 0.*ONE )
c
c
c  Next, ensure that polynomial functions do not deviate beyond the
c  range [<al(k)>,<ar(k)>]
c
      do k = 1,NZ

        dpc = cvert(k) / dz2(ixy,k)
        if( (ar(k)-dpc)*(dpc-al(k)) .le. 0. ) then
          al(k) = dpc
          ar(k) = dpc
        endif

        if( (ar(k)-al(k))*( dpc - 0.5*(al(k)+ar(k)) )
     $      .gt. 1./6.*(ar(k)-al(k))**2 )
     $     al(k) = 3.*dpc - 2.*ar(k)

        if( (ar(k)-al(k))*( dpc - 0.5*(al(k)+ar(k)) )
     $      .lt. -1./6.*(ar(k)-al(k))**2 )
     $     ar(k) = 3.*dpc - 2.*al(k)

      enddo
c
c
c  Calculate fluxes across each layer boundary
c
      do k = 1,NZ

        dpc = cvert(k) / dz2(ixy,k)
        dela(k) = ar(k) - al(k)
        a6(k) = 6. * ( dpc - 0.5*(ar(k)+al(k)) )

      enddo

      do k = 1,NZ-1

        com2  = ( dz2(ixy,k) + dz2(ixy,k+1) ) / 2.
        x = vtrans(k+1)*dtime/dz2(ixy,k)
        xpos = abs(x)
c
c  Upward transport rate
c
        if( vtrans(k+1) .gt. 0. )then

          if( x .lt. 1. .and. cvert(k) .ne. 0. )then
            vertadvu(k+1) = ( vtrans(k+1) * com2 ) *
     $                 ( ( ar(k) - 0.5*dela(k)*x +
     $                 (x/2. - x**2/3.)*a6(k) ) / cvert(k) )
c
c  If Courant # > 1, use upwind advection
c
          else
            vertadvu(k+1) = vtrans(k+1)
          endif
c
c  Downward transport rate
c
        elseif( vtrans(k+1) .lt. 0. )then

          if( x .gt. -1. .and. cvert(k+1) .ne. 0. )then
            vertadvd(k+1) = ( -vtrans(k+1) * com2 ) *
     $               ( ( al(k+1) + 0.5*dela(k+1)*xpos +
     $               ( xpos/2. - xpos**2/3.)*a6(k+1) ) / cvert(k+1) )
          else
            vertadvd(k+1) = -vtrans(k+1)
          endif

        endif

      enddo    ! k = 1,NZ-1
c
c
c  Lower boundary transport rates:  If I_FIXED_CONC boundary
c  condtion is selected, then use concentration assumed just beyond
c  the lowest layer edge to calculate the transport rate across
c  the bottom boundary of the model.
c
      if( ibbnd .eq. I_FIXED_CONC ) then

        com2  = ( dz2(ixy,1) + dz2(ixy,itwo) ) / 2.
        x = vtrans(1)*dtime/dz2(ixy,1)
        xpos = abs(x)
        cvert0 = cvert_bbnd
        if( vtrans(1) .gt. 0. )then

          if( x .lt. 1. .and. cvert0 .ne. 0. )then
            vertadvu(1) = vtrans(1)/cvert0*com2
     $                 * ( ar(1) - 0.5*dela(1)*x +
     $                 (x/2. - x**2/3.)*a6(1) )
          else
            vertadvu(1) = vtrans(1)
          endif

        elseif( vtrans(1) .lt. 0. )then

          if( x .gt. -1. .and. cvert(1) .ne. 0. )then
            vertadvd(1) = -vtrans(1)/
     $                 cvert(1)*com2
     $                 * ( al(1) + 0.5*dela(1)*xpos +
     $                 (xpos/2. - xpos**2/3.)*a6(1) )
          else
            vertadvd(1) = -vtrans(1)
          endif

        endif

      endif
c
c
c  Upper boundary transport rates
c
      if( itbnd .eq. I_FIXED_CONC ) then

        com2  = ( dz2(ixy,NZ) + dz2(ixy,nzm1) ) / 2.
        x = vtrans(NZ+1)*dtime/dz2(ixy,NZ)
        xpos = abs(x)
        cvertnzp1 = cvert_tbnd

        if( vtrans(NZ+1) .gt. 0. )then

          if( x .lt. 1. .and. cvert(NZ) .ne. 0. )then
            vertadvu(NZ+1) = vtrans(NZ+1)/cvert(NZ)*com2
     $               * ( ar(NZ) - 0.5*dela(NZ)*x +
     $               (x/2. - x**2/3.)*a6(NZ) )
          else
            vertadvu(NZ+1) = vtrans(NZ+1)
          endif

        elseif( vtrans(NZ+1) .lt. 0. )then

          if( x .gt. -1. .and. cvertnzp1 .ne. 0. )then
            vertadvd(NZ+1) = -vtrans(NZ+1)/
     $               cvertnzp1*com2
     $               * ( al(NZ) + 0.5*dela(NZ)*xpos +
     $               (xpos/2. - xpos**2/3.)*a6(NZ) )
          else
            vertadvd(NZ+1) = -vtrans(NZ+1)
          endif

        endif

      endif
!       if (ibin .eq. 30) then
!!	   open(unit=99, file='vertadv1.txt', status='unknown') !EJL
!!	   write(99,*) 'vtrans', vtrans
!!	   write(99,*) 'vertadvd', vertadvd
!!       write(99,*) 'vertadvu', vertadvu
!!	   close(unit=99)
!       endif
c
c
c  Return to caller with vertical transport rates.
c
      return
      end
