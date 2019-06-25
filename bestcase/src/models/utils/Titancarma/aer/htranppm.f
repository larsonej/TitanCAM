       subroutine htranppm(idir,nadv)
c
c
c  @(#) htranppm.f  Jensen  Jan-1997
c  This routine calculates horizontal advection rates using
c  Piecewise Polynomial Method [Colela and Woodard, J. Comp. Phys.,
c  54, 174-201, 1984]
c  Currently only valid for uniform grid and periodic B.C.s
c
c  Argument list input:
c    idir
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
      dimension dela(NXORNY), delma(NXORNY), aju(NXORNY),
     $  ar(NXORNY), al(NXORNY), a6(NXORNY), cold(NXORNY),
     $  flux(NXORNY), htransbnd(NXORNY+1),
     $  fluxu(NXORNY+1), fluxd(NXORNY+1),
     $  ctempu(NXORNY), ctempl(NXORNY),
     $  bl(NXORNY), dl(NXORNY), ul(NXORNY),
     $          el(NXORNY), fl(NXORNY)
c
c
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter htranppm'
c
c
c  <idir> = IDIRX: x-direction; <idir>=IDIRY: y-direction
c
c      if( idir .eq. IDIRX ) then
c        nadv = NX
c      end if
c      if( idir .eq. IDIRY ) then
c        nadv = NY
c      end if
c
c
c  Store old values of chor
c
      do i = 1,nadv
        cold(i) = chor(i)
      enddo
c
c
c  Initialize fluxes to zero
c
      do i = 1,nadv
        flux(i) = 0.
      enddo
c
c
c  First, use cubic fits to estimate concentration values at layer
c  boundaries
c
      do i = 2,nadv-1

        dpc = chor(i) / dhor(i)
        dpc1 = chor(i+1) / dhor(i+1)
        dpcm1 = chor(i-1) / dhor(i-1)
        ratt1 = dhor(i) /
     $    ( dhor(i-1) + dhor(i) + dhor(i+1) )
        ratt2 = ( 2.*dhor(i-1) + dhor(i) ) /
     $          ( dhor(i+1) + dhor(i) )
        ratt3 = ( 2.*dhor(i+1) + dhor(i) ) /
     $          ( dhor(i-1) + dhor(i) )

        dela(i) = ratt1 *
     $           ( ratt2*(dpc1-dpc) + ratt3*(dpc-dpcm1) )

        if( (dpc1-dpc)*(dpc-dpcm1) .gt. 0. .and. dela(i) .ne. 0. ) then
          delma(i) = min( abs(dela(i)), 2.*abs(dpc-dpc1),
     $         2.*abs(dpc-dpcm1) ) * abs(dela(i))/dela(i)
        else
          delma(i) = 0.
        endif

      enddo     ! i = 2,nadv-1

      do i = 2,nadv-2
        dpc = chor(i) / dhor(i)
        dpc1 = chor(i+1) / dhor(i+1)
        dpcm1 = chor(i-1) / dhor(i-1)
        rat1 = dhor(i) /
     $          ( dhor(i) + dhor(i+1) )
        rat2 = 2. * dhor(i+1) * dhor(i) /
     $         ( dhor(i) + dhor(i+1) )
        rat3 = ( dhor(i-1) + dhor(i) ) /
     $         ( 2.*dhor(i) + dhor(i+1) )
        rat4 = ( dhor(i+2) + dhor(i+1) ) /
     $         ( 2.*dhor(i+1) + dhor(i) )
        den1 = dhor(i-1) + dhor(i) +
     $         dhor(i+1) + dhor(i+2)
c
c
c  <aju(i)> is the estimate for concentration (dn/dhor) at layer
c  boundary <i>+1/2.
c
        aju(i) = dpc + rat1*(dpc1-dpc) + 1./den1 *
     $           ( rat2*(rat3-rat4)*(dpc1-dpc) -
     $           dhor(i)*rat3*delma(i+1) +
     $           dhor(i+1)*rat4*delma(i) )

      enddo     ! i = 2,nadv-2
c
c
c  For x-direction, use periodic boundary conditions
c
      if( idir .eq. IDIRX ) then

        i = 1
        dpc = chor(i) / dhor(i)
        dpc1 = chor(i+1) / dhor(i+1)
        dpcm1 = chor(nadv) / dhor(nadv)
        ratt1 = dhor(i) /
     $    ( dhor(nadv) + dhor(i) + dhor(i+1) )
        ratt2 = ( 2.*dhor(nadv) + dhor(i) ) /
     $          ( dhor(i+1) + dhor(i) )
        ratt3 = ( 2.*dhor(i+1) + dhor(i) ) /
     $          ( dhor(nadv) + dhor(i) )

        dela(i) = ratt1 *
     $           ( ratt2*(dpc1-dpc) + ratt3*(dpc-dpcm1) )

        if( (dpc1-dpc)*(dpc-dpcm1) .gt. 0. .and. dela(i) .ne. 0. ) then
          delma(i) = min( abs(dela(i)), 2.*abs(dpc-dpc1),
     $         2.*abs(dpc-dpcm1) ) * abs(dela(i))/dela(i)
        else
          delma(i) = 0.
        endif

        i = nadv
        dpc = chor(i) / dhor(i)
        dpc1 = chor(1) / dhor(1)
        dpcm1 = chor(i-1) / dhor(i-1)
        ratt1 = dhor(i) /
     $    ( dhor(i-1) + dhor(i) + dhor(1) )
        ratt2 = ( 2.*dhor(i-1) + dhor(i) ) /
     $          ( dhor(1) + dhor(i) )
        ratt3 = ( 2.*dhor(1) + dhor(i) ) /
     $          ( dhor(i-1) + dhor(i) )

        dela(i) = ratt1 *
     $           ( ratt2*(dpc1-dpc) + ratt3*(dpc-dpcm1) )

        if( (dpc1-dpc)*(dpc-dpcm1) .gt. 0. .and. dela(i) .ne. 0. ) then
          delma(i) = min( abs(dela(i)), 2.*abs(dpc-dpc1),
     $         2.*abs(dpc-dpcm1) ) * abs(dela(i))/dela(i)
        else
          delma(i) = 0.
        endif

        i = 1
        dpc = chor(i) / dhor(i)
        dpc1 = chor(i+1) / dhor(i+1)
        dpcm1 = chor(nadv) / dhor(nadv)
        rat1 = dhor(i) /
     $          ( dhor(i) + dhor(i+1) )
        rat2 = 2. * dhor(i+1) * dhor(i) /
     $         ( dhor(i) + dhor(i+1) )
        rat3 = ( dhor(nadv) + dhor(i) ) /
     $         ( 2.*dhor(i) + dhor(i+1) )
        rat4 = ( dhor(i+2) + dhor(i+1) ) /
     $         ( 2.*dhor(i+1) + dhor(i) )
        den1 = dhor(nadv) + dhor(i) +
     $         dhor(i+1) + dhor(i+2)
        aju(i) = dpc + rat1*(dpc1-dpc) + 1./den1 *
     $           ( rat2*(rat3-rat4)*(dpc1-dpc) -
     $           dhor(i)*rat3*delma(i+1) +
     $           dhor(i+1)*rat4*delma(i) )

        i = nadv
        dpc = chor(i) / dhor(i)
        dpc1 = chor(1) / dhor(1)
        dpcm1 = chor(i-1) / dhor(i-1)
        rat1 = dhor(i) /
     $          ( dhor(i) + dhor(1) )
        rat2 = 2. * dhor(1) * dhor(i) /
     $         ( dhor(i) + dhor(1) )
        rat3 = ( dhor(i-1) + dhor(i) ) /
     $         ( 2.*dhor(i) + dhor(1) )
        rat4 = ( dhor(2) + dhor(1) ) /
     $         ( 2.*dhor(1) + dhor(i) )
        den1 = dhor(i-1) + dhor(i) +
     $         dhor(1) + dhor(2)
        aju(i) = dpc + rat1*(dpc1-dpc) + 1./den1 *
     $           ( rat2*(rat3-rat4)*(dpc1-dpc) -
     $           dhor(i)*rat3*delma(1) +
     $           dhor(1)*rat4*delma(i) )

        i = nadv-1
        dpc = chor(i) / dhor(i)
        dpc1 = chor(nadv) / dhor(nadv)
        dpcm1 = chor(i-1) / dhor(i-1)
        rat1 = dhor(i) / ( dhor(i) + dhor(nadv) )
        rat2 = 2. * dhor(nadv) * dhor(i) /
     $         ( dhor(i) + dhor(nadv) )
        rat3 = ( dhor(i-1) + dhor(i) ) /
     $         ( 2.*dhor(i) + dhor(nadv) )
        rat4 = ( dhor(1) + dhor(nadv) ) /
     $         ( 2.*dhor(nadv) + dhor(i) )
        den1 = dhor(i-1) + dhor(i) +
     $         dhor(nadv) + dhor(1)
        aju(i) = dpc + rat1*(dpc1-dpc) + 1./den1 *
     $           ( rat2*(rat3-rat4)*(dpc1-dpc) -
     $           dhor(i)*rat3*delma(nadv) +
     $           dhor(nadv)*rat4*delma(i) )

c
c
c  For y-direction, extrapolate gradient across boundary
c
      elseif( idir .eq. IDIRY ) then

        i = 1
        dpc = chor(i) / dhor(i)
        dpc1 = chor(i+1) / dhor(i+1)
        dpcm1 = dpc
        ratt1 = dhor(i) /
     $    ( dhor(i) + dhor(i) + dhor(i+1) )
        ratt2 = ( 2.*dhor(i) + dhor(i) ) /
     $          ( dhor(i+1) + dhor(i) )
        ratt3 = ( 2.*dhor(i+1) + dhor(i) ) /
     $          ( dhor(i) + dhor(i) )
        dela(i) = ratt1 *
     $           ( ratt2*(dpc1-dpc) + ratt3*(dpc-dpcm1) )
        if( (dpc1-dpc)*(dpc-dpcm1) .gt. 0. .and. dela(i) .ne. 0. ) then
          delma(i) = min( abs(dela(i)), 2.*abs(dpc-dpc1),
     $         2.*abs(dpc-dpcm1) ) * abs(dela(i))/dela(i)
        else
          delma(i) = 0.
        endif

        i = nadv
        dpc = chor(i) / dhor(i)
        dpc1 = dpc
        dpcm1 = chor(i-1) / dhor(i-1)
        ratt1 = dhor(i) /
     $    ( dhor(i-1) + dhor(i) + dhor(i) )
        ratt2 = ( 2.*dhor(i-1) + dhor(i) ) /
     $          ( dhor(i) + dhor(i) )
        ratt3 = ( 2.*dhor(i) + dhor(i) ) /
     $          ( dhor(i-1) + dhor(i) )

        dela(i) = ratt1 *
     $           ( ratt2*(dpc1-dpc) + ratt3*(dpc-dpcm1) )

        if( (dpc1-dpc)*(dpc-dpcm1) .gt. 0. .and. dela(i) .ne. 0. ) then
          delma(i) = min( abs(dela(i)), 2.*abs(dpc-dpc1),
     $         2.*abs(dpc-dpcm1) ) * abs(dela(i))/dela(i)
        else
          delma(i) = 0.
        endif

        i = 1
        dpc = chor(i) / dhor(i)
        dpc1 = chor(i+1) / dhor(i+1)
        dpcm1 = dpc
        rat1 = dhor(i) /
     $          ( dhor(i) + dhor(i+1) )
        rat2 = 2. * dhor(i+1) * dhor(i) /
     $         ( dhor(i) + dhor(i+1) )
        rat3 = ( dhor(i) + dhor(i) ) /
     $         ( 2.*dhor(i) + dhor(i+1) )
        rat4 = ( dhor(i+2) + dhor(i+1) ) /
     $         ( 2.*dhor(i+1) + dhor(i) )
        den1 = dhor(i) + dhor(i) +
     $         dhor(i+1) + dhor(i+2)
        aju(i) = dpc + rat1*(dpc1-dpc) + 1./den1 *
     $           ( rat2*(rat3-rat4)*(dpc1-dpc) -
     $           dhor(i)*rat3*delma(i+1) +
     $           dhor(i+1)*rat4*delma(i) )

        i = nadv
        dpc = chor(i) / dhor(i)
        dpc1 = dpc
        dpcm1 = chor(i-1) / dhor(i-1)
        rat1 = dhor(i) /
     $          ( dhor(i) + dhor(i) )
        rat2 = 2. * dhor(i) * dhor(i) /
     $         ( dhor(i) + dhor(i) )
        rat3 = ( dhor(i-1) + dhor(i) ) /
     $         ( 2.*dhor(i) + dhor(i) )
        rat4 = ( dhor(i) + dhor(i) ) /
     $         ( 2.*dhor(i) + dhor(i) )
        den1 = dhor(i-1) + dhor(i) +
     $         dhor(i) + dhor(i)
        aju(i) = dpc + rat1*(dpc1-dpc) + 1./den1 *
     $           ( rat2*(rat3-rat4)*(dpc1-dpc) -
     $           dhor(i)*rat3*delma(i) +
     $           dhor(i)*rat4*delma(i) )

        i = nadv-1
        dpc = chor(i) / dhor(i)
        dpc1 = chor(i+1) / dhor(i+1)
        dpcm1 = chor(i-1) / dhor(i-1)
        rat1 = dhor(i) /
     $          ( dhor(i) + dhor(i+1) )
        rat2 = 2. * dhor(i+1) * dhor(i) /
     $         ( dhor(i) + dhor(i+1) )
        rat3 = ( dhor(i-1) + dhor(i) ) /
     $         ( 2.*dhor(i) + dhor(i+1) )
        rat4 = ( dhor(i+1) + dhor(i+1) ) /
     $         ( 2.*dhor(i+1) + dhor(i) )
        den1 = dhor(i-1) + dhor(i) +
     $         dhor(i+1) + dhor(i+1)
        aju(i) = dpc + rat1*(dpc1-dpc) + 1./den1 *
     $           ( rat2*(rat3-rat4)*(dpc1-dpc) -
     $           dhor(i)*rat3*delma(i+1) +
     $           dhor(i+1)*rat4*delma(i) )

      endif
c
c
c  Now construct polynomial functions in each layer
c
      do i = 2,nadv

        al(i) = aju(i-1)
        ar(i) = aju(i)

      enddo

      al(1) = aju(nadv)
      ar(1) = aju(1)
c
c
c  Next, ensure that polynomial functions do not deviate beyond the
c  range [<al(i)>,<ar(i)>]
c
      do i = 1,nadv

        dpc = chor(i) / dhor(i)
        if( (ar(i)-dpc)*(dpc-al(i)) .le. 0. ) then
          al(i) = dpc
          ar(i) = dpc
        endif

        term1 = (ar(i)-al(i)) * ( dpc - 0.5*(al(i)+ar(i)) )
        term2 = 1./6. * (ar(i)-al(i))**2
        if( term1 .gt. term2 )
     $     al(i) = 3.*dpc - 2.*ar(i)

        if( term1 .lt. -term2 )
     $     ar(i) = 3.*dpc - 2.*al(i)

      enddo
c
c
c  Construct polynomials
c
      do i = 1,nadv

        dpc = chor(i) / dhor(i)
        dela(i) = ar(i) - al(i)
        a6(i) = 6. * ( dpc - 0.5*(ar(i)+al(i)) )

      enddo
c
c  Calculate wind speeds at boundaries by averaging
c
      do i = 1,nadv-1
        htransbnd(i+1) = 0.5 * ( htrans(i) + htrans(i+1) )
      enddo
      htransbnd(1) = 0.5 * ( htrans(nadv) + htrans(1) )
      htransbnd(nadv+1) = 0.5 * ( htrans(nadv) + htrans(1) )

      do i = 1,nadv+1
        fluxu(i) = 0
        fluxd(i) = 0
      enddo
      do i = 1,nadv-1

        x = htransbnd(i+1)*dtime/dhor(i)
        xpos = abs(x)
c
c  Upward transport rate across i+1/2 boundary
c
        if( htransbnd(i+1) .gt. 0. )then

          if( x .lt. 1. )then
            fluxu(i+1) = ( htransbnd(i+1) * dhor(i) ) *
     $                 ( ( ar(i) - 0.5*dela(i)*x +
     $                 (x/2. - x**2/3.)*a6(i) ) / chor(i) )
          else
            fluxu(i+1) = htransbnd(i+1)
          endif
c
c  Downward transport rate across i+1/2 boundary
c
        elseif( htransbnd(i+1) .lt. 0. )then

          if( x .gt. -1. )then
            fluxd(i+1) = ( -htransbnd(i+1) * dhor(i) ) *
     $               ( ( al(i+1) + 0.5*dela(i+1)*xpos +
     $               ( xpos/2. - xpos**2/3.)*a6(i+1) ) / chor(i+1) )
          else
            fluxd(i+1) = -htransbnd(i+1)
          endif

        endif

      enddo    ! i = 1,nadv-1

      if( idir .eq. IDIRX ) then

        i = 0
        x = htransbnd(i+1)*dtime/dhor(nadv)
        xpos = abs(x)
        if( htransbnd(i+1) .gt. 0. )then
          if( x .lt. 1. )then
            fluxu(i+1) = ( htransbnd(i+1) * dhor(nadv) ) *
     $               ( ( ar(nadv) - 0.5*dela(nadv)*x +
     $               (x/2. - x**2/3.)*a6(nadv) ) / chor(nadv) )
          else
            fluxu(i+1) = htransbnd(i+1)
          endif
        elseif( htransbnd(i+1) .lt. 0. )then
          if( x .gt. -1. )then
            fluxd(i+1) = ( -htransbnd(i+1) * dhor(nadv) ) *
     $               ( ( al(i+1) + 0.5*dela(i+1)*xpos +
     $               ( xpos/2. - xpos**2/3.)*a6(i+1) ) / chor(i+1) )
          else
            fluxd(i+1) = -htransbnd(i+1)
          endif
        endif

        i = nadv
        x = htransbnd(i+1)*dtime/dhor(nadv)
        xpos = abs(x)
        if( htransbnd(i+1) .gt. 0. )then
          if( x .lt. 1. )then
            fluxu(i+1) = ( htransbnd(i+1) * dhor(i) ) *
     $                 ( ( ar(i) - 0.5*dela(i)*x +
     $                 (x/2. - x**2/3.)*a6(i) ) / chor(i) )
          else
            fluxu(i+1) = htransbnd(i+1)
          endif
        elseif( htransbnd(i+1) .lt. 0. )then
          if( x .gt. -1. )then
            fluxd(i+1) = ( -htransbnd(i+1) * dhor(nadv) ) *
     $               ( ( al(1) + 0.5*dela(1)*xpos +
     $               ( xpos/2. - xpos**2/3.)*a6(1) ) / chor(1) )
          else
            fluxd(i+1) = -htransbnd(i+1)
          endif
        endif

      endif

      if( idir .eq. IDIRY ) then

        fluxu(1) = 0.
        fluxd(1) = 0.
        fluxu(nadv+1) = 0.
        fluxd(nadv+1) = 0.
        fluxu(nadv) = 0.
        fluxd(nadv) = 0.

      endif

      uc = 0.
      do i = 1,nadv
        cour = dhor(i)/dtime - ( fluxu(i+1) + fluxd(i) )
        if( cour .lt. 0. ) uc = 1.
      enddo

      do i = 2,nadv
        ctempl(i) = chor(i-1)
        ctempu(i-1) = chor(i)
      enddo
      ctempl(1) = chor(nadv)
      ctempu(nadv) = chor(1)

      do i = 1,nadv
        al(i) = uc * fluxd(i+1)
        bl(i) = -( uc*(fluxd(i)+fluxu(i+1)) + dhor(i)/dtime )
        ul(i) = uc * fluxu(i)
        dl(i) = chor(i) * ( (1.-uc)*(fluxd(i)+fluxu(i+1)) -
     $                      dhor(i)/dtime ) - (1.-uc) *
     $          ( fluxu(i)*ctempl(i) + fluxd(i+1)*ctempu(i) )
      enddo

      el(1) = dl(1)/bl(1)
      fl(1) = al(1)/bl(1)
      do i = 2,nadv
        denom = bl(i) - ul(i) * fl(i-1)
        el(i) = ( dl(i) - ul(i)*el(i-1) ) / denom
        fl(i) = al(i) / denom
      enddo

      chor(nadv) = el(nadv)
      do i = nadv-1,1,-1
        chor(i) = el(i) - fl(i)*chor(i+1)
      enddo
c
c
c  Return to caller with new concentrations after horizontal transport.
c
      return
      end
