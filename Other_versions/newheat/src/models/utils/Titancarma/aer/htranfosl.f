       subroutine htranfosl(idir)
c
c
c  @(#) htranfosl.f  Jensen  Jan-1997
c  This routine calculates horizontal advection rates using
c  First-order semi-Lagrangian technique [Bates and McDonald,
c  Mon. Wea. Rev., 110, 1831-1842, 1982].
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
      dimension cold(NXORNY)
c
c
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter htranfosl'
c
c
c  <idir> = IDIRX: x-direction; <idir>=IDIRY: y-direction
c
      if( idir .eq. IDIRX ) then
        nadv = NX
      end if
      if( idir .eq. IDIRY ) then
        nadv = NY
      end if
c
c
c  Store old values of chor
c
      do i = 1,nadv
        cold(i) = chor(i)
      enddo
c
c
c  Interpolate to estimate concentration at departure point
c  based on concentrations at neighboring grid points.
c
      do i = 2,nadv-1

        dep = abs( htrans(i) ) * dtime
        frac = dep/dhor(i)
        if( htrans(i) .ge. 0. ) then
          chor(i) = (1.-frac)*cold(i) + frac*cold(i-1)
        else
          chor(i) = (1.-frac)*cold(i) + frac*cold(i+1)
        endif

      enddo
c
c
c  Periodic boundary conditions
c
      i = 1
      dep = abs( htrans(i) ) * dtime
      frac = dep/dhor(i)
      if( htrans(i) .ge. 0. ) then
        chor(i) = (1.-frac)*cold(i) + frac*cold(nadv)
      else
        chor(i) = (1.-frac)*cold(i) + frac*cold(i+1)
      endif

      i = nadv
      dep = abs( htrans(i) ) * dtime
      frac = dep/dhor(i)
      if( htrans(i) .ge. 0. ) then
        chor(i) = (1.-frac)*cold(i) + frac*cold(i-1)
      else
        chor(i) = (1.-frac)*cold(i) + frac*cold(1)
      endif
c
c
c  Return to caller with new concentrations after horizontal transport.
c
      return
      end
