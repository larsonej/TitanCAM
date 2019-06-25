       subroutine hordif(idir)
c
c
c  @(#) hordif.f  Jensen  Jan-1997
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter hordif'
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
c  Calculate new concentrations using simple first-order
c  difference approximation for diffusion.
c
      do i = 2,nadv-1

        chor(i) = cold(i) + hdiff(i)*
     $            ( cold(i+1) + cold(i-1) - 2.*cold(i) ) /
     $            dhor(i)**2

      enddo
c
c
c  Periodic boundary conditions
c
      i = 1
      chor(i) = cold(i) + hdiff(i)*
     $          ( cold(i+1) + cold(nadv) - 2.*cold(i) ) /
     $          dhor(i)**2

      i = nadv
      chor(i) = cold(i) + hdiff(i)*
     $          ( cold(1) + cold(i-1) - 2.*cold(i) ) /
     $          dhor(i)**2
c
c
c  Return to caller with new concentrations after horizontal diffusion.
c
      return
      end
