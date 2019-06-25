       subroutine supersat
c
c
c  @(#) supersat.f  Ackerman  Dec-1995
c  This routine evaluates supersaturations <supsatl> and <supsati> for all gases.
c
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial point's indices.
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
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter supersat'
c
c-------------------------------------------------------------------------------
c
c
c   Calculate vapor pressures.
c
      call vaporp
c
c
c   Loop over all gases
c
      do igas = 1,NGAS
c
c
c   Define gas constant for this gas
c
        rvap = RGAS/gwtmol(igas)

        gc_cgs = gc3(ixyz,igas) / (zmet3(ixyz)*xmet3(ixyz)*ymet3(ixyz))

        supsatl3(ixyz,igas) = ( gc_cgs * rvap * t3(ixyz) -
     $     pvapl3(ixyz,igas) ) / pvapl3(ixyz,igas)

        supsati3(ixyz,igas) = ( gc_cgs * rvap * t3(ixyz) -
     $     pvapi3(ixyz,igas) ) / pvapi3(ixyz,igas)
 
c
c    Add a perturbation:  Increase s at 10 km by a factor of 5
c
c       if(itime .eq. ibtime + 50)        !start perturbation 50 timesteps into run
c    $    supsati3(6,igas) = supsati3(6,igas)*5.

      enddo
c
c
c  Return to caller with supersaturations evaluated.
c
      return
      end
