      subroutine smallconc(ibin,ie)
c
c
c  @(#) smallconc.f  Ackerman  Oct-1997
c  This routine ensures limits all particle concentrations in a grid box
c  to SMALL_PC.  In bins where this limitation results in the particle 
c  concentration changing, the core mass fraction and second moment fraction 
c  are set to <FIX_COREF>. 
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
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter smallconc'
   
      ig = igelem(ie)
      ip = ienconc(ig)

c     if( pc3(ixyz,ibin,ip) .le. small_val(ibin,ip) .or.
c    $    pc3(ixyz,ibin,ie) .le. small_val(ibin,ie) )then
c         pc3(ixyz,ibin,ie) = small_val(ibin,ie)

      if( pc3(ixyz,ibin,ip) .lt. 0. .or.
     $    pc3(ixyz,ibin,ie) .lt. 0. )then
          write(LUNOPRT,'(/,a)') 'Negative particle concentration'
          pc3(ixyz,ibin,ie) = 0.

      endif
 
      return
      end
