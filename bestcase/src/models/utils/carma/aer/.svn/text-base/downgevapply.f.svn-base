       subroutine downgevapply
c
c
c  @(#) downgevapply.f  Ackerman  Dec-1995
c  This routine applies evaporation and nucleation production terms to
c  particle concentrations.
c
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
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
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter downgevapply'
c
c
c  Visit each radius bin for each element to compute particle production 
c  due to evaporation and element transfer processes for which the source
c  element number is greater than the target element number
c
      do ielem = 1,NELEM
        do ibin = 1,NBIN

            pc3(ixyz,ibin,ielem) = pc3(ixyz,ibin,ielem) +
     $                         dtime * ( evappe(ibin,ielem) +
     $                                   rnucpe(ibin,ielem) )

        enddo
      enddo
c
c
c  Return to caller with evaporation and down-grid element transfer
c  production terms applied to particle concentrations.
c
      return
      end
