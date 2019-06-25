      subroutine molfrac
c
c
c  @(#) molfrac.f  Barth  Jul-2007
c  Calculate mole fractions in mixed droplets
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c  Define formats
c
    1 format(a,i4,'  pc(gcore): ',1pe11.4,'  df: ',e11.4,
     $       '  binmass: ',e11.4,'  f: ',e11.4,'  fmax: ',
     $       e11.4,'  volmass: ',e11.4,'  cmf: ',e11.4,
     $       '  pc(cloud): ',e11.4)

c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter molfrac'
   
      do igrp = 1,NGROUP

       iepart = ienconc(igrp)
       imaingas = igrowgas(iepart)

      !Only need to do calculations if cloud group
       if( imaingas .ne. 0 ) then

        do igas=1,NGAS
          xfrac(igas) = ZERO
        enddo

        do ibin = 1,NBIN

         coretot = ZERO
         binmass = pc3(ixyz,ibin,iepart) * rmass(ibin,igrp)

       !Find involatile coremass fraction
         do ic = 1,ncore(igrp)
          iecore = icorelem(ic,igrp)
          if( itype(iecore) .eq. I_COREMASS ) 
     $     coretot = coretot + pc3(ixyz,ibin,iecore) 
         enddo !core elements
         coreavg = coretot / pc3(ixyz,ibin,iepart)
         coreavg = min( rmass(ibin,igrp),coreavg )
         cmf(ibin,igrp) = coreavg / rmass(ibin,igrp)

       !Total volatile mass of particles in ibin
         volmass = binmass * (ONE - cmf(ibin,igrp))

       !Sum mole fractions for each secondary condensate
         do ic = 1,ncore(igrp)
          iecore = icorelem(ic,igrp)
          if( itype(iecore) .eq. I_GROWCORE ) then
             igas = igrowgas(iecore)
             xfrac(igas) = xfrac(igas) + pc3(ixyz,ibin,iecore)/volmass 
          endif
         enddo !core elements

        enddo !bins

       !Get average mole fraction for the layer
        xfrac(igas) = xfrac(igas)/real(NBINS)

        
       endif !imaingas ne 0
      enddo !group
c
c
c  Return to caller with growcore fractions fixed
c
      return
      end
