      subroutine evap_ingrp(ibin,ig,ip)
c
c
c  @(#) evap_ingrp.f  Ackerman  Aug-2001
c  This routine calculates particle source terms <evappe> of droplets
c  evaporating within a particle group.
c
c  Distinct evaporation of cores has not been treated.
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter evap_ingrp'
c
c
c-------------------------------------------------------------------------------
c
c
c  The smallest bin cannot be a source to smaller bins in same group
c
      if( ibin .eq. 1 )then
        return
      endif
c
c
c  Evaluate evaporation source term <evappe> for all elements in group
c
      do isub = 1, nelemg(ig)
        ie = ip + isub - 1
        evappe(ibin-1,ie) = evappe(ibin-1,ie) +
     $     pc3(ixyz,ibin,ie)*evaplg(ibin,ig)
      enddo

      return
      end
