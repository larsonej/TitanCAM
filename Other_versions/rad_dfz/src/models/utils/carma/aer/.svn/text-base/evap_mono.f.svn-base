      subroutine evap_mono(ibin,ig,iavg,ieto,igto)
c
c
c  @(#) evap_mono.f  Ackerman  Aug-2001
c  This routine calculates particle source terms <evappe> due to total
c  evaporation from bin <ibin> group <ig> into a monodisperse
c  distribution.
c
c  Distinct evaporation of cores has not been treated.
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Local declarations
c
      logical conserve_mass
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter evap_mono'
c
c
c-------------------------------------------------------------------------------
c
c
c  Define option to conserve mass or number when a choice must be made
c  during monodisperse total evaporation beyond CN grid -- should be done in setupaer()
c
      conserve_mass = .true.
c
c
c  Set automatic flag for total evaporation used in gasexchange()
c
      totevap(ibin,ig) = .true.
c
c
c  Possibly put all of core mass into largest, smallest, or
c  smallest nucelated CN bin 
c
      if( too_big .or. too_small .or. nuc_small )then

        if( too_big )then
          jbin = NBIN
        elseif( too_small )then
          jbin = 1
        else
          jbin = inucmin(igto)
        endif

        if( conserve_mass )then
          factor = coreavg/rmass(jbin,igto)
        else
          factor = ONE
        endif
c
c  First the CN number concentration element
c 
        evappe(jbin,ieto) = evappe(jbin,ieto) + factor*evdrop
c
c  Now the CN cores
c
        do ic = 2, ncore(ig)
          iecore = icorelem(ic,ig)
          ie2cn  = ievp2elem(iecore)
          evappe(jbin,ie2cn) = evappe(jbin,ie2cn) +
     $       factor*evcore(ic)*rmass(jbin,igto)
        enddo

      else
c
c
c  Partition core mass between two CN bins, conserving total core mass
c  and number.  The number will be subdivided into bins <iavg> and <iavg>-1.
c
       if( iavg .le. 1 .or. iavg .gt. NBIN )then
         print*,' stop in evap_mono: bad iavg = ', iavg
         call endcarma
       endif

       fracmass = ( rmass(iavg,igto) - coreavg ) /
     $    diffmass(iavg,igto,iavg-1,igto)
c      fracmass = max( 0., min( ONE, fracmass ) )
c
c  First the CN number concentration element
c 
       evappe(iavg-1,ieto) = evappe(iavg-1,ieto) + evdrop*fracmass
       evappe(iavg,ieto) = evappe(iavg,ieto) + evdrop*( ONE - fracmass )
c
c  Now the cores
c
       do ic = 2, ncore(ig)
         iecore = icorelem(ic,ig)
         ie2cn  = ievp2elem(iecore)
         evappe(iavg-1,ie2cn) = evappe(iavg-1,ie2cn) +
     $      rmass(iavg-1,igto)*evcore(ic)*fracmass
         evappe(iavg,ie2cn) = evappe(iavg,ie2cn) +
     $      rmass(iavg,igto)*evcore(ic)*( ONE - fracmass )
       enddo
 
      endif

      return
      end
