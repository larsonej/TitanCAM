      subroutine evap_poly(ibin,ig,iavg,ieto,igto)
c
c
c  @(#) evap_poly.f  Ackerman  Aug-2001
c  This routine calculates particle source terms <evappe> due to
c  total evaporation into a polydisperse CN distribution by assuming
c  that the pdf of core mass is log-normal skewed by mass raised to
c  the -3/2 power (which guarantees average core mass from pdf is the
c  same as average core mass).
c
c  Distinct evaporation of cores has not been treated.
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Local declarations
c
      dimension prob(NBIN)
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter evap_poly'
c
c
c-------------------------------------------------------------------------------
c
c
c  Treat total evaporation from a polydisperse core mass distribution:
c  assume a log-normal CN size distribution and conserve number and mass as
c  described by Turco (NASA Technical Paper 1362).  Don't put anything
c  in CN bins with index less than <inucmin>.
c  
c  Set automatic flag for total evaporation used in gasexchange()
c
      totevap(ibin,ig) = .true.
c
c
c  Calculate number <rn_norms,rn_norml> and mass <rm_norms,rm_norml>
c  normalization factors for cores smaller and larger than <rmass(m,igto)>.
c
      rn_norms = 0.
      rn_norml = 0.
      rm_norms = 0.
      rm_norml = 0.
      kount_s = 0
      kount_l = 0

      do ito = inucmin(igto), NBIN

        rmassto = rmass(ito,igto)           
        dmto = dm(ito,igto)
c
c  <prob> is probability that core mass is in CN bin <ito>.
c
        if( coreavg .gt. 0. .and. coresig .gt. 0. )then
          expon = -log( rmassto/coreavg )**2 / ( 2.*coresig )
          expon = max(-POWMAX, expon)
        else
          expon = 0.
        endif
        prob(ito) = rmassto**(-1.5) * exp( expon )

c       if( ito .le. iavg )then
        if( ito .lt. iavg )then     ! Kevin M
          rn_norms = rn_norms + prob(ito)*dmto
          rm_norms = rm_norms + prob(ito)*dmto*rmassto
          kount_s = kount_s + 1
        else
          rn_norml = rn_norml + prob(ito)*dmto
          rm_norml = rm_norml + prob(ito)*dmto*rmassto
          kount_l = kount_l + 1
        endif

      enddo
c   
c
c  Calculate mass weighting factors <weights,weightl> for small and
c  large cores.
c
      if( kount_s .eq. 0 )then
        weightl = ONE
      elseif( kount_l .eq. 0 )then
        weightl = 0.
      else
        rm_norms = rm_norms/rn_norms
        rm_norml = rm_norml/rn_norml
        weightl = (coreavg - rm_norms) / (rm_norml - rm_norms)
        if( weightl .gt. ALMOST_ONE )then
          weightl = ONE
        elseif( weightl .lt. ALMOST_ZERO )then
          weightl = 0.
        endif
      endif
      weights = ONE - weightl
c
c
c Renormalize probability distribution function and evaluate the CN
c evaporation source term <evappe>.
c
      do ito = inucmin(igto), NBIN

c       if( ito .le. iavg )then
        if( ito .lt. iavg )then      ! Kevin M
          prob(ito) = prob(ito)*weights/rn_norms
        else
          prob(ito) = prob(ito)*weightl/rn_norml
        endif
c
c  First the CN number concentration element
c 
        evappe(ito,ieto) = evappe(ito,ieto) +
     $     evdrop*prob(ito)*dm(ito,igto)

c
c  Now the CN core elements
c
        do ic = 2, ncore(ig)
          iecore = icorelem(ic,ig)
          ie2cn  = ievp2elem(iecore)
          evappe(ito,ie2cn) = evappe(ito,ie2cn) +
     $       rmass(ito,igto)*evcore(ic)*prob(ito)*dm(ito,igto)
        enddo

      enddo
 
      return
      end
