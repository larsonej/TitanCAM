      subroutine growevapl
c
c
c  @(#) growevapl.f  Ackerman  Dec-1995
c  This routine evaluate particle loss rates due to condensational
c  growth and evaporation for all condensing gases.
c
c  The loss rates for each group are <growlg> and <evaplg>.
c
c  Units are [s^-1].
c
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c
c
c  Argument list input:
c
c  Argument list output:
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Local declarations
c 
      dimension otherm(NELEM), ieother(NELEM),
     $  dela(NBIN), delma(NBIN), aju(NBIN),
     $  ar(NBIN), al(NBIN), a6(NBIN),dHv(NGAS),
     $  alpha(NBIN,NGAS),beta(NBIN,NGAS),gamma(NBIN,NGAS),
     $  dmdt(NBIN),gromat(NGAS,NGAS+1),xx(NGAS+1),
     $  fracmol(NELEM),elemmol(NELEM),
     $  dmdt_grow(NBIN),dmdt_evap(NBIN),rnormfrac(NBIN,NELEM)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter growevapl'
c
c
c  Calculate particle fluxes due to growth and evaporation for all
c  particle groups.
c
      do igroup = 1,NGROUP

       iepart = ienconc(igroup)     ! element of particle number concentration 
       imaingas = igrowgas(iepart)      ! condensing gas
       small = rmass(1,igroup) * SMALL_PC

       if( imaingas .ne. 0 )then
c
c
c  Bypass calculation if few particles are present 
c
        if( pconmax(ixyz,igroup) .gt. FEW_PC )then

         do ibin = 1,NBIN-1
c
c
c  Consider growth of average particle at radius <rup(ibin)>.
c 
c  Treat solute effect first: <asol> is solute factor.
c
c  Only need to treat solute effect if <nelemg(igroup)> > 1
c
          if( nelemg(igroup) .le. 1 )then

           argsol = 0.

          else
c
c  <condm> is mass concentration of condensed gas <igas> in particle.
c  <nother> is number of other elements in group having mass.
c  <otherm> are mass concentrations of other elements in particle group.
c  <othermtot> is total mass concentrations of other elements in particle.
c
           nother = 0
           othermtot = 0.
c
c
c  <ieoth_rel> is relative element number of other element in group.
c
           do ieoth_rel  = 2,nelemg(igroup)       
c
c
c  <ieoth_abs> is absolute element number of other element.
c
             ieoth_abs = iepart + ieoth_rel - 1    

             if( itype(ieoth_abs) .eq. I_COREMASS )then
              nother = nother + 1
              ieother(nother) = ieoth_abs
              otherm(nother) = pc3(ixyz,ibin,ieoth_abs)
              othermtot = othermtot + otherm(nother)
             endif

           enddo

           condm = rmass(ibin,igroup)*pc3(ixyz,ibin,iepart) - othermtot

           if( condm .le. 0. )then
c
c
c  Zero mass for the condensate -- <asol> is a small value << 1
c
             argsol = 1e6     

           else
c
c
c  Sum over masses of other elements in group for argument of solute factor.
c
             argsol = 0.
             do jother = 1,nother
               isol = isolelem(ieother(jother))
               argsol = argsol 
     $                + sol_ions(isol)*otherm(jother)/solwtmol(isol)
             enddo 
             argsol = argsol*gwtmol(igas)/condm

           endif 

          endif    ! nelemg(igroup) > 1
c
c
c  <akas> is combined kelvin (curvature) and solute factors.
c
c  Ignore solute factor for ice particles.
c
          if( is_grp_ice(igroup) )then
            expon = akelvini(iz,igas) / rup(ibin,igroup)
          else
            expon = akelvin(iz,igas)  / rup(ibin,igroup) - argsol 
          endif
          expon = max(-POWMAX, expon)
          akas = exp( expon )
c 
c 
c  Trick for removing haze droplets from droplet bins:
c  allows haze droplets to exist under supersaturated conditions;
c  when below supersaturation, haze droplets will evaporate.
c
          if( (.not. is_grp_ice(igroup)) .and. (akas .lt. 1.) .and.
     $        (supsatl3(ixyz,igas) .lt. 0.) ) akas = 1.
c
c
c  <dmdt> is growth rate in mass space [g/s].
c
          g0 =  gro(iz,ibin+1,igroup)
          g1 = gro1(iz,ibin+1,igroup)
          g2 = gro2(iz,igroup)
 
          if( is_grp_ice(igroup) )then
            ss = supsati3(ixyz,igas)
            pvap = pvapi3(ixyz,igas)
          else
            ss = supsatl3(ixyz,igas)
            pvap = pvapl3(ixyz,igas)
          endif
            
          dmdt(ibin) = pvap * ( ss + ONE - akas - 
     $                 qrad3(ixyz,ibin+1,igroup) * g1 * g2 ) *
     $                 g0 / ( 1. + g0 * g1 * pvap )

         enddo     ! ibin = 1,NBIN-1
c
c
c  Now calculate condensation/evaporation production and loss rates.
c  Use Piecewise Polynomial Method [Colela and Woodard, J. Comp. Phys.,
c  54, 174-201, 1984]
c
c  First, use cubic fits to estimate concentration values at bin
c  boundaries
c
         do ibin = 2,NBIN-1

          dpc = pc3(ixyz,ibin,iepart) / dm(ibin,igroup)
          dpc1 = pc3(ixyz,ibin+1,iepart) / dm(ibin+1,igroup)
          dpcm1 = pc3(ixyz,ibin-1,iepart) / dm(ibin-1,igroup)
          ratt1 = pratt(1,ibin,igroup)
          ratt2 = pratt(2,ibin,igroup)
          ratt3 = pratt(3,ibin,igroup)
          dela(ibin) = ratt1 *
     $             ( ratt2*(dpc1-dpc) + ratt3*(dpc-dpcm1) )
          delma(ibin) = 0.

          if( (dpc1-dpc)*(dpc-dpcm1) .gt. 0. )
     $       delma(ibin) = min( abs(dela(ibin)), 2.*abs(dpc-dpc1),
     $            2.*abs(dpc-dpcm1) ) * abs(dela(ibin))/dela(ibin)

         enddo     ! ibin = 2,NBIN-2

         do ibin = 2,NBIN-2

          dpc = pc3(ixyz,ibin,iepart) / dm(ibin,igroup)
          dpc1 = pc3(ixyz,ibin+1,iepart) / dm(ibin+1,igroup)
          dpcm1 = pc3(ixyz,ibin-1,iepart) / dm(ibin-1,igroup)
          rat1 = prat(1,ibin,igroup)
          rat2 = prat(2,ibin,igroup)
          rat3 = prat(3,ibin,igroup)
          rat4 = prat(4,ibin,igroup)
          den1 = pden1(ibin,igroup)
c
c
c  <aju(ibin)> is the estimate for concentration (dn/dm) at bin
c  boundary <ibin>+1/2.
c
          aju(ibin) = dpc + rat1*(dpc1-dpc) + 1./den1 *
     $             ( rat2*(rat3-rat4)*(dpc1-dpc) -
     $             dm(ibin,igroup)*rat3*delma(ibin+1) +
     $             dm(ibin+1,igroup)*rat4*delma(ibin) )

         enddo     ! ibin = 2,NBIN-2
c
c
c  Now construct polynomial functions in each bin
c
         do ibin = 3,NBIN-2

          al(ibin) = aju(ibin-1)
          ar(ibin) = aju(ibin)

         enddo
c
c
c  Use linear functions in first two and last two bins
c
         if( NBIN .gt. 1 )then
         ar(2) = aju(2)
         al(2) = pc3(ixyz,1,iepart)/dm(1,igroup) +
     $           palr(1,igroup) *
     $           (pc3(ixyz,2,iepart)/dm(2,igroup)-
     $           pc3(ixyz,1,iepart)/dm(1,igroup))
         ar(1) = al(2)
         al(1) = pc3(ixyz,1,iepart)/dm(1,igroup) +
     $           palr(2,igroup) *
     $           (pc3(ixyz,2,iepart)/dm(2,igroup)-
     $           pc3(ixyz,1,iepart)/dm(1,igroup))

         al(NBIN-1) = aju(NBIN-2)
         ar(NBIN-1) = pc3(ixyz,NBIN-1,iepart)/dm(NBIN-1,igroup) +
     $           palr(3,igroup) *
     $           (pc3(ixyz,NBIN,iepart)/dm(NBIN,igroup)-
     $           pc3(ixyz,NBIN-1,iepart)/dm(NBIN-1,igroup))
         al(NBIN) = ar(NBIN-1)
         ar(NBIN) = pc3(ixyz,NBIN-1,iepart)/dm(NBIN-1,igroup) +
     $           palr(4,igroup) *
     $           (pc3(ixyz,NBIN,iepart)/dm(NBIN,igroup)-
     $           pc3(ixyz,NBIN-1,iepart)/dm(NBIN-1,igroup))
         endif
c
c
c  Next, ensure that polynomial functions do not deviate beyond the
c  range [<al(ibin)>,<ar(ibin)>]
c
         do ibin = 1,NBIN

          dpc = pc3(ixyz,ibin,iepart) / dm(ibin,igroup)

          if( (ar(ibin)-dpc)*(dpc-al(ibin)) .le. 0. )then
            al(ibin) = dpc
            ar(ibin) = dpc
          endif

          test1 = (ar(ibin)-al(ibin))*(dpc - 0.5*(al(ibin)+ar(ibin)))
          test2 = 1./6.*(ar(ibin)-al(ibin))**2

          if( test1 .gt. test2 )then
             al(ibin) = 3.*dpc - 2.*ar(ibin)
          elseif( test1 .lt. -test2 )then
             ar(ibin) = 3.*dpc - 2.*al(ibin)
          endif

         enddo
c
c
c  Lastly, calculate fluxes across each bin boundary.
c
c  Use upwind advection when courant number > 1.
c  
         do ibin = 1,NBIN

          dpc = pc3(ixyz,ibin,iepart) / dm(ibin,igroup)
          dela(ibin) = ar(ibin) - al(ibin)
          a6(ibin) = 6. * ( dpc - 0.5*(ar(ibin)+al(ibin)) )

         enddo

         do ibin = 1,NBIN-1

          if( dmdt(ibin) .gt. 0 .and.
     $        pc3(ixyz,ibin,iepart) .gt. SMALL_PC )then

           x = dmdt(ibin)*dtime/dm(ibin,igroup)

           if( x .lt. 1. )then
            growlg(ibin,igroup) = dmdt(ibin)/pc3(ixyz,ibin,iepart)
     $                 * ( ar(ibin) - 0.5*dela(ibin)*x +
     $                 (x/2. - x**2/3.)*a6(ibin) )
           else
            growlg(ibin,igroup) = dmdt(ibin) / dm(ibin,igroup)
           endif

          elseif( dmdt(ibin) .lt. 0 .and.
     $        pc3(ixyz,ibin+1,iepart) .gt. SMALL_PC )then

           x = -dmdt(ibin)*dtime/dm(ibin+1,igroup)

           if( x .lt. 1. )then
            evaplg(ibin+1,igroup) = -dmdt(ibin)/
     $                 pc3(ixyz,ibin+1,iepart)
     $                 * ( al(ibin+1) + 0.5*dela(ibin+1)*x +
     $                 (x/2. - (x**2)/3.)*a6(ibin+1) )
           else
            evaplg(ibin+1,igroup) = -dmdt(ibin) / dm(ibin+1,igroup)
           endif
c
c
c  Boundary conditions: for evaporation out of first bin (with cores), 
c  use evaporation rate from second bin.
c  
           if( ibin .eq. 1 .and. ncore(igroup) .gt. 0 )then
             evaplg(1,igroup) = -dmdt(1) / dm(1,igroup)
           endif

          endif
c
c
c  Limit growth rates to "reasonable" values
c
          growlg(ibin,igroup) = min( growlg(ibin,igroup), 1e10*ONE )
          evaplg(ibin+1,igroup) = min( evaplg(ibin+1,igroup),
     $                                 1e10*ONE )

         enddo    ! ibin = 1,NBIN-1

        endif     ! (pconmax .gt. FEW_PC)
       endif      ! (igas = igrowgas(ielem)) .ne. 0 
      enddo       ! igroup = 1,NGROUP
c
c
c  Return to caller with particle loss rates for growth and evaporation
c  evaluated.
c
      return
      end
