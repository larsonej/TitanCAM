       subroutine setupgkern
c
c
c  @(#) setupgkern.f  Ackerman  Dec-1995
c  This routine defines radius-dependent but time-independent parameters
c  used to calculate condensational growth of particles.  Growth rates
c  are calculated at bin boundaries: the parameters calculated here 
c  ( <gro>, <gro1>, <gro2>, and <akelvin> ) 
c  are defined at lower bin boundaries through the growth rate expression
c  (for one particle) used in growevapl.f:
c
c    dm = gro*pvap*( S + 1 - Ak*As - gro1*gro2*qrad )
c    --   -------------------------------------------
c    dt               1 + gro*gro1*pvap
c
c   where 
c
c   S    = supersaturation
c   Ak   = exp(akelvin/r)
c   As   = exp(-sol_ions * solute_mass/solwtmol * gwtmol/condensate_mass)
c   pvap = saturation vapor pressure [dyne cm**-2]
c   qrad = radiative energy absorbed
c
c  This routine requires that vertical profiles of temperature <T>,
c  and pressure <p> are defined (initatm.f must be called before this).
c
c  This routine also requires that particle Reynolds' numbers are
c  defined (setupvfall.f must be called before this).
c
c  The vertical profile with ix = iy = 1 is used.
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
c  Local declarations
c    None.
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupgkern'
c
c-------------------------------------------------------------------------------
c
c
c  Specify radius-independent parameters.
c
c  Sticking coefficient for growth (also called mass accomodation coefficient)
c
      gstickl = 1.0
      gsticki = 0.93
c
c
c  Thermal accommodation coefficient
c
      tstick = 1.0
c
c-------------------------------------------------------------------------------
c
c
c   Loop over aerosol groups only (no radius, gas, or spatial dependence).
c   
c
      do igroup = 1, NGROUP
c
c   Use gstickl or gsticki, depending on whether group is ice or not
c
        if( is_grp_ice(igroup) ) then
          gstick = gsticki
        else
          gstick = gstickl
        endif
c
c   Non-spherical corrections (need a reference for these)
c
        if( ishape(igroup) .eq. 1 )then
c
c   Spheres
c
          cor = 1.
          phish = 1.

        else 

          if( ishape(igroup) .eq. 2 )then
c
c  Hexagons
c
            phish = 6./PI*tan(PI/6.)*( eshape(igroup) + 0.5 )
     $         * ( PI / ( 9.*eshape(igroup)*tan(PI/6.) ) )**(ONE*2./3.)

          else if( ishape(igroup) .eq. 3 )then
c
c  Spheroids
c
            phish = ( eshape(igroup) + 0.5 )
     $         * ( 2. / ( 3.*eshape(igroup) ) )**(ONE*2./3.)

          endif

          if( eshape(igroup) .lt. 1. )then
c
c  Oblate spheroids
c
            esh1 = 1. / eshape(igroup)
            a1 = sqrt(esh1**2 - 1.)
            cor = a1 / asin( a1 / esh1 ) / esh1**(ONE*2./3.)

          else 
c
c  Prolate spheroids
c
            a1 = sqrt( eshape(igroup)**2 - 1. )
            cor = a1 / log( eshape(igroup) + a1 ) 
     $          / eshape(igroup)**(ONE/3.)

          endif
        endif
c
c  Evaluate growth terms only for particle elements that grow.
c EJL 3-27-13 Added loop
       do ielemg = 1,nelemg(igroup)

        ielem = ienconc(igroup) - 1 + ielemg    ! particle number concentration element
        igas = igrowgas(ielem)     ! condensing gas is <igas>

        if( igas .ne. 0 )then
c
c  Loop over vertical layers (use column ix = iy = 1)
c
          ix = 1
          iy = 1

          do k = 1, NZ
c
c  Radius-independent parameters for condensing gas
c
c  This is <rhoa> in cgs units.
c
           rhoa_cgs = rhoa(ix,iy,k) /
     $                (xmet(ix,iy,k)*ymet(ix,iy,k)*zmet(ix,iy,k))
c
c EJL 2-12-13 added gasnames
c
          if (gasname(igas) .eq. 'meth' )then
 
            TcTPPD = 190.7 *ONE
            T1TPPD = 100.0 *ONE
            sig1TPPD = 16.33 *ONE
 
          else if (gasname(igas) .eq. 'eth' )then   

            TcTPPD = 305.42 *ONE
            T1TPPD = 153.2 *ONE
            sig1TPPD = 21.157 *ONE

          else
c
c   Condensing gas is not yet configured.
c
            write(LUNOPRT,'(/,a)') 'invalid igas in setupgkern.f'
            stop 1

          endif
c
c  <surfctwa> is surface tension of liquid-air interface
c  from Thermodynamic & Physical Property Data (Yaws 1992)
c
          surfctwa(k,igas) = sig1TPPD * ( (TcTPPD - t3(k)) /
     $                              (TcTPPD - T1TPPD) )**(11./9.)
c
c
c  <surfctia> is surface tension of ice-air interface
c  from Guez et al. 1992, Planet. Space Sci. 45(6) p.611-625
c
c         Latent heat of sublimation 
          rlhs = rlhm(k,igas) + rlhe(k,igas) 

          surfctia(k,igas) = (rlhs / rlhe(k,igas) )**2 * 
     $                                    surfctwa(k,igas)

c
c  <surfctiw> is surface tension of liquid-ice interface
c  (Antonoff's rule)

          surfctiw(k,igas) = surfctia(k,igas) - surfctwa(k,igas)
c
c EJL 2-12-13 changed below
c  <akelvin> is argument of exponential in kelvin curvature term.
c
            akelvin(k,igas) = 2.*gwtmol(igas)*surfctwa(k,igas) 
     $                      / ( t(ix,iy,k)*rhoelem(ielem)*RGAS )  !RHO_W
 
            akelvini(k,igas) = 2.*gwtmol(igas)*surfctia(k,igas) 
     $                      / ( t(ix,iy,k)*rhoelem(ielem)*RGAS )
 
c
c
c   Molecular free path of condensing gas 
c
          freep  = 3.*diffus(k,igas)
     $           * sqrt( ( PI*gwtmol(igas) ) / ( 8.*RGAS*t(ix,iy,k) ) )
c
c   Thermal free path of condensing gas
c
          freept = freep*thcond(k) /
     $             ( diffus(k,igas) * rhoa_cgs
     $           * ( CP - RGAS/( 2.*WTMOL_AIR ) ) )
c
c   Latent heat of condensing gas 
c
          if( is_grp_ice(igroup) )then
            rlh = rlhe(k,igas) + rlhm(k,igas)
          else
            rlh = rlhe(k,igas)
          endif
c
c    Radius-dependent parameters 
c
          do i = 1, NBIN

           br = rlow(i,igroup)     ! particle bin Boundary Radius
c
c  These are Knudsen numbers
c
           rknudn  = freep /br
           rknudnt = freept/br
c
c  These are "lambdas" used in correction for gas kinetic effects.
c
           rlam  = ( 1.33*rknudn  + 0.71 ) / ( rknudn  + 1. ) 
     $           + ( 4.*( 1. - gstick ) ) / ( 3.*gstick )

           rlamt = ( 1.33*rknudnt + 0.71 ) / ( rknudnt + 1. ) 
     $           + ( 4.*( 1. - tstick ) ) / ( 3.*tstick )
c
c  Diffusion coefficient and thermal conductivity modified for
c  free molecular limit and for particle shape.
c
           diffus1 = diffus(k,igas)*cor / ( 1. + rlam*rknudn*cor/phish )
           thcond1 = thcond(k)*cor / ( 1. + rlamt*rknudnt*cor/phish )
c
c  Reynolds' number based on particle shape <reyn_shape>
c
           if( ishape(igroup) .eq. 1 )then
            reyn_shape = re(k,i,igroup)

           else if( eshape(igroup) .lt. 1. )then
             reyn_shape = re(k,i,igroup) * ( 1. + 2.*eshape(igroup) )

           else
            reyn_shape = re(k,i,igroup) * PI*( 1.+2.*eshape(igroup) )
     $                 / ( 2.*( 1. + eshape(igroup) ) )
           endif
c
c   Particle Schmidt number
c   
           schn = rmu(k) / ( rhoa_cgs * diffus1 )
c
c   Prandtl number
c
           prnum = rmu(k)*CP/thcond1
c
c  Ventilation factors <fv> and <ft> from Pruppacher and Klett
c  
           x1 = schn **(ONE/3.) * sqrt( reyn_shape )
           x2 = prnum**(ONE/3.) * sqrt( reyn_shape )

           if( is_grp_ice(igroup) )then
c
c   Ice crystals
c
            if( x1 .le. 1. )then
              fv = 1.   + 0.14*x1**2
            else
              fv = 0.86 + 0.28*x1
            endif

            if( x2 .le. 1. )then
              ft(k,i,igroup) = 1.   + 0.14*x2**2
            else
              ft(k,i,igroup) = 0.86 + 0.28*x2
            endif

           else
c
c   Liquid water drops
c
            if( x1 .le. 1.4  )then
              fv = 1.   + 0.108*x1**2
            else
              fv = 0.78 + 0.308*x1
            endif

            if( x2 .le. 1.4 )then
              ft(k,i,igroup) = 1.   + 0.108*x2**2
            else
              ft(k,i,igroup) = 0.78 + 0.308*x2
            endif

           endif
c
c  Growth kernel for particle without radiation or heat conduction at
c  radius lower boundary [g cm^3 / erg / s]
c
           gro(k,i,igroup) = 4.*PI*br
     $                     * diffus1*fv*gwtmol(igas)
     $                     / ( BK*t(ix,iy,k)*AVG )
c
c  Coefficient for conduction term in growth kernel [s/g]
cEJL 2-12-13 changed igroup to ielem and below
           gro1(k,i,ielem) = gwtmol(igas)*rlh**2
     $         / ( RGAS*t(ix,iy,k)**2*ft(k,i,igroup)*thcond1 )
     $         / ( 4.*PI*br )
c
c  Coefficient for radiation term in growth kernel [g/erg]
c  (note: no radial dependence).
c
           if( i .eq. 1 )then
             gro2(k,ielem) = 1. / rlh
c             gro3(k,ielem) = t3(k)*gro2(k,ielem)
           endif
         
          enddo   ! i=1,NBIN
         enddo    ! k=1,NZ
c
c
c  Compute radius-dependent terms used in PPM advection scheme
c
          do i = 2,NBIN-1
            pratt(1,i,igroup) = dm(i,igroup) /
     $            ( dm(i-1,igroup) + dm(i,igroup) + dm(i+1,igroup) )
            pratt(2,i,igroup) = ( 2.*dm(i-1,igroup) + dm(i,igroup) ) /
     $            ( dm(i+1,igroup) + dm(i,igroup) )
            pratt(3,i,igroup) = ( 2.*dm(i+1,igroup) + dm(i,igroup) ) /
     $            ( dm(i-1,igroup) + dm(i,igroup) )
          enddo

          do i = 2,NBIN-2
            prat(1,i,igroup) = dm(i,igroup) /
     $              ( dm(i,igroup) + dm(i+1,igroup) )
            prat(2,i,igroup) = 2. * dm(i+1,igroup) * dm(i,igroup) /
     $             ( dm(i,igroup) + dm(i+1,igroup) )
            prat(3,i,igroup) = ( dm(i-1,igroup) + dm(i,igroup) ) /
     $             ( 2.*dm(i,igroup) + dm(i+1,igroup) )
            prat(4,i,igroup) = ( dm(i+2,igroup) + dm(i+1,igroup) ) /
     $             ( 2.*dm(i+1,igroup) + dm(i,igroup) )
            pden1(i,igroup) = dm(i-1,igroup) + dm(i,igroup) +
     $             dm(i+1,igroup) + dm(i+2,igroup)
          enddo

          if( NBIN .gt. 1 )then
            palr(1,igroup) =
     $           (rmassup(1,igroup)-rmass(1,igroup)) /
     $           (rmass(2,igroup)-rmass(1,igroup))
            palr(2,igroup) =
     $           (rmassup(1,igroup)/rmrat(igroup)-rmass(1,igroup)) /
     $           (rmass(2,igroup)-rmass(1,igroup))
            palr(3,igroup) =
     $           (rmassup(NBIN-1,igroup)-rmass(NBIN-1,igroup))
     $           / (rmass(NBIN,igroup)-rmass(NBIN-1,igroup))
            palr(4,igroup) =
     $           (rmassup(NBIN,igroup)-rmass(NBIN-1,igroup))
     $           / (rmass(NBIN,igroup)-rmass(NBIN-1,igroup))

          endif

        endif     ! igas ne 0
       enddo      ! ielemg = 1,nelem(igrp)
      enddo       ! igroup=1,NGROUP
c
c
c
c  Return to caller with time-independent particle growth 
c  parameters initialized.
c
      return
      end
