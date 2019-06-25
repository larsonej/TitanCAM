      subroutine coremelt
c
c
c  @(#) coremelt.f  Jensen  Jan-2000
c  This routine evaluates core mass loss rates due to melting in
c  mixed particles <evaplg>:
c  
c  Calculations are done at one spatial point per call.
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
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter coremelt'
c
c
c  Atmospheric density in cgs units:
c
      xyzmet = xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz)
      rhoa_cgs = rhoa3(ixyz) /  xyzmet
c
c
c  Loop over particle groups.
c
      do ielem = 1,NELEM
c
c
c  This calculation is only required for volatile core particle
c  types
c
       if( itype(ielem) .eq. I_VOLCORE ) then

        igroup = igelem( ielem )             ! particle group
        iepart = ienconc( igroup )           ! particle number density element
c
c
c  Bypass calculation if few particles are present 
c
        if( pconmax(ixyz,igroup) .gt. FEW_PC )then
c
c
c  Loop over particle bins.  Loop from largest to smallest for 
c  evaluation of index of smallest bin nucleated during time step <inucstep>.
c
         do ibin =NBIN,1,-1
c
c
c  Calculate the involatile core mass (g/cm^3)
c
          rmass_invcore = 0.
          if( ncore(igroup) .gt. 1 ) then
            do jcore = 1, ncore(igroup)
              iecore = icorelem(jcore,igroup)
              if( itype(iecore) .eq. I_COREMASS ) then
                rmass_invcore = rmass_invcore + pc3(ixyz,ibin,iecore)
              endif
            enddo
          endif
c
c
c  Calculate latent heat due to growth/evaporation of droplets in this bin
c  using <growlg>/<evaplg> calculated in growevapl.f; units: erg/s
c
          rlh_evap = ( growlg(ibin,igroup) - evaplg(ibin,igroup) ) *
     $                 rmass(ibin,igroup) * rlhm(iz,1)
c
c
c  Calculate the equilibrium temperature at the droplet surface by equating
c  the heat transfer throught the liquid layer to (heat conduction in air +
c  latent heat of droplet evaporation).
c
          thcond_h2o = 5.60e4 + (t3(ixyz)-T0)*1.65e2  ! From CRC values

          core_mass_frac = pc3(ixyz,ibin,ielem) /
     $            ( pc3(ixyz,ibin,iepart) *rmass(ibin,igroup) )

          if( core_mass_frac .lt. ONE-ALMOST_ZERO ) then

           r_ice = r(ibin,igroup) * core_mass_frac**(1./3.)
           rmass_ice = rmass(ibin,igroup) * core_mass_frac

           term1 = thcond_h2o * r_ice
           term2 = thcond(iz) * ft(iz,ibin,igroup) *
     $             (r(ibin,igroup)-r_ice)
           term3 = rlh_evap * (r(ibin,igroup)-r_ice)
           t_drop = ( term1*T0 + term2*t3(ixyz) + term3 ) /
     $            ( term1 + term2 )
c
c
c  Calculate core mass melting rate (g/s) by equating heat transfer
c  through the liquid layer to latent heat of melting.  This expression
c  works both for core melting at T > T0 and core growth at T < T0
c
           if( abs(t_drop-T0) .gt. 1.e-8 ) then
             dmi_dt = 4.*PI*r(ibin,igroup)*r_ice*thcond_h2o/
     $              (r(ibin,igroup)-r_ice) * (T0-t_drop) / rlhm(iz,1)
           else
             dmi_dt = 0.
           endif
c
c
c  Decrease ice core mass based on melt rate.  If <dmi_dt> > 0, then
c  use an explicit expression, otherwise use an implicit expression.
c  Do not let the total core mass fraction exceed 1.
c
           pc_premelt = pc3(ixyz,ibin,ielem)
           if( dmi_dt .gt. 0. ) then
             ppd = dmi_dt * pc3(ixyz,ibin,iepart)   ! g/cm^3/s
             rmax_coremass = pc3(ixyz,ibin,iepart) *
     $                       rmass(ibin,igroup) - rmass_invcore
             pc3(ixyz,ibin,ielem) = min( rmax_coremass,
     $                       pc3(ixyz,ibin,ielem) + ppd*dtime )
           else
             pls = -dmi_dt / rmass_ice
             pc3(ixyz,ibin,ielem) = pc3(ixyz,ibin,ielem) /
     $                            ( ONE + pls*dtime )
           endif
c
c
c  Calculate latent heat associated with core melting
c
           rmice_change = pc3(ixyz,ibin,ielem) - pc_premelt
           rlheat = rlheat + rmice_change * rlhm(iz,1) /
     $                       ( CP * rhoa_cgs )

          endif

         enddo     ! end of do ibin loop
        endif    ! end of pconmax .gt. FEW_PC
       endif      ! end of itype(ielem) .eq. I_VOLCORE
      enddo       ! end of ielem = 1, NELEM
c
c
c  Return to caller with particle loss rates due to nucleation evaluated.
c
      return
      end
