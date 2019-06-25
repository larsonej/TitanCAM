       subroutine tsolve
c
c  @(#) tsolve.f  Ackerman Oct-1997
c
c  This routine calculates new potential temperature concentration 
c  (and updates temperature) due to microphysical and radiative forcings.
c  The equation solved (the first law of thermodynamics in flux form) is
c
c    d(rhostar theta)     rhostar theta       d(qv)    1  dF
c    ---------------  = - ------------- * ( L ----- + --- -- )
c          dt                 Cp T              dt    rho dz
c
c  where
c    rhostar = scaled air density
c    theta   = potential temperature
c    t       = time
c    Cp      = specific heat (at constant pressure) of air
c    T       = air temperature
c    qv      = water vapor mixing ratio
c    L       = latent heat
c    F       = net radiative flux
c    z       = unscaled altitude
c   
c  Modified  Sep-1997  (McKie)
c  To calculate at one spatial point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter tsolve'
c
c-------------------------------------------------------------------------------
c
c
c  Solve for the new <ptc> due to latent heat exchange and radiative heating.
c  Latent and radiative heating rates are in units of [deg_K/s].
c
      ptc3(ixyz) = ptc3(ixyz) * 
     $   ( 1. + dtime / t3(ixyz) * ( rlheat + radheat3(ixyz) ) )
c
c
c  Update <t>, using latest <rhoa> from mass continuity and <p> from hydrostasy
c 
      pt = ptc3(ixyz) / rhoa3(ixyz)
      t3(ixyz) = pt * ( p3(ixyz) / PREF )**RKAPPA 
c
c
c  Return to caller with new temperature.
c
      return
      end
