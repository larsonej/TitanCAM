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
c
c  Define formats
c
c   1 format(1pe20.15,3x,f4.1,3x,f10.6,3x,1pe13.6)
    1 format(7(1pe13.6,3x))
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
c  Define potential temp to sensible temp conversion pressure ratio exponent for air.
c
      dtdt = 0.e-2
c
c
c  Solve for the new <ptc> due to latent heat exchange and radiative heating.
c  Latent and radiative heating rates are in units of [deg_K/s].
c
c     if(supsati3(ixyz,1) .gt. 0.) then 
c     write(*,1) time/3600.d0/24.d0,zl3(ixyz)/1.d5,
c    $           t3(ixyz),100.*(supsati3(ixyz,1)+1.),
c    $           ptc3(ixyz),rlheat,radheat3(ixyz)
c     endif

      ptc3(ixyz) = ptc3(ixyz) * 
     $   ( 1. + dtime / t3(ixyz) * ( rlheat + radheat3(ixyz) +
     $     dtdt )  )

c
c
c  Update <t>, using latest <rhoa> from mass continuity and <p> from hydrostasy
c 
      pt = ptc3(ixyz) / rhoa3(ixyz)
      t3(ixyz) = pt * ( p3(ixyz) / PREF )**RKAPPA 

c     if(zl3(ixyz) .eq. 58.d5) then 
c     if(supsati3(ixyz,1) .gt. 0.) then 
c     write(*,1) time/3600.d0/24.d0,zl3(ixyz)/1.d5,
c    $           t3(ixyz),100.*(supsati3(ixyz,1)+1.),
c    $           ptc3(ixyz),rhoa3(ixyz),p3(ixyz)
c     endif
 
c  Update <t> from sinusoidal oscillation (NOTE: This ignores any changes from ptc calc)
c    amplitude = 1K, period = 
c
      if( time .ge. tpstart ) then
       if(zl3(ixyz) .ge. 24.d5 .and. zl3(ixyz) .lt. 28.d5) then
c        t3(ixyz) = tinit(ixyz) - 3.d0*sin(1.04d-5 * (time-tpstart))   !p = 1 week
         t3(ixyz) = tinit(ixyz) - 3.d0*sin(2.42d-6 * (time-tpstart))   !p = 1 month
c        t3(ixyz) = tinit(ixyz) + 1.d0*sin(1.21d-6 * time)   !p = 2 months
c        t3(ixyz) = tinit(ixyz) + 1.d0*sin(1.99d-7 * time)   !p = 1 year
c        t3(ixyz) = tinit(ixyz) + 1.d0*sin(1.99d-8 * time)   !p = 10 years
c        phist = 1.8d3  ! half hour
         if( 1.d0+supsati3(ixyz,2) .gt. 0.98 ) then
           phist = 3.6d3  ! hour
         else
           phist = 8.64d4 / 2.
         endif
c        phist = 7.2d3  ! 2 hours
cc       if( t3(ixyz) .gt. tlast ) t3(ixyz) = tlast
cc       tlast = t3(ixyz)
       endif
cc    else
cc     if(zl3(ixyz) .ge. 24.d5 .and. zl3(ixyz) .lt. 26.d5)
cc   $    tlast = t3(ixyz)
      endif
c
c
c  Return to caller with new temperature.
c
      return
      end
