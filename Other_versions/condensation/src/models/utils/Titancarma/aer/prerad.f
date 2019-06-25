       subroutine prerad
c
c
c  @(#) prerad.f  Ackerman  Oct-1997
c  This routine loads arrays into the radiation interface common block
c  for a single column.
c  Note that vertical index in radiative transfer model domain is reversed
c  for cartesian coordinates.
c
c  Indices <ix> and <iy> are passed through global common block.
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter prerad'
c
c
c  Load profiles of temperature [K], water vapor [g/cm^2], and
c  aerosol particle concentrations [#/cm^2]
c
      igas = 1

      do iz = 1,NZ

        xymet = xmet(ix,iy,iz)*ymet(ix,iy,iz)
c
c
c  Reverse the vertical index when in cartesian coordinates
c
        if( igridv .eq. I_CART )then
           jz = NZ + 1 - iz
        else
           jz = iz
        endif

        t_aerad(jz) = t(ix,iy,iz)
        qv_aerad(jz) = gc(ix,iy,iz,igas) / rhoa(ix,iy,iz)

        do igroup = 1,NGROUP
          iep = ienconc(igroup)
          do ibin = 1,NBIN
            pc_aerad(jz,ibin,igroup) = pc(ix,iy,iz,ibin,iep) *
     $                                 dz(ix,iy,iz) / xymet
          enddo
        enddo
      enddo
c
c
c  Compute <u0> = cos( solar_zenith_angle ) 
c
      if( isolar_zen .eq. I_FIXED )then

        u0 = u0_fixed

      elseif( isolar_zen .eq. I_DIURNAL )then
c
c
c  <sun_angle> is solar hour angle from noon [rad]
c
        sun_angle = PI + ( time + rad_start )*2.*PI/SCDAY

        u0 = zsin(ix,iy) + zcos(ix,iy)*cos(sun_angle)
        u0 = max( 0.*ONE, u0 )

      endif
c
c
c  Load <u0> into interface common block
c
      u0_aerad = u0
c
c      
c  Turn on solar radiation only when sun is above horizon
c
      if( do_solar .and. ( u0 .gt. 0. ) )then
        isl_aerad = 1
      else
        isl_aerad = 0
      endif
c
c
c  Set heating rates to zero
c
      call zerorad
c
c
c  Return to caller with arrays loaded into radiation interface common block
c
      return
      end
