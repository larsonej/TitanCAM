       subroutine hydrostat
c
c
c  @(#) hydrostat.f  Ackerman  Jul-1997
c
c    This routine updates pressure by hydrostatic integration.
c    In sigma coordinates, it also updates vertical metric scale
c    factor <zmet>, temperature <t>, and scaled air density <rhoa>.
c
c  Argument list input:
c
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
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter hydrostat'
c
c
c-------------------------------------------------------------------------------
c
      return
c
c  Scan the spatial points
c
      do ixy = 1,NXY
c
c
c  cartesian coordinates
c
        if( igridv .eq. I_CART )then
c
c
c  <ptop> is pressure at top of layer
c
          iz = NZ
          ptop = p_top2(ixy)
          dp = dz2(ixy,iz)*GRAV*rhoa2(ixy,iz)
          p2(ixy,iz) = ptop * sqrt( 1. + dp/ptop )
   
          do iz = NZ-1,1,-1
            ptop = ptop + dp
            dp = dz2(ixy,iz)*GRAV*rhoa2(ixy,iz)
            p2(ixy,iz) = ptop * sqrt( 1. + dp/ptop )
          enddo
          p_surf2(ixy) = ptop + dp
c
c
c  Sigma coordinates: first integrate for total column mass, then
c  update pressures and scaled air density, then temperatures, and then
c  air density (in cgs units) and vertical metric scale factor.
c
        else if( igridv .eq. I_SIG )then

          pstar = 0.
          do iz = 1,NZ
            xymet = xmet2(ixy,iz) * ymet2(ixy,iz)
            pstar = pstar + dz2(ixy,iz) * GRAV * rhoa2(ixy,iz) / xymet
          enddo
          p_surf2(ixy) = p_top2(ixy) + pstar

          do iz = 1,NZ
            xymet = xmet2(ixy,iz) * ymet2(ixy,iz)
            p2(ixy,iz) = p_top2(ixy) + zc2(ixy,iz) * pstar
            rhoa2(ixy,iz) = pstar / GRAV * xymet
            t2(ixy,iz) = ptc2(ixy,iz) / rhoa2(ixy,iz) * 
     $                   ( p2(ixy,iz) / PREF )**(R_AIR/CP)
            rhoa_cgs = p2(ixy,iz) / ( R_AIR * t2(ixy,iz) )
            zmet2(ixy,iz) = pstar / ( GRAV * rhoa_cgs )
          enddo

        endif
 
      enddo
c
c
c  Return to caller with pressures hydrostatically balanced.
c
      return
      end
