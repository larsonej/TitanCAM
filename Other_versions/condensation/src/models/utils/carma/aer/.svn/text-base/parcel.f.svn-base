       subroutine parcel
c
c
c  @(#) parcel.f  Ackerman  Sep-1997
c
c  This routine calculates new altitudes, pressures, air densities,
c  and temperatures due to adiabatic forcing to a Lagrangian parcel.
c  The particle losses due to sedimentation are also applied
c  (there are no sources due to sedimentation for parcel model).
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
c   Define number of iterations for evaluating potential temperature and
c   temperature.
c
      parameter( N_ITERATE = 3 )
c
c
c   Define whether or not to allow sedimentation
c
      logical DO_SED
      parameter( DO_SED = .false. )

c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter parcel'
c
c-------------------------------------------------------------------------------
c

      ixy = 1
c
c
c  Apply sedimentation losses to particles (only when layer thickness
c  is non-zero).
c
      if( DO_SED )then
        do ielem = 1,NELEM
          ig = igelem(ielem)
          do ibin = 1,NBIN
            do iz = 1,NZ
              if ( dz2(ixy,iz) .gt. 0. )then
                pc2(ixy,iz,ibin,ielem) = pc2(ixy,iz,ibin,ielem) /
     $              ( 1. + dtime * vf(iz,ibin,ig) / dz2(ixy,iz) )
              endif
            enddo
          enddo
        enddo
      endif
c
c
c  Update altitudes, then iteratively update <p>, <t>, <rhoa>, and <ptc>.
c
      zl2(ixy,NZP1) = zl2(ixy,NZP1) + w2(ixy,NZP1)*dtime

      do iz = 1,NZ

        delz = w2(ixy,iz)*dtime
        zc2(ixy,iz) = zc2(ixy,iz) + delz
        zl2(ixy,iz) = zl2(ixy,iz) + delz
 
        pt = ptc2(ixy,iz) / rhoa2(ixy,iz)
        p_start = p2(ixy,iz)
 
        do i = 1,N_ITERATE
          p2(ixy,iz) = p_start - rhoa2(ixy,iz)*GRAV*delz
          t2(ixy,iz) = pt * ( p2(ixy,iz) / PREF )**RKAPPA 
          rhoa2(ixy,iz) = p2(ixy,iz) / ( R_AIR * t2(ixy,iz) )
        enddo
 
        pt = t2(ixy,iz) * ( PREF / p2(ixy,iz) )**RKAPPA 
        ptc2(ixy,iz) = pt * rhoa2(ixy,iz)

      enddo
c
c
c  Change particle and gas concentrations to account for adiabatic
c  expansion of the parcel
c
      do ixyz = 1, NXYZ
        do igas = 1, NGAS
          gc3(ixyz,igas) = gc3(ixyz,igas) *
     $                     rhoa3(ixyz)/rhoaold(ixyz)
        enddo
        do ibin = 1,NBIN
          do ielem = 1,NELEM
            pc3(ixyz,ibin,ielem) = pc3(ixyz,ibin,ielem) *
     $                             rhoa3(ixyz)/rhoaold(ixyz)
          enddo
        enddo
      enddo
c
c
c  Return to caller with new temperature.
c
      return
      end
