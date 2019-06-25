       subroutine rhopart
c
c
c  @(#) rhopart.f  Jensen  Oct-1995
c  This routine calculates new average particle densities.
c
c  The particle mass density can change at each time-step due to
c  changes in the core mass fraction.
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
c  Local declarations
c
      dimension rmshell(NXYZ,NBIN), vcore(NXYZ,NBIN)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter rhopart'
c
c
c  Calculate average particle mass density for each group
c
      do igroup = 1,NGROUP
c
c
c  Only need to do this if particle group has at least 1 core
c
       if( ncore(igroup) .ge. 1 )then
c
c
c  Define particle # concentration element index for current group
c
        iepart = ienconc(igroup)     ! element of particle number concentration 
c
c  Calculate volume of cores and the mass of shell material
c  <vcore> is the volume of core material and <rmshell> is the
c  mass of shell material.
c
        do ibin = 1,NBIN
          vcore(ixyz,ibin) = 0.
          rmshell(ixyz,ibin) = rmass(ibin,igroup) *
     $                         pc3(ixyz,ibin,iepart)
        enddo

        do jcore = 1,ncore(igroup)

          iecore = icorelem(jcore,igroup)    ! core element

          do ibin = 1,NBIN

            vcore(ixyz,ibin) = vcore(ixyz,ibin) + 
     $                pc3(ixyz,ibin,iecore) / rhoelem(iecore)
            rmshell(ixyz,ibin) = rmshell(ixyz,ibin) -
     $                           pc3(ixyz,ibin,iecore)

          enddo

        enddo
c
c
c  Calculate average density
c
        do ibin = 1,NBIN

          if( pc3(ixyz,ibin,iepart) .gt. SMALL_PC ) then
            rhop3(ixyz,ibin,igroup) = rmass(ibin,igroup) *
     $        pc3(ixyz,ibin,iepart) /
     $        ( vcore(ixyz,ibin) + rmshell(ixyz,ibin)/rhoelem(iepart) )
          else
            rhop3(ixyz,ibin,igroup) = rhoelem(iepart)
          endif

        enddo

       endif

      enddo
c
c
c  Return to caller with new particle number densities.
c
      return
      end
