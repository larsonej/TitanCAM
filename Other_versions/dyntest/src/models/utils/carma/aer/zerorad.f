      subroutine zerorad
c
c
c  @(#) zerorad.f  Ackerman  Oct-1997
c  This routine zeroes the radiative transfer heating rates.
c
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
c
c
c  Include global constants and variables.
c
      include 'globaer.h'
c
c
c  Announce entry to this routine.
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter zerorad'
c
c
c   Loop over all grid points.
c
      do ixyz = 1, NXYZ
c
c
c   Heating rates for air.
c
        radheat3(ixyz) = 0.
c
c
c   Heating rates for particles
c
        do igroup = 1,NGROUP
          iep = ienconc(igroup)
          do ibin = 1,NBIN
            qrad3(ixyz,ibin,igroup) = 0.
          enddo
        enddo

      enddo
c
c
c   Albedos and optical depths 
c
      do ixy = 1, NXY

        alb_toai2(ixy) = 0.
        alb_tomi2(ixy) = 0.

        do isol = 1,NSOL
          alb_toa2(ixy,isol) = 0.
          opd2(ixy,isol) = 0.
        enddo

      enddo
c
c
c  Return to caller with radiative heating rates etc zeroed.
c
      return
      end
