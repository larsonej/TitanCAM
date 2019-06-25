       subroutine postrad
c
c
c  @(#) postrad.f  Ackerman  Oct-1997
c  This routine gets output from the radiation interface common block.
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
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter postrad'
c
c
c   Horizontal spatial index
c
      ixy = NX*( iy - 1 ) + ix 
c
c
c   Heating rates for air and particles.
c
      do iz = 1,NZ

        ixyz = NXY*( iz - 1 ) + ixy
c
c
c  Reverse the vertical index when in cartesian coordinates
c
        if( igridv .eq. I_CART )then
           jz = NZ + 1 - iz
        else
           jz = iz
        endif

        if( isl_aerad .eq. 1 ) then
         radheat3(ixyz) = heats_aerad(jz) + heati_aerad(jz)
        else
         radheat3(ixyz) = heati_aerad(jz)
        endif

        do igroup = 1,NGROUP
         iep = ienconc(igroup)
         do ibin = 1,NBIN
          if( isl_aerad .eq. 1 ) then
           qrad3(ixyz,ibin,igroup) = qrad_aerad(ibin,jz,igroup)
          else
           qrad3(ixyz,ibin,igroup) = 0.
          endif
         enddo
        enddo
      enddo
c
c
c   Fluxes
c
      do iz = 1,NZ_RAD

        if( igridv .eq. I_CART )then
           jz = NZ_RAD + 1 - iz
        else
           jz = iz
        endif

        fsl_up2(ixy,iz) = fsl_up_aerad(jz)
        fsl_dn2(ixy,iz) = fsl_dn_aerad(jz)
        fir_up2(ixy,iz) = fir_up_aerad(jz)
        fir_dn2(ixy,iz) = fir_dn_aerad(jz)

      enddo
c
c
c   Albedos and optical depths 
c
      alb_toai2(ixy) = alb_toai_aerad
      alb_tomi2(ixy) = alb_tomi_aerad

      do isol = 1,NSOL
        alb_toa2(ixy,isol) = alb_toa_aerad(isol)
        opd2(ixy,isol) = opd_aerad(isol)
      enddo
c
c
c   Switch wavelength bins 11 and 12)
c
      alb_temp = alb_toa2(ixy,11)
      alb_toa2(ixy,11) = alb_toa2(ixy,12)
      alb_toa2(ixy,12) = alb_temp

      opdtemp = opd2(ixy,11)
      opd2(ixy,11) = opd2(ixy,12)
      opd2(ixy,12) = opdtemp
c
c
c  Return to caller with output copied from the radiation interface 
c  common block
c
      return
      end
