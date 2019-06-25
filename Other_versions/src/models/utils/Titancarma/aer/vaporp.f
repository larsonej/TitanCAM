       subroutine vaporp
c
c
c  @(#) vaporp.f  Ackerman  Dec-1995
c  This routine calculates the vapor pressure for all gases 
c  over the entire spatial grid:
c
c  <pvapl> and <pvapi> are vapor pressures in units of [dyne/cm^2]
c
c  Uses temperature <t> as input.
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
c  Include constants for vapor pressure equations
c
      include 'vpconstants.h'
c
c  Local declarations
c
      dimension pvl_pure(NGAS)
c
c  Define formats
c
    1 format('T = ',1pe12.3,a,i6,a,1pe11.3)
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter vaporp'
c
c-------------------------------------------------------------------------------
c
c  Loop over all gases.
c
      do igas = 1, NGAS ! Calculate liq/ice vapor pressure for pure components
c
c  Calculate vapor pressures based on which gas(es) are defined in <setupaer>
c
        if( gasname(igas) .eq. 'methane' ) then
c
c  Saturation vapor pressure of methane over its liquid
c  (from Moses et al., Icarus, 99, 318-346, (1992) )
c
          pvapl3(ixyz,igas) = vplA_CH4 - vplB_CH4/t3(ixyz) + 
     $            vplC_CH4/t3(ixyz)**2 - vplD_CH4/t3(ixyz)**3
 
          pvapl3(ixyz,igas) = Ratm2cgs * 10**pvapl3(ixyz,igas)
c  Possibly adjust vapor pressure if simulating mixing with N2
c  (comment this out when not running CH4+N2 test)
          if(zl3(ixyz) .lt. 40.d5) 
     $     pvapl3(ixyz,igas) = VP_CH4_adjust(ixyz)*pvapl3(ixyz,igas)
c
c  Saturation vapor pressure of methane over its ice
c  (from Moses et al., Icarus, 99, 318-346, (1992) )
c 
          pvapi3(ixyz,igas) = vpiA_CH4 - vpiB_CH4/t3(ixyz) - 
     $                        vpiC_CH4/t3(ixyz)**2 +  
     $           vpiD_CH4/t3(ixyz)**3 - vpiE_CH4/t3(ixyz)**4
 
          pvapi3(ixyz,igas) = Ratm2cgs * 10**pvapi3(ixyz,igas)
 
        else if( gasname(igas) .eq. 'ethane') then
c
c  Saturation vapor pressure of ethane over its ice
c  (from Moses et al., Icarus, 99, 318-346, (1992) )
 
           pvapi3(ixyz,igas) = vpiA_C2H6 - vpiB_C2H6/
     $                                  (t3(ixyz) - vpiC_C2H6)
     $                                  
 
           pvapi3(ixyz,igas) = RmmHg2cgs * 10**pvapi3(ixyz,igas)
c
c  Saturation vapor pressure over liquid ethane
c  (from Handbook of Vapor Pressure, vol.1, 1994)
c
           pvapl3(ixyz,igas) = vplA_C2H6 + vplB_C2H6/t3(ixyz) +
     $             vplC_C2H6*dlog10(t3(ixyz)) + vplD_C2H6*t3(ixyz) +
     $                          vplE_C2H6*t3(ixyz)**2

           pvapl3(ixyz,igas) = RmmHg2cgs * 10**pvapl3(ixyz,igas)
c
c  Saturation vapor pressure over liquid nitrogen
c  (from Graves et al. 2007+)
c
        else if( gasname(igas) .eq. 'nitrogen') then
           pvapl3(ixyz,igas) = vplA_N2 * 
     $                          10**(vplB_N2 - vplC_N2/t3(ixyz))
           pvapl3(ixyz,igas) = RPa2cgs * 10**pvapl3(ixyz,igas)
c
c  Report unexpected gas index
        else
            write(LUNOPRT,'(/,a,a)') 'invalid <igas> in vaporp.f',
     $                             gasname(igas)
            call endcarma
        endif

        pvl_pure(igas) = pvapl3(ixyz,igas)
 
      enddo
cc
cc Now loop over gases again to adjust liquid vapor pressures for mixing
cc
cc    do igas = 1, NGAS  ! Calculate activity coefficients
cc      gam(igas) = 0.
cc      do i = 1, NGAS
cc        if(i.ne.igas) gam(igas) = gam(igas) + 
cc   $                    pvap_alpha(igas,i)*xfrac(i)**pvap_m(igas,i)
cc      enddo
cc      gam(igas) = exp( gam(igas)/t3(ixzy) )
cc    enddo
cc
cc    do igas = 1, NGAS
cc
cc Possibly adjust vapor pressure if simulating mixing with N2
cc (comment this out when not running CH4+N2 test)
cc        if(zl3(ixyz) .lt. 40.d5) 
cc   $     pvapl3(ixyz,igas) = VP_CH4_adjust(ixyz)*pvapl3(ixyz,igas)
cc
cc      pvapl3(ixyz,igas) = gam(igas)*xfrac(igas) * pvl_pure(igas)
cc      do i = 1, NGAS
cc       if( i .ne. igas ) pvapl3(ixyz,igas) = pvapl3(ixyz,igas) + 
cc   $                          gam(i)*xfrac(i) * pvl_pure(i)
cc      enddo
cc    enddo
c
c
c  Return to caller with vapor pressures evaluated.
c
      return
      end
