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
c
c  Loop over all gases.
c
      do igas = 1, NGAS
c
c
c  Calculate vapor pressures based on which gas(es) are defined in <setupaer>
c
        if( gasname(igas) .eq. 'methane' ) then
           ! add CH4 + N2

c
c  Saturation vapor pressure of methane over its liquid
c  (from Moses et al., Icarus, 99, 318-346, (1992) )
c
          pvapl3(ixyz,igas) = vplA_CH4 - vplB_CH4/t3(ixyz) + 
     $            vplC_CH4/t3(ixyz)**2 - vplD_CH4/t3(ixyz)**3
 
          pvapl3(ixyz,igas) = Ratm2cgs * 10**pvapl3(ixyz,igas)
c
c  Saturation vapor pressure of methane over its ice
c  (from Moses et al., Icarus, 99, 318-346, (1992) )
c 
          pvapi3(ixyz,igas) = vpiA_CH4 - vpiB_CH4/t3(ixyz) - 
     $                        vpiC_CH4/t3(ixyz)**2 +  
     $           vpiD_CH4/t3(ixyz)**3 - vpiE_CH4/t3(ixyz)**4
 
          pvapi3(ixyz,igas) = Ratm2cgs * 10**pvapi3(ixyz,igas)
c 
c
c
          if( gasname(igas) .eq. 'CH4 + N2') then
c
c  Saturation vapor pressure for a binary mixture of N2 and CH4
c  (from Thompson et al., Icarus, 97, 187-199, (1992) )

           if(X_CH4(ixyz) .gt. 0.) then
 
            Rmks = RGAS / RJ2cgs

c  First find vapor pressure for each pure gas

c    j: 1 = N2, 2 = CH4
            do j = 1,2     

              a0 = ONE - Ptrip(j) / (Pcrit(j) - Ptrip(j) )
              a2 = b1IS(j) / Rmks / Ttrip(j)
              a1 = -(a0 - ONE) * exp(a2 - b0IS(j)/Rmks )
              a3 = (Tcrit(j) - Ttrip(j)) / Ttrip(j)

c      Dimensionless temperature
              Tdim = (t3(ixyz) - Ttrip(j)) / (Tcrit(j) - Ttrip(j))

c      Asymptotic form of vapor pressure near the triple point
              p_tp = a0 + a1*(a3*Tdim + ONE)**(b0IS(j)/Rmks) *
     $              exp( (-a2 + b0IS(j)/Rmks)/(a3*Tdim + ONE) )


              a5 = -0.11599104d0 + 0.29506258d0*a4IS(j)**2 - 
     $              0.00021222d0*a4IS(j)**5 
              a6 = -0.01546028d0 + 0.08978160d0*a4IS(j)**2 - 
     $              0.05322199d0*a4IS(j)**3
              a7 =  0.05725757d0 - 0.06817687d0*a4IS(j)    + 
     $              0.00047188d0*a4IS(j)**5

c      Asymptotic form of vapor pressure near the critical point
              p_cr = 2.d0 - a4IS(j)*(ONE - Tdim) + 
     $            a5*(ONE - Tdim)**(1.8d0) +
     $            a6*(ONE - Tdim)**3 + a7*(ONE - Tdim)**4

              rN = 87.d0 * Ttrip(j) / Tcrit(j)

c      Dimensionless pressure
              pdim = ( p_tp**rN + p_cr**rN )**(ONE/rN) - ONE

c         Saturation vapor pressure of pure gas 
c         from Iglesias-Silva et al., AIChE Journal, 33, 1550-1556, (1987) )

              psatl(j) = pdim * (Pcrit(j) - Ptrip(j)) 
     $                           + Ptrip(j)
c             print *,zl3(ixyz),j,psatl(j)

            enddo !pure gas component

c  Now combine the CH4 and N2 vapor pressures found above

c      Temperature dependence terms
            a = TZGa0 + TZGa1*t3(ixyz)**-1 + TZGa2*t3(ixyz)**-2
            b = TZGb0 + TZGb1*t3(ixyz)**-1
            c = TZGc0 + TZGc1*t3(ixyz)**-1

c  Vapor pressure equation for the mixture
            pvapl3(ixyz,igas) = (X_N2(ixyz) / phi_N2(ixyz))*psatl(1) *
     $                      exp(X_CH4(ixyz)**2 * ( (a + 3.d0*b + 5.d0*c)
     $                             - 4.d0*(b + 4.d0*c)*X_CH4(ixyz)**2 
     $                             + 12.d0*c*X_CH4(ixyz)**2 ) )
     $                          +
     $                        (X_CH4(ixyz) / phi_CH4(ixyz))*psatl(2) *
     $                      exp(X_N2(ixyz)**2 * ( (a - 3.d0*b + 5.d0*c)
     $                             + 4.d0*(b + 4.d0*c)*X_N2(ixyz)**2 
     $                             + 12.d0*c*X_N2(ixyz)**2 ) )

            pvapl3(ixyz,igas) = RMB2CGS**2 * pvapl3(ixyz,igas) 


          endif !nonzero methane mole fraction
         endif !gasname = CH4 + N2
c
c
c
        else if( gasname(igas) .eq. 'ethane') then
c
c
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
     $                          vplC_C2H6*dlog10(t3(ixyz)) + vplD_C2H6*t3(ixyz) +
     $                          vplE_C2H6*t3(ixyz)**2

           pvapl3(ixyz,igas) = RmmHg2cgs * 10**pvapl3(ixyz,igas)
c
c
c
c
c
c  Report unexpected gas index
c
        else
            write(LUNOPRT,'(/,a)') 'invalid <igas> in vaporp.f'
            stop 1
        endif
 
      enddo
c
c
c  Return to caller with vapor pressures evaluated.
c
      return
      end
