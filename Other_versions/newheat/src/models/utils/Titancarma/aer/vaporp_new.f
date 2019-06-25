
c    Fugacity coefficients
      data phi_N2/   0.895, 0.907, 0.915, 0.924, 0.932, 0.938, 0.944, 
     $ 0.949, 0.954, 0.958, 0.962, 0.965, 0.967, 0.969, 0.971, 0.973, 
     $ 0.974, 0.975, 0.976, 0.976, 0.000, 0.000, 0.000, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /

      data phi_CH4/  0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.996, 
     $ 0.997, 0.997, 0.998, 0.998, 0.998, 0.999, 0.999, 0.999, 0.999, 
     $ 0.999, 0.999, 0.999, 0.999, 0.000, 0.000, 0.000, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /

c    Saturation mole fractions
      data X_N2/     0.179, 0.204, 0.214, 0.232, 0.250, 0.264, 0.272, 
     $ 0.282, 0.292, 0.300, 0.295, 0.289, 0.286, 0.274, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /

      data X_CH4/    0.092, 0.073, 0.064, 0.054, 0.046, 0.040, 0.035, 
     $ 0.031, 0.028, 0.024, 0.023, 0.022, 0.020, 0.020, 0.020, 0.020, 
     $ 0.020, 0.020, 0.020, 0.020, 0.000, 0.000, 0.000, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 
     $ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /

c      Critical and triple point values for N2 and CH4
       data Pc/ 34.002, 45.955 /   ! bar
       data Pt/ 0.1252, 0.1170 /   ! bar
       data Tc/ 126.20, 190.53 /   ! Kelvin
       data Tt/ 63.15,   90.68 /   ! Kelvin

       data a4IS/ 3.065972, 3.159023 /
       data b0IS/ -23.52451, -19.36816 /
       data b1IS/ 6103.604, 8799.140 /

c  Saturation vapor pressure for a binary mixture of N2 and CH4
c  (from Thompson et al., Icarus, 97, 187-199, (1992) )

c  First find vapor pressure for each pure gas

c    icomp: 1 = N2, 2 = CH4
       do icomp = 1,2     

          a0 = ONE - Pt(icomp) / (Pc(icomp) - Pt(icomp) )
          a2 = b1IS / RGAS / Tt(icomp)
          a1 = -(a0 - ONE) * exp(a2 - b0IS/RGAS )
          a3 = (Tc(icomp) - Tt(icomp)) / Tt(icomp)

c      Dimensionless temperature
          tdim = (t3(ixyz) - Tt(icomp) / (Tc(icomp) - Tt(icomp)

c      Asymptotic form of vapor pressure near the triple point
          p_tp = a0 + a1*(a3*tdim(ixyz) + ONE)**(b0IS/RGAS) *
     $            exp( (-a2 + b0IS/RGAS)/(a3*tdim(ixyz) + ONE) )


          a5 = -0.11599104d0 + 0.29506258d0*a4IS(icomp)**2 - 
     $              0.00021222d0*a4IS(icomp)**5 
          a6 = -0.01546028d0 + 0.08978160d0*a4IS(icomp)**2 - 
     $              0.05322199d0*a4IS(icomp)**3
          a7 =  0.05725757d0 - 0.06817687d0*a4IS(icomp)    + 
     $              0.00047188d0*a4IS(icomp)**5

c      Asymptotic form of vapor pressure near the critical point
          p_cr = 2.d0 - a4IS(icomp)*(ONE - tdim(ixyz)) + 
     $            a5*(ONE - tdim(ixyz))**(1.8d0) +
     $            a6*(ONE - tdim(ixyz))**3 + a7*(ONE - tdim(ixyz))**4

          rN = 87.d0 * Tt(icomp) / Tc(icomp)

c      Dimensionless pressure
          pdim = ( p_tp**rN + p_cr**rN )**(ONE/rN) - ONE

c         Saturation vapor pressure of pure gas 
c         from Iglesias-Silva et al., AIChE Journal, 33, 1550-1556, (1987) )

          psatl(icomp) = pdim * (Pc(icomp) - Pt(icomp)) 
     $                           + Pt(icomp)

       enddo !icomp

c  Now combine the CH4 and N2 vapor pressures found above

          a0 = 0.8096
          a1 = -52.07
          a2 = 5443.0
          b0 = -0.0829
          b1 = 9.34
          c0 = 0.0720
          c1 = -6.27

c      Temperature dependence terms
          a = a0 + a1*t3(ixyz)**-1 + a2*t3(ixyz)**-2
          b = b0 + b1*t3(ixyz)**-1
          c = c0 + c1*t3(ixyz)**-1

c  Vapor pressure equation for the mixture
          pvapl3(ixyz,igas) = (X_N2 / phi_N2)*psatl(1) *
     $                      exp(X_CH4**2 * ( (a + 3.d0*b + 5.d0*c)
     $                             - 4.d0*(b + 4.d0*c)*X_CH4**2 
     $                             + 12.d0*c*X_CH4**2 ) )
     $                          +
     $                        (X_CH4 / phi_CH4)*psatl(2) *
     $                      exp(X_N2**2 * ( (a - 3.d0*b + 5.d0*c)
     $                             + 4.d0*(b + 4.d0*c)*X_N2**2 
     $                             + 12.d0*c*X_N2**2 ) )
