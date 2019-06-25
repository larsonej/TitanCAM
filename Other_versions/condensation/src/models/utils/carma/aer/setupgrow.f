       subroutine setupgrow
c
c
c  @(#) setupgrow.f  Ackerman  Dec-1995
c  This routine defines time-independent parameters used to calculate
c  condensational growth/evaporation.
c
c  The parameters defined for each gas are 
c
c    gwtmol:   molecular weight [g/mol]
c    diffus:   diffusivity      [cm^2/s]
c    rlhe  :   latent heat of evaporation [cm^2/s^2]
c    rlhm  :   latent heat of melting [cm^2/s^2]
c
c  Time-independent parameters that depend on particle radius are
c  defined in setupgkern.f.
c
c  This routine requires that vertical profiles of temperature <T>,
c  and pressure <p> are defined (i.e., initatm.f must be called before this).
c  The vertical profile with ix = iy = 1 is used.
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
c  Include constants for vapor pressure equations
c
      include 'vpconstants.h'
c
c
c  Define formats
c
    1 format(a,':  ',12i6)
    2 format(a,':  ',i6)
    3 format(/' id  gwtmol   gasname',(/,i3,3x,f5.1,3x,a))
    5 format(/,'Particle growth mapping arrays (setupgrow):')
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupgrow'
c
c-----Check that values are valid------------------------------------------
c
      do ielem = 1, NELEM

        if( ielem .gt. 1 .and.
     $      itype(ielem) .eq. I_CORE2MOM )then 

          if( igrowgas(ielem) .ne. igrowgas(ielem-1) )then
            write(LUNOPRT,'(/,a)') 'igrowgas array is inconsistent'
            call endcarma
          endif

        endif

        if( igrowgas(ielem) .gt. NGAS )then
          write(LUNOPRT,'(/,a)') 'component of igrowgas > NGAS'
          call endcarma
        endif

      enddo
c
c--------------------------------------------------------------------------
c  
c  Define parameters with weak time-dependence to be used in
c  growth equation.
c
      do k = 1, NZ
        do igas = 1, NGAS


          if( gasname(igas) .eq. 'meth' .or. 
     $        gasname(igas) .eq. 'CH4 + N2' )   then  

c
c   Diffusivity at reference temperature and pressure from
c   Lorenz   Planet. Space Sci.  v.41, p.647-655 (1993)
c   units are [cm^2/s].
c
            diffus(k,igas) =  0.196 * ( p(1,1,k) / PREF )
     $              * ( t(1,1,k) / 273. )**1.75
c
c
c   Latent heat of evaporation for methane (from vapor pressure eqn)
c
            rlhe(k,igas) = RGAS/gwtmol(igas) *
     $                 (vplB_CH4 - 2.d0*vplC_CH4/t3(k) + 
     $                     3.d0*vplD_CH4/t3(k)**2)*dlog(1.d1)
c
c    Latent heat of sublimation (from vapor pressure eqn)
c
             rlhs = RGAS/gwtmol(igas) * 
     $                   (vpiB_CH4 + 2.d0*vpiC_CH4/t3(k) - 
     $                       3.d0*vpiD_CH4/t3(k)**2 + 
     $                           4.d0*vpiE_CH4/t3(k)**3)*dlog(1.d1)
c
c
c
            if( gasname(igas) .eq. 'CH4 + N2' ) then
              if( X_CH4(k) .gt. 0. ) then
 
c             Excess enthalpy (TZG Eqn. 8)
                aprime = TZGa1/t3(k) + 2.d0*TZGa2/t3(k)**2
                bprime = TZGb1/t3(k)
                cprime = TZGc1/t3(k)

                Hexcess = RGAS * t3(k) * X_N2(k) * X_CH4(k) *
     $                   ( aprime + bprime*(X_N2(k) - X_CH4(k)) +
     $                     cprime*(X_N2(k) - X_CH4(k))**2 )

                Hc_N2(k)  = RJ2cgs * Hc_N2(k)
                Hc_CH4(k) = RJ2cgs * Hc_CH4(k)

                deltaHc = X_N2(k)*Hc_N2(k) + X_CH4(k)*Hc_CH4(k) 
     $                      + Hexcess

                rlhe(k,igas) = deltaHc * gwtmol(igas)

              endif
            endif
c
c
c
          else if( gasname(igas) .eq. 'eth' ) then  
c
c
c
c   Diffusivity at reference temperature and pressure from
c   Boyd et al. J.Chem.Phys.  v.19, p.548 (1951)
c   units are [cm^2/s].
c
c
            diffus(k,igas) =  0.148 * ( 1.01325e+6 / p(1,1,k) )
     $              * ( t(1,1,k) / 298.15 )**1.94
 
c
c   Latent heat of evaporation for ethane (from vapor pressure eqn)
c
            rlhe(k,igas) = RGAS/gwtmol(igas) *
     $               ( vplC_C2H6*t(1,1,k) + (-vplB_C2H6 + 
     $                        vplD_C2H6*t(1,1,k)**2 +
     $                        2.*vplE_C2H6*t(1,1,k)**3)*dlog(1.d1) )


c    Latent heat of sublimation (from vapor pressure eqn)
            rlhs = vpiB_C2H6*dlog(1.d1)*RGAS/gwtmol(igas)*t3(k)**2 /
     $                               (t3(k) - vpiC_C2H6)**2

c
c
c  Report unexpected gas index
c
          else
            write(LUNOPRT,'(/,a)') 'invalid <igas> in vaporp.f'
            call endcarma
          endif
c
c
c   Latent heat of ice melting 
c   units are [cm^2/s^2]
c
          rlhm(k,igas) =  rlhs - rlhe(k,igas)
 
        enddo
      enddo

c Set vapor pressure adjustment factors for CH4+N2
c  Coefficients for polynomial fit to methane mole fractions
c   (gamma_ch4 * (1.-Xn2)  from Thompson92 table IV )
c  or ice vapor pressure above 30 km
      a0 = 0.76903175
      a1 = -8.1721842e-06
      a2 = 1.9970681e-08
      a3 = 1.7292750e-11
      do k=1,20
         y = p3(k)/1000. !Need pressure in mbar
         VP_CH4_adjust(k) = a0 + a1*y + a2*y**2 + a3*y**3
c EJL 4-8-13 comment out this write
c         write(*,*) k,VP_CH4_adjust(k)
      enddo
c
c
c--------------------------------------------------------------------------
c
c
c  Report some initialization values
c
      if (do_print_setup) then
      write(LUNOPRT,5)
      write(LUNOPRT,2) 'NGAS    ',NGAS
      write(LUNOPRT,1) 'igrowgas',(igrowgas(i),i=1,NELEM)
      write(LUNOPRT,3) (i,gwtmol(i),gasname(i),i=1,NGAS)
      endif
c
c
c  Return to caller with particle growth mapping arrays and time-dependent
c  parameters initialized.
c
      return
      end
